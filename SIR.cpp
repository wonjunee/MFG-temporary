#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstring>
// #include "poisson_solver.h"
#include "poisson_solver_3d.h"
#include "helper.h"
#include "initializer.h"
#include "method.h"

using namespace std;

// This is SIR model without moving
class SIR_model{
public:
    int n1;
    int n2;
    int nt;
    double dx;
    double dy;
    double dt;

    double* tau;
    double* sigma;
    double h;

    int max_iteration;
    double tolerance;

    double* mx[3];
    double* my[3];
    double* phi[3];
    double* phitmps[3];
    double* rhotmps[3];

    double* convarr;

    double* phiT;
    double* phiTtmp; 

    double* xi;
    double* xitmp;
    double TOTAL_SUM;

    double* karr;

    double energy;
    double previous_energy;
    double previous_dual;

    double m;

    double M0;

    double c0;
    double c1;
    double c2;

    double* etalist;
    double* alphalist;

    int convN;
    double conv_sum;
    double conv_r;

    // For SIR model

    double beta;
    double gamma;
    double var;
    double var_gamma;

    // ------------------------

    SIR_model(){
    }

    SIR_model(int n1, int n2, int nt, double dx, double dy, double dt, 
           double tau, double sigma, int max_iteration, double tolerance, 
           double c0, double c1, double c2, 
           double alpha1, double alpha2, double alpha3,
           double* etalist, double var){

        this->n1=n1;
        this->n2=n2;
        this->nt=nt;
        this->dx=dx;
        this->dy=dy;
        this->dt=dt;
        
        this->max_iteration=max_iteration;
        this->tolerance=tolerance;

        this->tau  =new double[3];
        this->sigma=new double[3];

        this->tau[0] = tau;
        this->tau[1] = tau;
        this->tau[2] = tau;

        this->sigma[0] = sigma;
        this->sigma[1] = sigma;
        this->sigma[2] = sigma;

        this->c0=c0;
        this->c1=c1;
        this->c2=c2;

        this->etalist = etalist; // this contains the eta values for viscosity terms
        this->var     = var;

        alphalist = new double[3];
        alphalist[0] = alpha1;
        alphalist[1] = alpha2;
        alphalist[2] = alpha3;

        karr = new double[3];
        karr[0] = 1e-2;
        karr[1] = 1e-2;
        karr[2] = 1e-2;

        xi    = new double[nt];
        xitmp = new double[nt];

        phiT     = new double[n1*n2];
        phiTtmp  = new double[n1*n2];

        convarr  = new double[n1*n2];

        M0 = 0.5;
        conv_r = 0.2;
        convN  = conv_r*n1;
        convN  = 2 * convN + 1;
        // convN  = 1;
        cout << "convN: " <<convN << endl;

        var_gamma = 0.01; // DON'T TOUCH THIS. This one is for regularization.
    }


    void run(double* rho[], double* const* f, int skip=1){

        printf("beta: %f\tgamma: %f", beta, gamma);

        for(int n=0;n<nt-1;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    rho[0][(n+1)*n1*n2+i*n1+j] = rho[0][(n)*n1*n2+i*n1+j] + 1.0/nt * ( - beta * rho[0][n*n1*n2+i*n1+j] * rho[1][n*n1*n2+i*n1+j]);
                    rho[1][(n+1)*n1*n2+i*n1+j] = rho[1][(n)*n1*n2+i*n1+j] + 1.0/nt * ( beta * rho[0][n*n1*n2+i*n1+j] * rho[1][n*n1*n2+i*n1+j] - gamma * rho[1][n*n1*n2+i*n1+j]);
                    rho[2][(n+1)*n1*n2+i*n1+j] = rho[2][(n)*n1*n2+i*n1+j] + 1.0/nt * (gamma * rho[1][n*n1*n2+i*n1+j]);
                }
            }
        }

        create_csv_file(rho[0],"./data/rho0.csv",n1,n2,nt);
        create_csv_file(rho[1],"./data/rho1.csv",n1,n2,nt);
        create_csv_file(rho[2],"./data/rho2.csv",n1,n2,nt);

        cout<<"The method is done!!"<<endl;
    }

}; // Method class

void remove_mass_obstacle(double** rho,double* obstacle, int n1, int n2, int nt){
    for(int n=0;n<nt;++n){
        for(int i=0;i<n1*n2;++i){
            if(obstacle[i] > 0){
                for(int k=0;k<4;++k) rho[k][n*n1*n2+i] = 0;
            }
        }
    }
}

int main(int argc, char **argv)
{

    if(argc!=12){
        cout << "Need to do the following : " << endl;
        cout << argv[0] <<  " [n1] [n2] [nt] [tau] [sigma] [tolerance] [max iteration] [eta] [skip] [beta] [gamma]" << endl;
        exit(1);
    }

    // Parameters for Grids
    
    int n1 = stoi(argv[1]);  // x panels
    int n2 = stoi(argv[2]);  // y panels
    int nt = stoi(argv[3]);  // dt
    double tau  =stod(argv[4]);
    double sigma=stod(argv[5]);
    double tolerance = stod(argv[6]);
    int max_iteration=stoi(argv[7]);
    double eta = stod(argv[8]);
    int skip=stoi(argv[9]);
    double beta  = stod(argv[10]);
    double gamma = stod(argv[11]);

    double base=0;

    // coefficients for velocities
    double alphalist[4] = {10.0, 10.0, 10.0, 0.005};
    // coefficients for SIR model
    double var = 0.005;


    double* rho[4];

    for(int i=0;i<4;++i) rho[i] = new double[n1*n2*nt];

    double* obstacle = new double[n1*n2];
    double etalist[] = {eta, eta, eta};  

    /*
        Initialize rho0, rho1, rho2
    */

    Initializer init(n1,n2,nt,base);

    init.intialize_rho0(rho[0]);
    init.intialize_rho1(rho[1]);
    init.intialize_rho2(rho[2]);
    init.intialize_rho3(rho[3]);

    for(int i=0;i<n1*n2;++i) obstacle[i] = 0;

    // for(int i=0;i<n2;++i){
    //     for(int j=0;j<n1;++j){
    //         double x = (j+0.5)/n1;
    //         double y = (i+0.5)/n2;
    //         if(fabs(y-0.2) < 0.125 && fabs(x-0.35)<0.02) obstacle[i*n1+j] = 99999999;
    //         if(fabs(y-0.2) < 0.125 && fabs(x-0.65)<0.02) obstacle[i*n1+j] = 99999999;
    //         if(fabs(y-0.8) < 0.125 && fabs(x-0.35)<0.02) obstacle[i*n1+j] = 99999999;
    //         if(fabs(y-0.8) < 0.125 && fabs(x-0.65)<0.02) obstacle[i*n1+j] = 99999999;
    //     }
    // }

    create_bin_file(obstacle, n1*n2, "./data/obstacle.csv");

    remove_mass_obstacle(rho,obstacle,n1,n2,nt);

    // initialize the method
    Method method(n1, n2, nt, tau, sigma, max_iteration, tolerance, alphalist, etalist, var);
    // SIR_model method(n1, n2, nt, dx, dy, dt, tau, sigma, max_iteration, tolerance, c0, c1, c2, alpha1, alpha2, alpha3, etalist, var);

    // coefficients for SIR system
    method.beta  = beta;
    method.gamma = gamma;

    cout<<"XXX G-Prox PDHG XXX"<<endl;
    cout<<endl;
    cout<<"n1 : "<<n1<<" n2 : "<<n2<<" nt : "<<nt<<" base : "<<base<<endl;
    cout<<fixed;
    cout<<"tau           : "<<tau<<endl;
    cout<<"sigma         : "<<sigma<<endl;
    cout<<"max_iteration : "<<max_iteration<<endl;
    cout<<"tolerance     : "<<scientific<<tolerance<<endl;
    cout<<"eta           : "<<scientific<<eta<<endl;
    cout<<"beta          : "<<scientific<<method.beta<<endl;
    cout<<"gamma         : "<<scientific<<method.gamma<<endl;
    cout<<"var           : "<<scientific<<method.var<<endl;

    create_csv_file_for_parameters(n1,n2,nt,0,0,0,method.beta,method.gamma,var);

    cout<<"\nXXX Starting Iterations XXX"<<endl;

    clock_t t;
    t = clock();

    method.run(rho,obstacle,skip);   

    t = clock() - t;

    printf ("CPU time for Iterations: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);

    for(int k=0;k<4;++k){
        delete[] rho[k];
    }
    delete[] obstacle;
}