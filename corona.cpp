#include <iostream>
#include <fftw3.h>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstring>
// #include "poisson_solver.h"
// #include "poisson_solver_3d.h"
#include "helper.h"
#include "initializer.h"
#include "method.h"

using namespace std;

int main(int argc, char **argv)
{

    if(argc!=10){
        cout << "Need to do the following : " << endl;
        cout << "./main.exe [n1] [n2] [nt] [tau] [sigma] [tolerance] [max iteration] [eta] [skip]" << endl;
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

    double base=0;

    double dx=1.0/n1;
    double dy=1.0/n2;
    double dt=1.0/nt;

    // coefficients for the energy functions

    double c0=0.001;
    double c1=1.0;
    double c2=0.01;

    // coefficients for velocities

    double alpha1 = 1.0;
    double alpha2 = 2.0;
    double alpha3 = 1.0;

    // coefficients for SIR model
    double beta  = 0.9;
    double gamma = 0.05;

    double var = 0.1;

    double* rho[3];

    for(int i=0;i<3;++i) rho[i] = new double[n1*n2*nt];

    double** f_arr = new double*[3]; // f is an obstacle

    for(int k=0;k<3;++k) f_arr[k] = new double[n1*n2];

    double Clist[] = {eta, eta, eta};  

    /*
        Initialize rho0, rho1, rho2
    */

    Initializer init(n1,n2,nt,dx,dy,dt,base);

    init.intialize_rho0(rho[0]);
    init.intialize_rho1(rho[1]);
    init.intialize_rho2(rho[2]);
    // init.renormalize_all(rho);
    // init.intialize_f(f);
    for(int k=0;k<3;++k) for(int i=0;i<n1*n2;++i) f_arr[k][i] = 0;

    // for(int i=0;i<n2;++i){
    //     for(int j=0;j<n1;++j){
    //         double x = (j+0.5)/n1;
    //         double y = (i+0.5)/n2;

    //         f_arr[1][i*n1+j] = 10 * ((x-0.9)*(x-0.9) + (y-0.9)*(y-0.9));
    //     }
    // }

    // initialize the method
    Method method(n1, n2, nt, dx, dy, dt, tau, sigma, max_iteration, tolerance, c0, c1, c2, alpha1, alpha2, alpha3, Clist, var);

    // coefficients for SIR system
    method.h = 1;
    method.m = 2;
    method.beta  = beta;
    method.gamma = gamma;

    cout<<"XXX G-Prox PDHG XXX"<<endl;
    cout<<endl;
    cout<<"n1 : "<<n1<<" n2 : "<<n2<<" nt : "<<nt<<" base : "<<base<<endl;
    cout<<"dx : "<<scientific<<dx<<endl;
    cout<<"dy : "<<scientific<<dy<<endl;
    cout<<"dt : "<<scientific<<dt<<endl;
    cout<<fixed;
    cout<<"tau           : "<<tau<<endl;
    cout<<"sigma         : "<<sigma<<endl;
    cout<<"max_iteration : "<<max_iteration<<endl;
    cout<<"tolerance     : "<<scientific<<tolerance<<endl;
    cout<<"eta           : "<<scientific<<eta<<endl;
    cout<<"beta          : "<<scientific<<method.beta<<endl;
    cout<<"gamma         : "<<scientific<<method.gamma<<endl;
    cout<<"var           : "<<scientific<<method.var<<endl;

    create_csv_file_for_parameters(n1,n2,nt,c0,c1,c2,method.beta,method.gamma,var);

    cout<<"\nXXX Starting Iterations XXX"<<endl;

    clock_t t;
    t = clock();

    string filename;    

    filename="./data/rho0-"+to_string(0)+".csv";
    create_bin_file(&rho[0][(nt-1)*n1*n2],n1*n2,filename);
    filename="./data/rho1-"+to_string(0)+".csv";
    create_bin_file(&rho[1][(nt-1)*n1*n2],n1*n2,filename);
    filename="./data/rho2-"+to_string(0)+".csv";
    create_bin_file(&rho[2][(nt-1)*n1*n2],n1*n2,filename);

    for(int iter=0; iter<1; ++iter)
    {
        method.run(rho,f_arr,skip);   

        filename="./data/rho0-"+to_string(iter+1)+".csv";
        create_bin_file(&rho[0][(nt-1)*n1*n2],n1*n2,filename);
        filename="./data/rho1-"+to_string(iter+1)+".csv";
        create_bin_file(&rho[1][(nt-1)*n1*n2],n1*n2,filename);
        filename="./data/rho2-"+to_string(iter+1)+".csv";
        create_bin_file(&rho[2][(nt-1)*n1*n2],n1*n2,filename);

        for(int n=0;n<nt-1;++n){
            for(int i=0;i<n1*n2;++i){
                rho[0][n*n1*n2+i] = rho[0][(nt-1)*n1*n2+i];
                rho[1][n*n1*n2+i] = rho[1][(nt-1)*n1*n2+i];
                rho[2][n*n1*n2+i] = rho[2][(nt-1)*n1*n2+i];
            }
        }
    }
    

    t = clock() - t;

    printf ("CPU time for Iterations: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);

    for(int k=0;k<3;++k){
        delete[] rho[k];
        delete[] f_arr[k];
    }
    delete[] f_arr;
}