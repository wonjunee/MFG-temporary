#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstring>
// #include "poisson_solver.h"
#include "poisson_solver_3d.h"
#include "helper.h"
#include "method.h"

using namespace std;

void initialize_rho(double* rho, int n1, int n2, int nt, double base=0){
    int n = 0;
    double sum0 = 0, sum1 = 0;
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            double x = (j+0.5)/n1;
            double y = (i+0.5)/n2;

            n = 0;
            rho[n*n1*n2+i*n1+j] = exp(-40*(pow(x-0.3,2)+pow(y-0.3,2))) + base;
            sum0 += rho[n*n1*n2+i*n1+j];

            n = nt-1;
            rho[n*n1*n2+i*n1+j] = exp(-40*(pow(x-0.7,2)+pow(y-0.7,2))) + base;
            sum1 += rho[n*n1*n2+i*n1+j];
        }
    }
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            n = 0;
            rho[n*n1*n2+i*n1+j] *= n1*n2/sum0;
            n = nt-1;
            rho[n*n1*n2+i*n1+j] *= n1*n2/sum1;
        }
    }
    for(int n=1;n<nt-1;++n){
        for(int i=0;i<n1*n2;++i){
            rho[n*n1*n2+i] = rho[i];
        }
    }
}

int main(int argc, char **argv)
{

    if(argc!=9){
        cout << "Need to do the following : " << endl;
        cout << argv[0] <<  " [n1] [n2] [nt] [tau] [sigma] [tolerance] [max iteration] [skip]" << endl;
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
    int skip=stoi(argv[8]);
    
    double base=0.01;

    double* rho = new double[n1*n2*nt];

    // Initialize rho

    initialize_rho(rho,n1,n2,nt,base);

    // initialize the method

    Method method(n1, n2, nt, tau, sigma, max_iteration, tolerance);

    cout<<"XXX G-Prox PDHG XXX"<<endl;
    cout<<endl;
    cout<<"n1 : "<<n1<<" n2 : "<<n2<<" nt : "<<nt<<" base : "<<base<<endl;
    cout<<fixed;
    cout<<"tau           : "<<tau<<endl;
    cout<<"sigma         : "<<sigma<<endl;
    cout<<"max_iteration : "<<max_iteration<<endl;
    cout<<"tolerance     : "<<scientific<<tolerance<<endl;

    create_csv_file_for_parameters(n1,n2,nt);

    cout<<"\nXXX Starting Iterations XXX"<<endl;

    clock_t t;
    t = clock();

    method.run(rho,skip);   

    t = clock() - t;

    printf ("CPU time for Iterations: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);

    delete[] rho;
}