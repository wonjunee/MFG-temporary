#include <iostream>
#include <fftw3.h>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <time.h>
#include <cstring>
// #include "poisson_solver.h"
#include "poisson_solver_3d.h"
#include "helper.h"

using namespace std;

void create_bin_file(const double* A, int size, string filename){
    ofstream out(filename, ios::out | ios::binary);
    if(!out) {
        cout << "Cannot open file.";
        return;
    }

    out.write((char *) A, size*sizeof(double));
    out.close();
}

void create_csv_file_for_parameters(int n1,int n2,int nt){
    ofstream outfile;
    outfile.open("./data/parameters.csv");
    outfile<<n1<<","<<n2<<","<<nt;
    outfile.close();
}

void create_csv_file(const double* A,string filename,int n1,int n2,int nt){
    ofstream outfile;
    outfile.open(filename);
    for(int i=0;i<n1*n2*nt;++i){
        outfile<<A[i]<<"\n";
    }
    outfile.close();
}

double gaussian(const double x, const double y, const double mux, const double muy, const double sigmax, const double sigmay){
    return exp(-0.5*(pow(x-mux,2)/pow(sigmax,2) + pow(y-muy,2)/pow(sigmay,2) ));
}

void intialize_rho1(double* rho,double base,int n1,int n2,int nt,double dx,double dy,double dt){

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){

            double x=1.0*(j+0.5)*dx;
            double y=1.0*(i+0.5)*dy;

            // if(pow(x-0.5,2) + pow(y-0.5,2) >= pow(0.1,2)){
            //     rho[i*n1+j] = 1;    
            // }else{
            //     rho[i*n1+j] = 0;
            // }
            rho[i*n1+j] = 2*exp(-20*pow(x-0.5,2)-20*pow(y-0.5,2));
        }
    }

    for(int n=1;n<nt;++n){
        for(int i=0;i<n1*n2;++i){
            rho[n*n1*n2+i]=rho[i];
        }
    }
}

void intialize_rho2(double* rho,double base,int n1,int n2,int nt,double dx,double dy,double dt){

    double sum=0;

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){

            double x=1.0*(j+0.5)*dx;
            double y=1.0*(i+0.5)*dy;

            // if(pow(x-0.5,2) + pow(y-0.5,2) < pow(0.1,2)){
            //     rho[i*n1+j] = 1;    
            // }else{
            //     rho[i*n1+j] = 0;
            // }

            // rho[i*n1+j] = exp(-60*pow(x-0.5,2)-60*pow(y-0.5,2));

            rho[i*n1+j] = fmax(0.01-pow(x-0.5,2)-pow(y-0.5,2),0);
            sum += rho[i*n1+j];
        }
    }

    for(int i=0;i<n1*n2;++i){
        rho[i] *= (n1*n2)/sum * 0.1;
    }

    for(int n=1;n<nt;++n){
        for(int i=0;i<n1*n2;++i){
            rho[n*n1*n2+i]=rho[i];
        }
    }
}

void intialize_rho3(double* rho,double base,int n1,int n2,int nt,double dx,double dy,double dt){

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){

            double x=1.0*(j+0.5)*dx;
            double y=1.0*(i+0.5)*dy;

            rho[i*n1+j] = 0;
        }
    }

    for(int n=1;n<nt;++n){
        for(int i=0;i<n1*n2;++i){
            rho[n*n1*n2+i]=rho[i];
        }
    }
}

class Method{
public:
    int n1;
    int n2;
    int nt;
    double dx;
    double dy;
    double dt;

    double tau;
    double sigma;
    double h;

    int max_iteration;
    double tolerance;

    double* rhotmps[3];

    double* mx[3];
    double* my[3];
    double* mxtmps[3];
    double* mytmps[3];
    double* phi[3];

    double energy;
    double previous_energy;
    double previous_dual;

    double m;

    double M0;

    double c1;
    double c2;
    double c3;

    double* Clist;

    double* alphalist;

    // For SIR model

    double beta;
    double gamma;

    // ------------------------

    Method(){
    }

    Method(int n1, int n2, int nt, double dx, double dy, double dt, 
           double tau, double sigma, int max_iteration, double tolerance, 
           double c1, double c2, double c3, 
           double alpha1, double alpha2, double alpha3,
           double* Clist){
        this->n1=n1;
        this->n2=n2;
        this->nt=nt;
        this->dx=dx;
        this->dy=dy;
        this->dt=dt;
        this->tau=tau;
        this->sigma=sigma;
        this->max_iteration=max_iteration;
        this->tolerance=tolerance;

        this->c1=c1;
        this->c2=c2;
        this->c3=c3;

        this->Clist = Clist;

        alphalist = new double[3];
        alphalist[0] = alpha1;
        alphalist[1] = alpha2;
        alphalist[2] = alpha3;

        M0 = 0.5;

        for(int i=0;i<3;++i){
            mx[i]      = new double[n1*n2*nt];
            my[i]      = new double[n1*n2*nt];
            rhotmps[i] = new double[n1*n2*nt];
            mxtmps[i]  = new double[n1*n2*nt];
            mytmps[i]  = new double[n1*n2*nt];
            phi[i]     = new double[n1*n2*nt];
        }    
    }


    void destroy_all(){
        for(int i=0;i<3;++i){
            delete[] rhotmps[i];
            delete[] mx[i];
            delete[] my[i];
            delete[] mxtmps[i];
            delete[] mytmps[i];
            delete[] phi[i];
        }

        delete[] alphalist;
        delete[] Clist;
    }

    void update_m(double* mx, double* my, const double* rho, const double* phi, const int ind){

        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double nablaxphi = 0;
                    double nablayphi = 0;

                    if(j<n1-1){
                        nablaxphi = (phi[n*n1*n2+i*n1+j+1]-phi[n*n1*n2+i*n1+j])/dx; 
                    }else{
                        nablaxphi = 0;
                    }

                    if(i<n2-1){
                        nablayphi = (phi[n*n1*n2+(i+1)*n1+j]-phi[n*n1*n2+i*n1+j])/dy;
                    }else{
                        nablayphi = 0;
                    }

                    double rhovalx=0;
                    double rhovaly=0;

                    if(j==n1-1){
                        rhovalx = 0.5*(rho[n*n1*n2+i*n1+j]+rho[n*n1*n2+i*n1+j]);    
                    }else{
                        rhovalx = 0.5*(rho[n*n1*n2+i*n1+j+1]+rho[n*n1*n2+i*n1+j]);
                    }

                    if(i==n2-1){
                        rhovaly = 0.5*(rho[n*n1*n2+i*n1+j]+rho[n*n1*n2+i*n1+j]);    
                    }else{
                        rhovaly = 0.5*(rho[n*n1*n2+(i+1)*n1+j]+rho[n*n1*n2+i*n1+j]);
                    }
                    
                    mx[n*n1*n2+i*n1+j] = h * rhovalx/(tau * alphalist[ind] + h * rhovalx) * (mx[n*n1*n2+i*n1+j] - tau * nablaxphi);
                    my[n*n1*n2+i*n1+j] = h * rhovaly/(tau * alphalist[ind] + h * rhovaly) * (my[n*n1*n2+i*n1+j] - tau * nablayphi);

                    if(j==n1-1){
                        mx[n*n1*n2+i*n1+j]=0;
                    }
                    if(i==n2-1){
                        my[n*n1*n2+i*n1+j]=0;                       
                    }
                }
            }   
        }
        
    }

    void update_rho1(double* rho1,const double* rho2, const double* rho3,const double* mx,const double* my){

        for(int n=1;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int ind = n*n1*n2+i*n1+j;

                    double mxvalue, myvalue;
                    if(j==0){
                        mxvalue = 0.5*(mx[n*n1*n2+i*n1+j]);
                    }else{
                        mxvalue = 0.5*(mx[n*n1*n2+i*n1+j]+mx[n*n1*n2+i*n1+j-1]);
                    }
                    if(i==0){
                        myvalue = 0.5*(my[n*n1*n2+i*n1+j]);
                    }else{
                        myvalue = 0.5*(my[n*n1*n2+i*n1+j]+my[n*n1*n2+(i-1)*n1+j]);
                    }
                    double mvalue = sqrt(mxvalue*mxvalue + myvalue*myvalue);
                    double Dtphi;

                    if(n==nt-1){
                        Dtphi=1.0/dt*(phi[0][n*n1*n2+i*n1+j]-phi[0][(n-1)*n1*n2+i*n1+j]);
                    }else{
                        Dtphi=1.0/dt*(phi[0][(n+1)*n1*n2+i*n1+j]-phi[0][(n)*n1*n2+i*n1+j]);
                    }

                    double Deltaphi = -1.0/(dx*dx) * (-phi[0][n*n1*n2+i*n1+(int) fmax(0,j-1)]+2*phi[0][n*n1*n2+i*n1+j]-phi[0][n*n1*n2+i*n1+(int) fmin(n1-1,j+1)])
                                      -1.0/(dy*dy) * (-phi[0][n*n1*n2+(int) fmax(0,i-1)*n1+j]+2*phi[0][n*n1*n2+i*n1+j]-phi[0][n*n1*n2 +(int) fmin(n2-1,i+1)*n1+j]);


                    double newrhovalue=cubic_solve(tau*Dtphi + tau*Clist[0]*Deltaphi - rho1[ind] + tau*beta*rho2[ind]*(phi[1][ind]-phi[0][ind])  + tau*beta*M0*rho1[ind], 0, -0.5/h*tau*alphalist[0]*mvalue*mvalue);
                    rho1[n*n1*n2+i*n1+j]=fmax(0,newrhovalue);
                }
            }
        }
    }


    void update_rho2(const double* rho1,double* rho2, const double* rho3,const double* mx,const double* my){

        for(int n=1;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int ind = n*n1*n2+i*n1+j;

                    double mxvalue, myvalue;
                    if(j==0){
                        mxvalue = 0.5*(mx[n*n1*n2+i*n1+j]);
                    }else{
                        mxvalue = 0.5*(mx[n*n1*n2+i*n1+j]+mx[n*n1*n2+i*n1+j-1]);
                    }
                    if(i==0){
                        myvalue = 0.5*(my[n*n1*n2+i*n1+j]);
                    }else{
                        myvalue = 0.5*(my[n*n1*n2+i*n1+j]+my[n*n1*n2+(i-1)*n1+j]);
                    }
                    double mvalue = sqrt(mxvalue*mxvalue + myvalue*myvalue);
                    double Dtphi;

                    if(n==nt-1){
                        Dtphi=1.0/dt*(phi[1][n*n1*n2+i*n1+j]-phi[1][(n-1)*n1*n2+i*n1+j]);
                    }else{
                        Dtphi=1.0/dt*(phi[1][(n+1)*n1*n2+i*n1+j]-phi[1][(n)*n1*n2+i*n1+j]);
                    }


                    double Deltaphi = -1.0/(dx*dx) * (-phi[1][n*n1*n2+i*n1+(int) fmax(0,j-1)]+2*phi[1][n*n1*n2+i*n1+j]-phi[1][n*n1*n2+i*n1+(int) fmin(n1-1,j+1)])
                                      -1.0/(dy*dy) * (-phi[1][n*n1*n2+(int) fmax(0,i-1)*n1+j]+2*phi[1][n*n1*n2+i*n1+j]-phi[1][n*n1*n2 +(int) fmin(n2-1,i+1)*n1+j]);


                    double newrhovalue=cubic_solve(tau*Dtphi + tau*Clist[1]*Deltaphi - rho2[ind] + tau*beta*rho1[ind]*(phi[1][ind]-phi[0][ind]) + tau*gamma*(phi[2][ind] - phi[1][ind]) + tau*beta*M0*rho2[ind], 0, -0.5/h*tau*alphalist[1]*mvalue*mvalue);
                    rho2[ind]=fmax(0,newrhovalue);
                }
            }
        }

        // ----- rho(1,x) = dGstar(phi(1,x)) ------

        double sum1=0;
        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                int ind = (nt-1)*n1*n2+i*n1+j;
                // if(phi[1][ind] > 0){
                //     rho2[ind] = 1.0/c2 * phi[1][ind];
                // }else{
                //     rho2[ind] = 0;
                // }

                double eval = 1.0/c2 * phi[1][ind];
                double m = 2;

                if(eval > 0){
                    rho2[ind] = exp(1.0/(m-1) * log(eval));
                }else{
                    rho2[ind] = 0;
                }

            }
        }
    }

    void update_rho3(const double* rho1, const double* rho2, double* rho3,const double* mx,const double* my){

        for(int n=1;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int ind = n*n1*n2+i*n1+j;

                    double mxvalue, myvalue;
                    if(j==0){
                        mxvalue = 0.5*(mx[n*n1*n2+i*n1+j]);
                    }else{
                        mxvalue = 0.5*(mx[n*n1*n2+i*n1+j]+mx[n*n1*n2+i*n1+j-1]);
                    }
                    if(i==0){
                        myvalue = 0.5*(my[n*n1*n2+i*n1+j]);
                    }else{
                        myvalue = 0.5*(my[n*n1*n2+i*n1+j]+my[n*n1*n2+(i-1)*n1+j]);
                    }
                    double mvalue = sqrt(mxvalue*mxvalue + myvalue*myvalue);
                    double Dtphi;

                    if(n==nt-1){
                        Dtphi=1.0/dt*(phi[2][n*n1*n2+i*n1+j]-phi[2][(n-1)*n1*n2+i*n1+j]);
                    }else{
                        Dtphi=1.0/dt*(phi[2][(n+1)*n1*n2+i*n1+j]-phi[2][(n)*n1*n2+i*n1+j]);
                    }


                    double Deltaphi = -1.0/(dx*dx) * (-phi[2][n*n1*n2+i*n1+(int) fmax(0,j-1)]+2*phi[2][n*n1*n2+i*n1+j]-phi[2][n*n1*n2+i*n1+(int) fmin(n1-1,j+1)])
                                      -1.0/(dy*dy) * (-phi[2][n*n1*n2+(int) fmax(0,i-1)*n1+j]+2*phi[2][n*n1*n2+i*n1+j]-phi[2][n*n1*n2 +(int) fmin(n2-1,i+1)*n1+j]);


                    double newrhovalue=cubic_solve(tau*Dtphi + tau*Clist[2]*Deltaphi - rho3[ind], 0, -0.5/h*tau*alphalist[2]*mvalue*mvalue);
                    rho3[ind]=fmax(0,newrhovalue);
                }
            }
        }
    }

    double calculate_grad_mx(const double* mxTmp, int n, int i, int j){
        double mxval;

        if(j==0){
            mxval = (mxTmp[n*n1*n2+i*n1+j])/dx;
        }else{
            mxval = (mxTmp[n*n1*n2+i*n1+j]-mxTmp[n*n1*n2+i*n1+j-1])/dx;
        }

        return mxval;
    }

    double calculate_grad_my(const double* myTmp, int n, int i, int j){
        double myval;

        if(i==0){
            myval = (myTmp[n*n1*n2+i*n1+j])/dy;
        }else{
            myval = (myTmp[n*n1*n2+i*n1+j]-myTmp[n*n1*n2+(i-1)*n1+j])/dy;
        }

        return myval;
    }

    void update_phi1(poisson_solver& fftps, double* const rho[], const double* mx, const double* my){

        // Update 0 < t < 1

        for(int n=0;n<nt;++n){  
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double dtrho;
                    if(n==0){
                        dtrho=1.0/dt*(rho[0][n1*n2+i*n1+j]-rho[0][i*n1+j]); 
                    }else{
                        dtrho=1.0/dt*(rho[0][(n)*n1*n2+i*n1+j]-rho[0][(n-1)*n1*n2+i*n1+j]); 
                    }

                    double nablamx=calculate_grad_mx(mx,n,i,j);
                    double nablamy=calculate_grad_my(my,n,i,j);

                    int ind = n*n1*n2+i*n1+j;

                    double Deltarho = -1.0/(dx*dx) * (-rho[0][n*n1*n2+i*n1+(int) fmax(0,j-1)]+2*rho[0][n*n1*n2+i*n1+j]-rho[0][n*n1*n2+i*n1+(int) fmin(n1-1,j+1)])
                                      -1.0/(dy*dy) * (-rho[0][n*n1*n2+(int) fmax(0,i-1)*n1+j]+2*rho[0][n*n1*n2+i*n1+j]-rho[0][n*n1*n2 +(int) fmin(n2-1,i+1)*n1+j]);

                    fftps.workspace[n*n1*n2+i*n1+j]=-(dtrho+nablamx+nablamy + beta*rho[0][ind]*rho[1][ind] - Clist[0]*Deltarho); 
                }
            }
        }

        fftps.perform_inverse_laplacian();
        for(int i=0;i<n1*n2*nt;++i){
            phi[0][i] += sigma*fftps.workspace[i];
        }
    }

    void update_phi2(poisson_solver& fftps, double* const rho[], const double* mx, const double* my){

        // Update 0 < t < 1

        for(int n=0;n<nt;++n){  
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double dtrho;
                    if(n==0){
                        dtrho=1.0/dt*(rho[1][n1*n2+i*n1+j]-rho[1][i*n1+j]); 
                    }else{
                        dtrho=1.0/dt*(rho[1][(n)*n1*n2+i*n1+j]-rho[1][(n-1)*n1*n2+i*n1+j]); 
                    }

                    double nablamx=calculate_grad_mx(mx,n,i,j);
                    double nablamy=calculate_grad_my(my,n,i,j);

                    int ind = n*n1*n2+i*n1+j;

                    double Deltarho = -1.0/(dx*dx) * (-rho[1][n*n1*n2+i*n1+(int) fmax(0,j-1)]+2*rho[1][n*n1*n2+i*n1+j]-rho[1][n*n1*n2+i*n1+(int) fmin(n1-1,j+1)])
                                      -1.0/(dy*dy) * (-rho[1][n*n1*n2+(int) fmax(0,i-1)*n1+j]+2*rho[1][n*n1*n2+i*n1+j]-rho[1][n*n1*n2 +(int) fmin(n2-1,i+1)*n1+j]);

                    fftps.workspace[n*n1*n2+i*n1+j]=-(dtrho+nablamx+nablamy - beta*rho[0][ind]*rho[1][ind] + gamma*rho[1][ind] - Clist[1]*Deltarho); 
                }
            }
        }

        fftps.perform_inverse_laplacian();
        for(int i=0;i<n1*n2*nt;++i){
            phi[1][i] += sigma*fftps.workspace[i];
        }
    }

    void update_phi3(poisson_solver& fftps, double* const rho[], const double* mx, const double* my){

        // Update 0 < t < 1

        for(int n=0;n<nt;++n){  
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double dtrho;
                    if(n==0){
                        dtrho=1.0/dt*(rho[2][n1*n2+i*n1+j]-rho[2][i*n1+j]); 
                    }else{
                        dtrho=1.0/dt*(rho[2][(n)*n1*n2+i*n1+j]-rho[2][(n-1)*n1*n2+i*n1+j]); 
                    }

                    double nablamx=calculate_grad_mx(mx,n,i,j);
                    double nablamy=calculate_grad_my(my,n,i,j);

                    int ind = n*n1*n2+i*n1+j;


                    double Deltarho = -1.0/(dx*dx) * (-rho[2][n*n1*n2+i*n1+(int) fmax(0,j-1)]+2*rho[2][n*n1*n2+i*n1+j]-rho[2][n*n1*n2+i*n1+(int) fmin(n1-1,j+1)])
                                      -1.0/(dy*dy) * (-rho[2][n*n1*n2+(int) fmax(0,i-1)*n1+j]+2*rho[2][n*n1*n2+i*n1+j]-rho[2][n*n1*n2 +(int) fmin(n2-1,i+1)*n1+j]);


                    fftps.workspace[n*n1*n2+i*n1+j]=-(dtrho+nablamx+nablamy - gamma*rho[1][ind] - Clist[2]*Deltarho); 
                }
            }
        }

        fftps.perform_inverse_laplacian();
        for(int i=0;i<n1*n2*nt;++i){
            phi[2][i] += sigma*fftps.workspace[i];
        }
    }

    double calculate_energy1(const double* rho,const double* mx,const double* my) const{
        double sum1=0;
        double sum2=0;

        for(int n=0;n<nt-1;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int ind = n*n1*n2+i*n1+j;
                    double mxval=mx[ind];
                    double myval=my[ind];

                    double mval=sqrt(mxval*mxval+myval*myval);

                    double rhoval=rho[ind];

                    if(rhoval>0){
                        sum1 += alphalist[0]*mval*mval/(2.0*rhoval);
                    }
                }
            }
        }

        return sum1*dx*dy*dt + sum2*dx*dy;
    }


    double calculate_energy2(const double* rho,const double* mx,const double* my) const{
        double sum1=0;
        double sum2=0;

        for(int n=0;n<nt-1;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int ind = n*n1*n2+i*n1+j;
                    double mxval=mx[ind];
                    double myval=my[ind];

                    double mval=sqrt(mxval*mxval+myval*myval);

                    double rhoval=rho[ind];

                    if(rhoval>0){
                        sum1 += alphalist[1]*mval*mval/(2.0*rhoval);
                    }
                }
            }
        }

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                double rhoval=rho[(nt-1)*n1*n2+i*n1+j];
                sum2+=c2/m*exp(m*log(rhoval));
            }
        }

        return sum1*dx*dy*dt + sum2*dx*dy;
    }

    double calculate_energy3(double* const rho,const double* mx,const double* my) const{
        double sum1=0;
        double sum2=0;

        for(int n=0;n<nt-1;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int ind = n*n1*n2+i*n1+j;
                    double mxval=mx[ind];
                    double myval=my[ind];

                    double mval=sqrt(mxval*mxval+myval*myval);

                    double rhoval=rho[ind];

                    if(rhoval>0){
                        sum1 += alphalist[2]*mval*mval/(2.0*rhoval);
                    }
                }
            }
        }

        return sum1*dx*dy*dt + sum2*dx*dy;
    }


    double calculate_energy(double* const rho[]) const{
        double e1=calculate_energy1(rho[0],mx[0],my[0]);
        double e2=calculate_energy2(rho[1],mx[1],my[1]);
        double e3=calculate_energy3(rho[2],mx[2],my[2]);

        return e1+e2+e3;
    }

    double calculate_dual(double* const rho[]) const{   

        double mprime = m/(m-1);

        double term0=0;
        double term1=0;
        double term2=0;

        double c2_value = pow(c2,1.0/(m-1));
        for(int i=0;i<n1*n2;++i){
            term0 += phi[0][i]*rho[0][i] - phi[0][(nt-1)*n1*n2+i]*rho[0][(nt-1)*n1*n2+i];
            term1 += phi[1][i]*rho[1][i] - 1.0/(mprime*c2_value) * pow(fmax(0, phi[1][(nt-1)*n1*n2+i]), mprime);
            term2 += phi[2][i]*rho[2][i] - phi[2][(nt-1)*n1*n2+i]*rho[2][(nt-1)*n1*n2+i];
        }

        double term3=0;
        for(int n=0;n<nt;++n){
            for(int i=0;i<n1*n2;++i){
                term3 += beta*rho[0][n*n1*n2+i]*rho[1][n*n1*n2+i] * (phi[1][n*n1*n2+i]-phi[0][n*n1*n2+i]);
            }
        }

        return (term0+term1+term2)/(n1*n2) + term3/(n1*n2*nt);
    }

    void display_log(const int iterPDHG, const double tau, const double sigma, const double energy, const double dual, const double error, const double dual_gap) const{
        cout<<"iter : "<< setw(10) << iterPDHG+1<< " tau : "<<tau <<" sigma : "<<sigma <<" energy : "<<energy <<" dual : " << dual << " error : " << error <<" dual_gap : " << dual_gap << endl;  
    }

    void run(poisson_solver& fftps, double* rho[], int skip=1){

        previous_energy=1;
        previous_dual=1;
        double error=1, dual_gap=1, energy=1, dual=0;
        int iterPDHG;

        double beta_1 = 1.2;
        double beta_2 = 0.9;

        for(iterPDHG=0; iterPDHG<max_iteration; ++iterPDHG){

            // get the data before updates
            for(int i=0;i<3;++i){
                memcpy(rhotmps[i],rho[i],n1*n2*nt*sizeof(double));
                memcpy(mxtmps[i],mx[i],n1*n2*nt*sizeof(double));
                memcpy(mytmps[i],my[i],n1*n2*nt*sizeof(double));    
            }

            update_rho1(rho[0],rho[1],rho[2],mx[0],my[0]);
            update_rho2(rho[0],rho[1],rho[2],mx[1],my[1]);
            update_rho3(rho[0],rho[1],rho[2],mx[2],my[2]);            

            // normalize rho
            
            for(int n=1;n<nt;++n){
                double sumtmp = 0;
                double sum = 0;
                for(int k=0;k<3;++k){
                    for(int i=0;i<n1*n2;++i){
                        sumtmp += rhotmps[k][n*n1*n2+i];
                        sum    += rho[k][n*n1*n2+i];
                    }
                }
                for(int k=0;k<3;++k){
                    for(int i=0;i<n1*n2;++i){
                        if(sum > 0){
                            rho[k][n*n1*n2+i] *= sumtmp/sum;    
                        }
                    }
                }
            }


            update_m(mx[0],my[0],rho[0],phi[0],0);
            update_m(mx[1],my[1],rho[1],phi[1],1);
            update_m(mx[2],my[2],rho[2],phi[2],2);

            for(int k=0;k<3;++k){
                for(int i=0;i<n1*n2*nt;++i){
                    rhotmps[k][i] = 2*rho[k][i] - rhotmps[k][i];
                    mxtmps[k][i]  = 2*mx[k][i]  - mxtmps[k][i];
                    mytmps[k][i]  = 2*my[k][i]  - mytmps[k][i];
                }
            }

            update_phi1(fftps,rhotmps,mxtmps[0],mytmps[0]); 
            update_phi2(fftps,rhotmps,mxtmps[1],mytmps[1]); 
            update_phi3(fftps,rhotmps,mxtmps[2],mytmps[2]); 

            energy=calculate_energy(rho);
            error=fabs((energy-previous_energy)/previous_energy);
            previous_energy=energy;

            dual  =calculate_dual(rho);
            double relative_dual_error=abs((dual-previous_dual)/previous_dual);
            previous_dual=dual;


            if(error > relative_dual_error*beta_1){
                tau   *= beta_2;
                sigma /= beta_2;
            }

            if(error < relative_dual_error/beta_1){
                tau   /= beta_2;
                sigma *= beta_2;
            }

            tau = fmax(1e-5,fmin(0.01, tau));
            sigma = fmax(1e-3,fmin(1, sigma));

            double dual_gap = energy-dual;

                // adaptive step size

            if((iterPDHG+1)%skip==0 || error<tolerance){

                display_log(iterPDHG,  tau,  sigma,  energy,  dual,  error,  dual_gap);

                create_csv_file(rho[0],"./data/rho1.csv",n1,n2,nt);
                create_csv_file(rho[1],"./data/rho2.csv",n1,n2,nt);
                create_csv_file(rho[2],"./data/rho3.csv",n1,n2,nt);

                if(error<tolerance) break;
            }
        }

        cout<<"The method is done!!"<<endl;
        cout<<"iter : "<< setw(10) << iterPDHG+1<< " energy : "<<energy<< " error : " << error <<" dual : " << dual <<" dual_gap : " << dual_gap << endl;  

        cout<<"Iterations     : "<<iterPDHG<<endl;
        cout<<"Energy         : "<<energy<<endl;
        cout<<"Relative Error : "<<error<<endl;
    }

}; // Method class

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
    
    create_csv_file_for_parameters(n1,n2,nt);

    double base=0;

    double dx=1.0/n1;
    double dy=1.0/n2;
    double dt=1.0/nt;

    double c1=0.1;
    double c2=0.5;
    double c3=0.01;

    double alpha1 = 2.0;
    double alpha2 = 1.0;
    double alpha3 = 10.0;

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

    double* rho[3];

    for(int i=0;i<3;++i){
        rho[i] = new double[n1*n2*nt];
    }

    clock_t t;
    t = clock();

    double Clist[] = {eta, eta, eta};
    poisson_solver fftps(n1,n2,nt,dx,dy,dt,eta);

    t = clock() - t;
    printf ("\nCPU time for setting up FFT: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  

    intialize_rho1(rho[0],base,n1,n2,nt,dx,dy,dt);
    intialize_rho2(rho[1],base,n1,n2,nt,dx,dy,dt);
    intialize_rho3(rho[2],base,n1,n2,nt,dx,dy,dt);

    Method method(n1, n2, nt, dx, dy, dt, tau, sigma, max_iteration, tolerance, c1, c2, c3, alpha1, alpha2, alpha3, Clist);

    cout<<"\nXXX Starting Iterations XXX"<<endl;

    t = clock();

    method.h = 1;
    method.beta  = 0.5;
    method.gamma = 0.2;

    method.m = 2;

    string filename;

    filename="./data/rho1-"+to_string(0)+".csv";
    create_bin_file(&rho[0][(nt-1)*n1*n2],n1*n2,filename);
    filename="./data/rho2-"+to_string(0)+".csv";
    create_bin_file(&rho[1][(nt-1)*n1*n2],n1*n2,filename);
    filename="./data/rho3-"+to_string(0)+".csv";
    create_bin_file(&rho[2][(nt-1)*n1*n2],n1*n2,filename);

    for(int iter=0; iter<1; ++iter)
    {
        method.run(fftps,rho,skip);   

        filename="./data/rho1-"+to_string(iter+1)+".csv";
        create_bin_file(&rho[0][(nt-1)*n1*n2],n1*n2,filename);
        filename="./data/rho2-"+to_string(iter+1)+".csv";
        create_bin_file(&rho[1][(nt-1)*n1*n2],n1*n2,filename);
        filename="./data/rho3-"+to_string(iter+1)+".csv";
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

    // create_csv_file(rho,"./data/rho.csv",n1,n2,nt);

    create_csv_file_for_parameters(n1,n2,nt);

    for(int k=0;k<3;++k){
        delete[] rho[k];
    }

    fftps.destroy_all_fftps();
    method.destroy_all();
}