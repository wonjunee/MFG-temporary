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

            rho[i*n1+j] = exp(-20*pow(x-0.6,2)-20*pow(y-0.6,2));
            sum += rho[i*n1+j];
        }
    }

    // for(int i=0;i<n1*n2;++i){
    //     rho[i] /= sum/(n1*n2);
    // }

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

void initialize_f_g(double* f,double* g,int n1,int n2,double dx,double dy){

    double sum1=0;
    double sum2=0;

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){

            double x=1.0*(j+0.5)*dx;
            double y=1.0*(i+0.5)*dy;

            // f[i*n1+j] = gaussian(x,y,0.3,0.5,0.05,0.05) + gaussian(x,y,0.7,0.5,0.05,0.05);
            f[i*n1+j] = gaussian(x,y,0.5,0.7,0.04,0.04) + gaussian(x,y,0.2,0.5,0.04,0.04) + gaussian(x,y,0.8,0.5,0.04,0.04) + gaussian(x,y,0.5,0.3,0.04,0.04);
            
            // g[i*n1+j] = 10.0 * (pow(x-0.9,2)/pow(0.05,2) + pow(y-0.9,2)/pow(0.05,2));
            // g[i*n1+j] = 5.0* (pow(x-0.9,2)/pow(0.05,2) + pow(y-0.1,2)/pow(0.05,2)) + 5.0* (pow(x-0.1,2)/pow(0.05,2) + pow(y-0.1,2)/pow(0.05,2));

            // f[i*n1+j] = gaussian(x,y,0.5,0.5,0.05,0.05);
            g[i*n1+j] = 0.5*gaussian(x,y,0.9,0.9,0.05,0.05);

            sum1+=f[i*n1+j];
            sum2+=g[i*n1+j];
        }
    }

    for(int i=0;i<n1*n2;++i){
        f[i]/=sum1*dx*dy*5;
        // g[i]/=sum2*dx*dy/1.0;
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

    double* mx1;
    double* my1;

    double* mx2;
    double* my2;

    double* mx3;
    double* my3;

    double* phi1;
    double* phi2;
    double* phi3;

    double* phi1tmp;
    double* phi2tmp;
    double* phi3tmp;

    double* f;
    double* g;

    double energy;
    double previous_energy;
    double previous_dual;

    double c1;
    double c2;
    double c3;

    double* alphalist;

    // For SIR model

    double beta;
    double gamma;

    double* px;
    double* py;

    Method(){
    }

    Method(int n1, int n2, int nt, double dx, double dy, double dt, double tau, double sigma, int max_iteration, double tolerance, double c1, double c2, double c3, double alpha1, double alpha2, double alpha3){
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

        alphalist = new double[3];
        alphalist[0] = alpha1;
        alphalist[1] = alpha2;
        alphalist[2] = alpha3;

        mx1=new double[n1*n2*nt];
        my1=new double[n1*n2*nt];

        mx2=new double[n1*n2*nt];
        my2=new double[n1*n2*nt];

        mx3=new double[n1*n2*nt];
        my3=new double[n1*n2*nt];

        
        phi1=new double[n1*n2*nt];
        phi2=new double[n1*n2*nt];
        phi3=new double[n1*n2*nt];

        phi1tmp=new double[n1*n2*nt];
        phi2tmp=new double[n1*n2*nt];
        phi3tmp=new double[n1*n2*nt];

        px=new double[101];
        py=new double[101];
    }

    void setup_f_g(const double* f,const double* g){
        memcpy(this->f,f,n1*n2*sizeof(double));
        memcpy(this->g,g,n1*n2*sizeof(double));
    }

    void destroy_all(){
        delete[] mx1;
        delete[] my1;

        delete[] mx2;
        delete[] my2;

        delete[] mx3;
        delete[] my3;

        delete[] phi1;
        delete[] phi2;
        delete[] phi3;

        delete[] phi1tmp;
        delete[] phi2tmp;
        delete[] phi3tmp;
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

    void update_rho1(double* rho1,const double* rho2, const double* rho3,const double* mx,const double* my,const double* f,const double* g){

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
                        Dtphi=1.0/dt*(phi1[n*n1*n2+i*n1+j]-phi1[(n-1)*n1*n2+i*n1+j]);
                    }else{
                        Dtphi=1.0/dt*(phi1[(n+1)*n1*n2+i*n1+j]-phi1[(n)*n1*n2+i*n1+j]);
                    }

                    double newrhovalue=cubic_solve(tau*Dtphi - rho1[ind] + tau*beta*rho2[ind]*(phi2[ind]-phi1[ind]), 0, -0.5/h*tau*alphalist[0]*mvalue*mvalue);
                    rho1[n*n1*n2+i*n1+j]=fmax(0,newrhovalue);
                }
            }
        }
    }


    void update_rho2(const double* rho1,double* rho2, const double* rho3,const double* mx,const double* my,const double* f,const double* g){

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
                        Dtphi=1.0/dt*(phi2[n*n1*n2+i*n1+j]-phi2[(n-1)*n1*n2+i*n1+j]);
                    }else{
                        Dtphi=1.0/dt*(phi2[(n+1)*n1*n2+i*n1+j]-phi2[(n)*n1*n2+i*n1+j]);
                    }

                    double newrhovalue=cubic_solve(tau*Dtphi - rho2[ind] + tau*beta*rho1[ind]*(phi2[ind]-phi1[ind]) + tau*gamma*(phi3[ind] - phi2[ind]), 0, -0.5/h*tau*alphalist[1]*mvalue*mvalue);
                    rho2[ind]=fmax(0,newrhovalue);
                }
            }
        }

        // ----- rho(1,x) = dGstar(phi(1,x)) ------

        // double sum1=0;
        // for(int i=0;i<n2;++i){
        //     for(int j=0;j<n1;++j){
        //         int ind = (nt-1)*n1*n2+i*n1+j;
        //         // if(phi2[ind] > 0){
        //         //     rho2[ind] = 1.0/c2 * phi2[ind];
        //         // }else{
        //         //     rho2[ind] = 0;
        //         // }

        //         rho2[ind] = exp(phi2[ind]/c2 - 1);
        //     }
        // }
    }

    void update_rho3(const double* rho1, const double* rho2, double* rho3,const double* mx,const double* my,const double* f,const double* g){

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
                        Dtphi=1.0/dt*(phi3[n*n1*n2+i*n1+j]-phi3[(n-1)*n1*n2+i*n1+j]);
                    }else{
                        Dtphi=1.0/dt*(phi3[(n+1)*n1*n2+i*n1+j]-phi3[(n)*n1*n2+i*n1+j]);
                    }

                    double newrhovalue=cubic_solve(tau*Dtphi - rho3[ind], 0, -0.5/h*tau*alphalist[2]*mvalue*mvalue);
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

    void update_phi1(poisson_solver& fftps, const double* rho1, const double* rho2, const double* rho3, const double* mx, const double* my, const double* g){

        // Update 0 < t < 1

        for(int n=0;n<nt;++n){  
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double dtrho;
                    if(n==0){
                        dtrho=1.0/dt*(rho1[n1*n2+i*n1+j]-rho1[i*n1+j]); 
                    }else{
                        dtrho=1.0/dt*(rho1[(n)*n1*n2+i*n1+j]-rho1[(n-1)*n1*n2+i*n1+j]); 
                    }

                    double nablamx=calculate_grad_mx(mx,n,i,j);
                    double nablamy=calculate_grad_my(my,n,i,j);

                    int ind = n*n1*n2+i*n1+j;
                    fftps.workspace[n*n1*n2+i*n1+j]=-(dtrho+nablamx+nablamy + beta*rho1[ind]*rho2[ind]); 
                }
            }
        }

        fftps.perform_inverse_laplacian();
        for(int i=0;i<n1*n2*nt;++i){
            phi1[i] += sigma*fftps.workspace[i];
        }
    }

    void update_phi2(poisson_solver& fftps, const double* rho1, const double* rho2, const double* rho3, const double* mx, const double* my, const double* g){

        // Update 0 < t < 1

        for(int n=0;n<nt;++n){  
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double dtrho;
                    if(n==0){
                        dtrho=1.0/dt*(rho2[n1*n2+i*n1+j]-rho2[i*n1+j]); 
                    }else{
                        dtrho=1.0/dt*(rho2[(n)*n1*n2+i*n1+j]-rho2[(n-1)*n1*n2+i*n1+j]); 
                    }

                    double nablamx=calculate_grad_mx(mx,n,i,j);
                    double nablamy=calculate_grad_my(my,n,i,j);

                    int ind = n*n1*n2+i*n1+j;
                    fftps.workspace[n*n1*n2+i*n1+j]=-(dtrho+nablamx+nablamy - beta*rho1[ind]*rho2[ind] + gamma*rho2[ind]); 
                }
            }
        }

        fftps.perform_inverse_laplacian();
        for(int i=0;i<n1*n2*nt;++i){
            phi2[i] += sigma*fftps.workspace[i];
        }
    }

    void update_phi3(poisson_solver& fftps, const double* rho1, const double* rho2, const double* rho3, const double* mx, const double* my, const double* g){

        // Update 0 < t < 1

        for(int n=0;n<nt;++n){  
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double dtrho;
                    if(n==0){
                        dtrho=1.0/dt*(rho3[n1*n2+i*n1+j]-rho3[i*n1+j]); 
                    }else{
                        dtrho=1.0/dt*(rho3[(n)*n1*n2+i*n1+j]-rho3[(n-1)*n1*n2+i*n1+j]); 
                    }

                    double nablamx=calculate_grad_mx(mx,n,i,j);
                    double nablamy=calculate_grad_my(my,n,i,j);

                    int ind = n*n1*n2+i*n1+j;
                    fftps.workspace[n*n1*n2+i*n1+j]=-(dtrho+nablamx+nablamy - gamma*rho2[ind]); 
                }
            }
        }

        fftps.perform_inverse_laplacian();
        for(int i=0;i<n1*n2*nt;++i){
            phi3[i] += sigma*fftps.workspace[i];
        }
    }

    double calculate_energy1(const double* rho,const double* mx,const double* my){
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


    double calculate_energy2(const double* rho,const double* mx,const double* my){
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

        // for(int i=0;i<n2;++i){
        //     for(int j=0;j<n1;++j){
        //         double rhoval=rho[(nt-1)*n1*n2+i*n1+j];
        //         // sum2+=0.5*c2*rhoval*rhoval;
        //         sum2+=c2*rhoval*log(rhoval);
        //     }
        // }

        return sum1*dx*dy*dt + sum2*dx*dy;
    }

    double calculate_energy3(const double* rho,const double* mx,const double* my){
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


    double calculate_energy(const double* rho1,const double* rho2,const double* rho3){
        double e1=calculate_energy1(rho1,mx1,my1);
        double e2=calculate_energy2(rho2,mx2,my2);
        double e3=calculate_energy3(rho3,mx3,my3);

        return e1+e2+e3;
    }

    double calculate_dual(const double* rho,const double* f,const double* g){   

        // TODO

        return 0;
    }

    void run(poisson_solver& fftps, double* rho1, double* rho2, double* rho3, const double* f,const double* g, int skip=1){

        previous_energy=1;
        previous_dual=1;
        double error, dual_gap, energy, dual;
        int iterPDHG;

        for(iterPDHG=0; iterPDHG<max_iteration; ++iterPDHG){

            double sum_rho     = 0;
            double sum_rho_new = 0;

            // for(int i=0;i<nt*n1*n2;++i){
            //     sum_rho += rho1[i];
            //     sum_rho += rho2[i];
            //     sum_rho += rho3[i];
            // }

            update_rho1(rho1,rho2,rho3,mx1,my1,f,g);
            update_rho2(rho1,rho2,rho3,mx2,my2,f,g);
            update_rho3(rho1,rho2,rho3,mx3,my3,f,g);

            // for(int i=0;i<nt*n1*n2;++i){
            //     sum_rho_new += rho1[i];
            //     sum_rho_new += rho2[i];
            //     sum_rho_new += rho3[i];
            // }

            // for(int i=0;i<nt*n1*n2;++i){
            //     rho1[i] *= sum_rho/sum_rho_new;
            //     rho2[i] *= sum_rho/sum_rho_new;
            //     rho3[i] *= sum_rho/sum_rho_new;
            // }



            update_m(mx1,my1,rho1,phi1,0);
            update_m(mx2,my2,rho2,phi2,1);
            update_m(mx3,my3,rho3,phi3,2);

            memcpy(phi1tmp, phi1, n1*n2*nt*sizeof(double));
            memcpy(phi2tmp, phi2, n1*n2*nt*sizeof(double));
            memcpy(phi3tmp, phi3, n1*n2*nt*sizeof(double));

            update_phi1(fftps,rho1,rho2,rho3,mx1,my1,g); 
            update_phi2(fftps,rho1,rho2,rho3,mx2,my2,g); 
            update_phi3(fftps,rho1,rho2,rho3,mx3,my3,g); 

            for(int n=0;n<nt;++n){
                for(int i=0;i<n1*n2;++i){
                    phi1[n*n1*n2+i]=2*phi1[n*n1*n2+i]-phi1tmp[n*n1*n2+i];   
                    phi2[n*n1*n2+i]=2*phi2[n*n1*n2+i]-phi2tmp[n*n1*n2+i];   
                    phi3[n*n1*n2+i]=2*phi3[n*n1*n2+i]-phi3tmp[n*n1*n2+i];   
                }
            }

            energy=calculate_energy(rho1,rho2,rho3);

            dual=calculate_dual(rho1,f,g);

            error=fabs((energy-previous_energy)/previous_energy);

            previous_energy=energy;

            dual_gap=fabs(energy-dual);

            if((iterPDHG+1)%skip==0){

                cout<<"iter : "<< setw(10) << iterPDHG+1<< " energy : "<<energy<< " error : " << error <<" dual : " << dual <<" dual_gap : " << dual_gap << endl;  

                create_csv_file(rho1,"./data/rho1.csv",n1,n2,nt);
                create_csv_file(rho2,"./data/rho2.csv",n1,n2,nt);
                create_csv_file(rho3,"./data/rho3.csv",n1,n2,nt);
            }
            
            if(error<tolerance){
                break;
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
        cout << "./main.exe [n1] [n2] [nt] [tau] [sigma] [tolerance] [max iteration] [skip]" << endl;
        exit(1);
    }

    // Parameters for Grids
    
    int n1 = stoi(argv[1]);  // x panels
    int n2 = stoi(argv[2]);  // y panels
    double dt = stod(argv[3]);  // dt
    double T = stod(argv[4]);  // Total Time
    double tau=stod(argv[5]);
    double sigma=stod(argv[6]);
    double tolerance = stod(argv[7]);
    int max_iteration=stoi(argv[8]);
    int skip=stoi(argv[9]);

    int nt=T/dt;
    
    create_csv_file_for_parameters(n1,n2,nt);

    double base=0;

    double dx=1.0/n1;
    double dy=1.0/n2;

    double c1=0.1;
    double c2=0.5;
    double c3=0.01;

    double alpha1 = 0.1;
    double alpha2 = 1.0;
    double alpha3 = 1.0;

    cout<<"XXX G-Prox PDHG XXX"<<endl;
    cout<<endl;
    cout<<"n1 : "<<n1<<" n2 : "<<n2<<" nt : "<<nt<<" base : "<<base<<endl;
    cout<<"dx : "<<scientific<<dx<<endl;
    cout<<"dy : "<<scientific<<dy<<endl;
    cout<<"dt : "<<scientific<<dt<<endl;
    cout<<"Total Time : "<<scientific<<T<<endl;
    cout<<fixed;
    cout<<"tau   : "<<tau<<endl;
    cout<<"sigma : "<<sigma<<endl;
    cout<<"max_iteration : "<<max_iteration<<endl;
    cout<<"tolerance     : "<<scientific<<tolerance<<endl;

    double* rho1=new double[n1*n2*nt];
    double* rho2=new double[n1*n2*nt];
    double* rho3=new double[n1*n2*nt];
    double* f=new double[n1*n2];
    double* g=new double[n1*n2];

    clock_t t;
    t = clock();
    poisson_solver fftps(n1,n2,nt,dx,dy,dt);
    t = clock() - t;
    printf ("\nCPU time for setting up FFT: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  

    intialize_rho1(rho1,base,n1,n2,nt,dx,dy,dt);
    intialize_rho2(rho2,base,n1,n2,nt,dx,dy,dt);
    intialize_rho3(rho3,base,n1,n2,nt,dx,dy,dt);
    initialize_f_g(f,g,n1,n2,dx,dy);


    Method method(n1, n2, nt, dx, dy, dt, tau, sigma, max_iteration, tolerance, c1, c2, c3, alpha1, alpha2, alpha3);

    cout<<"\nXXX Starting Iterations XXX"<<endl;

    t = clock();

    method.h = 1;
    method.beta  = 0.5;
    method.gamma = 0.;

    string filename;

    filename="./data/rho1-"+to_string(0)+".csv";
    create_bin_file(&rho1[(nt-1)*n1*n2],n1*n2,filename);
    filename="./data/rho2-"+to_string(0)+".csv";
    create_bin_file(&rho2[(nt-1)*n1*n2],n1*n2,filename);
    filename="./data/rho3-"+to_string(0)+".csv";
    create_bin_file(&rho3[(nt-1)*n1*n2],n1*n2,filename);

    for(int iter=0; iter<1; ++iter)
    {
        method.run(fftps,rho1,rho2,rho3,f,g,skip);   

        filename="./data/rho1-"+to_string(iter+1)+".csv";
        create_bin_file(&rho1[(nt-1)*n1*n2],n1*n2,filename);
        filename="./data/rho2-"+to_string(iter+1)+".csv";
        create_bin_file(&rho2[(nt-1)*n1*n2],n1*n2,filename);
        filename="./data/rho3-"+to_string(iter+1)+".csv";
        create_bin_file(&rho3[(nt-1)*n1*n2],n1*n2,filename);

        for(int n=0;n<nt-1;++n){
            for(int i=0;i<n1*n2;++i){
                rho1[n*n1*n2+i] = rho1[(nt-1)*n1*n2+i];
                rho2[n*n1*n2+i] = rho2[(nt-1)*n1*n2+i];
                rho3[n*n1*n2+i] = rho3[(nt-1)*n1*n2+i];
            }
        }
    }
    

    t = clock() - t;

    printf ("CPU time for Iterations: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);

    // create_csv_file(rho,"./data/rho.csv",n1,n2,nt);

    create_csv_file_for_parameters(n1,n2,nt);

    delete[] rho1;
    delete[] rho2;
    delete[] rho3;
    delete[] f;
    delete[] g;

    fftps.destroy_all_fftps();
    method.destroy_all();
}