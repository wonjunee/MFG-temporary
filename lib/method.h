#ifndef METHOD_H
#define METHOD_H
#include <iostream>
#include <fftw3.h>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstring>
// #include "poisson_solver.h"
#include "poisson_solver_3d.h"
#include "helper.h"
#include "initializer.h"

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

    double c0;
    double c1;
    double c2;

    double* Clist;

    double* alphalist;

    poisson_solver* fftps;

    int convN;
    double conv_sum;
    double conv_r;

    // For SIR model

    double beta;
    double gamma;

    // function pointers for terminal functions
    // double (*E0)(double);
    // double (*E1)(double);
    // double (*E2)(double);

    // double (*deltaE0)(const double* rho);
    // double (*deltaE1)(const double* rho);
    // double (*deltaE2)(const double* rho);

    // ------------------------

    Method(){
        for(int i=0;i<3;++i){
            rhotmps[i]=NULL;
            mx[i]=NULL;my[i]=NULL;
            mxtmps[i]=NULL;mytmps[i]=NULL;
            phi[i]=NULL;    
        }
        
        Clist=NULL;
        alphalist=NULL;
        fftps=NULL;
    }

    Method(int n1, int n2, int nt, double dx, double dy, double dt, 
           double tau, double sigma, int max_iteration, double tolerance, 
           double c0, double c1, double c2, 
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

        this->c0=c0;
        this->c1=c1;
        this->c2=c2;

        this->Clist = Clist; // this contains the eta values for viscosity terms

        alphalist = new double[3];
        alphalist[0] = alpha1;
        alphalist[1] = alpha2;
        alphalist[2] = alpha3;

        M0 = 0.5;
        // convN = n1/4;
        // convN = convN + (convN % 2) -1; // to make it odd
        conv_r = 0.1;
        convN  = 2*conv_r*n1+1;
        convN  = convN + (convN % 2) -1; // to make it odd
        cout << "convN : " <<convN << endl;

        conv_sum = 0;
        int i = n2/2;
        int j = n1/2;
        double xx = (j+0.5)/n1;
        double yy = (i+0.5)/n2;
        for(int i1=i-convN/2;i1<i+convN/2+1;++i1){
            for(int j1=j-convN/2;j1<j+convN/2+1;++j1){
                
                double x=(j1+0.5)/n1;
                double y=(i1+0.5)/n2;        

                if((x-xx)*(x-xx) + (y-yy)*(y-yy) <= conv_r*conv_r){
	                double eval = calculate_K_xy(x,xx,y,yy);
	                conv_sum += eval;
	            }
            }
        }

        for(int i=0;i<3;++i){
            mx[i]      = new double[n1*n2*nt];
            my[i]      = new double[n1*n2*nt];
            rhotmps[i] = new double[n1*n2*nt];
            mxtmps[i]  = new double[n1*n2*nt];
            mytmps[i]  = new double[n1*n2*nt];
            phi[i]     = new double[n1*n2*nt];
        }    

        clock_t t;
        t = clock();
            
        double eta = Clist[0];
        fftps = new poisson_solver(n1,n2,nt,dx,dy,dt,eta);

        t = clock() - t;
        printf ("\nCPU time for setting up FFT: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
    }


    ~Method(){
        for(int i=0;i<3;++i){
            delete[] rhotmps[i];
            delete[] mx[i];
            delete[] my[i];
            delete[] mxtmps[i];
            delete[] mytmps[i];
            delete[] phi[i];
        }
        delete[] alphalist;
        delete fftps;
    }

    /* 
        Set function pointers
    */

    // void set_E0(double (*f)(double)){
    //     E0 = f;
    // }
    // void set_E1(double (*f)(double)){
    //     E1 = f;
    // }
    // void set_E2(double (*f)(double)){
    //     E2 = f;
    // }
    // void set_deltaE0(double (*f)(const double*)){
    //     deltaE0 = f;
    // }
    // void set_deltaE1(double (*f)(const double*)){
    //     deltaE1 = f;
    // }
    // void set_deltaE2(double (*f)(const double*)){
    //     deltaE2 = f;
    // }

    /* 
        Convolution functions
    */

    double calculate_K_xy(double x, double xx, double y, double yy) const{
        double var = 0.2;
        return 1.0/(var*sqrt(2*M_PI))*exp(- ((x-xx)*(x-xx)+(y-yy)*(y-yy))/(2*(var*var)));
    }


    double calculate_convval(const double* rho, const int i, const int j) const{
        double convval = 0;

        double xx=(j+0.5)/n1;
        double yy=(i+0.5)/n2;

        for(int i1=i-convN/2;i1<i+convN/2+1;++i1){
            for(int j1=j-convN/2;j1<j+convN/2+1;++j1){
        // for(int i1=fmax(0,i-convN/2);i1<fmin(n2,i+convN/2+1);++i1){
        //     for(int j1=fmax(0,j-convN/2);j1<fmin(n1,j+convN/2+1);++j1){
                double x=(j1+0.5)/n1;
                double y=(i1+0.5)/n2;        
                if((x-xx)*(x-xx) + (y-yy)*(y-yy) <= conv_r*conv_r){
	                double eval = calculate_K_xy(x,xx,y,yy);

	                int ii = fmin(n2-1, fmax(0, i1));
	                int jj = fmin(n1-1, fmax(0, j1));

	                convval += eval * rho[ii*n1+jj];
	            }
            }
        }

        convval/=conv_sum;
        return convval;
    }

    double calculate_convval2(const double* rho, const double* phi, const int i, const int j) const{
        double convval = 0;

        double xx=(j+0.5)/n1;
        double yy=(i+0.5)/n2;

        for(int i1=i-convN/2;i1<i+convN/2+1;++i1){
            for(int j1=j-convN/2;j1<j+convN/2+1;++j1){
                
                double x=(j1+0.5)/n1;
                double y=(i1+0.5)/n2;        
                if((x-xx)*(x-xx) + (y-yy)*(y-yy) <= conv_r*conv_r){
	                double eval = calculate_K_xy(x,xx,y,yy) ;

	                int ii = fmin(n2-1, fmax(0, i1));
	                int jj = fmin(n1-1, fmax(0, j1));

	                convval += eval * rho[ii*n1+jj] * phi[ii*n1+jj];
                }
                
            }
        }
        convval/=conv_sum;
        return convval;
    }

    /*
        Updaate m0, m1, m2;
    */

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

    /* 
        Update rho
    */

    void calculate_rho_related(double& mvalue, double& Dtphi, double& Deltaphi, const int n, const int i, const int j, const double* mx, const double* my, const double* phi){

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
        mvalue = sqrt(mxvalue*mxvalue + myvalue*myvalue);

        // if(n==nt-1){
        //     Dtphi=1.0/dt*(phi[n*n1*n2+i*n1+j]-phi[(n-1)*n1*n2+i*n1+j]);
        // }else{
        //     Dtphi=1.0/dt*(phi[(n+1)*n1*n2+i*n1+j]-phi[(n)*n1*n2+i*n1+j]);
        // }

        // if(n==nt-1){
        // 	Dtphi=0.5/dt*(phi[(n)*n1*n2+i*n1+j]-phi[(n-1)*n1*n2+i*n1+j]);
        // }else{
        // 	Dtphi=0.5/dt*(phi[(n+1)*n1*n2+i*n1+j]-phi[(n-1)*n1*n2+i*n1+j]);	
        // }

        Dtphi=1.0/dt*(phi[n*n1*n2+i*n1+j]-phi[(n-1)*n1*n2+i*n1+j]);
        
        

        Deltaphi = -1.0/(dx*dx) * (-phi[n*n1*n2+i*n1+(int) fmax(0,j-1)]+2*phi[n*n1*n2+i*n1+j]-phi[n*n1*n2+i*n1+(int) fmin(n1-1,j+1)])
                   -1.0/(dy*dy) * (-phi[n*n1*n2+(int) fmax(0,i-1)*n1+j]+2*phi[n*n1*n2+i*n1+j]-phi[n*n1*n2 +(int) fmin(n2-1,i+1)*n1+j]);

    }
    double calculate_E0prime(const double* rho, const double* f, const int i, const int j) const{
        int ind = (nt-1)*n1*n2+i*n1+j;
        // rho log rho
        // double eval = 1.0/c0 * (phi[0][ind] - f[i*n1+j]);
        // double val = 0;
        // if(f[i*n1+j] >= 0) val = exp(eval);

        double eval = 1.0/c0 * (phi[0][ind] - f[i*n1+j]) + rho[i*n1+j];

        if(eval <= 0) eval = 0;

        double val = 0;
        if(f[i*n1+j] >= 0 ){
            val = eval * phi[0][ind] - 0.5 * c0 * (eval - rho[i*n1+j]) * (eval - rho[i*n1+j]);
        }
        return val;
    }
    double calculate_E1prime(const double* rho, const double* f, int i, int j) const{
        int ind = (nt-1)*n1*n2+i*n1+j;
        double eval = 1.0/c1 * (phi[1][ind] - f[i*n1+j]);
        double val = 0;
        // if(eval > 0 && f[i*n1+j] >= 0) val = 1.0/(2*c1_value) * pow(fmax(0, (phi[1][ind]-f[i])), 2);

        if(eval > 0 && f[i*n1+j] >= 0) val = 0.5 * eval * (phi[1][ind]-f[i]);

        return val;
    }
    double calculate_E2prime(const double* rho, const double* f, int i, int j) const{
        int ind = (nt-1)*n1*n2+i*n1+j;
        // rho log rho
        double eval = 1.0/c2 * (phi[2][ind] - f[i*n1+j]);
        double val  = 0;
        if(f[i*n1+j] >= 0) val = c2*exp(eval);

        // double eval = 1.0/c2 * (phi[2][ind] - f[i*n1+j]) + rho[i*n1+j];

        // if(eval <= 0) eval = 0;

        // double val = 0;
        // if(f[i*n1+j] >= 0 ){
        //     val = eval * phi[2][ind] - 0.5 * c2 * (eval - rho[i*n1+j]) * (eval - rho[i*n1+j]);
        // }

        return val;
    }


    double calculate_deltaE0prime(const double* rho, const double* f, int i, int j){
        int ind = (nt-1)*n1*n2+i*n1+j;
        // rho log rho
        // double eval = 1.0/c0 * (phi[0][ind] - f[i*n1+j]);
        // double val = 0;
        // if(f[i*n1+j] >= 0) val = exp(eval);

        double eval = 1.0/c0 * (phi[0][ind] - f[i*n1+j]) + rho[i*n1+j];
        double val = 0;
        if(c0>0 && eval > 0 && f[i*n1+j] >= 0) val = eval;
        return val;
    }
    double calculate_deltaE1prime(const double* rho, const double* f, int i, int j){
        int ind = (nt-1)*n1*n2+i*n1+j;
        double eval = 1.0/c1 * (phi[1][ind] - f[i*n1+j]);
        double val = 0;
        if(eval > 0 && f[i*n1+j] >= 0) val = eval;
        return val;
    }
    double calculate_deltaE2prime(const double* rho, const double* f, int i, int j){
        int ind = (nt-1)*n1*n2+i*n1+j;
        // rho log rho
        double eval = 1.0/c2 * (phi[2][ind] - f[i*n1+j]);
        double val = 0;
        if(f[i*n1+j] >= 0) val = exp(eval);

        // double eval = 1.0/c2 * (phi[2][ind] - f[i*n1+j]) + rho[i*n1+j];
        // double val = 0;
        // if(c2>0 && eval > 0 && f[i*n1+j] >= 0) val = eval;
        return val;
    }

    void update_rho0(double* rho0,const double* rho1, const double* rho2,const double* mx,const double* my,const double* f){


        for(int n=1;n<nt;++n){

            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int ind = n*n1*n2+i*n1+j;

                    double mvalue=0;
                    double Dtphi =0;
                    double Deltaphi=0;

                    calculate_rho_related(mvalue, Dtphi, Deltaphi, n, i, j, mx, my, phi[0]);

                    double convval  =  calculate_convval(&rho1[n*n1*n2],i,j);
                    double convval2 = calculate_convval2(&rho1[n*n1*n2],&phi[1][n*n1*n2],i,j);

                    // double newrhovalue=cubic_solve(tau*Dtphi + tau*Clist[0]*Deltaphi - rho0[ind] + tau*beta*rho1[ind]*(phi[1][ind]-phi[0][ind]), 0, -0.5/h*tau*alphalist[0]*mvalue*mvalue);
                    double newrhovalue=cubic_solve(tau*Dtphi + tau*Clist[0]*Deltaphi - rho0[ind] + tau*beta*(convval2 - phi[0][ind]*convval), 0, -0.5/h*tau*alphalist[0]*mvalue*mvalue);
                    rho0[n*n1*n2+i*n1+j]=fmax(0,newrhovalue);
                }
            }
        }



        // ----- rho(1,x) = dGstar(phi(1,x)) ------

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                int ind = (nt-1)*n1*n2+i*n1+j;
                double val = calculate_deltaE0prime(rho0, f, i, j);

                // rho0[ind] = 1.0/(1.0+tau*0.5)*(rho0[ind] + tau*0.5*val);  
                rho0[ind] += 0.05 * (val - rho0[ind]);
                // rho0[ind] = val;
            }
        }


    }

    void update_rho1(const double* rho0,double* rho1, const double* rho2,const double* mx,const double* my,const double* f){

        for(int n=1;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int ind = n*n1*n2+i*n1+j;

                    
                    double mvalue=0;
                    double Dtphi =0;
                    double Deltaphi=0;

                    calculate_rho_related(mvalue, Dtphi, Deltaphi, n, i, j, mx, my, phi[1]);

                    double convval  = calculate_convval(&rho0[n*n1*n2],i,j);
                    double convval2 = calculate_convval2(&rho0[n*n1*n2],&phi[0][n*n1*n2],i,j);

                    // double newrhovalue=cubic_solve(tau*Dtphi + tau*Clist[1]*Deltaphi - rho1[ind] + tau*beta*rho0[ind]*(phi[1][ind]-phi[0][ind]) + tau*gamma*(phi[2][ind] - phi[1][ind]), 0, -0.5/h*tau*alphalist[1]*mvalue*mvalue);
                    double newrhovalue=cubic_solve(tau*Dtphi + tau*Clist[1]*Deltaphi - rho1[ind] + tau*beta*(phi[1][ind]*convval - convval2) + tau*gamma*(phi[2][ind] - phi[1][ind]), 0, -0.5/h*tau*alphalist[1]*mvalue*mvalue);
                    rho1[ind]=fmax(0,newrhovalue);
                }
            }
        }

        // ----- rho(1,x) = dGstar(phi(1,x)) ------

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                int ind = (nt-1)*n1*n2+i*n1+j;
                double val = calculate_deltaE1prime(rho1, f, i, j);

                // rho1[ind] = 1.0/(1.0+tau)*(rho1[ind] + tau*0.5*val);  
                rho1[ind] += 0.05 * (val - rho1[ind]);
                // rho1[ind] = val;
            }
        }


    }

    void update_rho2(const double* rho0, const double* rho1, double* rho2,const double* mx,const double* my,const double* f){


        for(int n=1;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int ind = n*n1*n2+i*n1+j;

                    
                    double mvalue=0;
                    double Dtphi =0;
                    double Deltaphi=0;

                    calculate_rho_related(mvalue, Dtphi, Deltaphi, n, i, j, mx, my, phi[2]);

                    double newrhovalue=cubic_solve(tau*Dtphi + tau*Clist[2]*Deltaphi - rho2[ind], 0, -0.5/h*tau*alphalist[2]*mvalue*mvalue);
                    rho2[ind]=fmax(0,newrhovalue);
                }
            }
        }


        // ----- rho(1,x) = dGstar(phi(1,x)) ------

        if(c2 > 0){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int ind = (nt-1)*n1*n2+i*n1+j;
                    double val = calculate_deltaE2prime(rho2, f, i, j);

                    // rho2[ind] = 1.0/(1.0+tau*0.5)*(rho2[ind] + tau*0.5*val);  
                    rho2[ind] += 0.05 * (val - rho2[ind]);     
                    // rho2[ind] = val;
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



    double calculate_dtrho(const double* rho, const int n, const int i, const int j){
    	double dtrho=0;
    	if(n==nt-1){
            dtrho=0;
        }else{
            dtrho=1.0/dt*(rho[(n+1)*n1*n2+i*n1+j]-rho[(n)*n1*n2+i*n1+j]); 
        }

        // if(n==0){
        //     dtrho=0.5/dt*(rho[(n+1)*n1*n2+i*n1+j]-rho[(n)*n1*n2+i*n1+j]); 
        // }else if(n==nt-1){
        // 	dtrho=0.5/dt*(rho[(n)*n1*n2+i*n1+j]-rho[(n-1)*n1*n2+i*n1+j]); 
        // }else{
        //     dtrho=0.5/dt*(rho[(n+1)*n1*n2+i*n1+j]-rho[(n-1)*n1*n2+i*n1+j]); 
        // }
        return dtrho;
    }

    void update_phi0(double* const rho[], const double* mx, const double* my, const double* f){

        // Update 0 < t < 1

        for(int i=0;i<n2;++i){
        	for(int j=0;j<n1;++j){
        		int ind = (nt-1)*n1*n2+i*n1+j;
        		double val = c0 * (rho[0][ind] - rho[0][i*n1+j]);
        		phi[0][ind] += sigma/nt * (val - phi[0][ind]);
        	}
        }

        for(int n=0;n<nt;++n){  
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double dtrho = calculate_dtrho(rho[0], n, i, j);

                    double nablamx=calculate_grad_mx(mx,n,i,j);
                    double nablamy=calculate_grad_my(my,n,i,j);

                    int ind = n*n1*n2+i*n1+j;

                    double Deltarho = -1.0/(dx*dx) * (-rho[0][n*n1*n2+i*n1+(int) fmax(0,j-1)]+2*rho[0][n*n1*n2+i*n1+j]-rho[0][n*n1*n2+i*n1+(int) fmin(n1-1,j+1)])
                                      -1.0/(dy*dy) * (-rho[0][n*n1*n2+(int) fmax(0,i-1)*n1+j]+2*rho[0][n*n1*n2+i*n1+j]-rho[0][n*n1*n2 +(int) fmin(n2-1,i+1)*n1+j]);

                    double convval = calculate_convval(&rho[1][n*n1*n2],i,j);

                    // fftps->workspace[n*n1*n2+i*n1+j]=-(dtrho+nablamx+nablamy + beta*rho[0][ind]*rho[1][ind] - Clist[0]*Deltarho); 
                    fftps->workspace[n*n1*n2+i*n1+j]=-(dtrho+nablamx+nablamy + beta*rho[0][ind]*convval - Clist[0]*Deltarho); 
                }
            }
        }

        fftps->perform_inverse_laplacian(c0);
        for(int i=0;i<n1*n2*nt;++i){
            phi[0][i] += sigma*fftps->workspace[i];
        }
    }

    void update_phi1(double* const rho[], const double* mx, const double* my, const double* f){

        // Update 0 < t < 1

        for(int i=0;i<n2;++i){
        	for(int j=0;j<n1;++j){
        		int ind = (nt-1)*n1*n2+i*n1+j;
        		double val = c1 * rho[1][ind];
        		phi[1][ind] += sigma/nt * (val - phi[1][ind]);
        	}
        }

        for(int n=0;n<nt;++n){  
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double dtrho = calculate_dtrho(rho[1], n, i, j);


                    double nablamx=calculate_grad_mx(mx,n,i,j);
                    double nablamy=calculate_grad_my(my,n,i,j);

                    int ind = n*n1*n2+i*n1+j;

                    double Deltarho = -1.0/(dx*dx) * (-rho[1][n*n1*n2+i*n1+(int) fmax(0,j-1)]+2*rho[1][n*n1*n2+i*n1+j]-rho[1][n*n1*n2+i*n1+(int) fmin(n1-1,j+1)])
                                      -1.0/(dy*dy) * (-rho[1][n*n1*n2+(int) fmax(0,i-1)*n1+j]+2*rho[1][n*n1*n2+i*n1+j]-rho[1][n*n1*n2 +(int) fmin(n2-1,i+1)*n1+j]);

                    double convval = calculate_convval(&rho[0][n*n1*n2],i,j);

                    // fftps->workspace[n*n1*n2+i*n1+j]=-(dtrho+nablamx+nablamy - beta*rho[0][ind]*rho[1][ind] + gamma*rho[1][ind] - Clist[1]*Deltarho); 
                    fftps->workspace[n*n1*n2+i*n1+j]=-(dtrho+nablamx+nablamy - beta*convval*rho[1][ind] + gamma*rho[1][ind] - Clist[1]*Deltarho); 
                }
            }
        }

        fftps->perform_inverse_laplacian(c1);
        for(int i=0;i<n1*n2*nt;++i){
            phi[1][i] += sigma*fftps->workspace[i];
        }
    }

    void update_phi2(double* const rho[], const double* mx, const double* my, const double* f){

        // Update 0 < t < 1

        for(int i=0;i<n2;++i){
        	for(int j=0;j<n1;++j){
        		int ind = (nt-1)*n1*n2+i*n1+j;

        		double val = 0;
        		if(rho[2][ind] > 0){
        			val = c2*log(rho[2][ind]);
        		}
        		phi[2][ind] += sigma/nt * (val - phi[2][ind]);
        	}
        }    	

        // Update 0 < t < 1

        for(int n=0;n<nt;++n){  
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double dtrho = calculate_dtrho(rho[2], n, i, j);

                    double nablamx=calculate_grad_mx(mx,n,i,j);
                    double nablamy=calculate_grad_my(my,n,i,j);

                    int ind = n*n1*n2+i*n1+j;


                    double Deltarho = -1.0/(dx*dx) * (-rho[2][n*n1*n2+i*n1+(int) fmax(0,j-1)]+2*rho[2][n*n1*n2+i*n1+j]-rho[2][n*n1*n2+i*n1+(int) fmin(n1-1,j+1)])
                                      -1.0/(dy*dy) * (-rho[2][n*n1*n2+(int) fmax(0,i-1)*n1+j]+2*rho[2][n*n1*n2+i*n1+j]-rho[2][n*n1*n2 +(int) fmin(n2-1,i+1)*n1+j]);


                    fftps->workspace[n*n1*n2+i*n1+j]=-(dtrho + nablamx + nablamy - gamma*rho[1][ind] - Clist[2]*Deltarho); 
                }
            }
        }

        fftps->perform_inverse_laplacian(c2);
        for(int i=0;i<n1*n2*nt;++i){
            phi[2][i] += sigma*fftps->workspace[i];
        }
    }

    double calculate_energy_num(double* const rho[],double* const  mx[], double* const my[], const int num) const{
        double sum1=0;

        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int ind = n*n1*n2+i*n1+j;
                    double mxval=mx[num][ind];
                    double myval=my[num][ind];

                    double mval=sqrt(mxval*mxval+myval*myval);

                    double rhoval=rho[num][ind];

                    if(rhoval>0){
                        sum1 += alphalist[num]*mval*mval/(2.0*rhoval);
                    }
                }
            }
        }

        return sum1*dx*dy*dt;
    }


    double calculate_energy(double* const rho[]) const{
        double sum = 0;
        for(int i=0;i<3;++i){
            sum += calculate_energy_num(rho,mx,my,i);
        }

        double sum1 = 0;
        for(int i=0;i<n1*n2;++i){
            int ind = (nt-1)*n1*n2+i;
            // if(rho[0][ind] > 0) sum1 += c0*(rho[0][ind]*log(rho[0][ind]) - rho[0][ind]);
            if(rho[0][ind] >= 0) sum1 += c0/2.0*(rho[0][ind] - rho[0][i])*(rho[0][ind] - rho[0][i]);
            if(rho[1][ind] >= 0) sum1 += c1/m*rho[1][ind]*rho[1][ind];
            if(rho[2][ind] >  0) sum1 += c2*(rho[2][ind]*log(rho[2][ind]) - rho[2][ind]);
        }
        sum1 /= 1.0*n1*n2;
        
        return sum + sum1;
    }

    double calculate_dual(double* const rho[], const double* f) const{   

        double term0=0;
        double term1=0;
        double term2=0;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                if(f[i*n1+j]>=0){
                    // if(c0>0) term0 += phi[0][i]*rho[0][i] - c0*exp((phi[0][(nt-1)*n1*n2+i]-f[i])/c0);
                    // term1 += phi[1][i]*rho[1][i] - 1.0/(mprime*c1_value) * pow(fmax(0, (phi[1][(nt-1)*n1*n2+i]-f[i])), mprime);
                    // if(c2>0) term2 += phi[2][i]*rho[2][i] - c2*exp((phi[2][(nt-1)*n1*n2+i]-f[i])/c2);    

                    if(c0>0) term0 += phi[0][i]*rho[0][i] - calculate_E0prime(rho[0],f,i,j);
                    term1 += phi[1][i]*rho[1][i] - calculate_E1prime(rho[1],f,i,j);
                    if(c2>0) term2 += phi[2][i]*rho[2][i] - calculate_E2prime(rho[2],f,i,j);
                }    
            }
            
        }

        double term3=0;
        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double convval0 = calculate_convval(&rho[0][n*n1*n2],i,j);
                    double convval1 = calculate_convval(&rho[1][n*n1*n2],i,j);
                    term3 += beta * (rho[0][n*n1*n2+i*n1+j] * convval1 * phi[0][n*n1*n2+i*n1+j] - rho[1][n*n1*n2+i*n1+j] * convval0 * phi[1][n*n1*n2+i*n1+j]);
                }
            }
        }

        return (term0+term1+term2)/(n1*n2) + term3/(n1*n2*nt);
    }

    void display_log(const int iterPDHG, const double tau, const double sigma, const double energy, const double dual, const double error, const double dual_gap) const{
        cout<<"iter : "<< setw(8) << iterPDHG+1<< " tau : "<< setw(13) <<tau <<" sigma : " << setw(13)<<sigma <<" energy : " << setw(13)<<energy <<" dual : " << setw(13) << dual << " relative_error : " << setw(13) << error <<" dual_gap : " << setw(13) << dual_gap << endl;  
    }

    void run(double* rho[], const double* f, int skip=1){

        previous_energy=1;
        previous_dual=1;
        double error=1, dual_gap=1, energy=1, dual=0;
        int iterPDHG;

        double beta_1 = 1.2;
        double beta_2 = 1;

        for(iterPDHG=0; iterPDHG<max_iteration; ++iterPDHG){

            // get the data before updates
            for(int i=0;i<3;++i){
                memcpy(rhotmps[i],rho[i],n1*n2*nt*sizeof(double));
                memcpy(mxtmps[i],mx[i],n1*n2*nt*sizeof(double));
                memcpy(mytmps[i],my[i],n1*n2*nt*sizeof(double));    
            }

            update_rho0(rho[0],rho[1],rho[2],mx[0],my[0],f);
            update_rho1(rho[0],rho[1],rho[2],mx[1],my[1],f);
            update_rho2(rho[0],rho[1],rho[2],mx[2],my[2],f);            

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
            update_phi0(rhotmps,mxtmps[0],mytmps[0],f); 
            update_phi1(rhotmps,mxtmps[1],mytmps[1],f); 
            update_phi2(rhotmps,mxtmps[2],mytmps[2],f); 

            energy=calculate_energy(rho);
            error=fabs((energy-previous_energy)/previous_energy);
            previous_energy=energy;

            dual  =calculate_dual(rho,f);
            double relative_dual_error=fabs((dual-previous_dual)/previous_dual);
            previous_dual=dual;

            if(error > relative_dual_error*beta_1){
                tau   *= beta_2;
                // sigma /= beta_2;
            }

            if(error < relative_dual_error/beta_1){
                tau   /= beta_2;
                // sigma *= beta_2;
            }

            tau = fmax(1e-5,fmin(10, tau));
            sigma = fmax(1e-5,fmin(10, sigma));

            dual_gap = energy-dual;

                // adaptive step size

            if((iterPDHG+1)%skip==0){

                display_log(iterPDHG,  tau,  sigma,  energy,  dual,  error,  dual_gap);

                create_csv_file(rho[0],"./data/rho0.csv",n1,n2,nt);
                create_csv_file(rho[1],"./data/rho1.csv",n1,n2,nt);
                create_csv_file(rho[2],"./data/rho2.csv",n1,n2,nt);
            }
            if((iterPDHG>20 ) && (fabs(dual_gap)<tolerance)) break;
        }

        cout<<"The method is done!!"<<endl;
        display_log(iterPDHG,  tau,  sigma,  energy,  dual,  error,  dual_gap);

        cout<<"Iterations     : "<<iterPDHG<<endl;
        cout<<"Energy         : "<<energy<<endl;
        cout<<"Relative Error : "<<error<<endl;
    }

}; // Method class


#endif