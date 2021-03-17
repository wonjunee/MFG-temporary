#ifndef METHOD_H
#define METHOD_H
#include <iostream>

#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstring>
// #include <omp.h> // parallelization
#include <fftw3.h>
// #include "poisson_solver.h"
#include "poisson_solver_3d.h"
#include "helper.h"
#include "initializer.h"

class Method{
public:
    int n1;
    int n2;
    int nt;

    int im;
    int ip;
    int jm;
    int jp;

    double c_rho_V;

    double* tau;
    double* sigma;

    double v_rate;
    double v_eff;

    int max_iteration;
    double tolerance;

    double* mx[4];
    double* my[4];
    double* phi[4];
    double* phitmps[4];
    double* rhotmps[4];

    double* rhont0tmps[4];

    double* factory;

    double* convarr;

    double TOTAL_SUM;

    double c_rho_sum_; // constant for (rho_S + rho_I + rho_R)^2

    double energy;
    double previous_energy;
    double previous_dual;

    double m;

    double M0;

    double* etalist;
    double* alphalist;

    poisson_solver** fftps;
    poisson_solver_DST* fftpsDST;
    poisson_solver_2d* fftps2d;

    // For SIR model

    double beta;
    double gamma;
    double var;
    double var_gamma;

    // function pointers for terminal functions
    // double (*E0)(double);
    // double (*E1)(double);
    // double (*E2)(double);

    // double (*deltaE0)(const double* rho);
    // double (*deltaE1)(const double* rho);
    // double (*deltaE2)(const double* rho);

    // ------------------------

    Method(){
        for(int i=0;i<4;++i){
            phitmps[i]=nullptr;
            rhotmps[i]=nullptr;
            rhont0tmps[i]=nullptr;
            mx[i]=nullptr;
            my[i]=nullptr;
            phi[i]=nullptr;
        }
        
        etalist=nullptr;
        alphalist=nullptr;
        fftps=nullptr;
        fftps2d=nullptr;
    }

    Method(int n1, int n2, int nt, double tau, double sigma, int max_iteration, double tolerance, 
           double alphalist[4], double* etalist, double var){

        this->n1=n1;
        this->n2=n2;
        this->nt=nt;
        this->max_iteration=max_iteration;
        this->tolerance=tolerance;

        this->tau  =new double[4];
        this->sigma=new double[4];

        this->tau[0] = tau;
        this->tau[1] = tau;
        this->tau[2] = tau;
        this->tau[3] = tau;

        this->sigma[0] = sigma;
        this->sigma[1] = sigma;
        this->sigma[2] = sigma;
        this->sigma[3] = sigma;

        this->etalist = etalist; // this contains the eta values for viscosity terms
        this->var     = var;

        this->alphalist = new double[4];
        this->alphalist[0] = alphalist[0];
        this->alphalist[1] = alphalist[1];
        this->alphalist[2] = alphalist[2];
        this->alphalist[3] = alphalist[3];

        printf("Alpha values\n");
        for(int k=0;k<4;++k){
            printf("k: %d alpha: %f\n", k, this->alphalist[k]);
        }

        c_rho_sum_ = 0.4;
        c_rho_V    = 0.4;

        convarr  = new double[n1*n2];

        // factory is a function in space and time that represents the location of factories and the vaccines produced in factories.
        // It is a nonnegative function only positive on factory locations.
        factory = new double[n1*n2*nt];

        // defining factory
        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double x = (j+.5)/n1;
                    double y = (i+.5)/n2;
                    // factory[n*n1*n2+i*n1+j] = fmax(0, 0.8 * exp(-30*pow(x-0.5,2)-30*pow(y-0.2,2)) - 0.6);
                    factory[n*n1*n2+i*n1+j] = 0;
                }
            }
        }


        // Vaccine efficiency, i.e. the rate people can be vaccinated
        v_eff  = 0.5;
        // Vaccine rate, i.e. the number of vaccines the doctors can use per day.
        v_rate = 0.5;

        var_gamma = 0.01; // DON'T TOUCH THIS. This one is for regularization.

        for(int i=0;i<4;++i){
            mx[i]      = new double[n1*n2*nt];
            my[i]      = new double[n1*n2*nt];
            phi[i]     = new double[n1*n2*nt];
            phitmps[i] = new double[n1*n2*nt];
            rhotmps[i] = new double[n1*n2*nt];
            rhont0tmps[i] = new double[n1*n2];
        }    

        for(int k=0;k<4;++k){
            memset(mx[k],      0, n1*n2*nt*sizeof(double));
            memset(my[k],      0, n1*n2*nt*sizeof(double));
            memset(phi[k],     0, n1*n2*nt*sizeof(double));
            memset(phitmps[k], 0, n1*n2*nt*sizeof(double));
        }

        clock_t t;
        t = clock();
            
        double eta = etalist[0];
        fftps    = new poisson_solver*[4];
        fftps[0] = new poisson_solver(n1,n2,nt,eta);
        fftps[1] = new poisson_solver(n1,n2,nt,eta);
        fftps[2] = new poisson_solver(n1,n2,nt,eta);
        fftps[3] = new poisson_solver(n1,n2,nt,eta);
        fftpsDST = new poisson_solver_DST(n1,n2,nt,eta);
        fftps2d = new poisson_solver_2d(n1,n2,eta);

        t = clock() - t;
        printf ("\nCPU time for setting up FFT: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
    }


    ~Method(){
        for(int k=0;k<4;++k){
            delete[] mx[k];
            delete[] my[k];
            delete[] rhotmps[k];
            delete[] rhont0tmps[k];
            delete[] phitmps[k];
            delete[] phi[k];
            delete fftps[k];
        }
        delete[] alphalist;

        delete[] tau;
        delete[] sigma;
        delete[] fftps;
        delete[] factory;
        delete fftps2d;

        delete[] convarr;
        
    }

    void setup_indices(int& im, int& ip, int& jm, int& jp, const int i, const int j){
        im = fmax(0,i-1);
        ip = fmin(n2-1,i+1);

        jm = fmax(0,j-1);
        jp = fmin(n1-1,j+1);
    }

    void perform_upwind_scheme(double& muxp, double& muxm, double& muyp, double& muym, const double* phi, const int i, const int j){

        setup_indices(im,ip,jm,jp,i,j);

        muxp = 1.0 * n1 * (phi[i*n1+jp] - phi[i*n1+j]);
        muxm = 1.0 * n1 * (phi[i*n1+j] - phi[i*n1+jm]);
        muyp = 1.0 * n2 * (phi[ip*n1+j] - phi[i*n1+j]);
        muym = 1.0 * n2 * (phi[i*n1+j] - phi[im*n1+j]);
    }

    // centered difference
    void update_m(double* mx, double* my, const double* rho, const double* phi, const int SIR){

        int im,ip,jm,jp;

        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double nablaxphi = 0;
                    double nablayphi = 0;

                    setup_indices(im,ip,jm,jp,i,j);

                    // nablaxphi = 0.5*n1*(phi[n*n1*n2+i*n1+jp]-phi[n*n1*n2+i*n1+jm]);
                    // nablayphi = 0.5*n2*(phi[n*n1*n2+ip*n1+j]-phi[n*n1*n2+im*n1+j]);

                    nablaxphi = 1.0*n1*(phi[n*n1*n2+i*n1+jp]-phi[n*n1*n2+i*n1+j]);
                    nablayphi = 1.0*n2*(phi[n*n1*n2+ip*n1+j]-phi[n*n1*n2+i*n1+j]);
                    
                    double rhovalx=rho[n*n1*n2+i*n1+j];
                    double rhovaly=rho[n*n1*n2+i*n1+j];
                    
                    mx[n*n1*n2+i*n1+j] = rhovalx/(tau[SIR] * alphalist[SIR] + rhovalx) * (mx[n*n1*n2+i*n1+j] - tau[SIR] * nablaxphi);
                    my[n*n1*n2+i*n1+j] = rhovaly/(tau[SIR] * alphalist[SIR] + rhovaly) * (my[n*n1*n2+i*n1+j] - tau[SIR] * nablayphi);
                }
            }   
        }   
    }

    /* 
        Update rho
    */

    void calculate_rho_related(double& mvalue, double& Dtphi, double& Deltaphi, const int n, const int i, const int j, const double* mx, const double* my, const double* phi){

        double mxvalue, myvalue;
        int im,ip,jm,jp;
        setup_indices(im,ip,jm,jp,i,j);

        mxvalue = mx[n*n1*n2+i*n1+j];
        myvalue = my[n*n1*n2+i*n1+j];
        
        mvalue = sqrt(mxvalue*mxvalue + myvalue*myvalue);

        if(n>0) Dtphi=1.0*nt*(phi[n*n1*n2+i*n1+j]-phi[(n-1)*n1*n2+i*n1+j]);
        else    Dtphi=0;
        Deltaphi = - n1*n1 * (-phi[n*n1*n2+i*n1+jm]+2*phi[n*n1*n2+i*n1+j]-phi[n*n1*n2+i*n1+jp])
                   - n2*n2 * (-phi[n*n1*n2+im*n1+j]+2*phi[n*n1*n2+i*n1+j]-phi[n*n1*n2+ip*n1+j]);

    }
    double calculate_Eprime(const double* rho, const double* obstacle, const int i, const int j, const double constant=1) const{
        int ind = (nt-1)*n1*n2+i*n1+j;
        double eval = 1.0/constant * (phi[1][ind] - obstacle[i*n1+j]);
        return fmax(0, eval);
    }
    
    double calculate_deltaEprime(const double* phi, const double* obstacle, const double constant, int i, int j){
        double eval = 1.0/constant * (phi[i*n1+j] - obstacle[i*n1+j]);
        if(obstacle[i*n1+j] > 0) return 0;
        return fmax(eval,0);
    }

    void update_rho0(double* rho0,const double* rho1, const double* rho2, const double* rho3,const double* mx,const double* my,const double* obstacle){

    	double newrhovalue = -1;
		double tauval = tau[0];
		double betaval  = beta;

        memcpy(&rho0[0], rhont0tmps[0], n1*n2*sizeof(double));

		for(int n=1;n<nt;++n){

            // calculate convolution
            fftps2d->perform_convolution(&rho1[n*n1*n2],var);
            for(int i=0;i<n1*n2;++i){
                convarr[i] = fftps2d->workspace[i];
            }

            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int ind = n*n1*n2+i*n1+j;

                    double mvalue=0;
                    double Dtphi =0;
                    double Deltaphi=0;

                    calculate_rho_related(mvalue, Dtphi, Deltaphi, n, i, j, mx, my, phi[0]);
                    double convval  =  convarr[i*n1+j];
                    double aval = 1.0/(tauval*c_rho_sum_+1) * 
                                (
                                    tauval * (Dtphi + etalist[0]*Deltaphi + c_rho_sum_*rho1[ind] + c_rho_sum_*rho2[ind]) 
                                    - rho0[ind] 
                                    + tauval * (
                                                betaval*convval*(phi[1][ind] - phi[0][ind])
                                              + v_eff * rho3[ind] * (phi[2][ind] - phi[0][ind])
                                              - v_rate * phi[3][ind] * rho3[ind]
                                              + obstacle[i*n1+j]
                                              )
                                );
                    double cval = -0.5*tauval/(tauval*c_rho_sum_+1)*alphalist[0]*mvalue*mvalue;

                    if(n==0)         aval += 1.0/(tauval*c_rho_sum_+1) * (tauval * phi[0][ind] * nt);
                    // else if(n==nt-1) aval -= 1.0/(tauval*c_rho_sum_+1) * (tauval * phi[0][ind] * nt);
                	newrhovalue=cubic_solve(aval, 0, cval);
                    rho0[n*n1*n2+i*n1+j] = fmin(1,fmax(0,newrhovalue));
                }
            }
            fftps2d->solve_heat_equation(&rho0[n*n1*n2], etalist[0]);
		}
    }

    void update_rho1(const double* rho0,double* rho1, const double* rho2,const double* mx,const double* my,const double* obstacle){

    	double tauval = tau[1];

        memcpy(&rho1[0], rhont0tmps[1], n1*n2*sizeof(double));

        for(int n=1;n<nt;++n){

        	double gammaval = gamma;
        	double betaval  = beta;

            // calculate convolution
            fftps2d->perform_convolution(&rho0[n*n1*n2],&phi[0][n*n1*n2],&phi[1][n*n1*n2],var);
            for(int i=0;i<n1*n2;++i){
                convarr[i] = fftps2d->workspace[i];
            }

            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int ind = n*n1*n2+i*n1+j;

                    
                    double mvalue=0;
                    double Dtphi =0;
                    double Deltaphi=0;

                    calculate_rho_related(mvalue, Dtphi, Deltaphi, n, i, j, mx, my, phi[1]);
                    double convval3 = convarr[i*n1+j];
                    double convval4 = (phi[2][ind] - phi[1][ind]);

                    double newrhovalue = 0;
                    double aval=0,cval=0;

                        aval = 1.0/(tauval*c_rho_sum_+1) * 
                                ( 
                                    tauval * (
                                        Dtphi + etalist[1]*Deltaphi + c_rho_sum_*rho0[ind] + c_rho_sum_*rho2[ind]
                                              + betaval * convval3 + gammaval * convval4
                                              + obstacle[i*n1+j]
                                             )
                                    - rho1[ind] 
                                );

                        if(n==0) aval += 1.0/(tauval*c_rho_sum_+1) * (tauval * phi[1][ind] * nt);
                        cval = -0.5*tauval/(tauval*c_rho_sum_+1)*alphalist[1]*mvalue*mvalue;


                    newrhovalue=cubic_solve(aval, 0, cval);                    
                    rho1[ind]=fmin(1,fmax(0,newrhovalue));
                }
            }
            fftps2d->solve_heat_equation(&rho1[n*n1*n2], etalist[1]);
        }
	}

    void update_rho2(const double* rho0, const double* rho1, double* rho2,const double* mx,const double* my,const double* obstacle){

        memcpy(&rho2[0], rhont0tmps[2], n1*n2*sizeof(double));

        double tauval = tau[2];

        for(int n=1;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int ind = n*n1*n2+i*n1+j;

                    double mvalue=0;
                    double Dtphi =0;
                    double Deltaphi=0;

                    calculate_rho_related(mvalue, Dtphi, Deltaphi, n, i, j, mx, my, phi[2]);

                    double aval  = 1.0/(tauval*c_rho_sum_+1) * (tauval * (Dtphi + etalist[2]*Deltaphi + c_rho_sum_*rho0[ind] + c_rho_sum_*rho1[ind] + + obstacle[i*n1+j]) - rho2[ind]);

                    if(n==0)         aval += 1.0/(tauval*c_rho_sum_+1) * (tauval * phi[2][ind] * nt);
                    else if(n==nt-1) aval -= 1.0/(tauval*c_rho_sum_+1) * (tauval * phi[2][ind] * nt);

                    double newrhovalue=cubic_solve(aval, 0, -0.5*tauval/(tauval*c_rho_sum_+1)*alphalist[2]*mvalue*mvalue);
                    rho2[ind]=fmin(1,fmax(0,newrhovalue));
                }
            }
            fftps2d->solve_heat_equation(&rho2[n*n1*n2], etalist[2]);
        }
    }

    void update_rho3(const double* rho0, const double* rho1, const double* rho2, double* rho3,const double* mx,const double* my,const double* obstacle){

        double tauval = tau[3];
        double betaval  = beta;

        memcpy(&rho3[0], rhont0tmps[3], n1*n2*sizeof(double));

        for(int n=1;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int ind = n*n1*n2+i*n1+j;

                    double mvalue=0;
                    double Dtphi =0;
                    double Deltaphi=0;

                    calculate_rho_related(mvalue, Dtphi, Deltaphi, n, i, j, mx, my, phi[3]);

                    double aval=0,cval=0;

                        aval =  ( 
                                    1.0/(c_rho_V + 1.0/tauval) * (
                                                    Dtphi 
                                                  // + etalist[3]*Deltaphi
                                                  - v_rate * rho0[ind] *  phi[3][ind]
                                                  + v_eff  * rho0[ind] * (phi[2][ind] - phi[0][ind])
                                                  - 1.0/tauval * rho3[ind]
                                                  + obstacle[i*n1+j]
                                             )
                                );
                        if(n==0)         aval += 1.0/(c_rho_V + 1.0/tauval) * phi[3][ind] * nt;
                        // else if(n==nt-1) aval -= tauval * phi[3][ind] * nt;

                        cval = - 0.5*alphalist[3]*mvalue*mvalue * 1.0/(c_rho_V + 1.0/tauval);

                    double newrhovalue=cubic_solve(aval, 0, cval); 
                    rho3[ind]=fmin(1,fmax(0,newrhovalue));
                }
            }
            // fftps2d->solve_heat_equation(&rho3[n*n1*n2], etalist[3]);
        }
    }

    double calculate_grad_mx(const double* mxTmp, int n, int i, int j){
        double mxval; int im,ip,jm,jp;
        setup_indices(im,ip,jm,jp,i,j);

        // mxval = 0.5*n1*(mxTmp[n*n1*n2+i*n1+jp]-mxTmp[n*n1*n2+i*n1+jm]);
        mxval = n1*(mxTmp[n*n1*n2+i*n1+j]-mxTmp[n*n1*n2+i*n1+jm]);

        return mxval;
    }

    double calculate_grad_my(const double* myTmp, int n, int i, int j){
        double myval; int im,ip,jm,jp;
        setup_indices(im,ip,jm,jp,i,j);
        // myval = 0.5*n2*(myTmp[n*n1*n2+ip*n1+j]-myTmp[n*n1*n2+im*n1+j]);        
        myval = n2*(myTmp[n*n1*n2+i*n1+j]-myTmp[n*n1*n2+im*n1+j]);

        return myval;
    }



    double calculate_dtrho(const double* rho, const int n, const int i, const int j){
    	double dtrho=0;

    	if(n==nt-1){
            dtrho=0;
        }else{
            dtrho=1.0*nt*(rho[(n+1)*n1*n2+i*n1+j]-rho[(n)*n1*n2+i*n1+j]); 
        }
        return dtrho;
    }

    double calculate_Delta_value(const double* rho, const int n, const int i, const int j){
        return - n1*n1 * (-rho[n*n1*n2+i*n1+(int) fmax(0,j-1)]+2*rho[n*n1*n2+i*n1+j]-rho[n*n1*n2+i*n1+(int) fmin(n1-1,j+1)])
               - n2*n2 * (-rho[n*n1*n2+(int) fmax(0,i-1)*n1+j]+2*rho[n*n1*n2+i*n1+j]-rho[n*n1*n2 +(int) fmin(n2-1,i+1)*n1+j]);
    }

    double update_phi_all(double* const rho[], double* const mx[], double* const my[], const double* obstacle){

    int n, i, j, ind;
    double dtrho0, nablamx0, nablamy0, dtrho1, nablamx1, nablamy1, dtrho2, nablamx2, nablamy2, Deltarho0, Deltarho1, Deltarho2, convval, convval_gamma;
    double dtrho3, nablamx3, nablamy3, Deltarho3;

        for(int n=0;n<nt;++n){  

            // calculate convolution
            fftps2d->perform_convolution(&rho[1][n*n1*n2],var);
            for(int i=0;i<n1*n2;++i){
                convarr[i] = fftps2d->workspace[i];
            }

            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    ind = n*n1*n2+i*n1+j;

                    dtrho0 = calculate_dtrho(rho[0], n, i, j);
                    nablamx0=calculate_grad_mx(mx[0],n,i,j);
                    nablamy0=calculate_grad_my(my[0],n,i,j);
                    Deltarho0 = calculate_Delta_value(rho[0],n,i,j);

                    dtrho1 = calculate_dtrho(rho[1], n, i, j);
                    nablamx1=calculate_grad_mx(mx[1],n,i,j);
                    nablamy1=calculate_grad_my(my[1],n,i,j);
                    Deltarho1 = calculate_Delta_value(rho[1],n,i,j);

                    dtrho2 = calculate_dtrho(rho[2], n, i, j);
                    nablamx2=calculate_grad_mx(mx[2],n,i,j);
                    nablamy2=calculate_grad_my(my[2],n,i,j);
                    Deltarho2 = calculate_Delta_value(rho[2],n,i,j);

                    dtrho3 = calculate_dtrho(rho[3], n, i, j);
                    nablamx3=calculate_grad_mx(mx[3],n, i, j);
                    nablamy3=calculate_grad_my(my[3],n, i, j);
                    Deltarho3 = calculate_Delta_value(rho[3],n,i,j);

                    convval = convarr[i*n1+j];

                    convval_gamma = rho[1][ind];

                    fftps[0]->u[ind]  = - (dtrho0 + nablamx0 + nablamy0 + beta*rho[0][ind]*convval + v_eff*rho[0][ind]*rho[3][ind]  - etalist[0]*Deltarho0 ); 
                    fftps[1]->u[ind]  = - (dtrho1 + nablamx1 + nablamy1 - beta*rho[0][ind]*convval + gamma*convval_gamma            - etalist[1]*Deltarho1 );
                    fftps[2]->u[ind]  = - (dtrho2 + nablamx2 + nablamy2 - gamma*convval_gamma      - v_eff*rho[0][ind]*rho[3][ind]  - etalist[2]*Deltarho2 ); 
                    fftps[3]->u[ind]  = - (dtrho3 + nablamx3 + nablamy3 - factory[ind]             + v_rate*rho[0][ind]*rho[3][ind] - etalist[3]*Deltarho3*0 ); 
                }
            }
        }

        // phiT
        n = nt-1;
        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                ind = n*n1*n2+i*n1+j;

                fftps[0]->u[ind] -= calculate_deltaEprime(&phi[0][n*n1*n2], obstacle, 1.0, i, j) * nt;
                fftps[0]->u[ind] += rho[0][ind] * nt;

                fftps[1]->u[ind] -= calculate_deltaEprime(&phi[1][n*n1*n2], obstacle, 1.0, i, j) * nt;
                fftps[1]->u[ind] += rho[1][ind] * nt;

                fftps[3]->u[ind] -= calculate_deltaEprime(&phi[3][n*n1*n2], obstacle, 0.2, i, j) * nt;
                fftps[3]->u[ind] += rho[3][ind] * nt;
            }
        }

        double error = 0;

        fftps[0]->perform_inverse_laplacian(beta + v_eff, etalist[0]);
        fftps[1]->perform_inverse_laplacian(beta + gamma, etalist[1]);
        fftps[2]->perform_inverse_laplacian(gamma+ v_rate, etalist[2]);
        fftps[3]->perform_inverse_laplacian(v_rate, 0);

        for(int k=0;k<4;++k){
            for(int i=0;i<n1*n2*nt;++i){
                phi[k][i] += sigma[k]       * fftps[k]->workspace[i];
                error     += fftps[k]->u[i] * fftps[k]->workspace[i];
            }
        }
        return error/(1.0*n1*n2*nt);

    }

    double calculate_energy_num(double* const rho[],double* const  mx[], double* const my[], const int num) const{
        double sum1=0;

        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int ind = n*n1*n2+i*n1+j;

                    double mxval=0.5*mx[num][ind];
                    double myval=0.5*my[num][ind];

                    if(j>0) mxval += 0.5*mx[num][n*n1*n2+i*n1+j-1];
                    if(i>0) myval += 0.5*my[num][n*n1*n2+(i-1)*n1+j];

                    double mval=sqrt(mxval*mxval+myval*myval);

                    double rhoval=rho[num][ind];

                    if(rhoval>0){
                        sum1 += alphalist[num]*mval*mval/(2.0*rhoval);
                    }
                }
            }
        }

        return sum1/(n1*n2*nt);
    }


    double calculate_energy(double* const rho[], const double* obstacle) const{
        double sum = 0;
        for(int k=0;k<4;++k){
            sum += calculate_energy_num(rho,mx,my,k);
        }

        double sum1 = 0;
        for(int i=0;i<n1*n2;++i){
            int ind = (nt-1)*n1*n2+i;
            // if(rho[0][ind] > 0) sum1 += c0*(rho[0][ind]*log(rho[0][ind]) - rho[0][ind]);
            // if(rho[0][ind] >= 0) sum1 += c0/2.0*(rho[0][ind] - rho[0][i])*(rho[0][ind] - rho[0][i]) + rho[0][ind]*obstacle[0][i];
            if(rho[1][ind] >= 0) sum1 += 1.0/m*rho[1][ind]*rho[1][ind] 							    + rho[1][ind]*obstacle[i];
            // if(rho[2][ind] >  0) sum1 += c2*(rho[2][ind]*log(rho[2][ind]) - rho[2][ind]) 		    + rho[2][ind]*obstacle[2][i];
        }
        sum1 /= 1.0*n1*n2;

        double sum2 = 0;
        for(int n=0;n<nt;++n){
            for(int i=0;i<n1*n2;++i){
                double val = rho[0][n*n1*n2+i] + rho[1][n*n1*n2+i] + rho[2][n*n1*n2+i];
                sum2 += 0.5 * c_rho_sum_ * (val*val);
            }
        }
        sum2 /= 1.0*n1*n2*nt;
        
        return sum + sum1;
    }

    double calculate_dual(double* const rho[], const double* obstacle) const{   

        double term0=0;
        double term1=0;
        double term2=0;
        double term3=0;
        double term4=0;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                term0 += phi[0][i*n1+j]*rho[0][i*n1+j] - calculate_Eprime(rho[0],obstacle,i,j,1.0);
                term1 += phi[1][i*n1+j]*rho[1][i*n1+j] - calculate_Eprime(rho[1],obstacle,i,j,1.0);
                term2 += phi[2][i*n1+j]*rho[2][i*n1+j] - phi[2][(nt-1)*n1*n2+i*n1+j]*rho[2][(nt-1)*n1*n2+i*n1+j];
                term3 += phi[3][i*n1+j]*rho[3][i*n1+j] - calculate_Eprime(rho[3],obstacle,i,j,0.2);
            }
            
        }
        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int idx = n*n1*n2+i*n1+j;
                    term4 += beta   * rho[0][idx] * rho[1][idx] * (phi[0][idx] - phi[1][idx]);
                    term4 += gamma  * rho[1][idx] * rho[2][idx] * (phi[1][idx] - phi[2][idx]);
                }
            }
        }

        return (term0+term1+term2+term3)/(n1*n2) + term4/(n1*n2*nt);
    }

    void update_step_sizes(double& tau, double& sigma, const double error, const double relative_dual_error, const double beta_1, const double beta_2){
    	if(error > relative_dual_error*beta_1){
            tau   *= beta_2;
            sigma /= beta_2;
        }

        if(error < relative_dual_error/beta_1){
            tau   /= beta_2;
            sigma *= beta_2;
        }
    }

    void display_log(const int iterPDHG, const double tau, const double sigma, const double energy, const double dual, const double error, const double dual_gap, const double sanity_value) const{
        // cout<<"iter : "<< setw(8) << iterPDHG+1<< " tau : "<< setw(13) <<tau <<" sigma : " << setw(13)<<sigma <<" energy : " << setw(13)<<energy <<" dual : " << setw(13) << dual << " relative_error : " << setw(13) << error <<" dual_gap : " << setw(13) << dual_gap<<" sanity_value : " << setw(13) << sanity_value << endl;  
        printf("iter: %5d tau: %5.2f sigma: %5.2f energy: %10.4e dual: %10.4e rel error: %10.4e dual gap: %10.4e sanity value: %10.4e\n", iterPDHG+1, tau, sigma, energy, dual, error, dual_gap, sanity_value);
    }

    void run(double* rho[], const double* obstacle, int skip=1){

        for(int k=0;k<4;++k){
            memcpy(rhont0tmps[k], &rho[k][0], n1*n2*sizeof(double));
        }

    	TOTAL_SUM = 0;
    	for(int k=0;k<4;++k){
    		for(int i=0;i<n1*n2;++i){
	    		TOTAL_SUM += rho[k][i];
	    	}	
    	}
    	TOTAL_SUM /= (n1*n2);
        previous_energy=1;
        previous_dual=1;
        double error=1, dual_gap=1, energy=1, dual=0, sanity_value=1;
        int iterPDHG;

        double beta_1 = 1.5;
        double beta_2 = 0.9;

        double sanity_value_previous = 100;

        for(iterPDHG=0; iterPDHG<max_iteration; ++iterPDHG){


            // update phi
            for(int k=0;k<4;++k){
                memcpy(phitmps[k],phi[k],n1*n2*nt*sizeof(double));
            }
            sanity_value  = update_phi_all(rho,mx,my,obstacle);
            sanity_value_previous = sanity_value;

            for(int k=0;k<4;++k){
                for(int i=0;i<n1*n2*nt;++i){
                    phi[k][i] = 2*phi[k][i] - phitmps[k][i];
                }
            }

            // get the data before updates
            for(int k=0;k<4;++k){
                memcpy(rhotmps[k],rho[k],n1*n2*nt*sizeof(double));
            }

            int skip2 = 100;

            update_rho0(rho[0],rhotmps[1],rhotmps[2],rhotmps[3],mx[0],my[0],obstacle);
            fftpsDST->solve_heat_equation_with_bdry(rho[0], rhont0tmps[0], 1e-5);
            update_m(mx[0],my[0],rho[0],phi[0],0);    
            update_rho1(rhotmps[0],rho[1],rhotmps[2],mx[1],my[1],obstacle);
            fftpsDST->solve_heat_equation_with_bdry(rho[1], rhont0tmps[1], 1e-5);
            update_m(mx[1],my[1],rho[1],phi[1],1);
            update_rho2(rhotmps[0],rhotmps[1],rho[2],mx[2],my[2],obstacle);
            fftpsDST->solve_heat_equation_with_bdry(rho[2], rhont0tmps[2], 1e-5);
            update_m(mx[2],my[2],rho[2],phi[2],2);
            update_rho3(rhotmps[0],rhotmps[1],rhotmps[2],rho[3],mx[3],my[3],obstacle);
            fftpsDST->solve_heat_equation_with_bdry(rho[3], rhont0tmps[3], 1e-5);
            update_m(mx[3],my[3],rho[3],phi[3],3);        

            // CALCULATE ENERGY
            error=fabs((energy-previous_energy)/previous_energy);
            previous_energy=energy;

            dual  =calculate_dual(rho,obstacle);
            double relative_dual_error=fabs((dual-previous_dual)/previous_dual);
            previous_dual=dual;

            dual_gap = energy-dual;

            if((iterPDHG+1)%skip==0){

                display_log(iterPDHG,  tau[0],  sigma[0],  energy,  dual,  error,  dual_gap, sanity_value);
                create_bin_file(rho[0], n1*n2*nt, "./data/rho0.csv");
                create_bin_file(rho[1], n1*n2*nt, "./data/rho1.csv");
                create_bin_file(rho[2], n1*n2*nt, "./data/rho2.csv");
                create_bin_file(rho[3], n1*n2*nt, "./data/rho3.csv");
            }
            if((iterPDHG>20 ) && (fabs(dual_gap)<tolerance)) break;
        }

        cout<<"The method is done!!"<<endl;
        display_log(iterPDHG,  tau[0],  sigma[0],  energy,  dual,  error,  dual_gap, sanity_value);

        cout<<"Iterations     : "<<iterPDHG<<endl;
        cout<<"Energy         : "<<energy<<endl;
        cout<<"Relative Error : "<<error<<endl;
    }




}; // Method class


#endif