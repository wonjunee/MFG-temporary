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

    poisson_solver** fftps;
    poisson_solver_2d* fftps2d;

    int convN;
    double conv_sum;
    double conv_r;

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
        for(int i=0;i<3;++i){
            phitmps[i]=nullptr;
            rhotmps[i]=nullptr;
            mx[i]=nullptr;
            my[i]=nullptr;
            phi[i]=nullptr;
        }
        
        etalist=nullptr;
        alphalist=nullptr;
        fftps=nullptr;
        fftps2d=nullptr;
    }

    Method(int n1, int n2, int nt, double dx, double dy, double dt, 
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
        this->sigma[1] = sigma*0.6;
        this->sigma[2] = sigma*0.3;

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
        karr[0] = 5e-2;
        karr[1] = 5e-2;
        karr[2] = 5e-2;

        xi    = new double[nt];
        xitmp = new double[nt];

        phiT     = new double[n1*n2];
        phiTtmp  = new double[n1*n2];

        M0 = 0.5;
        // convN = n1/4; 
        // convN = convN + (convN % 2) -1; // to make it odd
        conv_r = 0.2;
        convN  = conv_r*n1;
        convN  = 2 * convN + 1;
        // convN  = 7;
        cout << "convN: " <<convN << endl;

        var_gamma = 0.1; // DON'T TOUCH THIS. This one is for regularization.

        for(int i=0;i<3;++i){
            mx[i]      = new double[n1*n2*nt];
            my[i]      = new double[n1*n2*nt];
            phi[i]     = new double[n1*n2*nt];
            phitmps[i] = new double[n1*n2*nt];
            rhotmps[i] = new double[n1*n2*nt];
        }    

        for(int k=0;k<3;++k){
        	for(int i=0;i<n1*n2*nt;++i){
        		mx[k][i]      = 0;
	            my[k][i]      = 0;
	            phi[k][i]     = 0;
	            phitmps[k][i] = 0;
        	}
        }

        clock_t t;
        t = clock();
            
        double eta = etalist[0];
        fftps    = new poisson_solver*[3];
        fftps[0] = new poisson_solver(n1,n2,nt,dx,dy,dt,eta);
        fftps[1] = new poisson_solver(n1,n2,nt,dx,dy,dt,eta);
        fftps[2] = new poisson_solver(n1,n2,nt,dx,dy,dt,eta);
        fftps2d = new poisson_solver_2d(n1,n2,dx,dy,eta);

        t = clock() - t;
        printf ("\nCPU time for setting up FFT: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
    }


    ~Method(){
        for(int k=0;k<3;++k){
            delete[] mx[k];
            delete[] my[k];
            delete[] rhotmps[k];
            delete[] phitmps[k];
            delete[] phi[k];
            delete fftps[k];
        }
        delete[] alphalist;
        delete[] xi;
        delete[] xitmp;
        delete[] phiT;
        delete[] phiTtmp;
        delete[] karr;

        delete[] tau;
        delete[] sigma;
        delete[] fftps;
        delete fftps2d;
        
    }

    /* 
        Convolution functions
    */

    double calculate_K_xy(double x, double xx, double y, double yy, double var) const{
        return 1.0/(var*sqrt(2*M_PI))*exp(- ((x-xx)*(x-xx)+(y-yy)*(y-yy))/(2*(var*var)));
    }


    double calculate_convval(const double* rho, const int i, const int j, const double var) const{
        double convval = 0;
        double conv_sum= 0;

        double xx=(j+0.5)/n1;
        double yy=(i+0.5)/n2;

        for(int i1=i-convN/2;i1<i+convN/2+1;++i1){
            for(int j1=j-convN/2;j1<j+convN/2+1;++j1){
        // for(int i1=fmax(0,i-convN/2);i1<fmin(n2,i+convN/2+1);++i1){
        //     for(int j1=fmax(0,j-convN/2);j1<fmin(n1,j+convN/2+1);++j1){
                double x=(j1+0.5)/n1;
                double y=(i1+0.5)/n2;        
                if((x-xx)*(x-xx) + (y-yy)*(y-yy) <= conv_r*conv_r){
	                double eval = calculate_K_xy(x,xx,y,yy,var);

	                int ii = fmin(n2-1, fmax(0, i1));
	                int jj = fmin(n1-1, fmax(0, j1));

	                convval += eval * rho[ii*n1+jj];
                    conv_sum+= eval;
	            }
            }
        }

        convval/=conv_sum;
        return convval;
    }

    double calculate_convval2(const double* rho, const double* phi, const int i, const int j, const double var) const{
        double convval = 0;
        double conv_sum= 0;

        double xx=(j+0.5)/n1;
        double yy=(i+0.5)/n2;

        for(int i1=i-convN/2;i1<i+convN/2+1;++i1){
            for(int j1=j-convN/2;j1<j+convN/2+1;++j1){
                
                double x=(j1+0.5)/n1;
                double y=(i1+0.5)/n2;        
                if((x-xx)*(x-xx) + (y-yy)*(y-yy) <= conv_r*conv_r){
	                double eval = calculate_K_xy(x,xx,y,yy,var);

	                int ii = fmin(n2-1, fmax(0, i1));
	                int jj = fmin(n1-1, fmax(0, j1));

	                convval += eval * rho[ii*n1+jj] * phi[ii*n1+jj];
                    conv_sum+= eval;
                }
                
            }
        }
        convval/=conv_sum;
        return convval;
    }

    double calculate_convval3(const double* rho, const double* phi0, const double* phi1, const int i, const int j, const double var) const{
        double convval = 0;
        double conv_sum= 0;

        double xx=(j+0.5)/n1;
        double yy=(i+0.5)/n2;

        for(int i1=i-convN/2;i1<i+convN/2+1;++i1){
            for(int j1=j-convN/2;j1<j+convN/2+1;++j1){
                
                double x=(j1+0.5)/n1;
                double y=(i1+0.5)/n2;        
                if((x-xx)*(x-xx) + (y-yy)*(y-yy) <= conv_r*conv_r){
                    double eval = calculate_K_xy(x,xx,y,yy,var);

                    int ii = fmin(n2-1, fmax(0, i1));
                    int jj = fmin(n1-1, fmax(0, j1));

                    convval += eval * rho[ii*n1+jj] * (phi1[ii*n1+jj] - phi0[ii*n1+jj]);
                    conv_sum+= eval;
                }
                
            }
        }
        convval/=conv_sum;
        return convval;
    }

    double calculate_convval4(const double* phi0, const double* phi1, const int i, const int j, const double var) const{
        double convval = 0;
        double conv_sum = 0;

        double xx=(j+0.5)/n1;
        double yy=(i+0.5)/n2;

        for(int i1=i-convN/2;i1<i+convN/2+1;++i1){
            for(int j1=j-convN/2;j1<j+convN/2+1;++j1){
                
                double x=(j1+0.5)/n1;
                double y=(i1+0.5)/n2;        
                if((x-xx)*(x-xx) + (y-yy)*(y-yy) <= conv_r*conv_r){
                    double eval = calculate_K_xy(x,xx,y,yy,var);

                    int ii = fmin(n2-1, fmax(0, i1));
                    int jj = fmin(n1-1, fmax(0, j1));

                    convval += eval * (phi1[ii*n1+jj] - phi0[ii*n1+jj]);
                    conv_sum+= eval;
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
                    
                    mx[n*n1*n2+i*n1+j] = h * rhovalx/(tau[ind] * alphalist[ind] + h * rhovalx) * (mx[n*n1*n2+i*n1+j] - tau[ind] * nablaxphi);
                    my[n*n1*n2+i*n1+j] = h * rhovaly/(tau[ind] * alphalist[ind] + h * rhovaly) * (my[n*n1*n2+i*n1+j] - tau[ind] * nablayphi);

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
        // if(n==nt-1){
        //     Dtphi=0.5/dt*(phi[n*n1*n2+i*n1+j]-phi[(n-1)*n1*n2+i*n1+j]);
        // }
        
        Deltaphi = - n1*n1 * (-phi[n*n1*n2+i*n1+(int) fmax(0,j-1)]+2*phi[n*n1*n2+i*n1+j]-phi[n*n1*n2+i*n1+(int) fmin(n1-1,j+1)])
                   - n2*n2 * (-phi[n*n1*n2+(int) fmax(0,i-1)*n1+j]+2*phi[n*n1*n2+i*n1+j]-phi[n*n1*n2 +(int) fmin(n2-1,i+1)*n1+j]);

    }
    double calculate_E0prime(const double* rho, const double* f, const int i, const int j) const{
        int ind = (nt-1)*n1*n2+i*n1+j;
        // rho log rho
        // double eval = 1.0/c0 * (phi[0][ind] - f[i*n1+j]);
        double val = 0;
        // if(f[i*n1+j] >= 0) val = exp(eval);

        // double eval = 1.0/c0 * (phi[0][ind] - f[i*n1+j]) + rho[i*n1+j];

        // if(eval <= 0) eval = 0;

        // double val = 0;
        // if(f[i*n1+j] >= 0 ){
        //     val = eval * phi[0][ind] - 0.5 * c0 * (eval - rho[i*n1+j]) * (eval - rho[i*n1+j]);
        // }

        val = phi[0][ind]*rho[ind];

        return val;
    }
    double calculate_E1prime(const double* rho, const double* f, int i, int j) const{
        int ind = (nt-1)*n1*n2+i*n1+j;
        double eval = 1.0/c1 * (phi[1][ind] - f[i*n1+j]);
        double val = 0;
        // if(eval > 0 && f[i*n1+j] >= 0) val = 1.0/(2*c1_value) * pow(fmax(0, (phi[1][ind]-f[i])), 2);

        // if(eval > 0 && f[i*n1+j] >= 0) val = 0.5 * eval * (phi[1][ind]-f[i]);

        val = phi[1][ind]*rho[ind];

        return val;
    }
    double calculate_E2prime(const double* rho, const double* f, int i, int j) const{
        int ind = (nt-1)*n1*n2+i*n1+j;
        // rho log rho
        double eval = 1.0/c2 * (phi[2][ind] - f[i*n1+j]);
        double val  = 0;

        // if(f[i*n1+j] >= 0) val = c2*exp(eval);

        val = phi[2][ind]*rho[ind];

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


        // ----- rho(1,x) = dGstar(phi(1,x)) ------

        // for(int i=0;i<n2;++i){
        //     for(int j=0;j<n1;++j){
        //         int ind = (nt-1)*n1*n2+i*n1+j;
        //         double val = calculate_deltaE0prime(rho0, f, i, j);

        //         // rho0[ind] = 1.0/(1.0+tau*0.5)*(rho0[ind] + tau*0.5*val);  
        //         rho0[ind] += 0.0001 * (val - rho0[ind]);
        //         rho0[ind] = fmax(rho0[ind],0);
        //         // rho0[ind] = val;
        //     }
        // }


    	double newrhovalue = -1;
		double tauval = tau[0];
		double betaval  = beta;

		for(int n=1;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int ind = n*n1*n2+i*n1+j;

                    double mvalue=0;
                    double Dtphi =0;
                    double Deltaphi=0;

                    calculate_rho_related(mvalue, Dtphi, Deltaphi, n, i, j, mx, my, phi[0]);

                    double convval  =  calculate_convval(&rho1[n*n1*n2],i,j,var);
                    // double convval2 = calculate_convval2(&rho1[n*n1*n2],&phi[1][n*n1*n2],i,j);
                    // double convval3 = calculate_convval3(&rho1[n*n1*n2],&phi[0][n*n1*n2],&phi[1][n*n1*n2],i,j);

                    // double newrhovalue=cubic_solve(tau*Dtphi + tau*etalist[0]*Deltaphi - rho0[ind] + tau*betaval*rho1[ind]*(phi[1][ind]-phi[0][ind]), 0, -0.5/h*tau*alphalist[0]*mvalue*mvalue);
                        
                    double aval = 1.0/(tauval*karr[0]+1) * 
                                (
                                    tauval * (Dtphi + etalist[0]*Deltaphi + xi[n] + karr[0]*rho1[ind] + karr[0]*rho2[ind]) 
                                    - rho0[ind] 
                                    + tauval*betaval*convval*(phi[1][ind] - phi[0][ind])
                                );
                    double cval = -0.5/h*tauval/(tauval*karr[0]+1)*alphalist[0]*mvalue*mvalue;
                	newrhovalue=cubic_solve(aval, 0, cval);

                    rho0[n*n1*n2+i*n1+j] = fmin(1,fmax(0,newrhovalue));
                    
                }
            }
		}



    }

    void update_rho1(const double* rho0,double* rho1, const double* rho2,const double* mx,const double* my,const double* f){

        // ----- rho(1,x) = dGstar(phi(1,x)) ------

    	double tauval = tau[1];
    	
  //   	for(int i=0;i<n2;++i){
  //           for(int j=0;j<n1;++j){
  //               fftps2d->u[i*n1+j] = phiT[i*n1+j];
  //           }
  //       }
  //       fftps2d->perform_inverse_laplacian(0);

		// for(int i=0;i<n1*n2;++i){
	 //    	double newrhovalue = rho1[(nt-1)*n1*n2+i] - tauval * fftps2d->workspace[i];
	 //        rho1[(nt-1)*n1*n2+i] =  fmax(0, newrhovalue);
	 //    }

        for(int i=0;i<n1*n2;++i){
            double newrhovalue = rho1[(nt-1)*n1*n2+i] - tauval * phiT[i];
            rho1[(nt-1)*n1*n2+i] =  fmax(0, newrhovalue);
        }

        for(int n=1;n<nt;++n){

        	double gammaval = gamma;
        	double betaval  = beta;

            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int ind = n*n1*n2+i*n1+j;

                    
                    double mvalue=0;
                    double Dtphi =0;
                    double Deltaphi=0;

                    calculate_rho_related(mvalue, Dtphi, Deltaphi, n, i, j, mx, my, phi[1]);

                    // double convval  = calculate_convval(&rho0[n*n1*n2],i,j);
                    double convval3 = calculate_convval3(&rho0[n*n1*n2],&phi[0][n*n1*n2],&phi[1][n*n1*n2],i,j,var);
                    double convval4 = calculate_convval4(&phi[1][n*n1*n2],&phi[2][n*n1*n2],i,j,var_gamma);

                    double newrhovalue = 0;
                    // double newrhovalue=cubic_solve(tau[1]*Dtphi + tau[1]*etalist[1]*Deltaphi - rho1[ind] + tau[1]*beta*rho0[ind]*(phi[1][ind]-phi[0][ind]) + tau[1]*gamma*(phi[2][ind] - phi[1][ind]), 0, -0.5/h*tau[1]*alphalist[1]*mvalue*mvalue);

                    double aval=0,cval=0;

                    // if(n==nt){
                    // 	// newrhovalue=cubic_solve(tauval * (Dtphi + etalist[1]*Deltaphi + xi[n] + phiT[i*n1+j]*nt + karr[0]*rho0[ind] + karr[0]*rho2[ind]) - rho1[ind] + tauval*betaval*(phi[1][ind] - phi[0][ind])*convval + tauval*gammaval*(phi[2][ind] - phi[1][ind]), 0, -0.5/h*tauval*alphalist[1]*mvalue*mvalue);
                    // }else{
                    //     aval = 1.0/(tauval*karr[1]+1) * 
                    //             ( 
                    //                 tauval * (Dtphi + etalist[1]*Deltaphi + xi[n] + karr[0]*rho0[ind] + karr[0]*rho2[ind]) 
                    //                 - rho1[ind] 
                    //                 + tauval*betaval*convval3
                    //                 + tauval*gammaval * (phi[2][ind] - phi[1][ind]) 
                    //             );
                    //     cval = -0.5/h*tauval/(tauval*karr[1]+1)*alphalist[1]*mvalue*mvalue;
                    // }

                        // aval = 1.0/(tauval*karr[1]+1) * 
                        //         ( 
                        //             tauval * (Dtphi + etalist[1]*Deltaphi + xi[n] + karr[0]*rho0[ind] + karr[0]*rho2[ind]) 
                        //             - rho1[ind] 
                        //             + tauval*betaval*convval3
                        //             + tauval*gammaval * (phi[2][ind] - phi[1][ind]) 
                        //         );
                        // cval = -0.5/h*tauval/(tauval*karr[1]+1)*alphalist[1]*mvalue*mvalue;

                        aval = 1.0/(tauval*karr[1]+1) * 
                                ( 
                                    tauval * (Dtphi + etalist[1]*Deltaphi + xi[n] + karr[0]*rho0[ind] + karr[0]*rho2[ind]) 
                                    - rho1[ind] 
                                    + tauval * betaval  * convval3
                                    + tauval * gammaval * convval4
                                );
                        cval = -0.5/h*tauval/(tauval*karr[1]+1)*alphalist[1]*mvalue*mvalue;


                    newrhovalue=cubic_solve(aval, 0, cval);
                    
                    rho1[ind]=fmin(1,fmax(0,newrhovalue));
                }
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

                    double newrhovalue=cubic_solve(1.0/(tau[2]*karr[2]+1) * (tau[2] * (Dtphi + etalist[2]*Deltaphi + xi[n] + karr[0]*rho0[ind] + karr[0]*rho1[ind]) - rho2[ind]), 0, -0.5/h*tau[2]/(tau[2]*karr[2]+1)*alphalist[2]*mvalue*mvalue);
                    rho2[ind]=fmin(1,fmax(0,newrhovalue));
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
        //     dtrho=0;
        // }else{
        //     dtrho=1.0/dt*(rho[(n)*n1*n2+i*n1+j]-rho[(n-1)*n1*n2+i*n1+j]); 
        // }

        // if(n==0){
        //     dtrho=0.5/dt*(rho[(n+1)*n1*n2+i*n1+j]-rho[(n)*n1*n2+i*n1+j]); 
        // }else if(n==nt-1){
        // 	dtrho=0.5/dt*(rho[(n)*n1*n2+i*n1+j]-rho[(n-1)*n1*n2+i*n1+j]); 
        // }else{
        //     dtrho=0.5/dt*(rho[(n+1)*n1*n2+i*n1+j]-rho[(n-1)*n1*n2+i*n1+j]); 
        // }
        return dtrho;
    }

    double calculate_Delta_value(const double* rho, const int n, const int i, const int j){
        return - n1*n1 * (-rho[n*n1*n2+i*n1+(int) fmax(0,j-1)]+2*rho[n*n1*n2+i*n1+j]-rho[n*n1*n2+i*n1+(int) fmin(n1-1,j+1)])
               - n2*n2 * (-rho[n*n1*n2+(int) fmax(0,i-1)*n1+j]+2*rho[n*n1*n2+i*n1+j]-rho[n*n1*n2 +(int) fmin(n2-1,i+1)*n1+j]);
    }

    double update_phi_all(double* const rho[], double* const mx[], double* const my[], double* const f[]){

        for(int n=0;n<nt;++n){  
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int ind = n*n1*n2+i*n1+j;

                    double dtrho0 = calculate_dtrho(rho[0], n, i, j);
                    double nablamx0=calculate_grad_mx(mx[0],n,i,j);
                    double nablamy0=calculate_grad_my(my[0],n,i,j);

                    double dtrho1 = calculate_dtrho(rho[1], n, i, j);
                    double nablamx1=calculate_grad_mx(mx[1],n,i,j);
                    double nablamy1=calculate_grad_my(my[1],n,i,j);

                    double dtrho2 = calculate_dtrho(rho[2], n, i, j);
                    double nablamx2=calculate_grad_mx(mx[2],n,i,j);
                    double nablamy2=calculate_grad_my(my[2],n,i,j);

                    double Deltarho0 = calculate_Delta_value(rho[0],n,i,j);
                    double Deltarho1 = calculate_Delta_value(rho[1],n,i,j);
                    double Deltarho2 = calculate_Delta_value(rho[2],n,i,j);

                    double convval = calculate_convval(&rho[1][n*n1*n2],i,j,var);
                    double convval_gamma = calculate_convval(&rho[1][n*n1*n2],i,j,var_gamma);

                    // fftps->u[n*n1*n2+i*n1+j]=-(dtrho+nablamx+nablamy + beta*rho[0][ind]*rho[1][ind] - etalist[0]*Deltarho); 
                    fftps[0]->u[ind]  = - (dtrho0 + nablamx0 + nablamy0 + beta*rho[0][ind]*convval - etalist[0]*Deltarho0); 

                    fftps[1]->u[ind]  = - (dtrho1 + nablamx1 + nablamy1 - beta*convval*rho[0][ind] + gamma*convval_gamma - etalist[1]*Deltarho1);
                    // fftps[1]->u[ind]  = - (dtrho1 + nablamx1 + nablamy1 - beta*convval*rho[0][ind] + gamma*rho[1][ind] - etalist[1]*Deltarho1);

                    fftps[2]->u[ind]  = - (dtrho2 + nablamx2 + nablamy2 - gamma*convval_gamma - etalist[2]*Deltarho2); 
                    // fftps[2]->u[ind]  = - (dtrho2 + nablamx2 + nablamy2 - gamma*rho[1][ind] - etalist[2]*Deltarho2); 
                }
            }
        }

        double error = 0;

        for(int k=0;k<3;++k){
            if(k==0){
                fftps[k]->perform_inverse_laplacian(0, etalist[0]);    
            }else if(k==1){
                fftps[k]->perform_inverse_laplacian(0, etalist[1]);    
            }else{
                fftps[k]->perform_inverse_laplacian(0, etalist[2]);    
            }
        
            for(int i=0;i<n1*n2*nt;++i){
                phi[k][i] += sigma[k]      * fftps[k]->workspace[i];
                error     += fftps[k]->u[i]* fftps[k]->workspace[i];
            }
        }
        return error/(1.0*n1*n2*nt);

    }

    void update_xi(double* const rho[]){
    	for(int n=1;n<nt;++n){
    		double rho_sum = 0;
    		for(int k=0;k<3;++k){
	    		for(int i=0;i<n1*n2;++i){
	    			rho_sum += rho[k][n*n1*n2+i];
	    		}
	    	}
	    	xi[n] += sigma[0] * (rho_sum/(n1*n2) - TOTAL_SUM);
    	}
    }

    void update_phiT(double * const rho[], double* const f[]){

    	double sigmaval = sigma[1];

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                double rhoval = rho[1][(nt-1)*n1*n2+i*n1+j];
                double DeltaphiT = calculate_Delta_value(phiT,0,i,j);

                fftps2d->u[i*n1+j] = (sigmaval * rhoval - DeltaphiT); 
            }
        }

        fftps2d->perform_inverse_laplacian_phiT(sigmaval/c1,phiT);
        for(int i=0;i<n1*n2;++i){
            phiT[i] = fftps2d->workspace[i];
        }

        // for(int i=0;i<n1*n2;++i){
        // 	double val = 1;
        // 	if(phiT[i]>0){
        // 		val += sigmaval/c1;
        // 	}
        // 	phiT[i] = 1.0/val * (sigmaval * rho[1][(nt-1)*n1*n2+i] + phiT[i]);
        // }

        // for(int i=0;i<n1*n2;++i){
        //     if(phiT[i] - f[1][i]>0){
        //         double val = sigmaval/c1 + 1;
        //         phiT[i] = 1.0/val * (sigmaval * rho[1][(nt-1)*n1*n2+i] + phiT[i] + sigmaval/c1 * f[1][i]);
        //     }else{
        //         phiT[i] =sigmaval * rho[1][(nt-1)*n1*n2+i] + phiT[i];
        //     }
            
        // }
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


    double calculate_energy(double* const rho[], double* const* f) const{
        double sum = 0;
        for(int k=0;k<3;++k){
            sum += calculate_energy_num(rho,mx,my,k);
        }

        double sum1 = 0;
        for(int i=0;i<n1*n2;++i){
            int ind = (nt-1)*n1*n2+i;
            // if(rho[0][ind] > 0) sum1 += c0*(rho[0][ind]*log(rho[0][ind]) - rho[0][ind]);
            // if(rho[0][ind] >= 0) sum1 += c0/2.0*(rho[0][ind] - rho[0][i])*(rho[0][ind] - rho[0][i]) + rho[0][ind]*f[0][i];
            if(rho[1][ind] >= 0) sum1 += c1/m*rho[1][ind]*rho[1][ind] 							    + rho[1][ind]*f[1][i];
            // if(rho[2][ind] >  0) sum1 += c2*(rho[2][ind]*log(rho[2][ind]) - rho[2][ind]) 		    + rho[2][ind]*f[2][i];
        }
        sum1 /= 1.0*n1*n2;

        double sum2 = 0;
        for(int n=0;n<nt;++n){
            for(int i=0;i<n1*n2;++i){
                double val = rho[0][n*n1*n2+i] + rho[1][n*n1*n2+i] + rho[2][n*n1*n2+i];
                sum2 += 0.5 * karr[0] * (val*val);
            }
        }
        sum2 /= 1.0*n1*n2*nt;
        
        return sum + sum1;
    }

    double calculate_dual(double* const rho[], double* const* f) const{   

        double term0=0;
        double term1=0;
        double term2=0;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                // if(c0>0) term0 += phi[0][i]*rho[0][i] - c0*exp((phi[0][(nt-1)*n1*n2+i]-f[i])/c0);
                // term1 += phi[1][i]*rho[1][i] - 1.0/(mprime*c1_value) * pow(fmax(0, (phi[1][(nt-1)*n1*n2+i]-f[i])), mprime);
                // if(c2>0) term2 += phi[2][i]*rho[2][i] - c2*exp((phi[2][(nt-1)*n1*n2+i]-f[i])/c2);    

                if(c0>0) term0 += phi[0][i*n1+j]*rho[0][i*n1+j] - calculate_E0prime(rho[0],f[0],i,j);
                if(c1>0) term1 += phi[1][i*n1+j]*rho[1][i*n1+j] - calculate_E1prime(rho[1],f[1],i,j);
                if(c2>0) term2 += phi[2][i*n1+j]*rho[2][i*n1+j] - calculate_E2prime(rho[2],f[2],i,j);
            }
            
        }

        double term3=0;
        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    // double convval0 = calculate_convval(&rho[0][n*n1*n2],i,j);
                    // double convval1 = calculate_convval(&rho[1][n*n1*n2],i,j);
                    // term3 += beta * (rho[0][n*n1*n2+i*n1+j] * convval1 * phi[0][n*n1*n2+i*n1+j] - rho[1][n*n1*n2+i*n1+j] * convval0 * phi[1][n*n1*n2+i*n1+j]);
                    term3 += beta * rho[0][n*n1*n2+i*n1+j] * rho[1][n*n1*n2+i*n1+j] * (phi[0][n*n1*n2+i*n1+j] - phi[1][n*n1*n2+i*n1+j]);
                }
            }
        }

        return (term0+term1+term2)/(n1*n2) + term3/(n1*n2*nt);
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
        cout<<"iter : "<< setw(8) << iterPDHG+1<< " tau : "<< setw(13) <<tau <<" sigma : " << setw(13)<<sigma <<" energy : " << setw(13)<<energy <<" dual : " << setw(13) << dual << " relative_error : " << setw(13) << error <<" dual_gap : " << setw(13) << dual_gap<<" sanity_value : " << setw(13) << sanity_value << endl;  
    }

    void run(double* rho[], double* const* f, int skip=1){

    	TOTAL_SUM = 0;
    	for(int k=0;k<3;++k){
    		for(int i=0;i<n1*n2;++i){
	    		TOTAL_SUM += rho[k][i];
	    	}	
    	}
    	TOTAL_SUM /= (n1*n2);
    	for(int n=0;n<nt;++n){
    		xi[n] = 0;
    	}

        previous_energy=1;
        previous_dual=1;
        double error=1, dual_gap=1, energy=1, dual=0, sanity_value=1;
        int iterPDHG;

        double beta_1 = 1.5;
        double beta_2 = 0.9;

        double sanity_value_previous = 100;

        for(iterPDHG=0; iterPDHG<max_iteration; ++iterPDHG){

            // get the data before updates


            update_m(mx[0],my[0],rho[0],phi[0],0);
            update_m(mx[1],my[1],rho[1],phi[1],1);
            update_m(mx[2],my[2],rho[2],phi[2],2);

            for(int k=0;k<3;++k){
                memcpy(rhotmps[k],rho[k],n1*n2*nt*sizeof(double));
            }

            update_rho0(rho[0],rhotmps[1],rhotmps[2],mx[0],my[0],f[0]);
            update_rho1(rhotmps[0],rho[1],rhotmps[2],mx[1],my[1],f[1]);

            if(iterPDHG > 0){
                update_rho2(rhotmps[0],rhotmps[1],rho[2],mx[2],my[2],f[2]);    
            }

            /*
                    UPDATE PHI
            */

            for(int k=0;k<3;++k){
                memcpy(phitmps[k],phi[k],n1*n2*nt*sizeof(double));
            }
            // double error0 = update_phi0(rho,mx[0],my[0],f[0]); 
            // double error1 = update_phi1(rho,mx[1],my[1],f[1]);
            // double error2 = update_phi2(rho,mx[2],my[2],f[2]); 

            sanity_value  = update_phi_all(rho,mx,my,f);

            if(iterPDHG > 0)
            {
                if(sanity_value > sanity_value_previous){           
                    for(int k=0;k<3;++k){
                        tau[k]   *= 0.99;
                        sigma[k] *= 0.99;
                    }
                }

                if(sanity_value < 0.95 * sanity_value_previous){           
                    for(int k=0;k<2;++k){
                        tau[k]   *= 1.01;
                        sigma[k] *= 1.01;
                    }
                }
            }
            
            sanity_value_previous = sanity_value;

            for(int k=0;k<3;++k){
                for(int i=0;i<n1*n2*nt;++i){
                    phi[k][i] = 2*phi[k][i] - phitmps[k][i];
                }
            }
             
            
            /*
                UPDATE XI
            */
            memcpy(xitmp, xi, nt*sizeof(double));
            update_xi(rho);
            for(int n=0;n<nt;++n){
            	xi[n] = 2*xi[n] - xitmp[n];
            }


            /*
                UPDATE PSI
            */
            memcpy(phiTtmp, phiT, n1*n2*sizeof(double));
            update_phiT(rho,f);
            for(int i=0;i<n1*n2;++i){
            	phiT[i] = 2*phiT[i] - phiTtmp[i];
            }

            /*
                CALCULATE ENERGY
            */
            energy=calculate_energy(rho,f);
            error=fabs((energy-previous_energy)/previous_energy);
            previous_energy=energy;

            dual  =calculate_dual(rho,f);
            double relative_dual_error=fabs((dual-previous_dual)/previous_dual);
            previous_dual=dual;

            // update_step_sizes(tau[0], sigma[0], error, relative_dual_error, beta_1, beta_2);
            // update_step_sizes(tau[1], sigma[1], error, relative_dual_error, beta_1, beta_2);
            // update_step_sizes(tau[2], sigma[2], error, relative_dual_error, beta_1, beta_2);

            // for(int k=0;k<3;++k){
            // 	tau[k]   = fmin(1.0,fmax(1e-3,tau[k]));
            // 	sigma[k] = fmin(1.0,fmax(1e-3,sigma[k]));
            // }

            dual_gap = energy-dual;

            /*
                CALCULATE ERROR
            */

            // sanity_value = error;



            if((iterPDHG+1)%skip==0){

                display_log(iterPDHG,  tau[0],  sigma[0],  energy,  dual,  error,  dual_gap, sanity_value);

                create_csv_file(rho[0],"./data/rho0.csv",n1,n2,nt);
                create_csv_file(rho[1],"./data/rho1.csv",n1,n2,nt);
                create_csv_file(rho[2],"./data/rho2.csv",n1,n2,nt);
            }
            if((iterPDHG>20 ) && (fabs(dual_gap)<tolerance)) break;

            // if(fabs(dual_gap)<1e-2 && sanity_value<0.5) break;
        }

        cout<<"The method is done!!"<<endl;
        display_log(iterPDHG,  tau[0],  sigma[0],  energy,  dual,  error,  dual_gap, sanity_value);

        cout<<"Iterations     : "<<iterPDHG<<endl;
        cout<<"Energy         : "<<energy<<endl;
        cout<<"Relative Error : "<<error<<endl;
    }

}; // Method class


#endif