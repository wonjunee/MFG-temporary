#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H

#include <iostream>
#include <fftw3.h>

using namespace std;

class poisson_solver{
public:
    fftw_plan planIn;
    fftw_plan planOut;
    double *workspace;
    double *u;
    double *kernel;
    double *tridiagonalWorkspace;

    int n1;
    int n2;
    int nt;
    double dx;
    double dy;
    double dt;

    poisson_solver(int n1, int n2, int nt, double dx, double dy, double dt) {
    	this->n1=n1;
    	this->n2=n2;
    	this->nt=nt;
    	this->dx=dx;
    	this->dy=dy;
    	this->dt=dt;

        workspace =(double*) fftw_malloc(n1*n2*sizeof(double));

		planIn = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT10,FFTW_REDFT10, FFTW_MEASURE);
		planOut = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT01,FFTW_REDFT01, FFTW_MEASURE);

		u=new double[n1*n2*nt];
        kernel=new double[n1*n2];
        tridiagonalWorkspace=new double[n1*n2*nt];

        create_negative_laplacian_kernel_2d();
    }

    void create_negative_laplacian_kernel_2d(){
	    
	    for(int i=0;i<n2;i++){
	        for(int j=0;j<n1;j++){
	            double negativeLaplacian=2/(dx*dx)*(1-cos(M_PI*(j)*dx)) + 2/(dy*dy)*(1-cos(M_PI*(i)*dy));
	            kernel[i*n1+j]=1+negativeLaplacian;
	        }
	    }
	}

	void forward_tridiagonal_sweep(){

	    for(int i=0;i<n1*n2;i++){
	        tridiagonalWorkspace[i]=0.0/1.0; // c
	        u[i]= u[i]/1.0;	// d
	    }
	     
	    for(int k=1;k<nt-1;k++){
	        for(int i=0;i<n1*n2;i++){
	            double alpha=kernel[i]*dt*dt;
	            tridiagonalWorkspace[k*n1*n2+i]=-1/(2+alpha+tridiagonalWorkspace[(k-1)*n1*n2+i]);
	            u[k*n1*n2+i]=(u[k*n1*n2+i]+u[(k-1)*n1*n2+i])/(2+alpha+tridiagonalWorkspace[(k-1)*n1*n2+i]);
	        }
	    }
	    
	    for(int i=0;i<n1*n2;i++){
	        u[(nt-1)*n1*n2+i]=u[(nt-1)*n1*n2+i]/1.0;
	    }


	    // for(int i=0;i<n1*n2;i++){
	    //     double alpha=kernel[i]/(nt*nt);
	        
	    //     tridiagonalWorkspace[i]=-1/(1+alpha);
	    //     u[i]/=1+alpha;
	        
	    // }
	    
	    
	    // for(int k=1;k<nt-1;k++){
	    //     for(int i=0;i<n1*n2;i++){
	    //         double alpha=kernel[i]/(nt*nt);
	            
	    //         tridiagonalWorkspace[k*n1*n2+i]=-1/(2+alpha+tridiagonalWorkspace[(k-1)*n1*n2+i]);
	    //         u[k*n1*n2+i]=(u[k*n1*n2+i]+u[(k-1)*n1*n2+i])/(2+alpha+tridiagonalWorkspace[(k-1)*n1*n2+i]);
	            
	    //     }
	    // }
	    
	    // u[(nt-1)*n1*n2]=0;
	    
	    // for(int i=1;i<n1*n2;i++){
	    //     double alpha=kernel[i]/(nt*nt);
	    //     u[(nt-1)*n1*n2+i]=(u[(nt-1)*n1*n2+i]+u[(nt-2)*n1*n2+i])/(1+alpha+tridiagonalWorkspace[(nt-2)*n1*n2+i]);
	        
	    // }
	    
	}

	void backward_tridiagonal_sweep(){
	    for(int k=nt-2;k>=0;k--){
	        for(int i=0;i<n1*n2;i++){
	            u[k*n1*n2+i]=u[k*n1*n2+i]-tridiagonalWorkspace[k*n1*n2+i]*u[(k+1)*n1*n2+i];
	        }
	        
	    }
	}

	void perform_inverse_laplacian(){

		fftw_execute(planIn);

		// workspace[0]=0;

		for(int i=0;i<n1*n2;++i){
			workspace[i]/=4*(n1)*(n2)*(kernel[i]);
		}

		fftw_execute(planOut);
	}

	void get_fourier_coefficients(double* u){

		for(int i=0;i<n1*n2;++i){
			workspace[i]/=1.0/(dt*dt);	
		}

		fftw_execute(planIn);

		for(int i=0;i<n1*n2;++i){
			u[i]=workspace[i]/(4.0*n1*n2);
		}
	}
	void back_to_real_space(){
		fftw_execute(planOut);
	}

	void destroy_all_fftps(){
		delete[] u;
	    delete[] workspace;
	    delete[] kernel;
	    fftw_destroy_plan(planIn);
	    fftw_destroy_plan(planOut);
	}
};

#endif