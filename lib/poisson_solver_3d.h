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
    unique_ptr<double[]> u;
    unique_ptr<double[]> kernel;

    int n1;
    int n2;
    int nt;
    double eta;

    poisson_solver(){
        workspace=NULL;
    }
    poisson_solver(const int n1, const int n2, const int nt, const double eta=1) {
        this->n1=n1;
        this->n2=n2;
        this->nt=nt;

        this->eta=eta;

        workspace =(double*) fftw_malloc(n1*n2*nt*sizeof(double));

        planIn =fftw_plan_r2r_3d(nt, n2, n1, workspace, workspace, FFTW_REDFT10, FFTW_REDFT10, FFTW_REDFT10, FFTW_MEASURE);
        planOut=fftw_plan_r2r_3d(nt, n2, n1, workspace, workspace, FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_MEASURE);

        u      = make_unique<double[]>(n1*n2*nt);
        kernel = make_unique<double[]>(n1*n2*nt);

        create_negative_laplacian_kernel_2d();
    }

    ~poisson_solver(){
        fftw_free(workspace);
        fftw_destroy_plan(planIn);
        fftw_destroy_plan(planOut);
    }

    void create_negative_laplacian_kernel_2d(){
        
        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;i++){
                for(int j=0;j<n1;j++){
                    double xpart = 2*n1*n2*(1-cos(M_PI*(1.0*j)/n1));
                    double ypart = 2*n2*n2*(1-cos(M_PI*(1.0*i)/n2));
                    double tpart = 2*nt*nt*(1-cos(M_PI*(1.0*n)/nt));  // DCT
                    // double tpart = 2*nt*nt*(1-cos(M_PI*(1.0*n+1)/nt));// DST
                    double negativeLaplacian = xpart+ypart+tpart;
                    kernel[n*n1*n2+i*n1+j]   = negativeLaplacian;
                }
            }
        }
            
    }

    void perform_inverse_laplacian(){

        for(int i=0;i<n1*n2*nt;++i) workspace[i] = u[i];

        fftw_execute(planIn);

        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    double val = (0.1 + kernel[n*n1*n2+i*n1+j]);
                    if(val==0){
                        workspace[n*n1*n2+i*n1+j]=0;    
                    }else{
                        workspace[n*n1*n2+i*n1+j]/=8*(n1)*(n2)*(nt)*val;
                    }
                }
            }
        }

        fftw_execute(planOut);
    }

        // function for general rho. Boundary for t=0 only
    void solve_heat_equation(double* rho, double tau){
        memcpy(workspace,rho,n1*n2*nt*sizeof(double));
        fftw_execute(planIn);
        for(int i=0;i<n1*n2*nt;++i) workspace[i] /= 8.0*n1*n2*nt * (1 + tau * kernel[i]);
        fftw_execute(planOut);
        memcpy(rho,workspace,n1*n2*nt*sizeof(double));
    }
};


class poisson_solver_DST{
public:
    fftw_plan planIn;
    fftw_plan planOut;
    double *workspace;
    unique_ptr<double[]> u;
    unique_ptr<double[]> kernel;

    int n1;
    int n2;
    int nt;
    double eta;

    poisson_solver_DST(){
        workspace=NULL;
    }
    poisson_solver_DST(const int n1, const int n2, const int nt, const double eta=1) {
        this->n1=n1;
        this->n2=n2;
        this->nt=nt;

        this->eta=eta;

        workspace =(double*) fftw_malloc(n1*n2*nt*sizeof(double));

        planIn =fftw_plan_r2r_3d(nt, n2, n1, workspace, workspace, FFTW_RODFT10, FFTW_REDFT10, FFTW_REDFT10, FFTW_MEASURE);
        planOut=fftw_plan_r2r_3d(nt, n2, n1, workspace, workspace, FFTW_RODFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_MEASURE);

        u      = make_unique<double[]>(n1*n2*nt);
        kernel = make_unique<double[]>(n1*n2*nt);


        create_negative_laplacian_kernel_3d();
    }

    ~poisson_solver_DST(){
        fftw_free(workspace);
        fftw_destroy_plan(planIn);
        fftw_destroy_plan(planOut);
    }

    void create_negative_laplacian_kernel_3d(){
        
        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;i++){
                for(int j=0;j<n1;j++){
                    double xpart = 2*n1*n2*(1-cos(M_PI*(1.0*j)/n1));
                    double ypart = 2*n2*n2*(1-cos(M_PI*(1.0*i)/n2));
                    // double tpart = 2*nt*nt*(1-cos(M_PI*(1.0*n)/nt));  // DCT
                    double tpart = 2*nt*nt*(1-cos(M_PI*(1.0*n+1)/nt));// DST
                    double negativeLaplacian = tpart;
                    kernel[n*n1*n2+i*n1+j]   = negativeLaplacian;
                }
            }
        }
            
    }
    // function for general rho. Boundary for t=0 only
    void solve_heat_equation_with_bdry(double* rho, double* bdry, double tau){
        memcpy(workspace,rho,n1*n2*nt*sizeof(double));
        // boundary
        for(int i=0;i<n1*n2;++i) workspace[i]              += 2*tau*nt*nt*bdry[i];
        for(int i=0;i<n1*n2;++i) workspace[(nt-1)*n1*n2+i] += 2*tau*nt*nt*rho[(nt-1)*n1*n2+i];
        fftw_execute(planIn);
        for(int i=0;i<n1*n2*nt;++i) workspace[i] /= 8.0*n1*n2*nt * (1 + tau * kernel[i]);
        fftw_execute(planOut);
        memcpy(rho,workspace,n1*n2*nt*sizeof(double));
    }

    // function for general rho. Boundary for t=0 only
    void solve_heat_equation_with_bdry(double* rho, double tau){
        memcpy(workspace,rho,n1*n2*nt*sizeof(double));
        // boundary
        for(int i=0;i<n1*n2;++i) workspace[i]              += 2*tau*nt*nt*rho[i];
        for(int i=0;i<n1*n2;++i) workspace[(nt-1)*n1*n2+i] += 2*tau*nt*nt*rho[(nt-1)*n1*n2+i];
        fftw_execute(planIn);
        for(int i=0;i<n1*n2*nt;++i) workspace[i] /= 8.0*n1*n2*nt * (1 + tau * kernel[i]);
        fftw_execute(planOut);
        memcpy(rho,workspace,n1*n2*nt*sizeof(double));
    }
};


class poisson_solver_2d{
public:
    fftw_plan planIn;
    fftw_plan planOut;
    double *workspace;
    unique_ptr<double[]> u;
    unique_ptr<double[]> kernel;

    int n1;
    int n2;

    double eta;

    poisson_solver_2d(){
        workspace=NULL;
    }
    poisson_solver_2d(const int n1, const int n2, const double eta=1) {
        this->n1=n1;
        this->n2=n2;

        this->eta=eta;

        workspace =(double*) fftw_malloc(n1*n2*sizeof(double));

        // planIn = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT10,FFTW_REDFT10, FFTW_MEASURE);
        // planOut = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT01,FFTW_REDFT01, FFTW_MEASURE);

        planIn =fftw_plan_r2r_2d(n2, n1, workspace, workspace, FFTW_REDFT10, FFTW_REDFT10, FFTW_MEASURE);
        planOut=fftw_plan_r2r_2d(n2, n1, workspace, workspace, FFTW_REDFT01, FFTW_REDFT01, FFTW_MEASURE);

        u      = make_unique<double[]>(n1*n2);
    	kernel = make_unique<double[]>(n1*n2);

        create_negative_laplacian_kernel_2d();
    }

    ~poisson_solver_2d(){
        fftw_free(workspace);
        fftw_destroy_plan(planIn);
        fftw_destroy_plan(planOut);
    }

    void create_negative_laplacian_kernel_2d(){
        
        for(int i=0;i<n2;i++){
            for(int j=0;j<n1;j++){
                double xpart = 2*n1*n2*(1-cos(M_PI*(1.0*j)/n1));
                double ypart = 2*n2*n2*(1-cos(M_PI*(1.0*i)/n2));
                // double negativeLaplacian= xpart + ypart + eta * ( xpart*xpart + ypart*ypart ) + tpart;
                double negativeLaplacian= xpart + ypart;
                kernel[i*n1+j]=negativeLaplacian;
            }
        }
            
    }

    void perform_inverse_laplacian(const double c){

    	memcpy(workspace, u.get(), n1*n2*sizeof(double));

        fftw_execute(planIn);

        for(int i=0;i<n1*n2;++i){
            double val = c + kernel[i];
            if(val==0){
                workspace[i]=0; 
            }else{
                workspace[i]/=4*(n1)*(n2)*val;
            }
        }

        fftw_execute(planOut);
    }

    void perform_inverse_laplacian_phiT(const double c){

    	memcpy(workspace, u.get(), n1*n2*sizeof(double));

        fftw_execute(planIn);

        for(int i=0;i<n1*n2;++i){
            double val = c + kernel[i];
            if(val==0){
                workspace[i]=0; 
            }else{
                workspace[i]/=4*(n1)*(n2)*val;
            }
        }

        fftw_execute(planOut);
    }

    void perform_convolution(const double* rho, const double var){
        memcpy(workspace, rho, n1*n2*sizeof(double));
        fftw_execute(planIn);
        for(int i=0;i<n1*n2;++i) workspace[i] *= exp(-kernel[i]*var*var) / (4*n1*n2);
        fftw_execute(planOut);
    }


    void perform_convolution(const double* rho, const double* phi0, const double* phi1, const double var){
        for(int i=0;i<n1*n2;++i){
            workspace[i] = rho[i] * (phi1[i] - phi0[i]);
        }

        fftw_execute(planIn);

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                double x = (1.0*j)/n1;
                double y = (1.0*i)/n2;
                workspace[i*n1+j] *= exp(-kernel[i*n1+j]*var*var);
            }
        }

        for(int i=0;i<n1*n2;++i){
            workspace[i] /= 4*n1*n2;
        }

        fftw_execute(planOut);
    }

    void solve_heat_equation(double* rho, const double tau){
        memcpy(workspace,rho,n1*n2*sizeof(double));

        fftw_execute(planIn);
        
        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                workspace[i*n1+j]*=exp(-2*tau*kernel[i*n1+j]);
            }
        }

        for(int i=0;i<n1*n2;++i){
            workspace[i]/=4.0*(n1)*(n2);
        }

        fftw_execute(planOut);

        memcpy(rho,workspace,n1*n2*sizeof(double));
    }
};

#endif