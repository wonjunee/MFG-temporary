#ifndef METHOD_H
#define METHOD_H
#include <iostream>

#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstring>
// #include <omp.h> // parallelization
#include <fftw3.h>
#include <memory>
// #include "poisson_solver.h"
#include "poisson_solver_3d.h"
#include "helper.h"
#include <future>
#include <thread>
#include <vector>


#ifndef ASYNC
#define ASYNC 0
#endif

int K = ASYNC;

std::mutex mtx;

static const double LARGE_NUMBER = 99999999999;    

class Method{
public:
    // member variables
        int n1;
        int n2;
        int nt;

        int im;
        int ip;
        int jm;
        int jp;

        double tau;
        double sigma;

        double m1_coeff_;
        double m2_coeff_;

        int max_iteration;
        double tolerance;

        double energy;
        double previous_energy;
        double previous_dual;

        double M0;

        static unique_ptr<poisson_solver>     fftps;
        static unique_ptr<poisson_solver_DST> fftpsDST;
        static unique_ptr<poisson_solver_2d>  fftps2d;

        // dynamical arrays
        shared_ptr<double[]> mx;
        shared_ptr<double[]> my;
        shared_ptr<double[]> m2;
        shared_ptr<double[]> phi;
        shared_ptr<double[]> phitmp;
        shared_ptr<double[]> rhotmp;

    // default constructor
    Method(){}

    // constructor
    Method(const int n1, const int n2, const int nt, const double tau, const double sigma, const int max_iteration, const double tolerance){

        this->n1=n1;                        // grid size x-axis
        this->n2=n2;                        // grid size y-axis
        this->nt=nt;                        // grid size t-axis
        this->max_iteration=max_iteration;  // the maximum iteration of PDHG
        this->tolerance    =tolerance;      // the tolerance of the PDHG
        this->tau          =tau;            // the step size for the primal variables (rho, m)
        this->sigma        =sigma;          // the step size for the dual variable (phi)

        // initializing m, phi, rhotmp
        mx      = shared_ptr<double[]>(new double[n1*n2*nt]);
        my      = shared_ptr<double[]>(new double[n1*n2*nt]);
        m2      = shared_ptr<double[]>(new double[n1*n2*nt]);
        phi     = shared_ptr<double[]>(new double[n1*n2*nt]);
        phitmp  = shared_ptr<double[]>(new double[n1*n2*nt]);
        rhotmp  = shared_ptr<double[]>(new double[n1*n2*nt]);

        // initialize the arrays to be 0
        memset(mx.get(),      0, n1*n2*nt*sizeof(double));
        memset(my.get(),      0, n1*n2*nt*sizeof(double));
        memset(m2.get(),      0, n1*n2*nt*sizeof(double));
        memset(phi.get(),     0, n1*n2*nt*sizeof(double));
        memset(phitmp.get(),  0, n1*n2*nt*sizeof(double));
        memset(rhotmp.get(),  0, n1*n2*nt*sizeof(double));

        // timing the initialization of FFT
        clock_t t;
        t = clock();
            
        // initialize 3d fftps and 2d fftps
        fftps    = make_unique<poisson_solver>(n1,n2,nt);
        fftpsDST = make_unique<poisson_solver_DST>(n1,n2,nt);
        fftps2d  = make_unique<poisson_solver_2d>(n1,n2);
        
        t = clock() - t;
        printf ("\nCPU time for setting up FFT: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);

        // setup coefficients for m1 and m2
        m1_coeff_ = 1;
        m2_coeff_ = 1;
    }

    // Destructor
    ~Method(){     }

    inline void setup_indices(int& im, int& ip, int& jm, int& jp, const int i, const int j){
        im = fmax(0,i-1);
        ip = fmin(n2-1,i+1);

        jm = fmax(0,j-1);
        jp = fmin(n1-1,j+1);
    }

    inline void perform_upwind_scheme(double& muxp, double& muxm, double& muyp, double& muym, const shared_ptr<const double[]>& phi, const int i, const int j){

        setup_indices(im,ip,jm,jp,i,j);
        muxp = 1.0 * n1 * (phi[i*n1+jp] - phi[i*n1+j]);
        muxm = 1.0 * n1 * (phi[i*n1+j] - phi[i*n1+jm]);
        muyp = 1.0 * n2 * (phi[ip*n1+j] - phi[i*n1+j]);
        muym = 1.0 * n2 * (phi[i*n1+j] - phi[im*n1+j]);
    }

    // update mx and my

    void update_m_for_loop_per_nt(double* mx, double* my, const double* rho, const double* phi, const int n_start, const int n_end){
        // mtx.lock();
        int im,ip,jm,jp;
        for(int n=n_start;n<n_end;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int idx = n*n1*n2+i*n1+j;

                    setup_indices(im,ip,jm,jp,i,j);

                    double nablaxphi = 1.0*n1*(phi[n*n1*n2+i*n1+jp]-phi[n*n1*n2+i*n1+j]);
                    double nablayphi = 1.0*n2*(phi[n*n1*n2+ip*n1+j]-phi[n*n1*n2+i*n1+j]);
                    
                    double rhoval=rho[idx];
                    double mxval =mx[idx];
                    double myval =my[idx];

                    mx[idx] = (tau*rhoval)/(tau + rhoval) * (mxval/tau + nablaxphi);
                    my[idx] = (tau*rhoval)/(tau + rhoval) * (myval/tau + nablayphi);
                }
            }   
        }
        // mtx.unlock();
    }

    void update_m(shared_ptr<double[]>& mx, shared_ptr<double[]>& my, const shared_ptr<const double[]>& rho, const shared_ptr<const double[]>& phi){
#if ASYNC>0
        std::vector<std::future<void> > changes;
        for(int k=0;k<K;++k){
            changes.push_back(std::async(std::launch::async|std::launch::deferred, &Method::update_m_for_loop_per_nt, this, mx.get(), my.get(), rho.get(), phi.get(), k*nt/K, (k+1)*nt/K));
        }
        for(int k=0;k<K;++k){
            changes[k].get();
        }
#else
        update_m_for_loop_per_nt(mx.get(), my.get(), rho.get(), phi.get(), 0, nt);
#endif

    }

    void update_m2_for_loop_per_nt(double* m2, const double* rho, const double* phi, const int n_start, const int n_end){
        // mtx.lock();
        int im,ip,jm,jp;

        for(int n=n_start;n<n_end;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int idx = n*n1*n2+i*n1+j;
                    double rhoval=rho[idx];
                    double phival=phi[idx];
                    double m2val =m2 [idx];
                    m2[idx] = (tau*rhoval)/(tau + rhoval) * (m2val/tau + phival);
                }
            }   
        }   
        // mtx.unlock();
    }
    // update m2
    void update_m2(shared_ptr<double[]>& m2, const shared_ptr<const double[]>& rho, const shared_ptr<const double[]>& phi){
#if ASYNC>0
        std::vector<std::future<void> > changes;
        for(int k=0;k<K;++k){
            changes.push_back(std::async(std::launch::async|std::launch::deferred, &Method::update_m2_for_loop_per_nt, this, m2.get(), rho.get(), phi.get(), k*nt/K, (k+1)*nt/K));
        }
        for(int k=0;k<K;++k){
            changes[k].get();
        }
#else
        update_m2_for_loop_per_nt(m2.get(), rho.get(), phi.get(), 0, nt);
#endif


    }

    // update m, phi values that are needed to calculate rho
    inline void calculate_rho_related(double& mvalue, double& Dtphi, const int n, const int i, const int j, const double* mx, const double* my){

        double mxvalue, myvalue;
        int im,ip,jm,jp;
        setup_indices(im,ip,jm,jp,i,j);

        mxvalue = mx[n*n1*n2+i*n1+j];
        myvalue = my[n*n1*n2+i*n1+j];
        
        mvalue = sqrt(mxvalue*mxvalue + myvalue*myvalue);

        if(n>0) Dtphi=1.0*nt*(phi[n*n1*n2+i*n1+j]-phi[(n-1)*n1*n2+i*n1+j]);
        else    Dtphi=0;
    }

    // update rho
    void update_rho_for_loop_per_nt(double* rho, const double* rhotmp, const double* mx, const double* my, const double* m2, const int n_start, const int n_end){

        for(int n=n_start;n<n_end;++n){
            if(n==0 || n==nt-1) continue;
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int idx = n*n1*n2+i*n1+j;

                    // initialize mval and Dtphi
                    double mval  =0;
                    double Dtphi =0;

                    // calculate |m| and Dt phi
                    calculate_rho_related(mval, Dtphi, n, i, j, mx, my);

                    // get the rest of the files
                    double rhoval    = rho[idx];
                    double m2val     = m2[idx];
                    
                    // solving cubic polynomical: x^3 + a x^2 + b x + c = 0
                    double aval = - rhoval - tau * Dtphi;
                    double bval = 0.0;
                    double cval = - 0.5 * tau * (mval*mval + m2_coeff_ * m2val * m2val);

                    // find cubic solution using cubic_solve function from helper function.
                    double newrhovalue=cubic_solve(aval, bval, cval);

                    // update rho
                    rho[idx] = fmax(0,newrhovalue);
                }
            }
        }
    }

    // update rho
    void update_rho(const shared_ptr<double[]>& rho, const shared_ptr<const double[]>& rhotmp, const shared_ptr<const double[]>& mx, const shared_ptr<const double[]>& my, const shared_ptr<const double[]>& m2){

#if ASYNC>0
        std::vector<std::future<void> > changes;
        for(int k=0;k<K;++k){
            changes.push_back(std::async(std::launch::async|std::launch::deferred, &Method::update_rho_for_loop_per_nt, this, rho.get(), rhotmp.get(), mx.get(), my.get(), m2.get(), k*nt/K, (k+1)*nt/K));
        }
        for(int k=0;k<K;++k){
            changes[k].get();
        }
#else
        update_rho_for_loop_per_nt(rho.get(), rhotmp.get(), mx.get(), my.get(), m2.get(), 0, nt);
#endif
    }

    inline double calculate_grad_mx(const double* mxTmp, const int n, const int im, const int i, const int ip, const int jm, const int j, const int jp){
        return n1*(mxTmp[n*n1*n2+i*n1+j]-mxTmp[n*n1*n2+i*n1+jm]);
    }

    inline double calculate_grad_my(const double* myTmp, const int n, const int im, const int i, const int ip, const int jm, const int j, const int jp){
        return n2*(myTmp[n*n1*n2+i*n1+j]-myTmp[n*n1*n2+im*n1+j]);
    }

    inline double calculate_dtrho(const double* rho, const int n, const int i, const int j){
    	double dtrho=0;
    	if(n==nt-1) dtrho=0;
        else        dtrho=1.0*nt*(rho[(n+1)*n1*n2+i*n1+j]-rho[(n)*n1*n2+i*n1+j]); 
        return dtrho;
    }   

    void update_phi_for_loop_per_nt(const double* rho, const double* mx, const double* my, const double* m2, const int n_start, const int n_end){

        int n,ip,im,jp,jm,idx;

        for(int n=n_start;n<n_end;++n){  
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    idx = n*n1*n2+i*n1+j;

                    setup_indices(im,ip,jm,jp,i,j);

                    double dtrho  =calculate_dtrho(rho, n, i, j);
                    double nablamx=calculate_grad_mx(mx,n,im,i,ip,jm,j,jp);
                    double nablamy=calculate_grad_my(my,n,im,i,ip,jm,j,jp);
                    double m2val  =m2[idx];

                    fftps->u[idx]  = (dtrho + nablamx + nablamy - m2_coeff_ * m2val); 
                }
            }
        }
    }

    double update_phi_for_loop_fftps(const int n_start, const int n_end){
        double error = 0;
        for(int n=n_start;n<n_end;++n){
            for(int i=0;i<n1*n2;++i){
                int idx = n*n1*n2+i;
                double fftps_val = fftps->workspace[idx];
                phi[idx] += sigma * fftps_val;
                error    += fftps->u[idx] * fftps_val;
            }
        }
        return error;
    }

    double update_phi(const shared_ptr<const double[]>& rho, const shared_ptr<const double[]>& mx, const shared_ptr<const double[]>& my, const shared_ptr<const double[]>& m2){

#if ASYNC>0
        std::vector<std::future<void> > changes;
        for(int k=0;k<K;++k){
            changes.push_back(std::async(std::launch::async|std::launch::deferred, &Method::update_phi_for_loop_per_nt, this, rho.get(), mx.get(), my.get(), m2.get(), k*nt/K, (k+1)*nt/K));
        }
        for(int k=0;k<K;++k){
            changes[k].get();
        }
#else
        update_phi_for_loop_per_nt(rho.get(), mx.get(), my.get(), m2.get(), 0, nt);
#endif

        fftps->perform_inverse_laplacian();
        
        double error = 0;

#if ASYNC>0
        std::vector<std::future<double> > changes2;

        for(int k=0;k<K;++k){
            changes2.push_back(std::async(std::launch::async|std::launch::deferred, &Method::update_phi_for_loop_fftps, this, k*nt/K, (k+1)*nt/K));
        }
        for(int k=0;k<K;++k){
            error += changes2[k].get();
        }            
#else
        error = update_phi_for_loop_fftps(0, nt);
#endif
        
        return error/(1.0*n1*n2*nt);

    }


    double calculate_energy(const shared_ptr<const double[]>& rho, const shared_ptr<const double[]>& mx, const shared_ptr<const double[]>& my, const shared_ptr<const double[]>& m2){
        double sum = 0;

        // TODO
        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    // setup index
                    int idx = n*n1*n2+i*n1+j;

                    double mval=0;
                    double Dtphi =0;
                    calculate_rho_related(mval, Dtphi, n, i, j, mx.get(), my.get());

                    double rhoval = rho[idx];
                    double m2val  = m2 [idx];
                    if(rhoval < 1e-8){
                        sum += 0;
                    }else{
                        sum += 0.5 * (mval*mval + m2val*m2val) / rhoval;
                    }
                    
                }
            }
        }

        return sum / (1.0*n1*n2*nt);
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

    void display_log(const int iterPDHG, const double tau, const double sigma, const double energy, const double error, const double sanity_value) const{
        printf("iter: %5d tau: %5.2f sigma: %5.2f energy: %10.4e rel error: %10.4e dual error: %10.4e\n", iterPDHG+1, tau, sigma, energy, error, sanity_value);
    }

    void run(const shared_ptr<double[]>& rho, int skip=1){

        if(K > 0){
            std::cout << "--------------------\n";
            std::cout << "Running on " << K << " threads\n";
            std::cout << "--------------------\n";
        }

        previous_energy=1;
        previous_dual=1;
        double error=1, dual_gap=1, energy=1, dual=0, sanity_value=1;
        int iterPDHG;

        double beta_1 = 1.5;
        double beta_2 = 0.9;

        double sanity_value_previous = 100;

        double smooth_param = 1e-6;

        for(iterPDHG=0; iterPDHG<max_iteration; ++iterPDHG){


            // update phi
            memcpy(phitmp.get(),phi.get(),n1*n2*nt*sizeof(double));
            sanity_value  = update_phi(rho,mx,my,m2);
            sanity_value_previous = sanity_value;

#if ASYNC>0
            std::vector<std::future<void> > changes;

            for(int k=0;k<K;++k){
                changes.push_back(std::async(std::launch::async|std::launch::deferred, [&](const int n_start, const int n_end){
                                                                                                    for(int n=n_start; n<n_end;++n){
                                                                                                        for(int i=0;i<n1*n2;++i){
                                                                                                            int idx = n*n1*n2+i;
                                                                                                            phi[idx] = 2*phi[idx] - phitmp[idx];
                                                                                                        }
                                                                                                    }
                                                                                                }, 
                                                            k*nt/K, (k+1)*nt/K));
            }
            for(int k=0;k<K;++k){
                changes[k].get();
            }  
#else            
            for(int i=0;i<n1*n2*nt;++i){
                phi[i] = 2*phi[i] - phitmp[i];
            }
#endif

            // get the data before updates
            memcpy(rhotmp.get(),rho.get(),n1*n2*nt*sizeof(double));

            int skip2 = 100;

            update_rho(rho,rhotmp,mx,my,m2);
            update_m (mx,my,rhotmp,phi);
            // update_m2(m2,rhotmp,phi);
    
            
            // CALCULATE ENERGY
            energy = calculate_energy(rho, mx, my, m2);
            error=fabs((energy-previous_energy)/previous_energy);
            previous_energy=energy;

            if((iterPDHG+1)%skip==0){
                display_log(iterPDHG,  tau,  sigma,  energy,  error, sanity_value);
                create_bin_file(rho.get(), n1*n2*nt, "./data/rho.csv");
            }

            if((iterPDHG>20 ) && (fabs(error)<tolerance)) break;

            if(false){ // smoothing the density using the heat equation
                memcpy(rhotmp.get(),rho.get(),n1*n2*nt*sizeof(double));
                fftpsDST->solve_heat_equation_with_bdry(rho.get(), smooth_param);
                memcpy(&rho[0],    &rhotmp[0],    n1*n2*sizeof(double));
                memcpy(&rho[(nt-1)*n1*n2], &rhotmp[(nt-1)*n1*n2], n1*n2*sizeof(double));    
            }
            
        }

        cout<<"The method is done!!"<<endl;
        display_log(iterPDHG,  tau,  sigma,  energy,  error, sanity_value);

        cout<<"Iterations     : "<<iterPDHG<<endl;
        cout<<"Energy         : "<<energy<<endl;
        cout<<"Relative Error : "<<error<<endl;
    }
}; // Method class

// initializing static member variable
unique_ptr<poisson_solver>     Method::fftps;
unique_ptr<poisson_solver_DST> Method::fftpsDST;
unique_ptr<poisson_solver_2d>  Method::fftps2d;

#endif