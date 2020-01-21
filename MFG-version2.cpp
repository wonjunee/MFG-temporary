#include <iostream>
#include <fftw3.h>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <time.h>
#include "poisson_solver.h"

using namespace std;


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
		outfile<<A[i]<<",";
	}
	outfile.close();
}

void create_csv_file_rho(const double* rho,const double* mu,const double* nu,string filename,int n1,int n2,int nt){
	ofstream outfile;
	outfile.open(filename);
	for(int i=0;i<n1*n2;++i){
		outfile<<mu[i]<<",";
	}
	for(int i=0;i<n1*n2*(nt-2);++i){
		outfile<<rho[n1*n2+i]<<",";
	}
	for(int i=0;i<n1*n2;++i){
		outfile<<nu[i]<<",";
	}
	outfile.close();
}

double gaussian(const double x, const double y, const double mux, const double muy, const double sigmax, const double sigmay){
	return exp(-0.5*(pow(x-mux,2)/pow(sigmax,2) + pow(y-muy,2)/pow(sigmay,2) ));
}

void intialize_rho_top(double* rho,double base,int n1,int n2,int nt,double dx,double dy,double dt){

	double sum=0;

	double mux=0.5;
	double muy=0.8;
	double sigmax=0.05;
	double sigmay=0.05;

	for(int i=0;i<n2;++i){
		for(int j=0;j<n1;++j){

			double x=1.0*j*dx;
			double y=1.0*i*dy;

			rho[i*n1+j] = gaussian(x,y,mux,muy,sigmax,sigmay);
			sum += rho[i*n1+j];
		}
	}

	for(int i=0;i<n1*n2;++i){
		rho[i] /= sum*dx*dy;
	}

	for(int n=1;n<nt;++n){
		for(int i=0;i<n1*n2;++i){
			rho[n*n1*n2+i]=1;
		}
	}
}

void initialize_f_g(double* f,double* g,int n1,int n2,double dx,double dy){

	double sum1=0;
	double sum2=0;

	for(int i=0;i<n2;++i){
		for(int j=0;j<n1;++j){

			double x=1.0*j*dx;
			double y=1.0*i*dy;

			f[i*n1+j] = gaussian(x,y,0.5,0.5,0.05,0.05);
			g[i*n1+j] = gaussian(x,y,0.5,0.2,0.05,0.1);

			sum1+=f[i*n1+j];
			sum2+=g[i*n1+j];
		}
	}

	for(int i=0;i<n1*n2;++i){
		f[i]/=sum1*dx*dy*5;
		g[i]/=sum2*dx*dy;
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

	int max_iteration;
	double tolerance;

    double* rhoTmp;
    double* mx;
    double* my;
    double* mxTmp;
    double* myTmp;
    double* Phi;
    double* Psi;
    double* PsiTmp;

    double energy;
    double previous_energy;
    double previous_dual;

    double c1;
    double c2;

    Method(){
    	rhoTmp=NULL;
    	mx=NULL;
    	my=NULL;
    	mxTmp=NULL;
    	myTmp=NULL;
    	Phi=NULL;
    }

    Method(int n1, int n2, int nt, double dx, double dy, double dt, double tau, double sigma, int max_iteration, double tolerance, double c1, double c2){
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

    	rhoTmp=new double[n1*n2*nt];
	    mx=new double[n1*n2*nt];
	    my=new double[n1*n2*nt];
	    mxTmp=new double[n1*n2*nt];
	    myTmp=new double[n1*n2*nt];
	    Phi=new double[n1*n2*nt];
	    Psi=new double[n1*n2*nt];
	    PsiTmp=new double[n1*n2*nt];
    }

    void destroy_all(){
    	delete[] rhoTmp;
    	delete[] mx;
    	delete[] my;
    	delete[] mxTmp;
    	delete[] myTmp;
    	delete[] Phi;
    }

    void update_m(double* rho){

    	for(int n=0;n<nt;++n){
    		for(int i=0;i<n2;++i){
	    		for(int j=0;j<n1;++j){
	    			double nablaxPhi = 0;
	    			double nablayPhi = 0;

	    			if(j<n1-1){
	    				nablaxPhi = (mx[n*n1*n2+i*n1+j+1]-mx[n*n1*n2+i*n1+j])/dx;	
	    			}else{
	    				nablaxPhi = 0;
	    			}

	    			if(i<n2-1){
	    				nablayPhi = (my[n*n1*n2+(i+1)*n1+j]-my[n*n1*n2+i*n1+j])/dy;
	    			}else{
	    				nablayPhi = 0;
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
	    			
	    			
	    			mx[n*n1*n2+i*n1+j] -= tau * (mx[n*n1*n2+i*n1+j]/rhovalx + nablaxPhi);
	    			my[n*n1*n2+i*n1+j] -= tau * (my[n*n1*n2+i*n1+j]/rhovaly + nablayPhi);

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

    void update_rho(double* rho,const double* f,const double* g){

    	// newton's method

    	memcpy(rhoTmp, rho, n1*n2*nt*sizeof(double));

    	double tolerance_newton=1e-10;
    	int max_iteration_newton=100;
    	double error=1;
    	for(int iterNewton=0; iterNewton<max_iteration_newton; ++iterNewton){

    		error = 0;

    		for(int n=1;n<nt;++n){
    			for(int i=0;i<n2;++i){
    				for(int j=0;j<n1;++j){
    					double Fval=0;
    					double DFval=0;

    					double mval=sqrt(mx[n*n1*n2+i*n1+j]*mx[n*n1*n2+i*n1+j]+my[n*n1*n2+i*n1+j]*my[n*n1*n2+i*n1+j]);
    					double rhoval=rho[n*n1*n2+i*n1+j];
    					double fval=f[i*n1+j];
    					double dtPhi=0;

    					if(n==nt-1){
    						dtPhi=(Phi[(n)*n1*n2+i*n1+j]-Phi[(n-1)*n1*n2+i*n1+j])/(dt);	
    					}else{
    						dtPhi=(Phi[(n+1)*n1*n2+i*n1+j]-Phi[(n-1)*n1*n2+i*n1+j])/(2.0*dt);
    					}

    					if(n==nt-1){
    						Fval  = -pow(mval,2)/(2*pow(rhoval,2)) + c1*log(rhoval*exp(fval/c1)) + c1 + dtPhi + 1.0/tau*(rhoval-rhoTmp[n*n1*n2+i*n1+j]) + 1.0/dt * (Psi[i*n1+j] - Phi[(nt-1)*n1*n2+i*n1+j]);
    						DFval = pow(mval,2)/pow(rhoval,3) + c1/rhoval + 1.0/tau;
    					}else{
    						Fval  = -pow(mval,2)/(2*pow(rhoval,2)) + c1*log(rhoval*exp(fval/c1)) + c1 + dtPhi + 1.0/tau*(rhoval-rhoTmp[n*n1*n2+i*n1+j]);
    						DFval = pow(mval,2)/pow(rhoval,3) + c1/rhoval + 1.0/tau;
    					}

    					rho[n*n1*n2+i*n1+j] = fmax(0.,rhoval-0.001*Fval/DFval);

    					error += pow(Fval/DFval,2);
    				}
    			}

    			// if(n<nt-1){
    	// 		{
    	// 			double rhosum=0;
					// for(int i=0;i<n1*n2;++i){
					// 	rhosum+=rho[n*n1*n2+i];
					// }
					// for(int i=0;i<n1*n2;++i){
					// 	rho[n*n1*n2+i]/=rhosum*dx*dy;
					// }	
    	// 		}
    			
    		}


	  //   	for(int n=1;n<nt-1;++n){
			// 	double rhosum=0;
			// 	for(int i=0;i<n1*n2;++i){
			// 		rhosum+=rho[n*n1*n2+i];
			// 	}
			// 	for(int i=0;i<n1*n2;++i){
			// 		rho[n*n1*n2+i]/=rhosum*dx*dy;
			// 	}
			// }

    		if(sqrt(error*dx*dy*dt)<tolerance_newton){
    			break;
    		}
    	}


    }

    double calculate_grad_mx(double* mxTmp, int n, int i, int j){
		double mxval;

		if(j==0){
			mxval = (mxTmp[n*n1*n2+i*n1+j])/dx;
		}else{
			mxval = (mxTmp[n*n1*n2+i*n1+j]-mxTmp[n*n1*n2+i*n1+j-1])/dx;
		}

		return mxval;
    }

    double calculate_grad_my(double* myTmp, int n, int i, int j){
		double myval;

		if(i==0){
			myval = (myTmp[n*n1*n2+i*n1+j])/dy;
		}else{
			myval = (myTmp[n*n1*n2+i*n1+j]-myTmp[n*n1*n2+(i-1)*n1+j])/dy;
		}

		return myval;
    }

    void update_phi(poisson_solver& fftps){


	    // Update t = 0

		for(int i=0;i<n2;++i){
	    	for(int j=0;j<n1;++j){
	    		int n = 0;
		    	double dtrho=1.0/dt*(rhoTmp[(n+1)*n1*n2+i*n1+j]-rhoTmp[n*n1*n2+i*n1+j]);
				double nablamx=calculate_grad_mx(mxTmp,n,i,j);
				double nablamy=calculate_grad_my(myTmp,n,i,j);
				fftps.workspace[i*n1+j]=-(dtrho+nablamx+nablamy);	
	    	}
	    }

	    fftps.perform_inverse_laplacian();
	    for(int i=0;i<n1*n2;++i){
	    	Phi[i] += sigma*fftps.workspace[i];
	    }

	    // Update t = 1

	    for(int i=0;i<n2;++i){
	    	for(int j=0;j<n1;++j){
	    		int n = nt-1;
		    	double dtrho=1.0/dt*(rhoTmp[n*n1*n2+i*n1+j]-rhoTmp[(n-1)*n1*n2+i*n1+j]);
				double nablamx=calculate_grad_mx(mxTmp,n,i,j);
				double nablamy=calculate_grad_my(myTmp,n,i,j);
				fftps.workspace[i*n1+j]=-(dtrho+nablamx+nablamy);	

				fftps.workspace[i*n1+j]=-(dtrho+nablamx+nablamy);				
	    	}
	    }

	    fftps.perform_inverse_laplacian();
	    for(int i=0;i<n1*n2;++i){
	    	Phi[(nt-1)*n1*n2+i] += sigma*fftps.workspace[i];
	    }

		// Update 0 < t < 1

		for(int n=0;n<nt;++n){	
			for(int i=0;i<n2;++i){
				for(int j=0;j<n1;++j){
					if(n==0){
						fftps.workspace[i*n1+j]=0;
					}else if(n==nt-1){
						fftps.workspace[i*n1+j]=0;
					}else{
						double dtrho=0.5/dt*(rhoTmp[(n+1)*n1*n2+i*n1+j]-rhoTmp[(n-1)*n1*n2+i*n1+j]);
						double nablamx=calculate_grad_mx(mxTmp,n,i,j);
						double nablamy=calculate_grad_my(myTmp,n,i,j);

						fftps.workspace[i*n1+j]=-(dtrho+nablamx+nablamy);	
					}
				}
			}

			fftps.get_fourier_coefficients(&fftps.u[n*n1*n2]);
		}


	    fftps.forward_tridiagonal_sweep();    
	    fftps.backward_tridiagonal_sweep();
	    
	    for(int n=0;n<nt;n++){
	        for(int i=0;i<n1*n2;i++){
	            fftps.workspace[i]=fftps.u[n*n1*n2+i];
	        }
	        
	        fftps.back_to_real_space();
	        
	        for(int i=0;i<n1*n2;i++){
	            Phi[n*n1*n2+i]+=sigma*fftps.workspace[i];
	        }
	    } 
    }

    double calculate_norm_of_deltaGstar(const double* Psi,const double* g){
    	double sum=0;

    	for(int i=0;i<n2;++i){
    		for(int j=0;j<n1;++j){
    			sum += deltaGstar(Psi,i,j,g);
    		}
    	}

    	return sum*dx*dy;
    }

    double deltaGstar(const double* Psi,const int i,const int j,const double* g){
    	return exp(-Psi[i*n1+j]/c2)*g[i*n1+j];
    }

    void update_psi(poisson_solver& fftps,double* rhoTmp,const double* g){

    	double max_iteration_fp=200;

    	memcpy(PsiTmp,Psi,n1*n2*sizeof(double));

    	for(int iter_fp=0;iter_fp<max_iteration_fp;++iter_fp){

    		double norm_deltaGstar = calculate_norm_of_deltaGstar(Psi,g);

    		for(int i=0;i<n2;++i){
    			for(int j=0;j<n1;++j){
    				fftps.workspace[i*n1+j] = rhoTmp[(nt-1)*n1*n2+i*n1+j] - c2*deltaGstar(Psi,i,j,g)/norm_deltaGstar;
    			}
    		}

    		fftps.perform_inverse_laplacian();

    		for(int i=0;i<n2;++i){
    			for(int j=0;j<n1;++j){
    				Psi[i*n1+j] = PsiTmp[i*n1+j] + sigma * fftps.workspace[i*n1+j];
    			}
    		}
    	}
    }

    double calculate_energy(const double* rho,const double* f){
    	double sum=0;

    	for(int n=0;n<nt;++n){
    		for(int i=0;i<n2;++i){
    			for(int j=0;j<n1;++j){
    				double mxval=mx[n*n1*n2+i*n1+j];
    				double myval=my[n*n1*n2+i*n1+j];
    				double mval=sqrt(mxval*mxval+myval*myval);

    				double rhoval=rho[n*n1*n2+i*n1+j];

    				double fval=f[i*n1+j];
    				if(rhoval>0){
    					sum+=pow(mval,2)/(2.0*rhoval) + c1*rhoval*log(rhoval*exp(fval/c1));	
    				}
    			}
    		}
    	}

    	return sum*dx*dy*dt;
    }

    void run(poisson_solver& fftps, double* rho,const double* f,const double* g){

	    previous_energy=1;
	    previous_dual=1;
	    double error, dual_gap;
	    int iterPDHG;

    	for(iterPDHG=0; iterPDHG<max_iteration; ++iterPDHG){

    		memcpy(mxTmp,mx,n1*n2*nt*sizeof(double));
    		memcpy(myTmp,my,n1*n2*nt*sizeof(double));
    		update_m(rho);

    		memcpy(rhoTmp,rho,n1*n2*nt*sizeof(double));
    		update_rho(rho,f,g);

    		for(int i=0;i<nt*n1*n2;++i){
    			rhoTmp[i]=2*rho[i]-rhoTmp[i];
    			mxTmp[i]=2*mx[i]-mxTmp[i];
    			myTmp[i]=2*my[i]-myTmp[i];
    		}

    		update_phi(fftps);

    		update_psi(fftps,rho,g);

    		double energy=calculate_energy(rho,f);

    		error=fabs((energy-previous_energy)/previous_energy);

    		previous_energy=energy;

    		if(iterPDHG%10==0){
    			cout<<"iter : "<<iterPDHG<< " energy : "<<energy<< " error : " << error<<endl;	

    			create_csv_file(rho,"./data/rho.csv",n1,n2,nt);
    		}
    		
    		if(error<tolerance){
    			break;
    		}
    	}

    	cout<<"The method is done!!"<<endl;
    	cout<<"iterPDHG : "<<iterPDHG<< " energy : "<<energy<< " error : " << error<<endl;

	    cout<<"Iterations     : "<<iterPDHG<<endl;
	    cout<<"Energy         : "<<energy<<endl;
	    cout<<"Relative Error : "<<error<<endl;
    }

}; // Method class

int main(int argc, char **argv)
{

    if(argc!=8){
        cout << "Need to do the following : " << endl;
        cout << "./main.exe [n1] [n2] [nt] [tau] [sigma] [tolerance] [max iteration]" << endl;
        exit(1);
    }

    // Parameters for Grids
    
    int n1 = stoi(argv[1]);  // x panels
    int n2 = stoi(argv[2]);  // y panels
    int nt = stoi(argv[3]);  // t panels
    double tau=stod(argv[4]);
    double sigma=stod(argv[5]);
    double tolerance = stod(argv[6]);
    int max_iteration=stoi(argv[7]);
    
    create_csv_file_for_parameters(n1,n2,nt);

    double base=0;

    double dx=1.0/n1;
    double dy=1.0/n2;
    double dt=1.0/(nt-1.0);

    double c1=0.5;
    double c2=0.1;

    cout<<"XXX G-Prox PDHG XXX"<<endl;
    cout<<endl;
    cout<<"n1 : "<<n1<<" n2 : "<<n2<<" nt : "<<nt<<" base : "<<base<<endl;
    cout<<"dx : "<<scientific<<dx<<endl;
    cout<<"dy : "<<scientific<<dy<<endl;
    cout<<"dt : "<<scientific<<dt<<endl;
    cout<<fixed;
    cout<<"tau   : "<<tau<<endl;
    cout<<"sigma : "<<sigma<<endl;
    cout<<"max_iteration : "<<max_iteration<<endl;
    cout<<"tolerance     : "<<scientific<<tolerance<<endl;


    double* rho=new double[n1*n2*nt];
    double* f=new double[n1*n2];
    double* g=new double[n1*n2];

    clock_t t;
    t = clock();
    poisson_solver fftps(n1,n2,nt,dx,dy,dt);
    t = clock() - t;
    printf ("\nCPU time for setting up FFT: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  

	intialize_rho_top(rho,base,n1,n2,nt,dx,dy,dt);
	initialize_f_g(f,g,n1,n2,dx,dy);


	Method method(n1, n2, nt, dx, dy, dt, tau, sigma, max_iteration, tolerance, c1, c2);

    cout<<"\nXXX Starting Iterations XXX"<<endl;

    t = clock();

    method.run(fftps,rho,f,g);

    t = clock() - t;

    printf ("CPU time for Iterations: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);

    create_csv_file(rho,"./data/rho.csv",n1,n2,nt);

    create_csv_file_for_parameters(n1,n2,nt);
	create_csv_file(rho,"rho.csv",n1,n2,nt);

	delete[] rho;
	delete[] f;
	delete[] g;

    fftps.destroy_all_fftps();
}