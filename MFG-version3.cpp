#include <iostream>
#include <fftw3.h>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <time.h>
#include "poisson_solver.h"
#include "helper.h"

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
		outfile<<A[i]<<"\n";
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
			g[i*n1+j] = 10.0 * (pow(x-0.5,2)/pow(2.0,2) + pow(y-0.2,2)/pow(2.0,2));

			sum1+=f[i*n1+j];
			sum2+=g[i*n1+j];
		}
	}

	for(int i=0;i<n1*n2;++i){
		f[i]/=sum1*dx*dy*10;
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

    double* f;
    double* g;

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
    	Psi=NULL;
    	PsiTmp=NULL;
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

    void setup_f_g(const double* f,const double* g){
    	memcpy(this->f,f,n1*n2*sizeof(double));
    	memcpy(this->g,g,n1*n2*sizeof(double));
    }

    void destroy_all(){
    	delete[] rhoTmp;
    	delete[] mx;
    	delete[] my;
    	delete[] mxTmp;
    	delete[] myTmp;
    	delete[] Phi;
    	delete[] Psi;
    }

    void update_m(double* rho){

    	memcpy(mxTmp,mx,n1*n2*nt*sizeof(double));
		memcpy(myTmp,my,n1*n2*nt*sizeof(double));

    	for(int n=0;n<nt;++n){
    		for(int i=0;i<n2;++i){
	    		for(int j=0;j<n1;++j){
	    			double nablaxPhi = 0;
	    			double nablayPhi = 0;

	    			if(j<n1-1){
	    				nablaxPhi = (Phi[n*n1*n2+i*n1+j+1]-Phi[n*n1*n2+i*n1+j])/dx;	
	    			}else{
	    				nablaxPhi = 0;
	    			}

	    			if(i<n2-1){
	    				nablayPhi = (Phi[n*n1*n2+(i+1)*n1+j]-Phi[n*n1*n2+i*n1+j])/dy;
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
	    			
	    			mx[n*n1*n2+i*n1+j] = rhovalx/(tau+rhovalx) * (mx[n*n1*n2+i*n1+j] - tau * nablaxPhi);
	    			my[n*n1*n2+i*n1+j] = rhovaly/(tau+rhovaly) * (my[n*n1*n2+i*n1+j] - tau * nablayPhi);

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

    	memcpy(rhoTmp, rho, n1*n2*nt*sizeof(double));

    	for(int n=1;n<nt-1;++n){
    		for(int i=0;i<n2;++i){
    			for(int j=0;j<n1;++j){

    				double rhoTmpvalue=rhoTmp[n*n1*n2+i*n1+j];
    				double Psivalue=Psi[n*n1*n2+i*n1+j];
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
    				double mvalue=sqrt(pow(mxvalue,2)+pow(myvalue,2));
    				double DtPhi;

    				if(n==nt-1){
    					DtPhi=1.0/dt*(Phi[n*n1*n2+i*n1+j]-Phi[(n-1)*n1*n2+i*n1+j]);
    				}else{
    					// DtPhi=0.5/dt*(Phi[(n+1)*n1*n2+i*n1+j]-Phi[(n-1)*n1*n2+i*n1+j]);
                        DtPhi=1.0/dt*(Phi[(n+1)*n1*n2+i*n1+j]-Phi[(n)*n1*n2+i*n1+j]);
    				}


    				double newrhovalue=cubic_solve(-rhoTmpvalue+tau*DtPhi+tau*Psivalue, 0, -0.5*tau*mvalue*mvalue);
    				rho[n*n1*n2+i*n1+j]=fmax(0,newrhovalue);
    			}
    		}
    	}

    	for(int n=1;n<nt-1;++n){
    		double sum=0;
    		for(int i=0;i<n1*n2;++i){
    			sum+=rho[n*n1*n2+i];
    		}
    		for(int i=0;i<n1*n2;++i){
    			rho[n*n1*n2+i]/=sum*dx*dy;
    		}
    	}

    	// rho(1,x) = dGstar(phi(1,x))

    	double sum1=0;
    	for(int i=0;i<n2;++i){
    		for(int j=0;j<n1;++j){
    			// rho[(nt-1)*n1*n2+i*n1+j] = exp(Phi[(nt-1)*n1*n2+i*n1+j]/c2) * g[i*n1+j];
    			rho[(nt-1)*n1*n2+i*n1+j] = exp( (Phi[(nt-1)*n1*n2+i*n1+j] - g[i*n1+j] - c2)/c2);
    			sum1 += rho[(nt-1)*n1*n2+i*n1+j];
    		}
    	}
    	for(int i=0;i<n2;++i){
    		for(int j=0;j<n1;++j){
    			rho[(nt-1)*n1*n2+i*n1+j] /= sum1*dx*dy;
    		}
    	}
    }

    void update_psi(const double* rho, const double* f){

    	memcpy(PsiTmp, Psi, n1*n2*nt*sizeof(double));

    	int max_iteration_GD=5000;
    	double tolerance_GD=1e-12;
    	double alpha=0.5;

    	for(int iterGD=0;iterGD<max_iteration_GD;++iterGD){

    		double error=0;

    		for(int n=0;n<nt;++n){
    			for(int i=0;i<n2;++i){
    				for(int j=0;j<n1;++j){
    					double K1 = exp((Psi[n*n1*n2+i*n1+j]-f[i*n1+j])/c1 - 1);
    					double F  = rho[n*n1*n2+i*n1+j] - K1 - 1.0/sigma*(Psi[n*n1*n2+i*n1+j] - PsiTmp[n*n1*n2+i*n1+j]);
    					double DF = - 1.0/c1 * K1 - 1.0/sigma;
    					Psi[n*n1*n2+i*n1+j] -= alpha * F/DF;
    					error += pow(F/DF,2);
    				}
    			}
    		}

    		error = sqrt(error*dx*dy*dt);

    		if(error<tolerance_GD){
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

    void update_phi(poisson_solver& fftps,  const double* g){

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


	    	// for(int i=0;i<n2;++i){
		    // 	for(int j=0;j<n1;++j){
		    // 		Phi[(nt-1)*n1*n2+i*n1+j] = (c2 + c2*log(fmax(1e-10,rhoTmp[(nt-1)*n1*n2+i*n1+j])) + g[i*n1+j]);
		    // 	}
		    // }

	    for(int i=0;i<n2;++i){
	    	for(int j=0;j<n1;++j){
	    		int n = nt-1;
		    	double dtrho=1.0/dt*(rhoTmp[(n)*n1*n2+i*n1+j]-rhoTmp[(n-1)*n1*n2+i*n1+j]);
				double nablamx=calculate_grad_mx(mxTmp,n,i,j);
				double nablamy=calculate_grad_my(myTmp,n,i,j);
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
						// double dtrho=0.5/dt*(rhoTmp[(n+1)*n1*n2+i*n1+j]-rhoTmp[(n-1)*n1*n2+i*n1+j]);
                        double dtrho=1.0/dt*(rhoTmp[(n)*n1*n2+i*n1+j]-rhoTmp[(n-1)*n1*n2+i*n1+j]);
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
    					sum+=pow(mval,2)/(2.0*rhoval) + c1*rhoval*(log(rhoval)+fval);	
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

	    // Initialize Phi at t=1

		    // for(int i=0;i<n2;++i){
		    // 	for(int j=0;j<n1;++j){
		    // 		Phi[(nt-1)*n1*n2+i*n1+j] = c2 + c2*log(rho[(nt-1)*n1*n2+i*n1+j]+1) + g[i*n1+j];
		    // 	}
		    // }

    	for(iterPDHG=0; iterPDHG<max_iteration; ++iterPDHG){

    		update_rho(rho,f,g);
    		update_m(rho);
    		
    		for(int i=0;i<nt*n1*n2;++i){
    			mxTmp[i]=2*mx[i]-mxTmp[i];
    			myTmp[i]=2*my[i]-myTmp[i];
    		}

    		for(int n=1;n<nt-1;++n){
    			for(int i=0;i<n1*n2;++i){
    				rhoTmp[n*n1*n2+i]=2*rho[n*n1*n2+i]-rhoTmp[n*n1*n2+i];	
    			}
    		}

    		update_psi(rho,f);
    		update_phi(fftps,g);

    		double energy=calculate_energy(rho,f);

    		error=fabs((energy-previous_energy)/previous_energy);

    		previous_energy=energy;

    		if(iterPDHG%1==0){
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

    double c1=0.01;
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

    // create_csv_file(rho,"./data/rho.csv",n1,n2,nt);

    create_csv_file_for_parameters(n1,n2,nt);

	delete[] rho;
	delete[] f;
	delete[] g;

    fftps.destroy_all_fftps();
    method.destroy_all();
}