#ifndef INITIALIZER_H
#define INITIALIZER_H

#include <iostream>

class Initializer
{
public:
	int n1;
	int n2;
	int nt;
	double dx;
	double dy;
	double dt;

	double base;

	Initializer(){}

	Initializer(int n1, int n2, int nt, double dx, double dy, double dt, double base=0){
		this->n1=n1;
		this->n2=n2;
		this->nt=nt;
		this->dx=dx;
		this->dy=dy;
		this->dt=dt;
		this->base=base;
	}

	void intialize_rho0(double* rho){

		double sum=0;
		double vmax=0;
	    for(int i=0;i<n2;++i){
	        for(int j=0;j<n1;++j){

	            double x=1.0*(j+0.5)*dx;
	            double y=1.0*(i+0.5)*dy;

	            // if(pow(x-0.5,2) + pow(y-0.5,2) >= pow(0.1,2)){
	            //     rho[i*n1+j] = 1;    
	            // }else{
	            //     rho[i*n1+j] = 0;
	            // }
	            rho[i*n1+j] = exp(-20*pow(x-0.5,2)-20*pow(y-0.5,2)) + base;
	            // sum += rho[i*n1+j];
	            vmax = fmax(vmax, rho[i*n1+j]);
	        }
	    }

	    for(int i=0;i<n1*n2;++i){
	        // rho[i] *= (n1*n2)/sum;
	        rho[i] /= vmax * 2;
	    }


	    for(int n=1;n<nt;++n){
	        for(int i=0;i<n1*n2;++i){
	            rho[n*n1*n2+i]=rho[i];
	        }
	    }
	}

	void intialize_rho1(double* rho){

	    double sum=0;
	    double vmax = 0;

	    for(int i=0;i<n2;++i){
	        for(int j=0;j<n1;++j){

	            double x=1.0*(j+0.5)*dx;
	            double y=1.0*(i+0.5)*dy;

	            // if(pow(x-0.5,2) + pow(y-0.5,2) < pow(0.1,2)){
	            //     rho[i*n1+j] = 1;    
	            // }else{
	            //     rho[i*n1+j] = 0;
	            // }

	            // rho[i*n1+j] = exp(-60*pow(x-0.6,2)-60*pow(y-0.6,2)) + base;

	            rho[i*n1+j] = fmax(0.03-pow(x-0.5,2)-pow(y-0.5,2),0);
	            // sum += rho[i*n1+j];
	            vmax = fmax(vmax, rho[i*n1+j]);
	        }
	    }

	    for(int i=0;i<n1*n2;++i){
	        // rho[i] *= (n1*n2)/sum*0.5;
	        rho[i] /= vmax * 3;
	    }

	    for(int n=1;n<nt;++n){
	        for(int i=0;i<n1*n2;++i){
	            rho[n*n1*n2+i]=rho[i];
	        }
	    }
	}

	void intialize_rho2(double* rho){

	    for(int i=0;i<n2;++i){
	        for(int j=0;j<n1;++j){
	            rho[i*n1+j] = 1e-4;
	        }
	    }

	    for(int n=1;n<nt;++n){
	        for(int i=0;i<n1*n2;++i){
	            rho[n*n1*n2+i]=rho[i];
	        }
	    }
	}

	/*
		C: a set of obstacle
		f = 1 if x in C
		f = 0 if x not in C
	*/
	void intialize_f(double* f){

	    for(int i=0;i<n2;++i){
	    	for(int j=0;j<n1;++j){
	    		double x = (j+0.5)/n1;
	    		double y = (i+0.5)/n2;
	    		if(fabs(x-0.5)<0.3 && fabs(y-0.5)<0.3){
	    			f[i*n1+j] = 0;
	    		}else{
	    			f[i*n1+j] = -1;
	    		}
	    	}
	    }
	}


	void renormalize_all(double* rho[]){
		{
			double maxval = 0;
			for(int i=0;i<n1*n2;++i){
				double s = 0;
				for(int k=0;k<3;++k){
					s += rho[k][k];
				}
				maxval = fmax(maxval, s);
			}
			maxval = fmax(1,maxval);
			for(int k=0;k<3;++k){
				for(int i=0;i<n1*n2;++i){
					rho[k][i] /= maxval*1.5;
				}
			}
		}

		for(int n=1;n<nt;++n){
			for(int k=0;k<3;++k){
				for(int i=0;i<n1*n2;++i){
					rho[k][n*n1*n2+i] = rho[k][i];
				}
			}
		}
	}
};




#endif