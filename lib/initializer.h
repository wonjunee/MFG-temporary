#ifndef INITIALIZER_H
#define INITIALIZER_H

#include <iostream>

class Initializer
{
public:
	int n1;
	int n2;
	int nt;

	double base;

	Initializer(){}

	Initializer(int n1, int n2, int nt, double base=0){
		this->n1=n1;
		this->n2=n2;
		this->nt=nt;
		this->base=base;
	}

	void intialize_rho0(double* rho){

		double sum=0;
		double vmax=0;
	    for(int i=0;i<n2;++i){
	        for(int j=0;j<n1;++j){

	            double x=1.0*(j+0.5)/n1;
	            double y=1.0*(i+0.5)/n2;

	            // if(pow(x-0.7,2) + pow(y-0.7,2) < pow(0.2,2)){
	            //     rho[i*n1+j] = 0.4;    
	            // }else{
	            //     rho[i*n1+j] = 0;
	            // }

	            /*
	            	Exp1
	            */
	            rho[i*n1+j] =  5*fmax(0, 0.6 * exp(-35*pow(x-0.7,2)-35*pow(y-0.7,2)) - 0.50);
	            // rho[i*n1+j] += 2*fmax(0, 0.6 * exp(-35*pow(x-0.2,2)-35*pow(y-0.2,2)) - 0.30);


	            /*
	            	Exp2
	            */
	            // rho[i*n1+j] = 0.5;

	            /*
	            	nonsymmetric example  
	            */
	            // if(x>0.1 && x<0.4 && y>0.1 && y<0.4)
	            // {
	            // 	rho[i*n1+j] = 0.3;
	            // }
	            
	            // if(x>0.3 && x<0.7 && y>0.3 && y<0.7)
	            // {
	            // 	rho[i*n1+j] = 0.3;
	            // }

	            // if(x>0.5 && x<0.8 && y>0.65 && y<0.9)
	            // {
	            // 	rho[i*n1+j] = 0.3;
	            // }

	            /* Exp2 DON'T TOUCH */
	            // rho[i*n1+j]  = 0.45 * exp(-15*pow(x-0.3,2)-15*pow(y-0.3,2)); rho[i*n1+j] += 0.45 * exp(-30*pow(x-0.8,2)-30*pow(y-0.35,2)); rho[i*n1+j] += 0.45 * exp(-25*pow(x-0.5,2)-25*pow(y-0.75,2));


	        }
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

	            double x=1.0*(j+0.5)/n1;
	            double y=1.0*(i+0.5)/n2;

	            // if(pow(x-0.7,2) + pow(y-0.7,2) < pow(0.1,2)){
	            //     rho[i*n1+j] = 0.4;    
	            // }else{
	            //     rho[i*n1+j] = 0;
	            // }

	            // if(fabs(x-0.5)<0.1 && fabs(y-0.5)<0.1) rho[i*n1+j] = 1;
	            rho[i*n1+j] = 2*fmax(0, 0.6 * exp(-45*pow(x-0.2,2)-45*pow(y-0.7,2)) - 0.4);

	            /* exp1 */
	            // rho[i*n1+j]  = 15*fmax(0.1-pow(x-0.6,2)-pow(y-0.6,2),0);

	            /* exp2 */
	            // rho[i*n1+j] = 0.6 * exp(-35*pow(x-0.6,2)-35*pow(y-0.6,2));
	            /*
	            	Ring
	            */

	            // double r = sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));
	            // rho[i*n1+j] = 70*fmax(0, 0.005 - (r-0.25)*(r-0.25));


	            /*
	            	multiple squares
	            */
	            
	            // if(x>0.35 && x<0.55 && y>0.65 && y<0.85){
	            // 	rho[i*n1+j] = 0.4;
	            // }

	            // if(x>0.55 && x<0.8 && y>0.2 && y<0.55){
	            // 	rho[i*n1+j] = 0.4;
	            // }

	            // if(x>0.14 && x<0.28 && y>0.14 && y<0.28)
	            // {
	            // 	rho[i*n1+j] = 0.4;
	            // }

	            // if(x>0.22 && x<0.36 && y>0.22 && y<0.36)
	            // {
	            // 	rho[i*n1+j] = 0.4;
	            // }

            	/* Exp2 DON'T TOUCH */
	            // rho[i*n1+j]   = 10*fmax(0.04-pow(x-0.2,2)-pow(y-0.65,2),0); rho[i*n1+j]  += 12*fmax(0.03-pow(x-0.5,2)-pow(y-0.2,2),0); rho[i*n1+j]  += 12*fmax(0.03-pow(x-0.8,2)-pow(y-0.55,2),0);
	        }
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
	            rho[i*n1+j] = 0;
	        }
	    }

	    for(int n=1;n<nt;++n){
	        for(int i=0;i<n1*n2;++i){
	            rho[n*n1*n2+i]=rho[i];
	        }
	    }
	}


	// density for vaccine
	void intialize_rho3(double* rho){

	    for(int i=0;i<n2;++i){
	        for(int j=0;j<n1;++j){
	        	double x = (j+0.5)/n1;
	        	double y = (i+0.5)/n2;
	        	rho[i*n1+j] = 0*fmax(0, 0.8 * exp(-30*pow(x-0.3,2)-30*pow(y-0.3,2)) - 0.6);
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