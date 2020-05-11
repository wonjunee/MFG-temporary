
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

    double eta;

    poisson_solver(){
    	workspace=NULL; u=NULL; kernel=NULL; tridiagonalWorkspace=NULL;
    }
    poisson_solver(int n1, int n2, int nt, double dx, double dy, double dt, double eta=1) {
    	this->n1=n1;
    	this->n2=n2;
    	this->nt=nt;
    	this->dx=dx;
    	this->dy=dy;
    	this->dt=dt;

    	this->eta=eta;

        workspace =(double*) fftw_malloc(n1*n2*nt*sizeof(double));

		// planIn = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT10,FFTW_REDFT10, FFTW_MEASURE);
		// planOut = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT01,FFTW_REDFT01, FFTW_MEASURE);

		planIn=fftw_plan_r2r_3d(nt, n2, n1, workspace, workspace,
                                 FFTW_REDFT10, FFTW_REDFT10, FFTW_REDFT10,
                                 FFTW_MEASURE);
    	planOut=fftw_plan_r2r_3d(nt, n2, n1, workspace, workspace,
                                  FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01,
                                  FFTW_MEASURE);

		u=new double[n1*n2*nt];
        kernel=new double[n1*n2*nt];
        tridiagonalWorkspace=new double[n1*n2*nt];

        create_negative_laplacian_kernel_2d();
    }

    ~poisson_solver(){
		delete[] u;
	    delete[] kernel;
	    delete[] tridiagonalWorkspace;
	    fftw_free(workspace);
	    fftw_destroy_plan(planIn);
	    fftw_destroy_plan(planOut);
	}

    void create_negative_laplacian_kernel_2d(){
	    
	    for(int n=0;n<nt;++n){
	    	for(int i=0;i<n2;i++){
		        for(int j=0;j<n1;j++){
		            double xpart = 2/(dx*dx)*(1-cos(M_PI*(1.0*j)/n1));
		        	double ypart = 2/(dy*dy)*(1-cos(M_PI*(1.0*i)/n2));
		        	double tpart = 2/(dt*dt)*(1-cos(M_PI*(1.0*n)/nt));
		            // double negativeLaplacian= xpart + ypart + eta * ( xpart*xpart + ypart*ypart ) + tpart;
		            double negativeLaplacian= xpart + ypart + eta * ( xpart*xpart + ypart*ypart ) + tpart;
		            kernel[n*n1*n2+i*n1+j]=negativeLaplacian;
		        }
		    }
	    }
		    
	}

	void perform_inverse_laplacian(const double c){

		fftw_execute(planIn);

		for(int i=0;i<n1*n2*nt;++i){
			if(kernel[i]==0){
				workspace[i]=0;	
			}else{
				workspace[i]/=8*(n1)*(n2)*(nt)*(1.0+kernel[i]);
			}
			
		}

		fftw_execute(planOut);
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
	}

	void backward_tridiagonal_sweep(){
	    for(int k=nt-2;k>=0;k--){
	        for(int i=0;i<n1*n2;i++){
	            u[k*n1*n2+i]=u[k*n1*n2+i]-tridiagonalWorkspace[k*n1*n2+i]*u[(k+1)*n1*n2+i];
	        }
	        
	    }
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
};