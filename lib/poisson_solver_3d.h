
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
			double val = (1+kernel[i]);
			if(val==0){
				workspace[i]=0;	
			}else{
				workspace[i]/=8*(n1)*(n2)*(nt)*val;
			}
			
		}

		fftw_execute(planOut);
	}
};


class poisson_solver_2d{
public:
    fftw_plan planIn;
    fftw_plan planOut;
    double *workspace;
    double *u;
    double *kernel;

    int n1;
    int n2;
    double dx;
    double dy;

    double eta;

    poisson_solver_2d(){
    	workspace=NULL; u=NULL; kernel=NULL;
    }
    poisson_solver_2d(int n1, int n2, double dx, double dy, double eta=1) {
    	this->n1=n1;
    	this->n2=n2;
    	this->dx=dx;
    	this->dy=dy;

    	this->eta=eta;

        workspace =(double*) fftw_malloc(n1*n2*sizeof(double));

		// planIn = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT10,FFTW_REDFT10, FFTW_MEASURE);
		// planOut = fftw_plan_r2r_2d(n2,n1, workspace, workspace, FFTW_REDFT01,FFTW_REDFT01, FFTW_MEASURE);

		planIn=fftw_plan_r2r_2d(n2, n1, workspace, workspace,
                                 FFTW_REDFT10, FFTW_REDFT10,
                                 FFTW_MEASURE);
    	planOut=fftw_plan_r2r_2d(n2, n1, workspace, workspace,
                                 FFTW_REDFT01, FFTW_REDFT01,
                                 FFTW_MEASURE);

		u=new double[n1*n2];
        kernel=new double[n1*n2];

        create_negative_laplacian_kernel_2d();
    }

    ~poisson_solver_2d(){
		delete[] u;
	    delete[] kernel;
	    fftw_free(workspace);
	    fftw_destroy_plan(planIn);
	    fftw_destroy_plan(planOut);
	}

    void create_negative_laplacian_kernel_2d(){
	    
    	for(int i=0;i<n2;i++){
	        for(int j=0;j<n1;j++){
	            double xpart = 2/(dx*dx)*(1-cos(M_PI*(1.0*j)/n1));
	        	double ypart = 2/(dy*dy)*(1-cos(M_PI*(1.0*i)/n2));
	            // double negativeLaplacian= xpart + ypart + eta * ( xpart*xpart + ypart*ypart ) + tpart;
	            double negativeLaplacian= xpart + ypart + eta * ( xpart*xpart + ypart*ypart );
	            kernel[i*n1+j]=negativeLaplacian;
	        }
	    }
		    
	}

	void perform_inverse_laplacian(const double c){

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

	void perform_inverse_laplacian_phiT(const double c, const double* phiT){

		fftw_execute(planIn);


		for(int i=0;i<n1*n2;++i){
			double val = kernel[i];
			if(phiT[i]>0){
				val += c;
			}
			if(val==0){
				workspace[i]=0;	
			}else{
				workspace[i]/=4*(n1)*(n2)*val;
			}
		}

		fftw_execute(planOut);
	}
};