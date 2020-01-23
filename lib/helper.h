
double sign(double x){
    
    double s= (x>0) - (x<0);
    
    return s;
    
}

double real3rdRoot1=-.5;   // equals cos(2*M_PI/3);
double im3rdRoot1=0.86602540378;   //equals sin(2*M_PI/3);
double real3rdRoot2=-.5;  //  equals cos(4*M_PI/3)=real3rdRoot1;
double im3rdRoot2=-0.86602540378;  //equals sin(4*M_PI/3)=-im3rdRoot1;


double cubic_solve(double b, double c, double d){
    
    double b3over3=(b/3)*(b/3)*(b/3);
    
    double p=c-b*(b/3);
    double q=d+2*b3over3-b*(c/3);
    double solution=0;
    
    if(p==0){
        
        solution=-sign(q)*exp(log(fabs(q))/3.0);
        
    }else{
        double discrim=(q/2)*(q/2)+(p/3)*(p/3)*(p/3);
        
        double s=sqrt(fabs(discrim));
        
        if(discrim<0){
            
            double theta=atan2(s,-q/2);
            
            double x=s*s+q*q/4;
            double rc=exp(log(x)/6);
            
            double thetac=theta/3;
            
            double real=rc*cos(thetac);
            double im=rc*sin(thetac);
            
            double solution1=2*real;
            
            
            double solution2=2*(real*real3rdRoot1-im*im3rdRoot1);
            double solution3=2*(real*real3rdRoot2-im*im3rdRoot2);
            
            solution=fmax(solution1,fmax(solution2,solution3));
            
            
        }else if(discrim>0){
            
            double u3=-q/2+s;
            double v3=-q/2-s;
            
            double u=sign(u3)*exp(log(fabs(u3))/3);
            double v=sign(v3)*exp(log(fabs(v3))/3);
            
            solution=u+v;
            
        }else{
            solution=fmax(3*q/p, -3*q/(2*p));
            
        }
    }
    
    return solution-b/3;
    
}


