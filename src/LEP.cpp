#include <math.h>
#include <vector>
#include <iostream>
#include <R.h>
using namespace std;

double soft(double z,double lambda)
{
	double s;
	int tmp1,tmp2;

	if(z>=-lambda)
		tmp1=1;
	else
		tmp1=0;
	if(z>lambda)
		tmp2=1;
	else
		tmp2=0;
	int ind=tmp1+tmp2;

	if(ind==2)
		s=z-lambda;
	else if(ind==0)
		s=z+lambda;
	else
		s=0;

	return s;
}

extern "C"{
void gaulep(double *X,double *Y,double *lambda,double *kappa,
	double *tol,int *max_ite,int *nrowa,int *ncola,double *beta)
{
	    int i,j;
		//1.convert R type to C++ type:
	    double **x=new double*[*nrowa];
		for(i=0;i<*nrowa;i++)
			x[i]=new double[*ncola];
	    for (i=0; i<*nrowa; i++){
	        for (j=0; j<*ncola; j++){
	            x[i][j] = X[j * (*nrowa) + i];
	        }
	    }
		double *y=new double[*nrowa];
		for(i=0;i<*nrowa;i++)
			y[i]=Y[i];

		//2.initial something:
		int n=(*nrowa),p=(*ncola);//n <- dim(x)[1],p <- dim(x)[2]
		double dif_b=1;int ite=1;//dif_b <- 1,ite <- 1
		double *beta_old=new double[p];
		for(i=0;i<p;i++)
			beta_old[i]=0;//beta_old <- rep(0,p)
		double *beta_new=new double[p];
		for(i=0;i<p;i++)
			beta_new[i]=0;//beta_new <- rep(0,p)
		double *r=new double [n];
        double temp;
		for(i=0;i<n;i++)
		{
			temp=0;
		    for(j=0;j<p;j++)
				temp+=x[i][j]*beta_old[j];
			r[i]=y[i]-temp;
		}//r <- y - x%*%beta_old

		//3.coordinate descent:
		double vj,z,sum1,sum2;
		while( (dif_b > *tol) && (ite< *max_ite) )
		{
			for(j=0;j<p;j++)
			{   
				temp=0;
				for(i=0;i<n;i++)
					temp+=(x[i][j]*x[i][j]);
				vj=temp/n;//  vj <- mean(x[,j]^2)

				z=0;
				for(i=0;i<n;i++)
					z+=x[i][j]*r[i];
				z/=n;
				z+=vj*beta_old[j];//z <- mean(x[,j]*r) + vj*beta_old[j]

				beta_new[j]=soft(z,(*lambda)/(*kappa)*exp(-fabs(beta_old[j])/(*kappa)))/vj;
			    // beta_new[j] <- soft(z, lambda/kappa*exp(-abs(beta_old[j])/kappa)) / vj

				for(i=0;i<n;i++)
					r[i]-=(beta_new[j]-beta_old[j])*x[i][j];
				//r <- r - (beta_new[j] - beta_old[j])*x[,j]
			}

			sum1=0;sum2=0;
            for(i=0;i<p;i++){
				sum1+=(beta_new[i]-beta_old[i])*(beta_new[i]-beta_old[i]);
				sum2+=(beta_old[i])*(beta_old[i]);
			}
			sum2+=0.01;
			dif_b=sqrt(sum1/sum2);
			for(i=0;i<p;i++) beta_old[i]=beta_new[i];//beta_old<-beta_new:
			ite+=1;

		}

		//4.intercept estimate:
		for(i=1;i<p+1;i++)
			beta[i]=beta_new[i-1];//beta <- beta_new

		double yBar=0;
		for(i=0;i<n;i++)
			yBar+=y[i];
		yBar/=n;//  yBar <- mean(y)

		double *xBar=new double[p],s;
		for(j=0;j<p;j++)
		{
			s=0;
			for(i=0;i<n;i++)
				s+=x[i][j];
			s/=n;
			xBar[j]=s;
		}//xBar <- apply(x,2,mean)

		double beta0=0;
		for(i=0;i<p;i++)
			beta0+=xBar[i]*beta[i+1];
		beta0=yBar-beta0;
		beta[0]=beta0;

		//5.End:last to do:delete the space,it is very important!
		for(i=0;i<n;i++)
			delete []x[i];
		delete []x;delete []y;delete []beta_old;delete []beta_new;
		delete []r;delete []xBar;


}

void binlep(double *X,double *Y,double *lambda,double *kappa,
	double *tol,int *max_ite,int *nrowa,int *ncola,double *beta)
{
	    int i,j,l,jj;
		//1.convert R type to C++ type:
	    double **x2=new double*[*nrowa];
		for(i=0;i<*nrowa;i++)
			x2[i]=new double[*ncola];
	    for (i=0; i<*nrowa; i++){
	        for (j=0; j<*ncola; j++){
	            x2[i][j] = X[j * (*nrowa) + i];
	        }
	    }

		//x<-cbind(1,x):
		double **x=new double *[*nrowa];
		for(i=0;i<*nrowa;i++)
			x[i]=new double[*ncola+1];
		for(i=0;i<*nrowa;i++)
			x[i][0]=1;

		for(i=0;i<*nrowa;i++){
			for(j=1;j<*ncola+1;j++){
				x[i][j]=x2[i][j-1];
			}
		}

		double *y=new double[*nrowa];
		for(i=0;i<*nrowa;i++)
			y[i]=Y[i];

		//2.initial something:
		int n=(*nrowa),p=(*ncola+1);//n <- dim(x)[1],p <- dim(x)[2]
		double dif_b=1;int ite=1;//dif_b <- 1,ite <- 1
		double *beta_old=new double[p];
		for(i=0;i<p;i++)
			beta_old[i]=0;//beta_old <- rep(0,p)
		double *beta_new=new double[p];
		for(i=0;i<p;i++)
			beta_new[i]=0;//beta_new <- rep(0,p)
		double *pr=new double[n];
		for(i=0;i<n;i++)
			pr[i]=0.5;// pr <- 1/2
		double *r=new double [n];
		for(i=0;i<n;i++)
			r[i]=(y[i]-pr[i])/pr[i];

		//3.coordinate descent:
		double vj,temp,z,sum1,sum2;
		double *x_times_beta_new=new double[n];
		while( (dif_b>*tol) && (ite<*max_ite) ){
			//3.1:first deal with the intercept term
			temp=0;
			for(i=0;i<n;i++)
				temp+=(x[i][0]*x[i][0]*pr[i]);
			vj=temp/n;//vj <- sum(pr*x[,1]^2)/n

			z=0;
			for(i=0;i<n;i++)
				z+=x[i][0]*r[i]*pr[i];
			z/=n;
			z+=vj*beta_old[0];//z<-sum(x[,1]*pr*r)/n + vj*beta_old[1]

			beta_new[0]=z/vj;// beta_new[1] <- z/vj

			for(i=0;i<n;i++){
				temp=0;
				for(j=0;j<p;j++)
					temp+=x[i][j]*beta_new[j];
				x_times_beta_new[i]=temp;
			}
			for(i=0;i<n;i++){
				temp=exp(x_times_beta_new[i]);
				pr[i]=temp/(temp+1);
			}// pr <- exp(x%*%beta_new)/(1+ exp(x%*%beta_new))

			for(i=0;i<n;i++)
			    r[i]=(y[i]-pr[i])/pr[i];// r <- (y - pr)/pr

			//3.2:then with the penalized column.
			for(j=1;j<p;j++){
				temp=0;
			    for(i=0;i<n;i++)
				   temp+=(x[i][j]*x[i][j]*pr[i]);
			    vj=temp/n;//vj <- sum(pr*x[,j]^2)/n
			
			    z=0;
			    for(i=0;i<n;i++)
				    z+=x[i][j]*r[i]*pr[i];
			    z/=n;
			    z+=vj*beta_old[j];//z <- sum(x[,j]*pr*r)/n + vj*beta_old[j]

			    beta_new[j]=soft(z,(*lambda)/(*kappa)*exp(-fabs(beta_old[j])/(*kappa)))/vj;
			    // beta_new[j] <- soft(z, lambda/kappa*exp(-abs(beta_old[j])/kappa)) / vj
			
				for(i=0;i<n;i++){
				temp=0;
				for(jj=0;jj<p;jj++)
					temp+=x[i][jj]*beta_new[jj];
				x_times_beta_new[i]=temp;
			   }
			    for(l=0;l<n;l++){
				    temp=exp(x_times_beta_new[l]);
				     pr[l]=temp/(temp+1);
				}// pr <- exp(x%*%beta_new)/(1+ exp(x%*%beta_new))

				for(i=0;i<n;i++)
			        r[i]=(y[i]-pr[i])/pr[i];// r <- (y - pr)/pr
			}

			//3.3:
			sum1=0;sum2=0;
            for(i=0;i<p;i++){
				sum1+=(beta_new[i]-beta_old[i])*(beta_new[i]-beta_old[i]);
				sum2+=(beta_old[i])*(beta_old[i]);
			}
			sum2+=0.01;
			dif_b=sqrt(sum1/sum2);
			for(i=0;i<p;i++) beta_old[i]=beta_new[i];//beta_old<-beta_new:
			ite+=1;

		}

		//4. beta <- beta_new
		for(i=0;i<p;i++)
			beta[i]=beta_new[i];

		//5.End:last to do:delete the space,it is very important!
		for(i=0;i<n;i++){
			delete []x2[i]; delete []x[i];}
		delete []x2;delete []y;delete []beta_old;delete []beta_new;
		delete []x;delete []x_times_beta_new;
		delete []r;delete []pr;
		
}

}
