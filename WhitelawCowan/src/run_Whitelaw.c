#include <R.h>
#include <time.h>
#include <stdio.h>
#include <Rmath.h>

//required for ceiling 
double ceil( double num );

void run(int *Nr,int *connecs,int *inv_strength,int *max_neigh, double *k, double *sij,int *Nt,double *timestep, double *cij,double *alpha,int *v,double *rmax, double *rmin, double *sijmin,int *istart, int *iend, int *conn_bij, int *inv_str_bij,int *max_n_bij){

// set seed 
 srand48((long)time(NULL));

GetRNGstate();

 int i;
 int Nrv= *Nr;
 int Ntv=*Nt; 
 int j;
 int m;
 int vv=*v;
 int start=*istart;
 int end=*iend;
 double kv=*k;
double tj[*Nt];
double tj2[*Nt];
 double tj3[*Nt];
double *update1=malloc((*v)*sizeof(double));
 double rmaxv=*rmax;
 double rminv=*rmin;
 double sijminv=*sijmin; 
 double a;
 double h;
 double timestepv =*timestep;
double l;
int dummy1;
double dummy2;
int u;
int w;
 int wb;
int z=*max_neigh;
 int op=*max_n_bij;
//start and end are numbers of epocs

	 for(i=start;i<=end;i=i+1){
	  

memset(update1,0,sizeof(double)*vv);
  memset(tj3,0,sizeof(double)*Ntv);
	 
	       
		// iterate through epoc 
		for(u=0;u<Nrv;u=u+1){
		 		       
			w=inv_strength[u];

 		       	for(j=0; j<Ntv; j=j+1){
   				l=0;
   				for(m=0; m<w; m=m+1){
				  
	 		   	   l=l+sij[j*Nrv+connecs[u*z+m]];
				
   				}
  				tj[j]=2*kv*l/((double)w);	
			       			}
			
// tectal cell activity 
	   		 for (j=0; j<Ntv;j=j+1){
	      			 l=0;
				
				wb=inv_str_bij[j];
			       
	      			for(m=0; m<wb;m=m+1){
		       
					l=l+tj[conn_bij[j*op+m]];
	      			}
	      			tj2[j]=l/((double)wb);

				
	    		}

	
// step update	   
	   		 for( j=0; j<Ntv;j=j+1){
				dummy1=j*Nrv;
			       
				tj3[j]=tj3[j]+tj2[j];
				
				for(m=0; m<w; m=m+1){
					
				  update1[connecs[u*z+m]+dummy1]=update1[connecs[u*z+m]+dummy1]+tj2[j]*cij[connecs[u*z+m]+dummy1]/((double)w);
			  				
				}				
			}
//add random number - we don't do this so remove it (its outside of the epoc loop so can threshold at end)
//for(j=0; j<vv;j=j+1){
//a=unif_rand();
//update1[j]=update1[j]+rminv+(rmaxv-rminv)*a;
//}

	
		}
	     
		//update weights matrix
		for(j=0;j<Ntv;j=j+1){
		  tj3[j]=tj3[j]*alpha[0];
		}

	       
		for(j=0;j<Ntv;j=j+1){
		  dummy1=j*Nrv;
		  for(m=0;m<Nrv;m=m+1){
		   
		    sij[m+dummy1]=sij[m+dummy1]+timestepv*(update1[m+dummy1]-tj3[j]);
		    //threshold
		    if(sij[m+dummy1]<sijminv){
		      sij[m+dummy1]=0;
		    }
		  
 		  }
		}
	    
		
	    
//normalise   
	    	for(j=0; j<Ntv;j=j+1){
	     		 h=0;
			dummy1=j*Nrv;
	      		for(m=0; m<Nrv; m=m+1){
				h=h+sij[dummy1+m];
			}
			dummy2=Nrv/h;
			
	      		for(m=0; m<Nrv; m=m+1){
				sij[dummy1+m]=dummy2*sij[dummy1+m];
			}
	   	}

	    	for(j=0; j<Nrv;j=j+1){
	      		h=0;
	      		for(m=0;m<Ntv;m=m+1){
				h=h+sij[j+m*Nrv];
	      		}
			dummy2=Ntv/h;
	      		for(m=0;m<Ntv;m=m+1){
				sij[j+m*Nrv]=dummy2*sij[j+m*Nrv];
	      		}
	    	}
 	
 }
	

	
PutRNGstate();

free(update1);

}
