
//compile: gcc -o LANG main.c -lm

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include"forc.c"
#include"gaussian.c"

/*
This is the main method for bent electrostatics simulations
*/

double main (int argc, char* argv[])
{
	
//argv[1] denotes the run number we would like to do

FILE *file0;

char adr[50];

sprintf(adr, "./run_%s/info", argv[1]);
file0=fopen(adr,"w");
fprintf(file0," \n");
fprintf(file0," ----------------------------------------- \n");	
fprintf(file0," This is the LANG code. \n");
fprintf(file0," --- \n");
fprintf(file0," \n");
fclose(file0);

srand( (unsigned)time(NULL) + atoi(argv[1]) ) ; //set random seed to cpu time plus the argv1 value (so parallel jobs have different random seeds)

//srand(11);//TEST

int cont=0,rri;
double pi=acos(-1.);

/// PARAMETERS /// /// /// lengths in units of DIAMETER "d"
    double sign=0.; //field in z direction
    int Neq=50000,Nst=100,Nsamp=5000,Npic=500; //these are good values to use but we must test w small vals also
    int Nc=0;//counterions
    int Np=18,Nn=36; //particle numbers
    double chp=+2.,chn=-1.;//charges
    double ew = 80.; //dielectric constant of water
    double el = 2.; //dielectric constant of membrane
    int nqx = 2; //number of slab undulations in x
    int nqy = 0; //number of slab undulations in y
    double coup=2.;
    double sig=0.;
    double Lx=10.;//sqrt((double)(Nc)/2./sig);
    double slab=14.;
    double hslab = slab/2.; //half the slab distance
    double Lc=60.;
    double Lz=Lc;
    double Ly=Lx;
    double qx = 2.*pi*nqx/Lx;
    double qy = 2.*pi*nqy/Ly;
    double q = sqrt(qx*qx+qy*qy);
    double hamp = 0.5; if (nqx==0 && nqy==0) hamp = 0; //amplitude of height fluctuations; 0 for flat case
    double fric = 0.001;
    double dt = 0.005;
    fric = fric/dt;
    double dt2 = dt*dt;
    double con = fric*dt;
    double c0 = exp(-con);
    double c1 = ((float)1-c0)/con;
    double c2 = ((float)1 - c1)/con;
    double cx = ((c0*c2)/c1);
    double cy = ((float)1/con)*((float)1-(c0/c1));
    double desv_r = dt *sqrt((float)1/con*((float)2-(float)1/con*((float)3-(float)4*c0+(c0*c0))));
    double desv_v = sqrt((float)1-(c0*c0));
    double crv = dt/con/desv_r/desv_v*((float)1-c0)*((float)1-c0);
    double crv2 = sqrt((float)1-(crv*crv));
/// /// /// /// /// /// ///

double dist,Qt,Mz;
int next=10; //number of m vectors to consider
int Ntot=Neq+Nsamp*Nst;

double rcerf=Lx/2.;
double cutwca=pow(2.,1./6.);

double ka = 5./Lx ;
double kapi = 2.*ka/sqrt(pi) ;
double XA2 = coup*2.*pi/(Lx*Ly*Lz) ;
double XA3 = coup*2.*pi/(Lx*Ly) ;
double MZC = 2.*XA2 ;

/////////////////////////////

int num_k,nx,ny,nz,kk;
kvectors(Lx , Ly , Lz , &num_k); num_k=num_k/2;
double kvec[3][num_k];
int num_m,m,mx,my,mm;
mvectors(next , Lx , Ly , &num_m) ;
double mvec[2][num_m];
double kra,kr2=4.*pi*pi*16./Lx/Lx;
int kcont=0;
for (nx=-4;nx<=4;nx++){
    for (ny=-4;ny<=4;ny++){
        for (nz=-4*(int)(Lz/Lx);nz<=4*(int)(Lz/Lx);nz++){

            if (nx==0&ny==0&nz==0){
                goto se2 ;
            }

            kra=4.*pi*pi*(nx*nx/Lx/Lx + ny*ny/Ly/Ly + nz*nz/Lz/Lz) ;

            if (kra<=kr2) {
                
                 for (kk=0;kk<kcont;kk++){
                     if (kvec[0][kk]==-2.*pi*nx/Lx&kvec[1][kk]==-2.*pi*ny/Ly&kvec[2][kk]==-2.*pi*nz/Lz) goto se2 ;
                 }                
                
                kvec[0][kcont]=2.*pi*nx/Lx ;
                kvec[1][kcont]=2.*pi*ny/Ly ;
                kvec[2][kcont]=2.*pi*nz/Lz ;
                
                kcont +=1;
                
            }

se2:  continue ;

        }
    }
}

int mcont=0;
for (mx=-next;mx<=next;mx++){
    for (my=-next;my<=next;my++){

            if (mx==0&my==0){
                goto se3 ;
            }
                mvec[0][mcont]=2.*pi*mx/Lx ;
                mvec[1][mcont]=2.*pi*my/Ly ;
                mcont +=1;

se3:  continue ;

    }
}

////////////////////////////

double V=Lx*Ly*Lc;
int N,i,j,k,s;
N=Np+Nn+Nc;
file0=fopen(adr,"a");
fprintf(file0," Num particles     = %d \n",N);
fprintf(file0," Np                = %i \n",Np);
fprintf(file0," Nn                = %i \n",Nn);
fprintf(file0," Charge_p          = %g \n",chp);
fprintf(file0," Charge_n          = %g \n",chn);
fprintf(file0," [-] Coupling      = %g \n",coup);
fprintf(file0," Lateral size      = %g diam \n",Lx);
fprintf(file0," Lz size           = %g diam \n",Lz);
fprintf(file0," Lc size           = %g diam \n",Lz);
fprintf(file0," Separation        = %g diam \n",slab);
fprintf(file0," Number of kvecs   = %d \n",num_k);
fprintf(file0," Number of mvecs   = %d \n",num_m);
fprintf(file0," Undulations in x  = %d \n",nqx);
fprintf(file0," Undulations in y  = %d \n",nqy);
fprintf(file0," Height amplitude  = %g \n",hamp);
fprintf(file0," epsilon1          = %g \n",ew);
fprintf(file0," epsilon2          = %g \n",el);
fclose(file0);
sprintf(adr, "./run_%s/dir_simul_data/param", argv[1]);
file0=fopen(adr,"w");
fprintf(file0,"%g %g %g %g %g %g %g %d %d %d %d %d \n",Lx,Ly,Lz,Lc,slab,dt,fric,N,Np,Nn,Nc,Nsamp);
fclose(file0);
double dist2,xOld,yOld,zOld;
double x[N],y[N],z[N],ch[N],rad[N],chi_x[N],chi_y[N],chi_z[N];

double pre_m[5][num_m],mmod[num_m],kappa[num_m],f[20][num_m];
double pre_k[num_k],A[num_k],B[num_k],dx,dy,dz,ran;
double vx[N],vy[N],vz[N],vaux[N],vauy[N],vauz[N];
double fxo[N],fyo[N],fzo[N],fxn[N],fyn[N],fzn[N],dts2,Fx,Fy,Fz;

dt2=dt*dt/2. ; dts2=dt/2. ;

/// PARTICLES CHARGES ///
memset(ch,0, sizeof(ch));
for (i=0;i<Np;i++){
     ch[i]=chp;rad[i]=0.5;
}
for (i=Np;i<N;i++){
     ch[i]=chn;rad[i]=0.5;
}
////////////////////////

Qt=0. ; for (i=0;i<N;i++) Qt = Qt + ch[i] ;

int contp=0,contn=0;

/// PARTICLES POSITIONS AND VELOCITIES ///

for (i=0;i<N;i++){

 vx[i]=0.;
 vy[i]=0.;
 vz[i]=0.;
    
 ov: ;
 
 x[i]=(1.-2.*(double)rand()/RAND_MAX)*Lx/2.; 
 y[i]=(1.-2.*(double)rand()/RAND_MAX)*Ly/2.; 


 //curved surface - slab plus the variation in h
 double slablocneg = -slab/2. + hamp*cos(qx*x[i]+qy*y[i]) ;
 double slablocpos = slab/2. + hamp*cos(qx*x[i]+qy*y[i])  ;

 if (i<Np/2) z[i]=-Lc/2.+rad[i]+(double)rand()/RAND_MAX*(Lc/2.+slablocneg-2.*rad[i]);
 if (i>=Np/2&&i<Np) z[i]=slablocpos+rad[i]+(double)rand()/RAND_MAX*(Lc/2.-slablocpos-2.*rad[i]);
 
 if (i>=Np&&i<Np+Nn/2) z[i]=-Lc/2.+rad[i]+(double)rand()/RAND_MAX*(Lc/2.+slablocneg-2.*rad[i]);
 if (i>=Np+Nn/2) z[i]=slablocpos+rad[i]+(double)rand()/RAND_MAX*(Lc/2.-slablocpos-2.*rad[i]);

      /// CHECK OVERLAPS 
      for (j=0;j<i;j++)
      {
         dx=fabs(x[i]-x[j]);
         if (dx>Lx-dx) dx=Lx-dx ;
        
         dy=fabs(y[i]-y[j]);
         if (dy>Ly-dy) dy=Ly-dy ;
        
         dz=fabs(z[i]-z[j]);
        
         dist2=dx*dx+dy*dy+dz*dz;
         
         if (dist2<rad[i]+rad[j]) {
			 goto ov;
		 }
		 
      }///over
}

/////////////////////Info for Ewald calculations///

    memset(A,0., sizeof(A));
    memset(B,0., sizeof(B));
    double kq;

for (k=0;k<num_k;k++){

    for (i=0;i<N;i++){

        A[k] = A[k] + ch[i]*cos(kvec[0][k]*x[i] + kvec[1][k]*y[i] + kvec[2][k]*z[i]) ;
        B[k] = B[k] - ch[i]*sin(kvec[0][k]*x[i] + kvec[1][k]*y[i] + kvec[2][k]*z[i]) ;

    }

        kq = kvec[0][k]*kvec[0][k] + kvec[1][k]*kvec[1][k] + kvec[2][k]*kvec[2][k] ;
        pre_k[k] = 4.*XA2/kq*exp(-kq/4./ka/ka) ;

}

Mz=0. ; for (i=0;i<N;i++) Mz = Mz + ch[i]*z[i] ; //2d
// Mz=0.; //3d

//DEFINITIONS FOR POLARIZATION TERMS
for (m=0;m<num_m;m++){

    f[0][m] = 0.; f[1][m] = 0.; f[2][m] = 0.; f[3][m] = 0.; f[4][m] = 0.; f[5][m] = 0.; f[6][m] = 0.;
    f[7][m] = 0.; f[8][m] = 0.; f[9][m] = 0.; f[10][m]= 0.; f[11][m]= 0.; f[12][m]= 0.; f[13][m]= 0.;
    f[14][m]= 0.; f[15][m]= 0.; f[16][m]= 0.; f[17][m]= 0.; f[18][m]= 0.; f[19][m]= 0.; 
    
    mmod[m]=sqrt(mvec[0][m]*mvec[0][m]+mvec[1][m]*mvec[1][m]);
    
    kappa[m]=sqrt((mvec[0][m]+qx)*(mvec[0][m]+qx)+(mvec[1][m]+qy)*(mvec[1][m]+qy));
    
    for (i=0;i<N;i++){
	if (z[i] < 0) {
	    f[0][m] = f[0][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(+mmod[m]*(z[i])) ;
            f[1][m] = f[1][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(+mmod[m]*(z[i])) ;
            //only add to these if kappa != 0
            if (kappa[m] != 0) {
            f[4][m] = f[4][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(+kappa[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;         
            f[5][m] = f[5][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(+kappa[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
            f[16][m] = f[16][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-(mmod[m]+kappa[m])*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i])) ;
            f[17][m] = f[17][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-(mmod[m]+kappa[m])*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i])) ;
            }
            f[8][m] = f[8][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i])) ;
            f[9][m] = f[9][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i])) ;
            f[12][m] = f[12][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
            f[13][m] = f[13][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
            
	}
	if (z[i] > 0) {
	    f[2][m] = f[2][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-mmod[m]*(z[i])) ;
            f[3][m] = f[3][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-mmod[m]*(z[i])) ;
            //only add to these if kappa != 0
            if (kappa[m] != 0) {
            f[6][m] = f[6][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-kappa[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
            f[7][m] = f[7][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-kappa[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
            f[18][m] = f[18][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp((mmod[m]+kappa[m])*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i])) ;
            f[19][m] = f[19][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp((mmod[m]+kappa[m])*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i])) ;
            }
            f[10][m] = f[10][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i])) ;
            f[11][m] = f[11][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i])) ;
            f[14][m] = f[14][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
            f[15][m] = f[15][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
            
	}  
    }
    	//prefactors for our energies
    	/*
    	4/19/23: converted these to the new standard (based on particle on the right side of the membrane)
    	*/
        pre_m[0][m] = (2.*exp(mmod[m]*2.*hslab)*(exp(4.*mmod[m]*hslab)-1.)*(ew*ew-el*el)*pi)/(mmod[m]*ew*(exp(4.*hslab*mmod[m])*(el+ew)*(el+ew)-(el-ew)*(el-ew))) / Lx / Ly * coup * ew;
        
        pre_m[1][m] = (-2.*(exp(4.*mmod[m]*hslab)-1.)*(el-ew)*(el-ew)*pi)/(mmod[m]*ew*(exp(4.*hslab*mmod[m])*(el+ew)*(el+ew)-(el-ew)*(el-ew))) / Lx / Ly * coup * ew ;
        
        //only define these terms if kappa != 0
        if (kappa[m] != 0) {
        pre_m[2][m] = -(exp(-hslab*(mmod[m]+kappa[m]))*(el-ew)*pi*(2.*(-1.+exp(4.*hslab*mmod[m]))*(-1.+exp(4.*hslab*kappa[m]))*el*el*el*mmod[m]*kappa[m]+(-1.+exp(4.*hslab*mmod[m]))*(-1.+exp(4.*hslab*kappa[m]))*ew*ew*ew*(mmod[m]*mmod[m]+kappa[m]*kappa[m]-q*q)+2.*el*ew*ew*((-1.+exp(4.*hslab*(mmod[m]+kappa[m])))*mmod[m]*mmod[m]+(1.+exp(4.*hslab*mmod[m])+exp(4.*hslab*kappa[m])-4.*exp(2.*hslab*(mmod[m]+kappa[m]))+exp(4.*hslab*(mmod[m]+kappa[m])))*mmod[m]*kappa[m]+(-1.+exp(4.*hslab*(mmod[m]+kappa[m])))*(kappa[m]-q)*(kappa[m]+q))+el*el*ew*((1.+exp(4.*hslab*mmod[m])+exp(4.*hslab*kappa[m])-4.*exp(2.*hslab*(mmod[m]+kappa[m]))+exp(4.*hslab*(mmod[m]+kappa[m])))*mmod[m]*mmod[m]+4.*(-1.+exp(4.*hslab*(mmod[m]+kappa[m])))*mmod[m]*kappa[m]+(1.+exp(4.*hslab*mmod[m])+exp(4.*hslab*kappa[m])-4.*exp(2.*hslab*(mmod[m]+kappa[m]))+exp(4.*hslab*(mmod[m]+kappa[m])))*(kappa[m]-q)*(kappa[m]+q))))/(2.*ew*mmod[m]*kappa[m]*(2.*el*ew*cosh(2.*hslab*mmod[m])+(el*el+ew*ew)*sinh(2.*hslab*mmod[m]))*(2.*el*ew*cosh(2.*hslab*kappa[m])+(el*el+ew*ew)*sinh(2.*hslab*kappa[m])))/Lx/Ly*coup*ew ;
        
        pre_m[3][m] = -4.*exp(hslab*(mmod[m]+kappa[m]))*el*(el-ew)*pi*sinh(hslab*(mmod[m]-kappa[m]))*((-2.*el*mmod[m]*kappa[m]+ew*(mmod[m]*mmod[m]+kappa[m]*kappa[m]-q*q))*cosh(hslab*(mmod[m]+kappa[m]))+(-2.*ew*mmod[m]*kappa[m]+el*(mmod[m]*mmod[m]+kappa[m]*kappa[m]-q*q))*sinh(hslab*(mmod[m]+kappa[m])))/(mmod[m]*kappa[m]*(2.*el*ew*cosh(2.*hslab*mmod[m])+(el*el+ew*ew)*sinh(2.*hslab*mmod[m]))*(2.*el*ew*cosh(2.*hslab*kappa[m])+(el*el+ew*ew)*sinh(2.*hslab*kappa[m])))/Lx/Ly*coup*ew;
        }
        else {
            pre_m[2][m]=0; pre_m[3][m]=0;
        }
        pre_m[4][m] = (4.*exp(mmod[m]*2.*hslab)*(exp(4.*mmod[m]*hslab)-1.)*(ew*ew-el*el)*pi)/(ew*(exp(4.*hslab*mmod[m])*(el+ew)*(el+ew)-(el-ew)*(el-ew))) / Lx / Ly * coup * ew;
}


////////////////////////////////////////
    /// CALCULATE INITIAL FORCE ON PARTICLES
    memset(fxo,0., sizeof(fxo));
    memset(fyo,0., sizeof(fyo));
    memset(fzo,0., sizeof(fzo));
    ////////////////////////////////////////
    
sprintf(adr, "./run_%s/info", argv[1]);
file0=fopen(adr,"a");
fprintf(file0," \n");
fprintf(file0," Equilibration ... \n");
fprintf(file0," \n");
fclose(file0);

/// STARTING SIMULATION ///
double energy, kt;
for (s=1;s<=Ntot;s++)
{

    energy=0.;
    /// PARTICLE UPDATE ///
    for (i=0;i<N;i++)
    {
            
             for (k=0;k<num_k;k++){
                 A[k] = A[k] - ch[i]*cos(kvec[0][k]*x[i] + kvec[1][k]*y[i] + kvec[2][k]*z[i]) ;
                 B[k] = B[k] + ch[i]*sin(kvec[0][k]*x[i] + kvec[1][k]*y[i] + kvec[2][k]*z[i]) ;
             }
            
             for (m=0;m<num_m;m++){
                 if (z[i]<0) {
                     f[0][m] = f[0][m] - ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(mmod[m]*(z[i])) ;
                     f[1][m] = f[1][m] - ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(mmod[m]*(z[i])) ;
                     //only add to these if kappa != 0
            	     if (kappa[m] != 0) {
                     f[4][m] = f[4][m] - ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(kappa[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
                     f[5][m] = f[5][m] - ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(kappa[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
                     f[16][m] = f[16][m] - ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-(mmod[m]+kappa[m])*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i])) ;
            	     f[17][m] = f[17][m] - ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-(mmod[m]+kappa[m])*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i])) ;
                     }
                     f[8][m] = f[8][m] - ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i])) ;
            	     f[9][m] = f[9][m] - ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i])) ;
            	     f[12][m] = f[12][m] - ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
            	     f[13][m] = f[13][m] - ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
                 }
                 if (z[i] > 0) {
                     f[2][m] = f[2][m] - ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-mmod[m]*(z[i])) ;
                     f[3][m] = f[3][m] - ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-mmod[m]*(z[i])) ;
                     //only add to these if kappa != 0
            	     if (kappa[m] != 0) {
                     f[6][m] = f[6][m] - ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-kappa[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
                     f[7][m] = f[7][m] - ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-kappa[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
                     f[18][m] = f[18][m] - ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp((mmod[m]+kappa[m])*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i])) ;
            	     f[19][m] = f[19][m] - ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp((mmod[m]+kappa[m])*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i])) ;
                     }
                     f[10][m] = f[10][m] - ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i])) ;
            	     f[11][m] = f[11][m] - ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i])) ;
           	     f[14][m] = f[14][m] - ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
          	     f[15][m] = f[15][m] - ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
                 }  
             }           
            
             Mz = Mz - ch[i]*z[i] ;

            xOld = x[i];
            yOld = y[i];
            zOld = z[i];
            chi_x[i] = gausran(); 
            x[i] = x[i] + (c1*vx[i]*dt) + (c2*fxo[i]*dt2) + (chi_x[i]*desv_r);
            chi_y[i] = gausran();
            y[i] = y[i] + (c1*vy[i]*dt) + (c2*fyo[i]*dt2) + (chi_y[i]*desv_r);
            chi_z[i] = gausran();
            z[i] = z[i] + (c1*vz[i]*dt) + (c2*fzo[i]*dt2) + (chi_z[i]*desv_r);
            
            //curved surface - slab plus the variation in h
 	    double slablocneg = -slab/2. + hamp*cos(qx*x[i]+qy*y[i]) ;
 	    double slablocpos = slab/2. + hamp*cos(qx*x[i]+qy*y[i]) ;
            
            //bounceback
            if ((fabs(z[i]) > (Lc/2.-rad[i])) || ((z[i] < slablocpos+rad[i])&&(z[i] > slablocneg-rad[i]))){
                x[i] = xOld;
                y[i] = yOld;
                z[i] = zOld;
                
                if (i<Np){
                    vx[i] = -vx[i];
                    vy[i] = -vy[i];
                    vz[i] = -vz[i];
                }
                else{
                    vx[i] = -vx[i];
                    vy[i] = -vy[i];
                    vz[i] = -vz[i];
                }
            }

         /// checking periodicity
         if (fabs(x[i])>Lx/2.) x[i]=x[i]-Lx*fabs(x[i])/x[i] ;
         if (fabs(y[i])>Ly/2.) y[i]=y[i]-Ly*fabs(y[i])/y[i] ;

             for (k=0;k<num_k;k++){
                 A[k] = A[k] + ch[i]*cos(kvec[0][k]*x[i] + kvec[1][k]*y[i] + kvec[2][k]*z[i]) ;
                 B[k] = B[k] - ch[i]*sin(kvec[0][k]*x[i] + kvec[1][k]*y[i] + kvec[2][k]*z[i]) ;
             }

             for (m=0;m<num_m;m++){
                 if (z[i]<0) {
                     f[0][m] = f[0][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(mmod[m]*(z[i])) ;
                     f[1][m] = f[1][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(mmod[m]*(z[i])) ;
                     //only add to these if kappa != 0
            	     if (kappa[m] != 0) {
                     f[4][m] = f[4][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(kappa[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
                     f[5][m] = f[5][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(kappa[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
                     f[16][m] = f[16][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-(mmod[m]+kappa[m])*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i])) ;
            	     f[17][m] = f[17][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-(mmod[m]+kappa[m])*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i])) ;
                     }
                     f[8][m] = f[8][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i])) ;
            	     f[9][m] = f[9][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i])) ;
            	     f[12][m] = f[12][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
            	     f[13][m] = f[13][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])+mmod[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
                 }
                 if (z[i] > 0) {
                     f[2][m] = f[2][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-mmod[m]*(z[i])) ;
                     f[3][m] = f[3][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-mmod[m]*(z[i])) ;
                     //only add to these if kappa != 0
            	     if (kappa[m] != 0) {
                     f[6][m] = f[6][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-kappa[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
                     f[7][m] = f[7][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(-kappa[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
                     f[18][m] = f[18][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp((mmod[m]+kappa[m])*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i])) ;
            	     f[19][m] = f[19][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp((mmod[m]+kappa[m])*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i])) ;
                     }
                     f[10][m] = f[10][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i])) ;
            	     f[11][m] = f[11][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i])) ;
           	     f[14][m] = f[14][m] + ch[i]*cos(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
          	     f[15][m] = f[15][m] + ch[i]*sin(mvec[0][m]*x[i]+mvec[1][m]*y[i])*exp(2*mmod[m]*hamp*cos(qx*x[i]+qy*y[i])-mmod[m]*(z[i]))*hamp*cos(qx*x[i]+qy*y[i]) ;
                 }  
             }                         
		             
             Mz = Mz + ch[i]*z[i] ;

    }

    /// CALCULATE NEW FORCE USING xyz[i] AND NEW VELOCITIES VAUX ///
    memset(fxn, 0., sizeof(fxn));
    memset(fyn, 0., sizeof(fyn));
    memset(fzn, 0., sizeof(fzn));
    for (i=0;i<N;i++)
    {

        ///short range electrostatics
        for (j=i+1;j<N;j++){

            force_es(N , i , j , x , y , z , ch , rad ,coup , cutwca , rcerf , ka , kapi , Lx , Ly , Lz , &Fx, &Fy, &Fz, &energy);

            fxn[i]=fxn[i]+Fx ; fyn[i]=fyn[i]+Fy ; fzn[i]=fzn[i]+Fz ;

            fxn[j]=fxn[j]-Fx ; fyn[j]=fyn[j]-Fy ; fzn[j]=fzn[j]-Fz ;

        }
        ////////////
                
        ///long range electrostatics
        force_el(num_k , x[i] , y[i] , z[i] , ch[i] , pre_k , A , B , slab , kvec , &Fx, &Fy, &Fz, &energy);
        
        fxn[i]=fxn[i]+Fx ; fyn[i]=fyn[i]+Fy ; fzn[i]=fzn[i]+Fz ;
        
        fzn[i]=fzn[i] - ch[i]*MZC*( Mz - Qt*z[i] ) ; ///CORRECTION Mz Gz
        /////////////

        ///POLARIZATION
        force_pol(num_m , x[i] , y[i] , z[i] , ch[i] , pre_m , mmod , f , slab , mvec , kappa , qx , qy , hamp , &Fx, &Fy, &Fz); 

        fxn[i]=fxn[i]+Fx ; fyn[i]=fyn[i]+Fy ; fzn[i]=fzn[i]+Fz ;
        ///////////////
                     
        fzn[i]=fzn[i] + ch[i]*sign ;

    }
	
	//ADD POLARIZATION ENERGY

	//add energy from the Mz correction
    energy = energy + XA2*Mz*Mz;

//////////////////////////
    for(i=0; i<N; i++){
            
            vx[i] = c0*vx[i] + (c1-c2)*fxo[i]*dt + c2*fxn[i]*dt + desv_v*(crv*chi_x[i] + crv2*gausran());
            vy[i] = c0*vy[i] + (c1-c2)*fyo[i]*dt + c2*fyn[i]*dt + desv_v*(crv*chi_y[i] + crv2*gausran());
            vz[i] = c0*vz[i] + (c1-c2)*fzo[i]*dt + c2*fzn[i]*dt + desv_v*(crv*chi_z[i] + crv2*gausran());

            //prepara pra proxima rodada:
            fxo[i]=fxn[i] ; fyo[i]=fyn[i] ; fzo[i]=fzn[i] ;
            
    }

    /// AVERAGES ///
     
    if (s>Neq&(s-Neq)%Nst==0){
     
        cont=cont+1;

        /// SAVING SAMPLES //////////
          sprintf(adr, "./run_%s/dir_simul_data/sample_%d.xyz", argv[1], cont);
          file0=fopen(adr,"w");
          for (i=0;i<N;i++)
          {
              fprintf(file0,"%g %g %g %g %g %g %g %g %g %g\n",x[i],y[i],z[i],vx[i],vy[i],vz[i],fxn[i],fyn[i],fzn[i],ch[i]);
          }
          fclose(file0);
          
        ////////////////////////////  


    }
    
    if (s==Neq) {
    	sprintf(adr, "./run_%s/info", argv[1]);
        file0=fopen(adr,"a");
        fprintf(file0," Equilibrated. \n");
        fprintf(file0," \n");
        fprintf(file0," Running... \n");
        fprintf(file0," \n");
        fclose(file0);
    }

}
////////////////////////////

sprintf(adr, "./run_%s/info", argv[1]);
file0=fopen(adr,"a");
fprintf(file0," End. \n");
fprintf(file0," \n");
fprintf(file0," ----------------------------------------- \n");
fprintf(file0," \n");
fclose(file0);

}
