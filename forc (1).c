#include<stdio.h>
#include<math.h>

///ELECTROSTATIC SHORT RANGE
void force_es(int N , int i , int j , double x[N] , double y[N] , double z[N] , double ch[N] , double rad[N] , double coup , double cutwca , double rcerf , double ka , double kapi , double Lx , double Ly , double Lz , double *pFx, double *pFy, double *pFz, double *penergy)
{                                                                                                                                                                           

        double Fx=0.,Fy=0.,Fz=0.,Fzim=0.,energy=*penergy,dx,dy,dz,c1,dist,c5,rr,ce;
                  
            dx=x[i]-x[j];
            dx=dx-round(dx/Lx)*Lx;
            
            dy=y[i]-y[j];
            dy=dy-round(dy/Ly)*Ly;

            dz=z[i]-z[j];
//          dz=dz-round(dz/Lz)*Lz;
            
            dist=sqrt(dx*dx+dy*dy+dz*dz);
            
            //Lennard-Jones
            rr=rad[i]+rad[j];
            if (dist<cutwca*rr){
                c5=2.*( - 6./pow(dist/rr,7) + 12./pow(dist/rr,13));
		ce=( - 1./pow(dist/rr,6) + 1./pow(dist/rr,12));
                if (dist<=0.9*rr) {
			c5=2.*( - 6./pow(0.9,7) + 12./pow(0.9,13));
			ce=(-1./pow(0.9,6) + 1./pow(0.9,12));
		}
                Fx = Fx + c5*dx/dist ; Fy = Fy + c5*dy/dist ; Fz = Fz + c5*dz/dist ; energy = energy + ce;
            }

            //SHORT RANGE ELECTRIC FORCE
            
            if (dist<rcerf&&dist>0.9*rr){
                c1=ch[i]*ch[j]*coup*(erfc(ka*dist)/dist/dist/dist + kapi*exp(-ka*ka*dist*dist)/dist/dist) ;
                ce=1./2.*ch[i]*ch[j]*coup*(erfc(ka*dist)/dist)-kapi*coup/2.*ch[i]*ch[i] ;
		Fx = Fx + c1*dx ; Fy = Fy + c1*dy ; Fz = Fz + c1*dz ; energy = energy + ce ;
	    }
            
            *pFx=Fx ; *pFy=Fy ; *pFz=Fz; *penergy=energy;
            
}

///ELECTROSTATIC LONG RANGE
void force_el(int num_k , double x1 , double y1 , double z1 , double ch1 , double pre_k[num_k] , double A[num_k] , double B[num_k] , double slab , double kvec[3][num_k] , double *pFx, double *pFy, double *pFz, double *penergy)
{                                                                                                                                                                           

        double Fx=0.,Fy=0.,Fz=0.,energy=*penergy,aux,ce;
        int k;

        for (k=0;k<num_k;k++){
            aux=ch1*pre_k[k]*( A[k]*sin(kvec[0][k]*x1+kvec[1][k]*y1+kvec[2][k]*z1)+
            B[k]*cos(kvec[0][k]*x1+kvec[1][k]*y1+kvec[2][k]*z1) ) ;
	    ce=pre_k[k]/4.*( A[k]*A[k] + B[k]*B[k] ) ;
            Fx = Fx + kvec[0][k]*aux ; Fy = Fy + kvec[1][k]*aux ; Fz = Fz + kvec[2][k]*aux ; energy = energy + ce;
        }

        *pFx=Fx ; *pFy=Fy ; *pFz=Fz; *penergy=energy;

}

///POLARIZATION
void force_pol(int num_m , double x1 , double y1 , double z1 , double ch1 , double pre_m[5][num_m] , double mmod[num_m] , double f[20][num_m] , double slab , double mvec[2][num_m] , double kappa[num_m] , double qx , double qy , double hamp , double *pFx, double *pFy, double *pFz)
{                                                                                                                                                                           

        double Fx=0.,Fy=0.,Fz=0.,aux_x,aux_y,aux_z,e1,e2,ant1,ant2,kap1,kap2,q1,q2,eh1,eh2,ekap1,ekap2;
        double  df0[3], df1[3], df2[3], df3[3], df4[3], df5[3], df6[3], df7[3], df8[3], df9[3];
        double df10[3],df11[3],df12[3],df13[3],df14[3],df15[3],df16[3],df17[3],df18[3],df19[3]; 
        int m;

        for (m=0;m<num_m;m++){


        e1=exp(mmod[m]*(z1));
        e2=exp(-mmod[m]*(z1));
        
        kap1=exp(kappa[m]*z1);
        kap2=exp(-kappa[m]*z1);
           
        ant1=ch1*cos(mvec[0][m]*x1+mvec[1][m]*y1);
        ant2=ch1*sin(mvec[0][m]*x1+mvec[1][m]*y1);
        
        q1=hamp*cos(qx*x1+qy*y1);
        q2=hamp*sin(qx*x1+qy*y1);
        
        eh1=exp(2*mmod[m]*q1) ;
        eh2=exp(-2*mmod[m]*q1) ;
        ekap1=exp((mmod[m]+kappa[m])*q1) ;
        ekap2=exp(-(mmod[m]+kappa[m])*q1) ;
            
        df0[0]= -mvec[0][m]*ant2*e1 ;
        df0[1]= -mvec[1][m]*ant2*e1 ;
        df0[2]= mmod[m]*ant1*e1    ;
        
        df1[0]= mvec[0][m]*ant1*e1 ;
        df1[1]= mvec[1][m]*ant1*e1 ;
        df1[2]= mmod[m]*ant2*e1    ;
        
        df2[0]= -mvec[0][m]*ant2*e2 ;
        df2[1]= -mvec[1][m]*ant2*e2 ;
        df2[2]= -mmod[m]*ant1*e2   ;
        
        df3[0]= mvec[0][m]*ant1*e2 ;
        df3[1]= mvec[1][m]*ant1*e2 ;
        df3[2]= -mmod[m]*ant2*e2 ;
        
        //f4 = ant1*kap1*q1
        df4[0]= kap1*(q1*-mvec[0][m]*ant2+ant1*-qx*q2) ;
        df4[1]= kap1*(q1*-mvec[1][m]*ant2+ant1*-qy*q2) ;
        df4[2]= kappa[m]*ant1*kap1*q1 ;
        
        //f5 = ant2*kap1*q1
        df5[0]= kap1*(q1*mvec[0][m]*ant1+ant2*-qx*q2) ;
        df5[1]= kap1*(q1*mvec[1][m]*ant1+ant2*-qy*q2) ;
        df5[2]= kappa[m]*ant2*kap1*q1 ;
        
        //f6 = ant1*kap2*q1
        df6[0]= kap2*(q1*-mvec[0][m]*ant2+ant1*-qx*q2) ;
        df6[1]= kap2*(q1*-mvec[1][m]*ant2+ant1*-qy*q2) ;
        df6[2]= -kappa[m]*ant1*kap2*q1 ;
        
        //f7 = ant2*kap2*q1
        df7[0]= kap2*(q1*mvec[0][m]*ant1+ant2*-qx*q2) ;
        df7[1]= kap2*(q1*mvec[1][m]*ant1+ant2*-qy*q2) ;
        df7[2]= -kappa[m]*ant2*kap2*q1 ;
        
        //f8
        df8[0]= e1*(-mvec[0][m]*ant2+2*mmod[m]*q2*qx*ant1)*eh2 ;
        df8[1]= e1*(-mvec[1][m]*ant2+2*mmod[m]*q2*qy*ant1)*eh2 ;
        df8[2]= mmod[m]*e1*ant1*eh2 ;
        
        //f9
        df9[0]= e1*(mvec[0][m]*ant1+2*mmod[m]*q2*qx*ant2)*eh2 ;
        df9[1]= e1*(mvec[1][m]*ant1+2*mmod[m]*q2*qy*ant2)*eh2 ;
        df9[2]= mmod[m]*e1*ant2*eh2 ;
        
        //f10
        df10[0]=e2*(-mvec[0][m]*ant2-2*mmod[m]*q2*qx*ant1)*eh1 ;
        df10[1]=e2*(-mvec[1][m]*ant2-2*mmod[m]*q2*qy*ant1)*eh1 ;
        df10[2]=-mmod[m]*e2*ant1*eh1 ;
        
        //f11
        df11[0]=e2*(mvec[0][m]*ant1-2*mmod[m]*q2*qx*ant2)*eh1 ;
        df11[1]=e2*(mvec[1][m]*ant1-2*mmod[m]*q2*qy*ant2)*eh1 ;
        df11[2]=-mmod[m]*e2*ant2*eh1 ;
        
        //f12
        df12[0]= e1*(-mvec[0][m]*ant2*q1+2*mmod[m]*qx*q2*ant1*q1-ant1*qx*q2)*eh2 ;
        df12[1]= e1*(-mvec[1][m]*ant2*q1+2*mmod[m]*qy*q2*ant1*q1-ant1*qy*q2)*eh2 ;
        df12[2]= mmod[m]*e1*ant1*eh2*q1 ;
        
        //f13
        df13[0]= e1*(mvec[0][m]*ant1*q1+2*mmod[m]*qx*q2*ant2*q1-ant2*qx*q2)*eh2 ;
        df13[1]= e1*(mvec[1][m]*ant1*q1+2*mmod[m]*qy*q2*ant2*q1-ant2*qy*q2)*eh2 ;
        df13[2]= mmod[m]*e1*ant2*eh2*q1 ;
        
        //f14
        df14[0]= e2*(-mvec[0][m]*ant2*q1-2*mmod[m]*qx*q2*ant1*q1-ant1*qx*q2)*eh1 ;
        df14[1]= e2*(-mvec[1][m]*ant2*q1-2*mmod[m]*qy*q2*ant1*q1-ant1*qy*q2)*eh1 ;
        df14[2]= -mmod[m]*e2*ant1*eh1*q1 ;
        
        //f15
        df15[0]= e2*(mvec[0][m]*ant1*q1-2*mmod[m]*qx*q2*ant2*q1-ant2*qx*q2)*eh1 ;
        df15[1]= e2*(mvec[1][m]*ant1*q1-2*mmod[m]*qy*q2*ant2*q1-ant2*qy*q2)*eh1 ;
        df15[2]= -mmod[m]*e2*ant2*eh1*q1 ;
        
        //f16
        df16[0]= e1*(-mvec[0][m]*ant2+(mmod[m]+kappa[m])*q2*qx*ant1)*ekap2 ;
        df16[1]= e1*(-mvec[1][m]*ant2+(mmod[m]+kappa[m])*q2*qy*ant1)*ekap2 ;
        df16[2]= mmod[m]*e1*ant1*ekap2 ;
        
        //f17
        df17[0]= e1*(mvec[0][m]*ant1+(mmod[m]+kappa[m])*q2*qx*ant2)*ekap2 ;
        df17[1]= e1*(mvec[1][m]*ant1+(mmod[m]+kappa[m])*q2*qy*ant2)*ekap2 ;
        df17[2]= mmod[m]*e1*ant2*ekap2 ;
        
        //f18
        df18[0]= e2*(-mvec[0][m]*ant2-(mmod[m]+kappa[m])*q2*qx*ant1)*ekap1 ;
        df18[1]= e2*(-mvec[1][m]*ant2-(mmod[m]+kappa[m])*q2*qy*ant1)*ekap1 ;
        df18[2]= -mmod[m]*e2*ant1*ekap1 ;
        
        //f19
        df19[0]= e2*(mvec[0][m]*ant1-(mmod[m]+kappa[m])*q2*qx*ant2)*ekap1 ;
        df19[1]= e2*(mvec[1][m]*ant1-(mmod[m]+kappa[m])*q2*qy*ant2)*ekap1 ;
        df19[2]= -mmod[m]*e2*ant2*ekap1 ;
        
        
        ////////////////////

        // find force on particle on left side
        if (z1 < 0) {
        	
        	aux_x = -0.5 * ( 
        		pre_m[0][m] * (f[8][m]*df0[0]+f[0][m]*df8[0]+f[9][m]*df1[0]+f[1][m]*df9[0]) + 
        		pre_m[1][m] * (f[2][m]*df8[0]+f[3][m]*df9[0]+f[10][m]*df0[0]+f[11][m]*df1[0]) - 
        		pre_m[2][m] * (f[16][m]*df4[0]+f[17][m]*df5[0]+f[4][m]*df16[0]+f[5][m]*df17[0]) + 
        		pre_m[3][m] * (-f[6][m]*df16[0]-f[7][m]*df17[0]+f[18][m]*df4[0]+f[19][m]*df5[0]) + 
        		pre_m[4][m] * (f[12][m]*df0[0]+f[13][m]*df1[0]+f[0][m]*df12[0]+f[1][m]*df13[0])
        		);

        	aux_y = -0.5 * ( 
        		pre_m[0][m] * (f[8][m]*df0[1]+f[0][m]*df8[1]+f[9][m]*df1[1]+f[1][m]*df9[1]) + 
        		pre_m[1][m] * (f[2][m]*df8[1]+f[3][m]*df9[1]+f[10][m]*df0[1]+f[11][m]*df1[1]) - 
        		pre_m[2][m] * (f[16][m]*df4[1]+f[17][m]*df5[1]+f[4][m]*df16[1]+f[5][m]*df17[1]) + 
        		pre_m[3][m] * (-f[6][m]*df16[1]-f[7][m]*df17[1]+f[18][m]*df4[1]+f[19][m]*df5[1]) + 
        		pre_m[4][m] * (f[12][m]*df0[1]+f[13][m]*df1[1]+f[0][m]*df12[1]+f[1][m]*df13[1])
        		);

        	aux_z = -0.5 * ( 
        		pre_m[0][m] * (f[8][m]*df0[2]+f[0][m]*df8[2]+f[9][m]*df1[2]+f[1][m]*df9[2]) + 
        		pre_m[1][m] * (f[2][m]*df8[2]+f[3][m]*df9[2]+f[10][m]*df0[2]+f[11][m]*df1[2]) - 
        		pre_m[2][m] * (f[16][m]*df4[2]+f[17][m]*df5[2]+f[4][m]*df16[2]+f[5][m]*df17[2]) + 
        		pre_m[3][m] * (-f[6][m]*df16[2]-f[7][m]*df17[2]+f[18][m]*df4[2]+f[19][m]*df5[2]) + 
        		pre_m[4][m] * (f[12][m]*df0[2]+f[13][m]*df1[2]+f[0][m]*df12[2]+f[1][m]*df13[2])
        		);
        	
        } 
        
        //find force on particle on right side
        if (z1 > 0) {

        	aux_x = -0.5 * (
        		pre_m[0][m] * (f[10][m]*df2[0]+f[11][m]*df3[0]+f[2][m]*df10[0]+f[3][m]*df11[0]) +
        		pre_m[1][m] * (f[8][m]*df2[0]+f[9][m]*df3[0]+f[0][m]*df10[0]+f[1][m]*df11[0]) +
        		pre_m[2][m] * (f[18][m]*df6[0]+f[19][m]*df7[0]+f[6][m]*df18[0]+f[7][m]*df19[0]) -
        		pre_m[3][m] * (f[16][m]*df6[0]+f[17][m]*df7[0]-f[4][m]*df18[0]-f[5][m]*df19[0]) -
        		pre_m[4][m] * (f[14][m]*df2[0]+f[15][m]*df3[0]+f[2][m]*df14[0]+f[3][m]*df15[0])
        		) ;
        		
        	aux_y = -0.5 * (
        		pre_m[0][m] * (f[10][m]*df2[1]+f[11][m]*df3[1]+f[2][m]*df10[1]+f[3][m]*df11[1]) +
        		pre_m[1][m] * (f[8][m]*df2[1]+f[9][m]*df3[1]+f[0][m]*df10[1]+f[1][m]*df11[1]) +
        		pre_m[2][m] * (f[18][m]*df6[1]+f[19][m]*df7[1]+f[6][m]*df18[1]+f[7][m]*df19[1]) -
        		pre_m[3][m] * (f[16][m]*df6[1]+f[17][m]*df7[1]-f[4][m]*df18[1]-f[5][m]*df19[1]) -
        		pre_m[4][m] * (f[14][m]*df2[1]+f[15][m]*df3[1]+f[2][m]*df14[1]+f[3][m]*df15[1])
        		) ;

        	aux_z = -0.5 * (
        		pre_m[0][m] * (f[10][m]*df2[2]+f[11][m]*df3[2]+f[2][m]*df10[2]+f[3][m]*df11[2]) +
        		pre_m[1][m] * (f[8][m]*df2[2]+f[9][m]*df3[2]+f[0][m]*df10[2]+f[1][m]*df11[2]) +
        		pre_m[2][m] * (f[18][m]*df6[2]+f[19][m]*df7[2]+f[6][m]*df18[2]+f[7][m]*df19[2]) -
        		pre_m[3][m] * (f[16][m]*df6[2]+f[17][m]*df7[2]-f[4][m]*df18[2]-f[5][m]*df19[2]) -
        		pre_m[4][m] * (f[14][m]*df2[2]+f[15][m]*df3[2]+f[2][m]*df14[2]+f[3][m]*df15[2])
        		) ;
        }

        Fx = Fx + aux_x ; Fy = Fy + aux_y ; Fz = Fz + aux_z ;

        }

        *pFx=Fx ; *pFy=Fy ; *pFz=Fz;

}

///KVEC INIT
void kvectors(double Lx , double Ly , double Lz , int *pnum_k)
{                             

int nx,ny,nz,kcont=0;
double kra,pi=acos(-1.);
double kr2=4.*pi*pi*16./Lx/Lx;

for (nx=-4;nx<=4;nx++){
    for (ny=-4;ny<=4;ny++){
        for (nz=-4*(int)(Lz/Lx);nz<=4*(int)(Lz/Lx);nz++){

            if (nx==0&ny==0&nz==0){
                goto se ;
            }

            kra=4.*pi*pi*(nx*nx/Lx/Lx + ny*ny/Ly/Ly + nz*nz/Lz/Lz) ;

            if (kra<=kr2) {
                kcont +=1;
            }

se:  continue ;

        }
    }
}

*pnum_k=kcont;

}

///MVEC INIT
void mvectors(int next, double Lx , double Ly , int *pnum_m)
{

int mx,my,mcont=0;
double pi=acos(-1.);

for (mx=-next;mx<=next;mx++){
    for (my=-next;my<=next;my++){

            if (mx==0&my==0){
                goto se2 ;
            }

                mcont +=1;
            
se2:  continue ;

    }
}

*pnum_m=mcont;

}
