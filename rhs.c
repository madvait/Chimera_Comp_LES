#include "variables.h"

extern PetscInt wavy,isothermal;
extern PetscInt block_number, NumberOfBodies;
extern PetscInt immersed, invicid;
extern PetscInt TwoD;
extern PetscInt moveframe,rotateframe;
extern PetscInt  les, clark, rans,central,calc_sigma;
extern PetscInt  channel;
extern PetscReal CMx_c, CMy_c, CMz_c;
//extern PetscReal Mach;
PetscReal ks;
PetscReal CI= 0.06;   // for les

Cmpnts MINUS(Cmpnts v1,Cmpnts v2);
Cmpnts AVERAGE4(Cmpnts v1,Cmpnts v2, Cmpnts v3, Cmpnts v5);
Cmpnts Cross(Cmpnts v1,Cmpnts v2);

/* Initial conditions:
 use option -init1 to set it in the control.dat

1: 
2: 


*/


PetscErrorCode SetInitialCondition(UserCtx *user)
{

  DM		cda = user->cda, da=user->da, fda=user->fda;
  DMDALocalInfo	info=user->info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  //  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k;

  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  CompVars	***q;
  PetscReal     u=1153, rho=0.00237, et=1719557; //Mach number M=1.5
  PetscReal     ***p, Twall=1, pi=3.14159265;;
  PetscInt      InitialGuessOne=0;

  PetscOptionsGetInt(NULL, NULL, "-init1", &InitialGuessOne, NULL);
  PetscOptionsGetReal(NULL, NULL,"-Twall", &Twall, NULL);

  VecSet(user->Q, 0.0);

  DMDAVecGetArray(cda, user->Q, &q);
  //  user->Ma=Mach;
	 
  Vec           Coor;
  Cmpnts        ***coor;
  DMGetCoordinatesLocal(da, &Coor);
  DMDAVecGetArray(fda, Coor, &coor); 

  for (k=zs; k<ze; k++){
    for (j=ys; j<ye; j++){
      for (i=xs; i<xe; i++){      
	if (InitialGuessOne==1) {
	  // Lax initial condition i-dire
	  if (i==1 && j==1 && k==1)
	    PetscPrintf(PETSC_COMM_WORLD, "!!init 1 Lax i dire\n");
	  if (i<mx/2) {
	    q[k][j][i].rho=0.445;
	    q[k][j][i].rhoU=0.311;
	    q[k][j][i].rhoE=8.928;
	  } else {
	    q[k][j][i].rho=0.5;
	    q[k][j][i].rhoU=0.;
	    q[k][j][i].rhoE=1.4275;
	  }
	} else if (InitialGuessOne==3) {
	  // Lax initial condition j-dire
	  if (i==1 && j==1 && k==1)
	    PetscPrintf(PETSC_COMM_WORLD, "!!init 3 Lax j dire\n");
	  if (j<my/2) {
	    q[k][j][i].rho=0.445;
	    q[k][j][i].rhoV=0.311;
	    q[k][j][i].rhoE=8.928;
	  } else {
	    q[k][j][i].rho=0.5;
	    q[k][j][i].rhoV=0.;
	    q[k][j][i].rhoE=1.4275;
	  }
	} else if (InitialGuessOne==4) {
	  // Lax initial condition k-dire
	  if (i==1 && j==1 && k==1)
	    PetscPrintf(PETSC_COMM_WORLD, "!!init 4 Lax k dire\n");
	  if (k<mz/2) {
	    q[k][j][i].rho=0.445;
	    q[k][j][i].rhoW=0.311;
	    q[k][j][i].rhoE=8.928;
	  } else {
	    q[k][j][i].rho=0.5;
	    q[k][j][i].rhoW=0.;
	    q[k][j][i].rhoE=1.4275;
	  }
	} else if (InitialGuessOne==2) {
	  /* initialize the flow-field with a moving shock, 
	     Rankine–Hugoniot relations 
	     L,R: left,right of shock conditions
	     M: Shock Mach number
	     gamma: cp/cv specific heats
	     R: gas constant
	     P,T: pressure, Temp
	     zeta_s=M^2+kappa(M^2-1) pressure jump across the shock aka shock strenght
	     kappa=(gamma-1)/(gamma+1)

	     From known values of M, PR, TR using the below set of equations
	     find L values
	  */	  
	  PetscReal M=user->Ma, gamma=user->gamma, PR=1., rhoR=1.; //PR=0.5 MPa= 50,000 Pa
	  PetscReal kappa, zeta_s, PL, cR, uL, rhoL;
	  PetscInt  k_s=mz*0.25;
	  PetscOptionsGetInt(NULL, NULL, "-shock_location_k", &k_s, NULL);
	  PetscOptionsGetReal(NULL, NULL, "-shock_rhoR", &rhoR, NULL);

	  kappa=(gamma-1)/(gamma+1);
	  zeta_s=M*M+kappa*(M*M-1);
	  //	  rhoR=PR/(Rair*TR);
	  cR=sqrt(gamma*PR/rhoR);
	  
	  PL=zeta_s*PR;
	  rhoL=(1+zeta_s/kappa)/(1.0/kappa+zeta_s)*rhoR;
	  uL=1/gamma*(zeta_s-1)/M*cR;
	  // non-dim
	  /* PR=PR/(rhoL*uL*uL); */
	  /* PL=PL/(rhoL*uL*uL); */
	  /* rhoR=rhoR/rhoL; */
	  /* rhoL=1; */
	  /* uL=1; */
	  
	  if (i==1 && j==1 && k==1) {
	  PetscPrintf(PETSC_COMM_WORLD, "!!init 2 gamma %le kappa %le zata_s %le\n", gamma, kappa, zeta_s);
	  PetscPrintf(PETSC_COMM_WORLD, "!!init 2 uL %le rhoL %le PL %le\n", uL, rhoL, PL);
	  PetscPrintf(PETSC_COMM_WORLD, "!!init 2 cR %le rhoR %le PR %le\n", cR, rhoR, PR);
	  }
	  if (k<k_s) {	
	    q[k][j][i].rho=rhoL;
	    q[k][j][i].rhoW=rhoL*uL;
	    q[k][j][i].rhoE=PL/(gamma-1)+0.5*rhoL*uL*uL;
	  }  else {
	   // q[k][j][i].rho=rhoR;
	    q[k][j][i].rho=rhoR; 
	    q[k][j][i].rhoW=0.0;
	    q[k][j][i].rhoE=PR/(gamma-1);
	  }	
	
	} else if (InitialGuessOne==5) {
	   /* supersonic flow M=3.5  
	    c=1/M => gamma*p/rho = 1/M^2 => p=1/(gamma* M^2*) 
	    or
	    u = M; c=1. & rho=gamma ==> p=1.
	    
	   */
	  PetscReal M=5.0, gamma=1.4;
	  M=user->Ma;
	  q[k][j][i].rho=1.; //gamma
	  //if (isothermal) q[k][j][i].rho=Twall;
	  q[k][j][i].rhoW=1.0+0.1*(rand()%(1000)-500.)/500.0; //gamma*M; 
	  q[k][j][i].rhoU=0.0;	   
	  q[k][j][i].rhoV=0.0;
	  /* if (isothermal) */
	  /*   q[k][j][i].rhoE=0.5+Twall/(M*M*gamma*(gamma-1)); */
	  /* else */
	  q[k][j][i].rhoE=0.5+1./(M*M*gamma*(gamma-1));	//*q[k][j][i].rhoW*q[k][j][i].rhoW
	  if (i==0 && j==0 && k==0)  
	    PetscPrintf(PETSC_COMM_WORLD, "!!init 5 M %f  rho %f rhoW 1  rhoE %f isothermal %d Tw %f \n", M, q[0][0][0].rho, q[0][0][0].rhoE, isothermal, Twall);
	} else if (InitialGuessOne==55) {
	  /* supersonic flow M=3.5  
	    c=1/M => gamma*p/rho = 1/M^2 => p=1/(gamma* M^2*) 
	    Compared to init 5, it has zero velocity, i.e., initially stagnant fluid
	  */
	  PetscReal M=3.5, gamma=1.4;
	  M=user->Ma;
	  q[k][j][i].rho=1.;
	  if (isothermal) q[k][j][i].rho=Twall;
	  q[k][j][i].rhoE=1./(M*M*gamma*(gamma-1));
	  if (i==0 && j==0 && k==0)  
	    PetscPrintf(PETSC_COMM_WORLD, "!!init 55 M %f  rho 1 rhoW 0 rhoE %f isothermal %d Tw %f \n", user->Ma,q[k][j][i].rhoE, isothermal, Twall);
	} else if (InitialGuessOne==6) { 
	   /* Inviscid vortex 
	      u= 1 + u_v/u_inf exp(1-(r/R)^2)/2 r/R 
	      where r=(z-zo)k+(y-yo)j
	      xo,yo initial position of the center of the vortex
	      R radius of the vortex
	      u_v/u_inf velocity of the vortex
	      rho/rho_inf=[1- (r-1)/2 * (u_v/c_inf)^2 exp(1-(r/R)^2)]^(1/(gamma-1))
	      p/p_inf=[1- (r-1)/2 * (u_v/c_inf)^2 exp(1-(r/R)^2)]^(gamma/(gamma-1))
	      set p_inf=rho_inf=u_inf=1 here
	      c_inf=sqrt(gamma)

	      initialize in y and z directions, x is symmetry
	    */
	  PetscReal gamma=1.4, p;
	  PetscReal c_inf2=gamma, u_inf=1., u_v=0.8*u_inf;
	  PetscReal y, z, r2;
	  z=coor[k][j][i].z-CMz_c;
	  y=coor[k][j][i].y-CMy_c;
	  r2=z*z+y*y;
	  q[k][j][i].rho=pow(1. - (gamma-1)/2. * (u_v*u_v/c_inf2) * exp(1-r2), 1./(gamma-1));
	  p             =pow(1. - (gamma-1)/2. * (u_v*u_v/c_inf2) * exp(1-r2), gamma/(gamma-1));
	  q[k][j][i].rhoV=q[k][j][i].rho * (        u_v* exp((1-r2)/2.) * z);
	  q[k][j][i].rhoW=q[k][j][i].rho * (u_inf - u_v* exp((1-r2)/2.) * y);
	  q[k][j][i].rhoE=p/(gamma-1)+0.5*(q[k][j][i].rhoV*q[k][j][i].rhoV+
					   q[k][j][i].rhoW*q[k][j][i].rhoW)/q[k][j][i].rho;
	  if (i==1 && j==1 && k==0) {
	    PetscPrintf(PETSC_COMM_WORLD, "!!init 2 gamma %le u_inf %le E %le\n", gamma, u_inf, q[k][j][i].rhoE);
	  }
	} else if (InitialGuessOne==7) {
	  /* initialize the flow-field with a moving shock, 
	     Rankine–Hugoniot relations 
	     L,R: left,right of shock conditions
	     M: Sock Mach number
	     gamma: cp/cv specific heats
	     R: gas constant
	     P,T: pressure, Temp
	     zeta_s=M^2+kappa(M^2-1) pressure jump across the shock aka shock strenght
	     kappa=(gamma-1)/(gamma+1)

	     From known values of M, PL, uL using the below set of equations
	     find R values

	     Compte to InitialGuess 2, the L conditions are known.                                                                            
	  */	   
	  PetscReal M=user->Ma, gamma=user->gamma, PL=1., rhoL=1.; //PR=0.5 MPa= 50,000 Pa
	  PetscReal kappa, zeta_s, PR, cR, uL=1, rhoR;
	  PetscInt  k_s=mz/4;
	  PetscOptionsGetInt(NULL, NULL, "-shock_location_k", &k_s, NULL);
	  PetscOptionsGetReal(NULL, NULL, "-shock_rhoR", &rhoR, NULL);
	  
	  kappa=(gamma-1)/(gamma+1);
	  zeta_s=M*M+kappa*(M*M-1);
	  PR=PL/zeta_s;
	  rhoR=rhoL/(1+zeta_s/kappa)*(1.0/kappa+zeta_s);
	  //	  rhoR=PR/(Rair*TR);
	  cR=sqrt(gamma*PR/rhoR);
	  
	  //	  PL=zeta_s*PR;
	  //	  rhoL=(1+zeta_s/kappa)/(1.0/kappa+zeta_s)*rhoR;
	  uL=1/gamma*(zeta_s-1)/M*cR;
	  // non-dim
	  /* PR=PR/(rhoL*uL*uL); */
	  /* PL=PL/(rhoL*uL*uL); */
	  /* rhoR=rhoR/rhoL; */
	  /* rhoL=1; */
	  /* uL=1; */
	  
	  if (i==1 && j==1 && k==1) {
	  PetscPrintf(PETSC_COMM_WORLD, "!!init 7 gamma %le kappa %le zata_s %le\n", gamma, kappa, zeta_s);
	  PetscPrintf(PETSC_COMM_WORLD, "!!init 7 uL %le rhoL %le PL %le\n", uL, rhoL, PL);
	  PetscPrintf(PETSC_COMM_WORLD, "!!init 7 cR %le rhoR %le PR %le\n", cR, rhoR, PR);
	  }
	  if (k<k_s) {	
	    q[k][j][i].rho=rhoL;
	    q[k][j][i].rhoW=rhoL*uL;
	    q[k][j][i].rhoE=PL/(gamma-1)+0.5*rhoL*uL*uL;
	  }  else {
	    q[k][j][i].rho=rhoR;
	    q[k][j][i].rhoW=0.;
	    q[k][j][i].rhoE=PR/(gamma-1);
	  }
	} else if (InitialGuessOne==8) { // isothermal wall
	  q[k][j][i].rho  = 1.0;
	  q[k][j][i].rhoU = 0.0;
	  q[k][j][i].rhoV = 0.0;
	  q[k][j][i].rhoW = 0.0;
	  q[k][j][i].rhoE = 0.1/(user->gamma*user->Ma*user->Ma*(user->gamma-1.0));
	} else if (InitialGuessOne==9) {//Comp Taylor Green Vortex//
	  PetscReal p0=100, Ma=user->Ma, p;
	  
	  if (i==1 && j==1 && k==1) 
	    PetscPrintf(PETSC_COMM_WORLD, "!!init 9 comp TGV gamma %le Ma\n", user->gamma, Ma);  
	  p = p0+(cos(2*coor[k][j][i].x)+cos(2*coor[k][j][i].y))*(cos(2*coor[k][j][i].z)+2.)/16.;
	  q[k][j][i].rho  = user->gamma*p*Ma*Ma;
	  q[k][j][i].rhoU =  sin(coor[k][j][i].x)*cos(coor[k][j][i].y)*cos(coor[k][j][i].z);
	  q[k][j][i].rhoV = -cos(coor[k][j][i].x)*sin(coor[k][j][i].y)*cos(coor[k][j][i].z);
	  q[k][j][i].rhoW = 0.;//sin(coor[k][j][i].z)*cos(coor[k][j][i].y);
	  q[k][j][i].rhoE = p/(user->gamma-1.)
	    +0.5*(q[k][j][i].rhoU*q[k][j][i].rhoU+
		  q[k][j][i].rhoV*q[k][j][i].rhoV+
		  q[k][j][i].rhoW*q[k][j][i].rhoW)/q[k][j][i].rho;

	  //DMDAVecRestoreArray(da, user->P, &p);
	
	} else if (InitialGuessOne==10) {//Incomp Taylor Green Vortex//
	  //DMDAVecGetArray(da, user->P, &p);
	  if (i==1 && j==1 && k==1) 
	    PetscPrintf(PETSC_COMM_WORLD, "!!init 10 incp TGV gamma %le Pr \n", user->gamma, user->Pr);  
	  q[k][j][i].rho  = 1.0;
	  q[k][j][i].rhoU = 0.0;
	  q[k][j][i].rhoW = -cos(coor[k][j][i].z)*sin(coor[k][j][i].y);
	  q[k][j][i].rhoV =  sin(coor[k][j][i].z)*cos(coor[k][j][i].y);
	  q[k][j][i].rhoE = (-0.25*(cos(2*coor[k][j][i].z)+cos(2*coor[k][j][i].y)))/(user->gamma-1.)
	    +0.5*(q[k][j][i].rhoV*q[k][j][i].rhoV+
		  q[k][j][i].rhoW*q[k][j][i].rhoW)/q[k][j][i].rho;

	  //DMDAVecRestoreArray(da, user->P, &p);
	} else if (InitialGuessOne==52) {//cos  inlet for poisueille test
	  q[k][j][i].rho=1.;
	  q[k][j][i].rhoW=cos(pi/2*(coor[k][j][i].y-1.))*cos(pi/2*(coor[k][j][i].y-1.));
	  q[k][j][i].rhoE=0.5*q[k][j][i].rhoW*q[k][j][i].rhoW
	    +1./(user->Ma*user->Ma*user->gamma*(user->gamma-1.)); 
	  if (i==1 && j==20 && k==1) 
	    PetscPrintf(PETSC_COMM_WORLD, "!!init 52 vel %le rhoE \n", q[k][j][i].rhoW, q[k][j][i].rhoE);  
	}//
	else {	  
	  q[k][j][i].rho=rho;
	  q[k][j][i].rhoW=rho*u;
	  q[k][j][i].rhoE=rho*et;
	}
	
      }
    }
  }

  //////// Initialization from a specific u.dat
  if  (InitialGuessOne==15) {//M
    
    int rank,N=64;
    
    /*  Cmpnts***array=(Cmpnts***)malloc(N*sizeof(Cmpnts**));
	for (i=0;i<N;i++){
	array[i]=(Cmpnts**)malloc(N*sizeof(*Cmpnts));
	for (j=0;j<N;j++){
	array[i][j]=(Cmpnts*)malloc(N*sizeof(Cmpnts)); 
	}	
	}  
    */
    double *uarray=(double*)malloc(N*N*N*sizeof(double)); 
    double *varray=(double*)malloc(N*N*N*sizeof(double)); 
    double *warray=(double*)malloc(N*N*N*sizeof(double)); 
    
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if(!rank) { // root processor read in the data
      FILE *fd;
      PetscPrintf(PETSC_COMM_SELF, "READ ibmdata\n");
      char filen[80];    char string[128];
      
      sprintf(filen,"u.dat");
      
      fd = fopen(filen, "r"); if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open u.dat");
      
      if (fd) {
	//    fgets(string, 128, fd);
	//   fgets(string, 128, fd);
	//    fgets(string, 128, fd);
	
	for (k=0; k<N; k++) {
	  for(j=0; j<N; j++){
	    for(i=0;i<N;i++){
	      
	      fscanf(fd, "%le %le %le",  &uarray[k*N*N+j*N+i], &varray[k*N*N+j*N+i], &warray[k*N*N+j*N+i]);//, &t, &t, &t);
	      //			if (k==0 && j==0) 
	      //			   PetscPrintf(PETSC_COMM_WORLD, "ux %f \n",uarray[k*N*N+j*N+i]);
	      
	      
	    }
	  }
	}			
      }
      
      
    }
    MPI_Bcast(uarray, N*N*N, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(varray, N*N*N, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(warray, N*N*N, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
    
    
    for (k=zs; k<ze; k++){
      for (j=ys; j<ye; j++){
	for (i=xs; i<xe; i++){      
	  
	  q[k][j][i].rhoU=uarray[k*N*N+j*N+i];
	  q[k][j][i].rhoV=varray[k*N*N+j*N+i];
	  q[k][j][i].rhoW=warray[k*N*N+j*N+i];
	  
	  q[k][j][i].rho=1.0;
	  
	}
      }
    }	
  }
  /////////////////////////////	
  
  //Amir
  //  

  if(channel||wavy){
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	//	if (k==2 && j==15 )
	//	  PetscPrintf(PETSC_COMM_WORLD, "i  %d  rhs!!   %le   %le   %le\n",i,rhs[k][j][i].rhoU,rhs[k][j][i].rhoV,rhs[k][j][i].rhoW);
	if (k==2 && i==5 )
	  
	  PetscPrintf(PETSC_COMM_WORLD, "j  %d  fp!!   %le   %le   %le  %le   %le\n",j,coor[k][j][i].x,coor[k][j][i].y, q[k][j][i].rhoU, q[k][j][i].rhoV, q[k][j][i].rhoW);
	
	
      }	  
    }  
  }
  }

/////////////////////////////////// Initialization from a specific u.dat


  DMDAVecRestoreArray(cda, user->Q, &q);
  DMDAVecRestoreArray(fda, Coor, &coor);

  // PetscPrintf(PETSC_COMM_WORLD, "!!init  uL  rhoL  PL \n");

  DMGlobalToLocalBegin(cda, user->Q, INSERT_VALUES, user->lQ);
  DMGlobalToLocalEnd(cda, user->Q, INSERT_VALUES, user->lQ);
  // PetscPrintf(PETSC_COMM_WORLD, "!!init 2 uL  rhoL  PL \n");
  VecCopy(user->Q, user->Qold);
  VecCopy(user->lQ, user->lQold);

  // PetscPrintf(PETSC_COMM_WORLD, "!!init 2 uL  rhoL  PL \n");
  EqOfState(user);
  // PetscPrintf(PETSC_COMM_WORLD, "!!init 2 uL  rhoL  PL \n");
  return(0);

}

PetscErrorCode SetRhoSolution(PetscReal rho[3],PetscReal rho_f, PetscReal *sol)
{
  if (rho[0]>0) {
    *sol=rho[0];
    if (rho[1]>0 && fabs(rho[1]-rho_f)<fabs(*sol-rho_f)) *sol=rho[1];
    if (rho[2]>0 && fabs(rho[2]-rho_f)<fabs(*sol-rho_f)) *sol=rho[2];
  } else if (rho[1]>0) {
    *sol=rho[1];
    if (rho[2]>0 && fabs(rho[2]-rho_f)<fabs(*sol-rho_f)) *sol=rho[2];
  } else if (rho[2]>0) 
    *sol=rho[2];
  else {
    PetscPrintf(PETSC_COMM_WORLD, "!! All rho solutions negative !!! %f %f %f\n", rho[0],rho[1],rho[2]);
    return(1);
  }

  return(0);
}  	 

/* Boundary condition defination (array user->bctype[0-5]):
0:	interpolation/interface
-1:  wallfunction
1:	solid wall (not moving)
2:	moving solid wall (U=1)
3:   slip wall/symmetry
5:	Inlet
51:  Lax initial condition
55:  supersonic Inlet
52:  subsonic inlet NRBC cos profile for Poiseuille flow with isothermal wall case from Poinsot & Lele
4:	Outlet
44:  supersonic outlet 1
45:  supersonic outlet 2 (Neumann all var)
49:  Charcterstic BC (Hoffmann & Chiang voll II section 12.9.3.5, page 191)
6:   farfield
7:   periodic
8:   Characteristic BC NRBC subsonic
10:  Oulet Junction (inviscid vortex)
14:  Outlet with Interface
*/

PetscErrorCode FormBCS(UserCtx *user, PetscReal alfa)
{
  DM		cda = user->cda, da = user->da;
  DMDALocalInfo	info=user->info;
  PetscInt	xs, xe, ys, ye, zs, ze, lxs, lxe, lys, lye, lzs, lze, ls, le; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  //  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k, tmp, itr_tmp=1, twoD;
  
  PetscReal	***p,***lp,***nvert,***aj, Twall=1., AA,BB, CC[4], rho_sol[3];
  PetscReal     a_sound, L[5],Lambda[5], d[5], beta;
  PetscReal     pi=3.14159265;
  Cmpnts 	***zet,***eta;
  Vec           coords; 
  Cmpnts        ***coor;
  
  
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;

  PetscOptionsGetReal(NULL, NULL,"-Twall", &Twall, NULL);
  PetscOptionsGetInt(NULL, NULL, "-twoD", &twoD, NULL);
  PetscPrintf(PETSC_COMM_WORLD, "BCS isothermal wall %d Tw %f\n", isothermal, Twall);
   
  DMGetCoordinatesLocal(da, &coords);
  DMDAVecGetArray(user->fda, coords, &coor);
  
  CompVars	***q,***lq,***qold;

  DMDAVecGetArray(cda, user->Q, &q);
  DMDAVecGetArray(cda, user->Qold, &qold);
  DMDAVecGetArray(cda, user->lQ, &lq);

  double A,nx,ny,nz,Un;
  
  
  DMDAVecGetArray(da, user->lP, &lp);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(user->fda, user->lZet, &zet);
  DMDAVecGetArray(user->fda, user->lEta, &eta);
  DMDAVecGetArray(da, user->lAj, &aj);
  
  
  ////////////mass and flux correction for channel 
  
  double fluxtemp=0.0,Area=0.0,vol=0.0,vol0=0.0,vol0sum=0.0,rhotemp0=0.,rhotemp=0.0,ratio=0.,fluxsum,Areasum,rhosum,ratiorho,volsum;
 

  //	if (istage==2){
  if (channel || wavy) {
    // for (k=zs; k<ze; k++){
    //  }
    
    
    for(k=zs; k<ze; k++){ 
      for (j=ys; j<ye; j++) {
		for (i=xs; i<xe; i++) { 
			if (nvert[k][j][i]<1.5){	
				//fluxtemp+=q[k][j][i].rhoW/aj[k][j][i];
				vol+= 1./aj[k][j][i];					  // Total Volume of Fluid Cells
				Area+=zet[k][j][i].z;				      // 1/J ∂ζ/∂z
				rhotemp+=q[k][j][i].rho*(1./aj[k][j][i]); // mass per unit cell
				fluxtemp+=q[k][j][i].rhoW*zet[k][j][i].z; // 1/J rhoW ∂ζ/∂z
				
			}
			if (nvert[k][j][i]<1.5 && nvert[k][j][i]>0.5){ // Immersed Boundary Node
				vol0+= 1./aj[k][j][i];
				rhotemp0+=q[k][j][i].rho*(1./aj[k][j][i])/(mx-1)/(mz-1);
			}
	  
		}
      }
    }	
    MPI_Allreduce( &Area, &Areasum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
    MPI_Allreduce( &fluxtemp, &fluxsum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce( &rhotemp, &rhosum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
    MPI_Allreduce( &vol, &volsum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
    MPI_Allreduce( &vol0, &vol0sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
    
    // if (ti==0 || ti==tistart) 
    //  user->flux=fluxsum;
    ratio= (Areasum-fluxsum)/Areasum;
    ratiorho= (volsum-rhosum)/volsum;
    user->ratio=ratio;
    
    PetscPrintf(PETSC_COMM_WORLD, " flux %.15le  volume  %15le  rhosume %le   ratio  %.15le    ratiorho  %15le    \n",fluxsum,volsum,rhosum,ratio,ratiorho );
    
    rhotemp=0.0;	
    for(k=zs; k<ze; k++){ 
      for (j=ys; j<ye; j++) {
		for (i=xs; i<xe; i++) { 
	  		if (nvert[k][j][i]<1.5){	
	    	//fluxtemp+=q[k][j][i].rhoW/aj[k][j][i];
	    	rhotemp+=q[k][j][i].rho*(1./aj[k][j][i]);
	    	q[k][j][i].rho+=ratiorho;	    	    	    
	  		}
		}
      }
    }	
    
    MPI_Allreduce( &rhotemp, &rhosum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
    ratiorho= volsum-rhosum;
    PetscPrintf(PETSC_COMM_WORLD, "  ratiorhoUpdated  %15le  \n",ratiorho );
    
    
  } // if wavy
  //istage	}	
  
  
  // i BC
  for (k=zs; k<ze; k++){
    for (j=ys; j<ye; j++){      	 
      if (xs==0) {
		i=xs;
		if (user->bctype[0]==3) {
		q[k][j][i].rho=q[k][j][i+1].rho;
		q[k][j][i].rhoU=0.;
		q[k][j][i].rhoV=q[k][j][i+1].rhoV;
		q[k][j][i].rhoW=q[k][j][i+1].rhoW;
		q[k][j][i].rhoE=q[k][j][i+1].rhoE;
		} else 	if (user->bctype[0]==51) {
		q[k][j][i].rho=0.445;
		q[k][j][i].rhoU=0.311;
		q[k][j][i].rhoV=0.;
		q[k][j][i].rhoW=0.;
		q[k][j][i].rhoE=8.928;
		}
      }
      if (xe==mx) {
	  i=mx-1;
		if (user->bctype[1]==3) {
		q[k][j][i].rho=q[k][j][i-1].rho;
		q[k][j][i].rhoU=0.;
		q[k][j][i].rhoV=q[k][j][i-1].rhoV;
		q[k][j][i].rhoW=q[k][j][i-1].rhoW;
		q[k][j][i].rhoE=q[k][j][i-1].rhoE;
		} else 	if (user->bctype[0]==51) {
		q[k][j][i].rho=0.5;
		q[k][j][i].rhoU=0.;
		q[k][j][i].rhoV=0.;
		q[k][j][i].rhoW=0.;
		q[k][j][i].rhoE=1.4275;
		}
      }      
      
    }
  }
  
  // j BC
  for (k=zs; k<ze; k++){
    for (i=xs; i<xe; i++){                 
      if (ys==0) {
	j=ys;
	if (user->bctype[2]==3) {
	  q[k][j][i].rho=q[k][j+1][i].rho;
	  q[k][j][i].rhoU=q[k][j+1][i].rhoU;
	  q[k][j][i].rhoV=0.;
	  q[k][j][i].rhoW=q[k][j+1][i].rhoW;
	  q[k][j][i].rhoE=q[k][j+1][i].rhoE;
	} 
	else if (user->bctype[2]==10) { //inviscid vortex
	q[k][j][i].rho=q[k][j+1][i].rho;
	q[k][j][i].rhoU=0.;
	q[k][j][i].rhoV=q[k][j+1][i].rhoV;
	q[k][j][i].rhoW=q[k][j+1][i].rhoW;
	q[k][j][i].rhoE=q[k][j+1][i].rhoE;
	} else if (user->bctype[2]==1) { // No Slip y wall
	  // q[k][j][i].rho=q[k][j+1][i].rho;
	  // q[k][j][i].rho=2.0*q[k][j+1][i].rho-q[k][j+2][i].rho;			
	  q[k][j][i].rhoU=0.0;
	  q[k][j][i].rhoV=0.0;
	  q[k][j][i].rhoW=0.0;
	  if (isothermal==2 && alfa!=0 && nvert[k][j][i]<1.5){
	    /* // A=Tw/(gamma M^2) */
	    /*  AA=Twall/(user->gamma*user->Ma*user->Ma); */
	    /* // B=(T_f-Tw)/(gamma M^2)=p_f/rho_f-Tw/(gamma M^2) */
	    /*  BB=lp[k][j+1][i]/q[k][j+1][i].rho-Twall/(user->gamma*user->Ma*user->Ma); */
	    /* // The Cubic equation */
	    /*
	    CC[0]=-q[k][j+1][i].rho*AA; // -A*rho_f
	    CC[1]= AA-BB; //A-B
	    CC[2]= lp[k][j+1][i]; //p_f
	    CC[3]= -AA;
	    SolveCubic(CC, &rho_sol);
	    */
	    // The quadrtic equation
	    /*
	    CC[0]=1+BB/AA;
	    CC[1]=-(1/q[k][j+1][i].rho+lp[k][j+1][i]/AA);
	    CC[2]=1.;
	    SolveQuadratic(CC, &rho_sol);

	    SetRhoSolution(rho_sol,q[k][j+1][i].rho,&q[k][j][i].rho); 
	    */
	    // lp[k][j][i]= q[k][j][i].rho*Twall/(user->gamma*user->Ma*user->Ma);
	    // lp[k][j][i]=1./(user->gamma*user->Ma*user->Ma)*q[k][j][i].rho;
	  

	    //q[k][j][i].rho=lp[k][j][i]*(user->gamma*user->Ma*user->Ma)/Twall;
	    // q[k][j][i].rho=lp[k][j][i]*(user->gamma*user->Ma*user->Ma);

	    /* This is the no-slip isothermal wall BC based on characteristics 
	     refer to Poinsot and Lele https://doi.org/10.1016/0021-9991(92)90046-2 */
	    
	    // compute the characteristics lambda5=u_n+c where c^2=gamma p/rho speed of sound
	    //a_sound=sqrt(user->gamma*lp[k][j+1][i]/lq[k][j+1][i].rho);
	    a_sound=sqrt(Twall/(user->Ma*user->Ma));
	    beta=sqrt(eta[k][j][i].y*eta[k][j][i].y +
		      eta[k][j][i].x*eta[k][j][i].x +
		      eta[k][j][i].z*eta[k][j][i].z);
	    //  Lambda[0]=lq[k][j+1][i].rhoV/lq[k][j+1][i].rho - a_sound; // lambda1 = u - a
	    Lambda[0]= -a_sound*beta; // lambda1 = u - a where u=0 at the wall
	    // L1=lambda1*(dp/dn - rho c du_n/dn
	    if (k==0 || k==mz-1) 
	      L[0]=0.;
	    else
	      for (tmp=0; tmp<itr_tmp; tmp++) {
		L[0]=Lambda[0]*aj[k][j][i]* // aj is for beta
		  ((-0.5*lp[k][j+2][i]+2*lp[k][j+1][i]-1.5*lp[k][j][i]) - 
		   lq[k][j][i].rho*a_sound/beta*(eta[k][j][i].y* 
		   //((lp[k][j+1][i]-lp[k][j][i]) - lq[k][j][i].rho*a_sound*
		   //((lp[k][j+2][i]-lp[k][j+1][i]) - lq[k][j+1][i].rho*a_sound*
		   //(lq[k][j+2][i].rhoV/lq[k][j+2][i].rho - lq[k][j+1][i].rhoV/lq[k][j+1][i].rho));
		   //(lq[k][j+1][i].rhoV/lq[k][j+1][i].rho - lq[k][j][i].rhoV/lq[k][j][i].rho));
		   (-0.5*lq[k][j+2][i].rhoV/lq[k][j+2][i].rho+2.*lq[k][j+1][i].rhoV/lq[k][j+1][i].rho
		    -1.5*lq[k][j][i].rhoV/lq[k][j][i].rho) + eta[k][j][i].z*
		   (-0.5*lq[k][j+2][i].rhoW/lq[k][j+2][i].rho+2.*lq[k][j+1][i].rhoW/lq[k][j+1][i].rho
		    -1.5*lq[k][j][i].rhoW/lq[k][j][i].rho)));
		if (twoD!=1) L[0]-=Lambda[0]*aj[k][j][i]*lq[k][j][i].rho*a_sound/beta*eta[k][j][i].x* 
		   (-0.5*lq[k][j+2][i].rhoU/lq[k][j+2][i].rho+2.*lq[k][j+1][i].rhoU/lq[k][j+1][i].rho
		    -1.5*lq[k][j][i].rhoU/lq[k][j][i].rho);

		
		//(lq[k][j+1][i].rhoV/lq[k][j+1][i].rho - lq[k][j][i].rhoV/lq[k][j][i].rho));
		q[k][j][i].rho=qold[k][j][i].rho - user->dt*alfa*L[0]/a_sound/a_sound;
		lp[k][j][i]= q[k][j][i].rho*Twall/(user->gamma*user->Ma*user->Ma);
	      }
	  } else { //else if (isothermal==0) //advait!
	    lp[k][j][i]=lp[k][j+1][i];
		//q[k][j][i].rho=q[k][j+1][i].rho; // Just check the boundary conditions again
		q[k][j][i].rho = lp[k][j][i]*(user->gamma*user->Ma*user->Ma);
		q[k][j][i].rhoE = lp[k][j][i]/(user->gamma-1);

	  }

	//   q[k][j][i].rhoE=lp[k][j][i]/(user->gamma-1);
	}	
      }

      if (ye==my) {
	j=my-1;
	if (user->bctype[3]==3) {
	q[k][j][i].rho=q[k][j-1][i].rho;
	q[k][j][i].rhoU=q[k][j-1][i].rhoU;
	q[k][j][i].rhoV=0.;
	q[k][j][i].rhoW=q[k][j-1][i].rhoW;
	q[k][j][i].rhoE=q[k][j-1][i].rhoE;
	}  else  if (user->bctype[3]==10) { //inviscid vortex
	q[k][j][i].rho=q[k][j-1][i].rho;
	q[k][j][i].rhoU=0.;
	q[k][j][i].rhoV=q[k][j-1][i].rhoV;
	q[k][j][i].rhoW=q[k][j-1][i].rhoW;
	q[k][j][i].rhoE=q[k][j-1][i].rhoE;
	}    
	/* if (user->bctype[3]==3) {	   */
	/*   q[k][j][i].rho=q[k][j-1][i].rho; */
	/*   q[k][j][i].rhoE=q[k][j-1][i].rhoE; */
	/*   A=sqrt(eta[k][j][i].z*eta[k][j][i].z + */
	/* 	 eta[k][j][i].y*eta[k][j][i].y + */
	/* 	 eta[k][j][i].x*eta[k][j][i].x); */
	/*   nx=eta[k][j][i].x/A; */
	/*   ny=eta[k][j][i].y/A; */
	/*   nz=eta[k][j][i].z/A; */
	/*   Un=(q[k][j-1][i].rhoU*nx+q[k][j-1][i].rhoV*ny+q[k][j-1][i].rhoW*nz); */
	/*   q[k][j][i].rhoU = q[k][j-1][i].rhoU-Un*nx; */
	/*   q[k][j][i].rhoV = q[k][j-1][i].rhoV-Un*ny; */
	/*   q[k][j][i].rhoW = q[k][j-1][i].rhoW-Un*nz;	*/					 
	/* }   */
	else if (user->bctype[3]==1) {
	  //	q[k][j][i].rho=	2.*q[k][j-1][i].rho-q[k][j-2][i].rho;
	  q[k][j][i].rhoU=0.0;
	  q[k][j][i].rhoV=0.0;
	  q[k][j][i].rhoW=0.0;
	  //lp[k][j][i]=lp[k][j-1][i];	  
	  
	  if (isothermal == 2){
	    /* // A=Tw/(gamma M^2) */
	    /*  AA=Twall/(user->gamma*user->Ma*user->Ma); */
	    /* // B=(T_f-Tw)/(gamma M^2)=p_f/rho_f-Tw/(gamma M^2) */
	    /*  BB=lp[k][j+1][i]/q[k][j+1][i].rho-Twall/(user->gamma*user->Ma*user->Ma); */
	    /* // The Cubic equation */
	    /*
	    CC[0]=-q[k][j+1][i].rho*AA; // -A*rho_f
	    CC[1]= AA-BB; //A-B
	    CC[2]= lp[k][j+1][i]; //p_f
	    CC[3]= -AA;
	    SolveCubic(CC, &rho_sol);
	    */
	    // The quadrtic equation
	    /*
	    CC[0]=1+BB/AA;
	    CC[1]=-(1/q[k][j+1][i].rho+lp[k][j+1][i]/AA);
	    CC[2]=1.;
	    SolveQuadratic(CC, &rho_sol);

	    SetRhoSolution(rho_sol,q[k][j+1][i].rho,&q[k][j][i].rho); 
	    */
	    // lp[k][j][i]= q[k][j][i].rho*Twall/(user->gamma*user->Ma*user->Ma);
	    // lp[k][j][i]=1./(user->gamma*user->Ma*user->Ma)*q[k][j][i].rho;
	  

	    //q[k][j][i].rho=lp[k][j][i]*(user->gamma*user->Ma*user->Ma)/Twall;
	    // q[k][j][i].rho=lp[k][j][i]*(user->gamma*user->Ma*user->Ma);

	    /* This is the no-slip isothermal wall BC based on characteristics 
	     refer to Poinsot and Lele https://doi.org/10.1016/0021-9991(92)90046-2 */
	    
	    // compute the characteristics lambda5=u_n+a where a^2=gamma p/rho speed of sound
	    a_sound=sqrt(user->gamma*lp[k][j-1][i]/lq[k][j-1][i].rho);
	    Lambda[0]=lq[k][j-1][i].rhoV/lq[k][j-1][i].rho - a_sound; //u-a
	    // L1=lambda1*(dp/dn - rho a du_n/dn
	    if (k==0 || k==mz-1 || Lambda[0]>0.) 
	      L[0]=0.;
	    else
	      L[0]=Lambda[0]*aj[k][j-1][i]*eta[k][j-1][i].y*
		((lp[k][j-1][i]-lp[k][j-2][i]) - lq[k][j-1][i].rho*a_sound*
		 // (lq[k][j-1][i].rhoV/lq[k][j-1][i].rho - lq[k][j-2][i].rhoV/lq[k][j-2][i].rho));
		 (lq[k][j][i].rhoV/lq[k][j][i].rho - lq[k][j-1][i].rhoV/lq[k][j-1][i].rho));
	    q[k][j][i].rho=qold[k][j][i].rho - user->dt*alfa*L[0]/a_sound/a_sound;
	    lp[k][j][i]= q[k][j][i].rho*Twall/(user->gamma*user->Ma*user->Ma);
	  }
	  else if (isothermal==1){
		q[k][j][i].rho=lp[k][j][i]*(user->gamma*user->Ma*user->Ma);
		lp[k][j][i]=lp[k][j-1][i];
	  }else {	    
	    q[k][j][i].rho=q[k][j-1][i].rho;
	    lp[k][j][i]=lp[k][j-1][i];
	  }
	  q[k][j][i].rhoE=lp[k][j][i]/(user->gamma-1.0);
	  
	}     
      }      
      
    }
  }

  
  PetscReal  u=user->Ma, rho=user->gamma, et; //Mach 3.
    //u=1153, rho=0.00237, et=1719557; //Mach number M=1.5
  if (user->bctype[2]==1) ls=lys; else ls=ys;
  if (user->bctype[3]==1) le=lye; else le=ye;

  // k BC
  for (j=ls; j<le; j++){ 
    for (i=xs; i<xe; i++){      	 
      if (zs==0) {
	k=zs;
	if (user->bctype[4]==3) { // slip
	  q[k][j][i].rho=q[k+1][j][i].rho;
	  q[k][j][i].rhoU=q[k+1][j][i].rhoU;
	  q[k][j][i].rhoV=q[k+1][j][i].rhoV;
	  q[k][j][i].rhoW=0.;
	  q[k][j][i].rhoE=q[k+1][j][i].rhoE;
	} else	if (user->bctype[4]==52) { //subsonic NRBC - cos inlet for poisueille test
	  AA=cos(pi/2*(coor[k][j][i].y-1.))*cos(pi/2*(coor[k][j][i].y-1.)); // velocity cos inlet
	  a_sound=sqrt(user->gamma*lp[k+1][j][i]/lq[k+1][j][i].rho);
	  Lambda[0]=lq[k+1][j][i].rhoW/lq[k+1][j][i].rho-a_sound; //lambda1=u-c
	  //L1=lambda1(dp/dn-rho a du/dn) use FD for p and BD for u because of the sign in the front
	  L[0]=Lambda[0]*aj[k+1][j][i]*zet[k+1][j][i].z*
	    ((lp[k+2][j][i]-lp[k+1][j][i])-
	     lq[k+1][j][i].rho*a_sound*
	     (lq[k+1][j][i].rhoW/lq[k+1][j][i].rho-AA));
	     //(lq[k+2][j][i].rhoW/lq[k+2][j][i].rho-lq[k+1][j][i].rhoW/lq[k+1][j][i].rho));	
	  L[1]=(user->gamma-1.)*L[0];
	  L[4]=L[0];
	  d[0]=1./(a_sound*a_sound)*(L[1]+0.5*(L[4]+L[0]));
	  q[k][j][i].rho=qold[k][j][i].rho-user->dt*alfa*d[0];
	  lp[k][j][i]= q[k][j][i].rho/(user->gamma*user->Ma*user->Ma); //T_inf=1
	  //q[k][j][i].rho=1.;
	  q[k][j][i].rhoW=q[k][j][i].rho*AA;
	  q[k][j][i].rhoU=0.;
	  q[k][j][i].rhoV=0.;
	  q[k][j][i].rhoE=0.5*q[k][j][i].rhoW*q[k][j][i].rhoW/q[k][j][i].rho+lp[k][j][i]/(user->gamma-1);
	  //	    +1./(user->Ma*user->Ma*user->gamma*(user->gamma-1.)); 
	} else	if (user->bctype[4]==55) { //supersonic inlet
	  q[k][j][i].rho=rho;
	  q[k][j][i].rhoW=rho*u;
	  q[k][j][i].rhoE=0.5*rho*u*u+1./(user->gamma-1.); //et;
	} else if (user->bctype[4]==44) { //supersonic outlet
	  q[k][j][i].rho=q[k+1][j][i].rho;
	  //q[k][j][i].rho=1.0;
	  //q[k][j][i].rho=2*q[k-1][j][i].rho-q[k-2][j][i].rho;
	  q[k][j][i].rhoU=0.;
	  q[k][j][i].rhoV=0.;
	  //	q[k][j][i].rhoW=2*q[k-1][j][i].rhoW-q[k-2][j][i].rhoW;
	  //q[k][j][i].rhoE=2*q[k-1][j][i].rhoE-q[k-2][j][i].rhoE;
	  q[k][j][i].rhoW=q[k+1][j][i].rhoW;//1.0;
	  q[k][j][i].rhoE=q[k+1][j][i].rhoE;
	} else if (user->bctype[4]==45) { //supersonic inflow
	  q[k][j][i].rho=1.;
	  q[k][j][i].rhoW=1.;
	  q[k][j][i].rhoE=.5+1./(user->Ma*user->Ma*user->gamma*(user->gamma-1.));
	} else  if (user->bctype[4]==51) { //Lax 
	    q[k][j][i].rho=0.445;
	    q[k][j][i].rhoW=0.311;
	    q[k][j][i].rhoE=8.928;
	} else if (user->bctype[4]==10) { //inviscid vortex
	  q[k][j][i].rho=1.;
	  q[k][j][i].rhoU=0.;
	  q[k][j][i].rhoV=0.;
	  q[k][j][i].rhoW=1.;//0.5*sqrt(1.4);
	  q[k][j][i].rhoE=1./0.4+0.5;//0.125*1.4;
	}
	else if (user->bctype[4]==47) { 
	  
	  double pi=3.141592, yy,yplus;
	  double c1=0.03, c2=-0.0032, c3=0.02, c4=-0.01,yp1=1.0,yp2=1.5,yp3=2.0;
	  double w0 =0.1, w1=0.25, w2=0.125, w3=0.0625; 
	  double phi1=0.0, phi2=0.1, phi3=0.15, b0=0.5*pi, b1=0.75*pi, b2=0.5*pi, b3=0.25*pi;
	  
	  yy=coor[k][j][i].y;
	  yplus=yy*0.065*user->ren;
	  if (yy<8.0){
	    q[k][j][i].rhoW=pow(yy/8.0,0.1428);
	  }else{ 
	    q[k][j][i].rhoW=1.0 ;
	    q[k][j][i].rhoU= 0.0 ;
	    
	    
	    
	    //	}
	    //	    }
	    //mean velocity
	  }	
	  if (yy<8.0){
	    q[k][j][i].rhoW+=c1*yplus*exp(-yplus/12.5)*sin(w0*user->ti*user->dt)*cos(b0*coor[k][j][i].x);
	    
	    q[k][j][i].rhoV=c2*yplus*yplus*exp(-yplus*yplus/12.5/12.5)*sin(w0*user->ti*user->dt)*cos(b0*coor[k][j][i].x);
	  } else if (yy>=8.0) {
	    
	    q[k][j][i].rhoW+=c3*yy/yp1*exp(-yy/yp1)*sin(w1*user->ti*user->dt)*cos(b1*coor[k][j][i].x+phi1)+
	      c3*yy/yp2*exp(-yy/yp2)*sin(w2*user->ti*user->dt)*cos(b2*coor[k][j][i].x+phi2)+
	      c3*yy/yp3*exp(-yy/yp3)*sin(w3*user->ti*user->dt)*cos(b3*coor[k][j][i].x+phi3);
	    
	    q[k][j][i].rhoV= c4*yy*yy/yp1/yp1*exp(-yy*yy/yp1/yp1)*sin(w1*user->ti*user->dt)*cos(b1*coor[k][j][i].x+phi1)+
	      c4*yy*yy/yp2/yp2*exp(-yy*yy/yp2/yp2)*sin(w2*user->ti*user->dt)*cos(b2*coor[k][j][i].x+phi2)+
	      c4*yy*yy/yp3/yp3*exp(-yy*yy/yp3/yp3)*sin(w3*user->ti*user->dt)*cos(b3*coor[k][j][i].x+phi3);
	    
	  }			
	
	  q[k][j][i].rho=1.0;//+0.8*q[k+1][j][i].rho;
	  q[k][j][i].rhoE=0.5*q[k][j][i].rhoW*q[k][j][i].rhoW+1./(user->Ma*user->Ma*user->gamma*(user->gamma-1.));
	  //			q[k][j][i].rhoE=1./(user->Ma*user->Ma*user->gamma*(user->gamma-1.));
	  
	  //		q[k][j][i].rhoV=0.0;
	  	  	  
	} // bcs 47
      }
      
      if (ze==mz) {
	k=mz-1;
	if (user->bctype[5]==3) {
	  q[k][j][i].rho=q[k-1][j][i].rho;
	  q[k][j][i].rhoU=q[k-1][j][i].rhoU;
	  q[k][j][i].rhoV=q[k-1][j][i].rhoV;
	  q[k][j][i].rhoW=0.;
	  q[k][j][i].rhoE=q[k-1][j][i].rhoE;
	} else 	if (user->bctype[5]==44) { //supersonic outflow 1
	  q[k][j][i].rho=q[k-1][j][i].rho;
	  //q[k][j][i].rho=2*q[k-1][j][i].rho-q[k-2][j][i].rho;
	  q[k][j][i].rhoU=0.;
	  q[k][j][i].rhoV=0.;
	  //	q[k][j][i].rhoW=2*q[k-1][j][i].rhoW-q[k-2][j][i].rhoW;
	  //q[k][j][i].rhoE=2*q[k-1][j][i].rhoE-q[k-2][j][i].rhoE;
	  q[k][j][i].rhoW=q[k-1][j][i].rhoW;
	  q[k][j][i].rhoE=q[k-1][j][i].rhoE;
	} else if (user->bctype[5]==45) { // supersonic outflow
	  q[k][j][i].rho =q[k-1][j][i].rho;
	  q[k][j][i].rhoU=q[k-1][j][i].rhoU;
	  q[k][j][i].rhoV=q[k-1][j][i].rhoV;
	  q[k][j][i].rhoW=q[k-1][j][i].rhoW;
	  q[k][j][i].rhoE=q[k-1][j][i].rhoE;
	} else if (user->bctype[5]==49) { // characteristic BC
	  PetscReal Vn, Rplus, Rminus, lambda1, lambda2, s, pOverrho, pp;
	  // I am just considering s= p/rho^gamma
	  Vn=q[k-1][j][i].rhoW/q[k-1][j][i].rho;
	  a_sound=sqrt(user->gamma*lp[k-1][j][i]/q[k-1][j][i].rho);	  	  
	  lambda1=Vn + a_sound;
	  lambda2=Vn - a_sound;
	  if (Vn<0.) { //treated as inflow
	    if (lambda1<0) {//set everthing equal to free stream
	      q[k][j][i].rhoU=0.;
	      q[k][j][i].rhoV=0.;
	      Rplus = Vn + 2/user->Ma/(user->gamma-1.);	 // Rplus = 1. + 2/user->Ma/(user->gamma-1.);	      
	      Rminus= Vn - 2/user->Ma/(user->gamma-1.);  // Rminus= 1. - 2/user->Ma/(user->gamma-1.);           	   
	      s=1./(user->Ma*user->Ma*user->gamma);
	      // q[k][j][i].rho=1.;
	      // q[k][j][i].rhoW=1.;
	      //q[k][j][i].rhoE=.5+1./(user->Ma*user->Ma*user->gamma*(user->gamma-1.));
	    } else { //R+ is interpolated, but  R-, s, and Vt are set to free stream
	      Rplus = Vn + 2*a_sound /(user->gamma-1.);	      
	      Rminus= Vn - 2/user->Ma/(user->gamma-1.); // Rminus= 1. - 2/user->Ma/(user->gamma-1.);
	      q[k][j][i].rhoU=0.;
	      q[k][j][i].rhoV=0.;
	      s=1./(user->Ma*user->Ma*user->gamma);
	    }
	  } else { // Vn>=0 treated as outflow
	    
	    if (lambda2 >=  0) { // interpolate everything
	      Rplus = Vn + 2*a_sound /(user->gamma-1.);	      
	      Rminus= Vn - 2*a_sound /(user->gamma-1.);	      
	      q[k][j][i].rhoU=q[k-1][j][i].rhoU;
	      q[k][j][i].rhoV=q[k-1][j][i].rho;
	      s = lp[k-1][j][i]/pow(q[k-1][j][i].rho, user->gamma);  
	    } else { // lambda2<0 set R- to freesteam and the others are inerpolated
	      Rminus= Vn - 2/user->Ma/(user->gamma-1.); // Rminus= 1. - 2/user->Ma/(user->gamma-1.);
	      Rplus = Vn + 2*a_sound /(user->gamma-1.);	
	      q[k][j][i].rhoU=q[k-1][j][i].rhoU;
	      q[k][j][i].rhoV=q[k-1][j][i].rhoV;
	      s = lp[k-1][j][i]/pow(q[k-1][j][i].rho, user->gamma);  
	    }
	  }
	  // computing variables
	  Vn=0.5*(Rplus+Rminus);
	  a_sound=0.25*(user->gamma-1.)*(Rplus-Rminus);
	  pOverrho=a_sound*a_sound/(user->gamma);
	  // rho= (pOverrho / s) ^ 1/(gamma-2) //Note s= p/rho^gamma here
	  q[k][j][i].rho = pow(pOverrho/s, 1./(user->gamma-2.));
	  q[k][j][i].rhoW=q[k][j][i].rho*Vn;
	  // p = pOverrho * rho
	  pp = pOverrho * q[k][j][i].rho;
	  // ideal gas law: rhoE=1/2 rho|u|^2+p/(gamma-1)
	  q[k][j][i].rhoE= pp/(user->gamma-1.) +0.5*(q[k][j][i].rhoU*q[k][j][i].rhoU+
						     q[k][j][i].rhoV*q[k][j][i].rhoV+
						     q[k][j][i].rhoW*q[k][j][i].rhoW)*q[k][j][i].rho;
	
	} else 	if (user->bctype[5]==51) { //Lax
	  q[k][j][i].rho=0.5;
	  q[k][j][i].rhoU=0.;
	  q[k][j][i].rhoV=0.;
	  q[k][j][i].rhoW=0.;
	  q[k][j][i].rhoE=1.4275;
	} else 	if (user->bctype[5]==10) { //inviscid vortex
	  q[k][j][i].rho=1.;
	  q[k][j][i].rhoU=0.;
	  q[k][j][i].rhoV=0.;
	  q[k][j][i].rhoW=1.;
	  q[k][j][i].rhoE=1/0.4+0.5;
	/* q[k][j][i].rho=q[k-1][j][i].rho; */
	/* q[k][j][i].rhoU=0.; */
	/* q[k][j][i].rhoV=q[k-1][j][i].rho; */
	/* q[k][j][i].rhoW=q[k-1][j][i].rho; */
	/* q[k][j][i].rhoE=q[k-1][j][i].rho; */
	} else 	if (user->bctype[5]==8) { //NRBC characteristic BC subsonic
	  a_sound=sqrt(user->gamma*lp[k-1][j][i]/lq[k-1][j][i].rho);
	  AA=q[k-1][j][i].rhoW/q[k-1][j][i].rho;
	  Lambda[0]=AA-a_sound; // lambda1=u-a
	  Lambda[1]=AA; Lambda[2]=AA; Lambda[3]=AA; //u
	  Lambda[4]=AA+a_sound; // lambda5=u+a
	  L[0]=0;
	  L[1]=Lambda[1]*aj[k-1][j][i]*zet[k-1][j][i].z*
	    (a_sound*a_sound*(lq[k-1][j][i].rho-lq[k-2][j][i].rho)-
	     (lp[k][j][i]-lp[k-1][j][i]));
	  //(lp[k-1][j][i]-lp[k-2][j][i]));
	  if (TwoD!=1) L[2]=Lambda[2]*aj[k-1][j][i]*zet[k-1][j][i].z*
			 (lq[k-1][j][i].rhoU/q[k-1][j][i].rho
			  -lq[k-2][j][i].rhoU/q[k-1][j][i].rho); 
	  else L[2]=0.;
	  if (TwoD!=2) L[3]=Lambda[3]*aj[k-1][j][i]*zet[k-1][j][i].z*
			 (lq[k-1][j][i].rhoV/q[k-1][j][i].rho-
			  lq[k-2][j][i].rhoV/q[k-2][j][i].rho);
	  else L[3]=0.;
	  // L5=lambda5(dp/dn + rho a du/dn) both BD because lambda +
	  L[4]=Lambda[4]*aj[k-1][j][i]*zet[k-1][j][i].z*
	    ((lp[k-1][j][i]-lp[k-2][j][i]) + q[k-1][j][i].rho*a_sound*
	     (lq[k-1][j][i].rhoW/q[k-1][j][i].rho-
	      lq[k-2][j][i].rhoW/q[k-2][j][i].rho)); 
	  
	  d[0]=1./(a_sound*a_sound)*(L[1]+0.5*(L[4]+L[0]));
	  d[1]=0.5*(L[0]+L[4]);
	  d[2]=0.5/(q[k-1][j][i].rho*a_sound)*(L[4]-L[0]);
	  d[3]=L[2];
	  d[4]=L[3];

	  q[k][j][i].rho =qold[k][j][i].rho-user->dt*alfa*d[0];
	  q[k][j][i].rhoE=qold[k][j][i].rhoE-
	    user->dt*alfa*(0.5*d[0]*
			   (//lq[k-1][j][i].rhoU/q[k-1][j][i].rho*lq[k-1][j][i].rhoU/q[k-1][j][i].rho+
			    //lq[k-1][j][i].rhoV/q[k-1][j][i].rho*lq[k-1][j][i].rhoV/q[k-1][j][i].rho+
			    lq[k-1][j][i].rhoW/q[k-1][j][i].rho*lq[k-1][j][i].rhoW/q[k-1][j][i].rho)+
			   d[1]/(user->gamma-1.)+lq[k-1][j][i].rhoW*d[2]); //need to add other terms later
	  q[k][j][i].rhoW=qold[k][j][i].rhoW - user->dt*alfa*
	    (lq[k-1][j][i].rhoW/q[k-1][j][i].rho*d[0] + q[k-1][j][i].rho*d[2]); // need to add other terms later
	}
      }      
      
    }
  }
  
  /* DMDAVecRestoreArray(da, user->lP, &lp); */
  
  
  /* DMDAVecRestoreArray(cda, user->Q, &q); */
  
  
  /* DMGlobalToLocalBegin(cda, user->Q, INSERT_VALUES, user->lQ); */
  /* DMGlobalToLocalEnd(cda, user->Q, INSERT_VALUES, user->lQ); */
  
  /* //   */
  
  
  /* DMDAVecGetArray(cda, user->Q, &q); */
  /* DMDAVecGetArray(cda, user->lQ, &lq); */
  DMDAVecGetArray(da, user->P, &p);
  /* DMDAVecGetArray(da, user->lP, &lp); */
  
  if (user->bctype[0]==7 || user->bctype[1]==7){
    if (xs==0){
      i=xs;
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  if(k>0 && k<user->KM && j>0 && j<user->JM){
	    q[k][j][i]=lq[k][j][i-2];
	    p[k][j][i]=lp[k][j][i-2];
	    
	  }
	}
      }
    }
    
    if (xe==mx){
      i=mx-1;
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  if(k>0 && k<user->KM && j>0 && j<user->JM){
	    q[k][j][i]=lq[k][j][i+2];
	    p[k][j][i]=lp[k][j][i+2];
	    
	  }
	}
      }
    }
  }
  
  
  
  /* DMDAVecRestoreArray(cda, user->Q, &q); */
  /* DMDAVecRestoreArray(cda, user->lQ, &lq); */
  /* DMDAVecRestoreArray(da, user->P, &p); */
  /* DMDAVecRestoreArray(da, user->lP, &lp); */
  
  
  
  /* DMGlobalToLocalBegin(cda, user->Q, INSERT_VALUES, user->lQ); */
  /* DMGlobalToLocalEnd(cda, user->Q, INSERT_VALUES, user->lQ); */
  /* DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP); */
  /* DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP); */
  
  /* DMDAVecGetArray(cda, user->Q, &q); */
  /* DMDAVecGetArray(cda, user->lQ, &lq); */
  /* DMDAVecGetArray(da, user->P, &p); */
  /* DMDAVecGetArray(da, user->lP, &lp); */
  
  if (user->bctype[2]==7 || user->bctype[3]==7){
    if (ys==0){
      j=ys;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  if(k>0 && k<user->KM){
	    q[k][j][i]=lq[k][j-2][i];
	    p[k][j][i]=lp[k][j-2][i];
	    
	  }
	}
      }
    }
    
    if (ye==my){
      j=my-1;
      for (k=lzs; k<lze; k++) {
	for (i=xs; i<xe; i++) {
	  if(k>0 && k<user->KM ){
	    q[k][j][i]=lq[k][j+2][i];
	    p[k][j][i]=lp[k][j+2][i];
	    
	  }
	}
      }
    }
  }
  
  
  
  /* DMDAVecRestoreArray(cda, user->Q, &q); */
  /* DMDAVecRestoreArray(cda, user->lQ, &lq); */
  /* DMDAVecRestoreArray(da, user->P, &p); */
  /* DMDAVecRestoreArray(da, user->lP, &lp); */
  
  
  
  /* DMGlobalToLocalBegin(cda, user->Q, INSERT_VALUES, user->lQ); */
  /* DMGlobalToLocalEnd(cda, user->Q, INSERT_VALUES, user->lQ); */
  /* DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP); */
  /* DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP); */
  
  /* DMDAVecGetArray(cda, user->Q, &q); */
  /* DMDAVecGetArray(cda, user->lQ, &lq); */
  /* DMDAVecGetArray(da, user->P, &p); */
  /* DMDAVecGetArray(da, user->lP, &lp); */
  
  if (user->bctype[4]==7 || user->bctype[5]==7){
    if (zs==0){
      k=zs;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  // if(k>0 && k<user->KM && j>0 && j<user->JM){
	  q[k][j][i]=lq[k-2][j][i];
	  p[k][j][i]=lp[k-2][j][i];
	  
	}
      }
      //	}
    }
    
    if (ze==mz){
      k=mz-1;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  //		if(k>0 && k<user->KM && j>0 && j<user->JM){
	  q[k][j][i]=lq[k+2][j][i];
	  p[k][j][i]=lp[k+2][j][i];
	  
	  //		}
	}
      }
    }
  }
  
  /* PetscInt    twoD=0; */
  /* PetscOptionsGetInt(NULL, NULL, "-twoD", &twoD, NULL); */
  
  /* if (twoD==1) { */
  /*   for(k=zs; k<ze; k++){  */
  /*     for (j=ys; j<ye; j++) { */
  /* 	for (i=xs; i<xe; i++) {  */
  /* 	  q[k][j][i].rhoU = 0.; */
  /* 	} */
  /*     } */
  /*   } */
  /* } */
  /*
    double fluxtemp=0.0,Area=0.0,vol=0.0,rhotemp=0.0,ratio,fluxsum,Areasum,rhosum,ratiorho,volsum;
    if (channel || wavy){
    if (zs==0){
    k=zs;
    for (j=lys; j<lye; j++) {
    for (i=xs; i<xe; i++) { 
    if (nvert[k][j][i]<0.1){	
    fluxtemp+=sqrt(zet[k][j][i].z*zet[k][j][i].z)*q[k][j][i].rhoW/q[k][j][i].rho;
    Area+=sqrt(zet[k][j][i].z*zet[k][j][i].z);
    //  rhotemp+=sqrt(zet[k][j][i].z*zet[k][j][i].z)*q[k][j][i].rho;
    
    }}}}
    
    
    
    for(k=lzs; k<lze; k++){ 
    for (j=lys; j<lye; j++) {
    for (i=lxs; i<lxe; i++) { 
    if (nvert[k][j][i]<0.1){	
    //fluxtemp+=q[k][j][i].rhoW/aj[k][j][i];
    vol+= sqrt(zet[k][j][i].z*zet[k][j][i].z);
    rhotemp+=q[k][j][i].rho*sqrt(zet[k][j][i].z*zet[k][j][i].z);
    
    }}}}	
    
    
    MPI_Allreduce( &Area, &Areasum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
    MPI_Allreduce( &fluxtemp, &fluxsum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce( &rhotemp, &rhosum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
    MPI_Allreduce( &vol, &volsum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
    
    ratio= 1.0-fluxsum/Areasum;
    ratiorho= 1.0-rhosum/volsum;
    
    PetscPrintf(PETSC_COMM_WORLD, "Area  %.15le flux %.15le  volume  %15le   ratio  %.15le    ratiorho  %15le \n",Areasum,fluxsum,volsum,ratio,ratiorho );
    
    for (k=zs;k<ze;k++){
    for (j=lys; j<lye; j++) {
    for (i=xs; i<xe; i++) {
    if (nvert[k][j][i]<0.1){
    q[k][j][i].rho+=ratiorho;
    q[k][j][i].rhoW+=ratio*q[k][j][i].rho;
    q[k][j][i].rhoE+= ratio*q[k][j][i].rhoW;
    
    }}}}
    }
  */
  DMDAVecRestoreArray(da, user->lAj, &aj);
  
  DMDAVecRestoreArray(user->fda, user->lEta,  &eta);
  DMDAVecRestoreArray(user->fda, user->lZet, &zet);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  
  
  
  DMDAVecRestoreArray(cda, user->Q, &q);
  DMDAVecRestoreArray(cda, user->Qold, &qold);
  DMDAVecRestoreArray(cda, user->lQ, &lq);
  DMDAVecRestoreArray(da, user->P, &p);
  DMDAVecRestoreArray(da, user->lP, &lp);
  

  
  DMGlobalToLocalBegin(cda, user->Q, INSERT_VALUES, user->lQ);
  DMGlobalToLocalEnd(cda, user->Q, INSERT_VALUES, user->lQ);
  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
  
  
  
  DMDAVecRestoreArray(user->fda, coords, &coor);  
  
  return(0);
  
}

/* ==================================================================================             */
/*   RHS subroutines */
/* ==================================================================================             */

PetscErrorCode CalcUCat(UserCtx *user)
{
  DM		fda = user->fda, cda=user->cda;
  Cmpnts	***ucat;
  CompVars	***q;
  
  DMDALocalInfo	info=user->info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  
  PetscInt	i, j, k;
  PetscInt    twoD=0;
  PetscOptionsGetInt(NULL, NULL, "-twoD", &twoD, NULL);

  DMDAVecGetArray(cda, user->lQ, &q);
  
  /* calculate the velocities */
  DMDAVecGetArray(fda, user->Ucat,  &ucat);
  
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (twoD!=1) ucat[k][j][i].x=q[k][j][i].rhoU/q[k][j][i].rho; else ucat[k][j][i].x = 0.;
	ucat[k][j][i].y=q[k][j][i].rhoV/q[k][j][i].rho;
	ucat[k][j][i].z=q[k][j][i].rhoW/q[k][j][i].rho;
	//	if (i==1 && (j==0 ||j==1 || j==2) && (k==21 || k==22|| k==20))
	//  PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d u is %.15le v is %.15le  w is %.15le \n",i,j,k,ucat[k][j][i].x,ucat[k][j][i].y,ucat[k][j][i].z );
      }
    }
  }
  
  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  DMDAVecRestoreArray(cda, user->lQ, &q);
  
  DMGlobalToLocalBegin(fda,user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda,user->Ucat, INSERT_VALUES, user->lUcat);
  
  return(0);
}

PetscErrorCode formViscous(UserCtx *user, Vec Visc)
{

  Vec		Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;
  
  Cmpnts        ***ucat, ucat_h;
  
  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***alpha1, ***alpha2, ***alpha3,t,Ss;
  Vec           Alpha1,Alpha2,Alpha3;
  /* Cmpnts	***icsi, ***ieta, ***izet; */
  /* Cmpnts	***jcsi, ***jeta, ***jzet; */
  /* Cmpnts	***kcsi, ***keta, ***kzet; */
  
  PetscReal	***nvert,S2;
  
  DM		da = user->da, fda = user->fda, cda=user->cda;
  DMDALocalInfo	info=user->info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k;
  Vec		Fp1, Fp2, Fp3;
  CompVars	***fp1, ***fp2, ***fp3;
  CompVars	***visc, ***q;
  PetscReal	***aj,***bbeta, ***p;//, ***iaj, ***jaj, ***kaj;
  Cmpnts        ***lsigma;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  
  PetscReal     ajc;
  
  PetscReal	dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;
  PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
  PetscReal	g11, g21, g31;
  PetscReal	r11, r21, r31, r12, r22, r32, r13, r23, r33;
  PetscReal     div, rho_c, sigma_c; 
  PetscReal     dpdc,dpde,dpdz;
  PetscScalar	solid,innerblank;
  FaceStencil rho;
  PetscReal u[3];

  solid      = 0.5;
  innerblank = 7.;
  PetscReal   pr=user->Pr;
  PetscReal   gamma=user->gamma, C_sutherland=110.4, Tref=273.15, Tr;
  PetscReal   nu = 1./user->ren, nu_t=0;
  PetscInt    twoD=0, Sutherland;
  PetscOptionsGetInt(NULL, NULL, "-twoD", &twoD, NULL);
  PetscOptionsGetInt(NULL, NULL, "-visT", &Sutherland, NULL);
  PetscOptionsGetReal(NULL, NULL,"-Tref", &Tref, NULL);

  DMDAVecGetArray(cda, user->lQ, &q);
  //  DMDAVecGetArray(fda, Ucont, &ucont);
  DMDAVecGetArray(cda, Visc,  &visc);
  
  DMDAVecGetArray(fda, Csi, &csi);
  DMDAVecGetArray(fda, Eta, &eta);
  DMDAVecGetArray(fda, Zet, &zet);
  
  DMDAVecGetArray(da, user->lNvert, &nvert);
  
  VecDuplicate(user->lQ, &Fp1);
  VecDuplicate(user->lQ, &Fp2);
  VecDuplicate(user->lQ, &Fp3);

  DMDAVecGetArray(cda, Fp1, &fp1);
  DMDAVecGetArray(cda, Fp2, &fp2);
  DMDAVecGetArray(cda, Fp3, &fp3);
  
  if (les) {
  VecDuplicate(user->lUcat, &Alpha1);
  VecDuplicate(user->lUcat, &Alpha2);
  VecDuplicate(user->lUcat, &Alpha3);
  
  DMDAVecGetArray(fda, Alpha1, &alpha1);
  DMDAVecGetArray(fda, Alpha2, &alpha2);
  DMDAVecGetArray(fda, Alpha3, &alpha3);
  DMDAVecGetArray(fda, user->lSigma, &lsigma);
  }


  DMDAVecGetArray(da, user->lAj, &aj);
  
  //  DMDAGetLocalInfo(da, &info);
  
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  
  /* First we calculate the flux on cell surfaces. Stored on the upper integer
     node. For example, along i direction, the flux are stored at node 0:mx-2*/
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  
  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;
  

  
  PetscReal ***lnu_t, nut_norm;

  VecNorm(user->Nu_t,NORM_INFINITY,&nut_norm);
  PetscPrintf(PETSC_COMM_WORLD,"Max Nu_t Set to %f\n", &nut_norm);
  
  if(les || rans) {
    DMDAVecGetArray(da, user->lNu_t, &lnu_t);
  }
  
  /* DMDAVecGetArray(da, user->lIAj, &iaj); */

  DMDAVecGetArray(da, user->lP, &p);
  
  DMDAVecGetArray(fda, user->lUcat,  &ucat);
  /* The visc flux on each surface center is stored at previous integer node */
  // i direction
  if (twoD!=1) {
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs-1; i<lxe; i++) {
	// ∂ui/∂ξ --> cartesian velocity gradient in curvilinear
	dudc = ucat[k][j][i+1].x - ucat[k][j][i].x;
	dvdc = ucat[k][j][i+1].y - ucat[k][j][i].y;
	dwdc = ucat[k][j][i+1].z - ucat[k][j][i].z;
	dpdc = p[k][j][i+1]/q[k][j][i+1].rho - p[k][j][i]/q[k][j][i].rho;
	
	// ∂ui/∂η --> cartesian velocity gradient in curvilinear
	if ((nvert[k][j+1][i  ]> solid && nvert[k][j+1][i  ]<innerblank)  ||
	    (nvert[k][j+1][i+1]> solid && nvert[k][j+1][i+1]<innerblank)) {
	  dude = (ucat[k][j  ][i+1].x + ucat[k][j  ][i].x -
		  ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.5;
	  dvde = (ucat[k][j  ][i+1].y + ucat[k][j  ][i].y -
		  ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.5;
	  dwde = (ucat[k][j  ][i+1].z + ucat[k][j  ][i].z -
		  ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.5;
	  dpde = (p[k][j  ][i+1]/q[k][j][i+1].rho + p[k][j  ][i]/q[k][j][i].rho -
		  p[k][j-1][i+1]/q[k][j-1][i+1].rho - p[k][j-1][i]/q[k][j-1][i].rho) * 0.5;
	}
	else if  ((nvert[k][j-1][i  ]> solid && nvert[k][j-1][i  ]<innerblank)  ||
		  (nvert[k][j-1][i+1]> solid && nvert[k][j-1][i+1]<innerblank)) {
	  dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x -
		  ucat[k][j  ][i+1].x - ucat[k][j  ][i].x) * 0.5;
	  dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y -
		  ucat[k][j  ][i+1].y - ucat[k][j  ][i].y) * 0.5;
	  dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z -
		  ucat[k][j  ][i+1].z - ucat[k][j  ][i].z) * 0.5;
	  dpde = (p[k][j+1][i+1]/q[k][j+1][i+1].rho + p[k][j+1][i]/q[k][j+1][i].rho -
		  p[k][j  ][i+1]/q[k][j  ][i+1].rho - p[k][j  ][i]/q[k][j  ][i].rho) * 0.5;
	}
	else {
	  dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x -
		  ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.25;
	  dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y -
		  ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.25;
	  dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z -
		  ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.25;
	  dpde =  (p[k][j+1][i+1]/q[k][j+1][i+1].rho + p[k][j+1][i]/q[k][j+1][i].rho -
		   p[k][j-1][i+1]/q[k][j-1][i+1].rho - p[k][j-1][i]/q[k][j-1][i].rho) * 0.25; 
	  
	}

	// ∂ui/∂ζ --> cartesian velocity gradient in curvilinear	
	if ((nvert[k+1][j][i  ]> solid && nvert[k+1][j][i  ]<innerblank)||
	    (nvert[k+1][j][i+1]> solid && nvert[k+1][j][i+1]<innerblank)) {
	  dudz = (ucat[k  ][j][i+1].x + ucat[k  ][j][i].x -
		  ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.5;
	  dvdz = (ucat[k  ][j][i+1].y + ucat[k  ][j][i].y -
		  ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.5;
	  dwdz = (ucat[k  ][j][i+1].z + ucat[k  ][j][i].z -
		  ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.5;
	  dpdz = (ucat[k  ][j][i+1].z + ucat[k  ][j][i].z -
		  ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.5;
	  
	  dpdz = (p[k][j][i+1]/q[k][j][i+1].rho + p[k][j][i]/q[k][j][i].rho -
		  p[k-1][j][i+1]/q[k-1][j][i+1].rho - p[k][j][i]/q[k][j][i].rho) * 0.5;		
	}
	else if ((nvert[k-1][j][i  ]> solid && nvert[k-1][j][i  ]<innerblank) ||
		 (nvert[k-1][j][i+1]> solid && nvert[k-1][j][i+1]<innerblank)) {
	  
	  dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x -
		  ucat[k  ][j][i+1].x - ucat[k  ][j][i].x) * 0.5;
	  dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y -
		  ucat[k  ][j][i+1].y - ucat[k  ][j][i].y) * 0.5;
	  dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z -
		  ucat[k  ][j][i+1].z - ucat[k  ][j][i].z) * 0.5;
	  dpdz = (p[k+1][j][i+1]/q[k+1][j][i+1].rho + p[k+1][j][i]/q[k+1][j][i].rho -
		  p[k  ][j][i+1]/q[k  ][j][i+1].rho - p[k  ][j][i]/q[k  ][j][i].rho) * 0.5;
	}
	else {
	  dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x -
		  ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.25;
	  dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y -
		  ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.25;
	  dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z -
		  ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.25;
	  dpdz = (p[k+1][j][i+1]/q[k+1][j][i+1].rho + p[k+1][j][i]/q[k+1][j][i].rho -
		  p[k-1][j][i+1]/q[k-1][j][i+1].rho - p[k-1][j][i]/q[k-1][j][i].rho) * 0.25;
	  
	}
	
	// metrics at i+1/2 (All these quantities have 1/J multiplied to it)
	// ∂ξ/∂x, ∂ξ/∂y, ∂ξ/∂z
	csi0 = 0.5*(csi[k][j][i].x+csi[k][j][i+1].x);
	csi1 = 0.5*(csi[k][j][i].y+csi[k][j][i+1].y);
	csi2 = 0.5*(csi[k][j][i].z+csi[k][j][i+1].z);
	
	// ∂η/∂x, ∂η/∂y, ∂η/∂z
	eta0 = 0.5*(eta[k][j][i].x+eta[k][j][i+1].x);
	eta1 = 0.5*(eta[k][j][i].y+eta[k][j][i+1].y);
	eta2 = 0.5*(eta[k][j][i].z+eta[k][j][i+1].z);
	
	// ∂ζ/∂x, ∂ζ/∂y, ∂ζ/∂z
	zet0 = 0.5*(zet[k][j][i].x+zet[k][j][i+1].x);
	zet1 = 0.5*(zet[k][j][i].y+zet[k][j][i+1].y);
	zet2 = 0.5*(zet[k][j][i].z+zet[k][j][i+1].z);
	
	/* Calculation of (2 mu) ∂ξq / ∂xr * S_mr : m = 1-> rhoU,2->rhoV,3->rhoW 
		S_mr = (1/2) * ([∂ξp / ∂xm * ∂ur / ∂ξp] + [∂ξp / ∂xr * ∂um / ∂ξp])
		substitute and expand:
		|∂ξq/∂xr * ∂ξp/∂xm * ∂ur/∂ξp| + |∂ξq/∂xr * ∂ξp/∂xr * ∂um/∂ξp|
		Remember that q is the index for flux direction, for this loop, q==1
		This term:
		∂ξq/∂xr * ∂ξp/∂xr
		Are the metrics of the transformation;
		The other term:
		∂ξp/∂xm * ∂ur/∂ξp
		calculate these terms and multiply with respective values
	*/
	// ∂ξr / ∂xj * ∂ξm / ∂xj
	g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;
	
	// ∂u/∂ξi* ∂ξi/∂xj
	r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
	r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
	r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;
	
	r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
	r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
	r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;
	
	r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
	r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
	r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;
	
	ajc = 0.5*(aj[k][j][i]+aj[k][j][i+1]); 
	
	 /* note that viscosity mu is non-dim by velocity scale U=sqrt(P_R/rho_R) 
	    , L, and density scale rho_R. P_R, rho_R, L, and U are all scales (consts). 
	    There is no need to multiply by the density as var. */
	
	// for les-> the viscosity is normalized with q
	if(les)nu = 1./user->ren*0.5*(q[k][j][i].rho+q[k][j][i+1].rho);
	
	//	double nu = 1./user->ren*0.5*(q[k][j][i].rho+q[k][j][i+1].rho), nu_t=0;
	if (Sutherland) {
	  //Tr should be 1 for free stream to get exactly 1/Re
	  // nu=1/Re * (1+C/Tref)/(T/Tref+C/Tref)(T/Tref)^3/2
	  Tr=user->gamma*user->Ma*user->Ma*
	    (p[k][j][i]+p[k][j][i+1])/(q[k][j][i].rho+q[k][j][i+1].rho); 
	  if (Tr<0) {
	    PetscPrintf(PETSC_COMM_WORLD, "!! k Tr<0 Sutherland %d Tr %f\n", Sutherland, Tr);
	    nu=1./user->ren;
	  } else 
	    nu=1./user->ren*(1+C_sutherland/Tref)/(Tr+C_sutherland/Tref)*pow(Tr, 1.5);
	}

	/* This is S̃_pp = (∂ξ/∂x)(∂ũ/∂ξ) + (∂η/∂x)(∂ũ/∂η) + (∂ζ/∂x)(∂ũ/∂ζ)
    				 + (∂ξ/∂y)(∂ṽ/∂ξ) + (∂η/∂y)(∂ṽ/∂η) + (∂ζ/∂y)(∂ṽ/∂ζ)
     					+ (∂ξ/∂z)(∂ŵ/∂ξ) + (∂η/∂z)(∂ŵ/∂η) + (∂ζ/∂z)(∂ŵ/∂ζ)
		Which appears in every CompVar term once as divergence */
	div = 2./3.*(dudc * csi0 + dude * eta0 + dudz * zet0 +
		     dvdc * csi1 + dvde * eta1 + dvdz * zet1 +
		     dwdc * csi2 + dwde * eta2 + dwdz * zet2 );
	
	fp1[k][j][i].rhoU = (g11 * dudc + g21 * dude + g31 * dudz+ r11 * csi0 + r21 * csi1 + r31 * csi2 ) * ajc * (nu);
	fp1[k][j][i].rhoV = (g11 * dvdc + g21 * dvde + g31 * dvdz+ r12 * csi0 + r22 * csi1 + r32 * csi2 ) * ajc * (nu);
	fp1[k][j][i].rhoW = (g11 * dwdc + g21 * dwde + g31 * dwdz+ r13 * csi0 + r23 * csi1 + r33 * csi2 ) * ajc * (nu);
	
	// div= 2/3 div(u) 
	// compressiblity component of tau_ij= -2/3 div(u) delta_ij
	// fpk .i -= div zeta^k_x_i 1/aj (nu+nu_t)
	fp1[k][j][i].rhoU -= div * csi0 * ajc * nu ;
	fp1[k][j][i].rhoV -= div * csi1 * ajc * nu ;
	fp1[k][j][i].rhoW -= div * csi2 * ajc * nu ;
	
	if (i==0 && user->bctype[0]==1)  // on the boundary;
	  ucat_h=ucat[k][j][i]; 
	else if (i==mx-2 && user->bctype[1]==1)
	  ucat_h=ucat[k][j][i+1];
	else {
	  ucat_h.x=0.5*(ucat[k][j][i+1].x + ucat[k][j][i].x);
	  ucat_h.y=0.5*(ucat[k][j][i+1].y + ucat[k][j][i].y);
	  ucat_h.z=0.5*(ucat[k][j][i+1].z + ucat[k][j][i].z);
	}

	// Calculating Kr: for some reason this is not calculated properly in the old code!!
	fp1[k][j][i].rhoE  =  (fp1[k][j][i].rhoU * ucat_h.x +
			       fp1[k][j][i].rhoV * ucat_h.y +
			       fp1[k][j][i].rhoW * ucat_h.z ) 
	  + (g11*dpdc+g21*dpde+g31*dpdz) *ajc *(gamma*nu/pr/(gamma-1.0)); //changed from -0.5*( --- iman
	
	
	if( les || rans) {//(rans && ti>0) ) {
	  //nu_t = pow( 0.5 * ( sqrt(lnu_t[k][j][i]) + sqrt(lnu_t[k][j][i+1]) ), 2.0) * Sabs;
	
		// Advait 
	  nu_t =  0.5 * (lnu_t[k][j][i]*lsigma[k][j][i].x + lnu_t[k][j][i+1]*lsigma[k][j][i+1].x);
	  //nu_t = 0.5 * (lnu_t[k][j][i]*lsigma[k][j][i].x + lnu_t[k][j][i+1]*lsigma[k][j][i+1].x);// sigma to calculate nu_t?
	  rho_c = 0.5* (q[k][j][i].rho+q[k][j][i+1].rho);
	  sigma_c = 0.5*(lsigma[k][j][i].x + lsigma[k][j][i+1].x);
	  if ( (user->bctype[0]==1 && i==0) || (user->bctype[1]==1 && i==mx-2) ) nu_t=0;    
		// t is the same as -2mut(Srm -1/3 Spp delrm) stress, S2 is sqrt(SrmSmr) and Ss is the Subgrid Kinetic Energy	  
	  t.x=  (g11 * dudc + g21 * dude + g31 * dudz + 
		 r11 * csi0 + r21 * csi1 + r31 * csi2 -div *csi0) * ajc * (nu_t*rho_c); 
	  t.y=  (g11 * dvdc + g21 * dvde + g31 * dvdz + 
		 r12 * csi0 + r22 * csi1 + r32 * csi2 -div *csi1) * ajc * (nu_t*rho_c);
	  t.z=  (g11 * dwdc + g21 * dwde + g31 * dwdz + 
		 r13 * csi0 + r23 * csi1 + r33 * csi2 -div *csi2) * ajc * (nu_t*rho_c);
	  S2=.5*((r11-div)*(r11-div)+r12*r12+r13*r13+r21*r21+
		 (r22-div)*(r22-div)+r23*r23+r31*r31+r32*r32+(r33-div)*(r33-div));
	  double filter  = pow( 1./aj[k][j][i],1./3.);
	  
	  Ss.x=-  S2*pow(filter,2.0)*CI*2.0/3.0*csi0*rho_c*sigma_c ;
	  Ss.y=-  S2*pow(filter,2.0)*CI*2.0/3.0*csi1*rho_c*sigma_c ;
	  Ss.z=-  S2*pow(filter,2.0)*CI*2.0/3.0*csi2*rho_c*sigma_c ;	  	  
	  
	  fp1[k][j][i].rhoU += t.x+Ss.x;
	  fp1[k][j][i].rhoV += t.y+Ss.y;
	  fp1[k][j][i].rhoW += t.z+Ss.z;
	  	  
	  // Not exactly a part of viscous flux, it is just Calculated for RHS
	  // Subgrid stress Calculated Separately?
	  alpha1[k][j][i].x=(t.x+Ss.x);		
	  alpha1[k][j][i].y=(t.y+Ss.y);		
	  alpha1[k][j][i].z=(t.z+Ss.z);			  	 	  
	  
	}
	
      }
    }
  }
  } // if !twoD
  /* DMDAVecRestoreArray(da, user->lIAj, &iaj);   */
  
  
  // j direction
  /* DMDAVecGetArray(da, user->lJAj, &jaj); */
  for (k=lzs; k<lze; k++) {
    for (j=lys-1; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	
	if ((nvert[k][j  ][i+1]> solid && nvert[k][j  ][i+1]<innerblank)||
	    (nvert[k][j+1][i+1]> solid && nvert[k][j+1][i+1]<innerblank)) {
	  dudc = (ucat[k][j+1][i  ].x + ucat[k][j][i  ].x -
		  ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.5;
	  dvdc = (ucat[k][j+1][i  ].y + ucat[k][j][i  ].y -
		  ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.5;
	  dwdc = (ucat[k][j+1][i  ].z + ucat[k][j][i  ].z -
		  ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.5;
	  dpdc = (p[k][j+1][i  ]/q[k][j+1][i  ].rho + p[k][j][i  ]/q[k][j][i  ].rho -
		  p[k][j+1][i-1]/q[k][j+1][i-1  ].rho - p[k][j][i-1]/q[k][j][i-1].rho) * 0.5;
	}
	else if ((nvert[k][j  ][i-1]> solid && nvert[k][j  ][i-1]<innerblank) ||
		 (nvert[k][j+1][i-1]> solid && nvert[k][j+1][i-1]<innerblank)) {
	  dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x -
		  ucat[k][j+1][i  ].x - ucat[k][j][i  ].x) * 0.5;
	  dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y -
		  ucat[k][j+1][i  ].y - ucat[k][j][i  ].y) * 0.5;
	  dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z -
		  ucat[k][j+1][i  ].z - ucat[k][j][i  ].z) * 0.5;
	  dpdc = (p[k][j+1][i+1]/q[k][j+1][i+1].rho + p[k][j][i+1]/q[k][j][i+1].rho -
		  p[k][j+1][i  ]/q[k][j+1][i].rho - p[k][j][i  ]/q[k][j][i].rho) * 0.5;
	  
	}
	else {
	  dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x -
		  ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.25;
	  dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y -
		  ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.25;
	  dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z -
		  ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.25;
	  dpdc = (p[k][j+1][i+1]/q[k][j+1][i+1].rho + p[k][j][i+1]/q[k][j][i+1].rho -
		  p[k][j+1][i-1]/q[k][j+1][i-1].rho - p[k][j][i-1]/q[k][j][i-1].rho) * 0.25;
	}
	
	dude = ucat[k][j+1][i].x - ucat[k][j][i].x;
	dvde = ucat[k][j+1][i].y - ucat[k][j][i].y;
	dwde = ucat[k][j+1][i].z - ucat[k][j][i].z;
	dpde = p[k][j+1][i]/q[k][j+1][i].rho - p[k][j][i]/q[k][j][i].rho;
	
	if ((nvert[k+1][j  ][i]> solid && nvert[k+1][j  ][i]<innerblank)||
	    (nvert[k+1][j+1][i]> solid && nvert[k+1][j+1][i]<innerblank)) {
	  dudz = (ucat[k  ][j+1][i].x + ucat[k  ][j][i].x -
		  ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.5;
	  dvdz = (ucat[k  ][j+1][i].y + ucat[k  ][j][i].y -
		  ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.5;
	  dwdz = (ucat[k  ][j+1][i].z + ucat[k  ][j][i].z -
		  ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.5;
	  dpdz = (p[k  ][j+1][i]/q[k  ][j+1][i].rho + p[k  ][j][i]/q[k  ][j][i].rho -
		  p[k-1][j+1][i]/q[k-1][j+1][i].rho - p[k-1][j][i]/q[k-1][j][i].rho) * 0.5;
	}
	else if ((nvert[k-1][j  ][i]> solid && nvert[k-1][j  ][i]<innerblank)||
		 (nvert[k-1][j+1][i]> solid && nvert[k-1][j+1][i]<innerblank)) {
	  dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x -
		  ucat[k  ][j+1][i].x - ucat[k  ][j][i].x) * 0.5;
	  dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y -
		  ucat[k  ][j+1][i].y - ucat[k  ][j][i].y) * 0.5;
	  dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z -
		  ucat[k  ][j+1][i].z - ucat[k  ][j][i].z) * 0.5;
	  dpdz = (p[k+1][j+1][i]/q[k+1][j+1][i].rho + p[k+1][j][i]/q[k+1][j][i].rho -
		  p[k  ][j+1][i]/q[k  ][j+1][i].rho - p[k  ][j][i]/q[k  ][j][i].rho) * 0.5;
	  
	}
	else {
	  dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x -
		  ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.25;
	  dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y -
		  ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.25;
	  dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z -
		  ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.25;
	  dpdz = (p[k+1][j+1][i]/q[k+1][j+1][i].rho + p[k+1][j][i]/q[k+1][j][i].rho -
		  p[k-1][j+1][i]/q[k-1][j+1][i].rho - p[k-1][j][i]/q[k-1][j][i].rho) * 0.25;
	  
	}
	
	csi0 = 0.5*(csi[k][j][i].x+csi[k][j+1][i].x);
	csi1 = 0.5*(csi[k][j][i].y+csi[k][j+1][i].y); 
	csi2 = 0.5*(csi[k][j][i].z+csi[k][j+1][i].z);
	
	eta0 = 0.5*(eta[k][j][i].x+eta[k][j+1][i].x);
	eta1 = 0.5*(eta[k][j][i].y+eta[k][j+1][i].y); 
	eta2 = 0.5*(eta[k][j][i].z+eta[k][j+1][i].z);
	
	zet0 = 0.5*(zet[k][j][i].x+zet[k][j+1][i].x);
	zet1 = 0.5*(zet[k][j][i].y+zet[k][j+1][i].y); 
	zet2 = 0.5*(zet[k][j][i].z+zet[k][j+1][i].z);

	g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
	g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;
	
	r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
	r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
	r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;
	
	r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
	r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
	r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;
	
	r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
	r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
	r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;
	
	//	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d dvdc is %.15le dvde is %.15le dvdz is %.15le \n",i,j,k,dvdc,dvde,dvdz);
				//	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d dwdc is %.15le dwde is %.15le dwdz is %.15le \n",i,j,k,dwdc,dwde,dwdz);
				//	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d jcsi is %.15le jeta is %.15le jzet is %.15le \n",i,j,k,jcsi[k][j][i].z,jeta[k][j][i].z,jzet[k][j][i].z);
				//	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d r13 is %.15le r23 is %.15le r33 is %.15le \n",i,j,k,r13,r23,r33);

	ajc = 0.5*(aj[k][j][i]+aj[k][j+1][i]);
	
	//	double nu = 1./user->ren*0.5*(q[k][j][i].rho+q[k][j+1][i].rho), nu_t = 0;

	if(les)nu = 1./user->ren*0.5*(q[k][j][i].rho+q[k][j][i+1].rho);

	if (Sutherland) {
	  //Tr should be 1 for free stream to get exactly 1/Re
	  // nu=1/Re * (1+C/Tref)/(T/Tref+C/Tref)(T/Tref)^3/2
	  Tr=user->gamma*user->Ma*user->Ma*
	    (p[k][j][i]+p[k][j+1][i])/(q[k][j][i].rho+q[k][j+1][i].rho); 
	  if (Tr<0) {
	    PetscPrintf(PETSC_COMM_WORLD, "!!Tr<0 Sutherland %d Tr %f\n", Sutherland, Tr);
	    nu=1./user->ren;
	  } else 
	    nu=1./user->ren*(1+C_sutherland/Tref)/(Tr+C_sutherland/Tref)*pow(Tr, 1.5);
	}

	div = 2./3.*(dudc * csi0 + dude * eta0 + dudz * zet0 +
		     dvdc * csi1 + dvde * eta1 + dvdz * zet1 +
		     dwdc * csi2 + dwde * eta2 + dwdz * zet2 );
	
	fp2[k][j][i].rhoU = (g11 * dudc + g21 * dude + g31 * dudz+ 
			     r11 * eta0 + r21 * eta1 + r31 * eta2 ) * ajc * (nu);
	fp2[k][j][i].rhoV = (g11 * dvdc + g21 * dvde + g31 * dvdz+ 
			     r12 * eta0 + r22 * eta1 + r32 * eta2 ) * ajc * (nu);
	fp2[k][j][i].rhoW = (g11 * dwdc + g21 * dwde + g31 * dwdz+ 
			     r13 * eta0 + r23 * eta1 + r33 * eta2 ) * ajc * (nu);

	// div= 2/3 div(u) 
	// compressiblity component of tau_ij= -2/3 div(u) delta_ij
	// fpk .i -= div zeta^k_x_i 1/aj (nu+nu_t)
	fp2[k][j][i].rhoU -= div * eta0 * ajc * nu ;
	fp2[k][j][i].rhoV -= div * eta1 * ajc * nu ;
	fp2[k][j][i].rhoW -= div * eta2 * ajc * nu ;

	if (j==0 && user->bctype[2]==1)  // on the boundary;
	  ucat_h=ucat[k][j][i]; 
	else if (j==my-2 && user->bctype[3]==1)
	  ucat_h=ucat[k][j+1][i];
	else {
	  ucat_h.x=0.5*(ucat[k][j+1][i].x + ucat[k][j][i].x);
	  ucat_h.y=0.5*(ucat[k][j+1][i].y + ucat[k][j][i].y);
	  ucat_h.z=0.5*(ucat[k][j+1][i].z + ucat[k][j][i].z);
	}
	
	fp2[k][j][i].rhoE  =  (fp2[k][j][i].rhoU * ucat_h.x  +
			       fp2[k][j][i].rhoV * ucat_h.y  +
			       fp2[k][j][i].rhoW * ucat_h.z) +
	  (g11*dpdc+g21*dpde+g31*dpdz) *ajc *(gamma*nu/pr/(gamma-1.0));
	
	if( les || rans) {// (rans && ti>0) ) {
	  //nu_t = pow( 0.5 * ( sqrt(lnu_t[k][j][i]) + sqrt(lnu_t[k][j+1][i]) ), 2.0) * Sabs;
	 // j flux multiplied by q? whyy?? Advait 
	  //nu_t = 0.5 * (lnu_t[k][j][i]* q[k][j][i].rho+ lnu_t[k][j+1][i]*q[k][j+1][i].rho); 
	  nu_t = 0.5 * (lnu_t[k][j][i]*lsigma[k][j][i].y+ lnu_t[k][j+1][i]*lsigma[k][j+1][i].y);
	  rho_c = 0.5* (q[k][j][i].rho+q[k][j+1][i].rho); 
	  sigma_c=0.5*(lsigma[k][j][i].y + lsigma[k][j+1][i].y);
	  if ( (user->bctype[2]==1 && j==0) || (user->bctype[3]==1 && j==my-2) ) nu_t=0;
	  
	  
	  t.x=  (g11 * dudc + g21 * dude + g31 * dudz + 
		 r11 * eta0 + r21 * eta1 + r31 * eta2 -div *eta0) * ajc * (nu_t*rho_c); 
	  t.y=  (g11 * dvdc + g21 * dvde + g31 * dvdz + 
		 r12 * eta0 + r22 * eta1 + r32 * eta2 -div *eta1) * ajc * (nu_t*rho_c);
	  t.z=  (g11 * dwdc + g21 * dwde + g31 * dwdz + 
		 r13 * eta0 + r23 * eta1 + r33 * eta2 -div *eta2) * ajc * (nu_t*rho_c);
	  S2=.5*((r11-div)*(r11-div)+r12*r12+r13*r13+r21*r21+
		 (r22-div)*(r22-div)+r23*r23+r31*r31+r32*r32+(r33-div)*(r33-div));
	  double filter  = pow( 1./aj[k][j][i],1./3.);
	  
	  Ss.x=-  S2*pow(filter,2.0)*CI*2.0/3.0*eta0*rho_c*sigma_c ;
	  Ss.y=-  S2*pow(filter,2.0)*CI*2.0/3.0*eta1*rho_c*sigma_c ;
	  Ss.z=-  S2*pow(filter,2.0)*CI*2.0/3.0*eta2*rho_c*sigma_c ;	  	  
	  
	  fp2[k][j][i].rhoU += t.x+Ss.x;
	  fp2[k][j][i].rhoV += t.y+Ss.y;
	  fp2[k][j][i].rhoW += t.z+Ss.z;	  
	  
	  alpha2[k][j][i].x=(t.x+Ss.x);		
	  alpha2[k][j][i].y=(t.y+Ss.y);		
	  alpha2[k][j][i].z=(t.z+Ss.z);			  	  
	}					
	

					
      }
    }
  }

  /* DMDAVecRestoreArray(da, user->lJAj, &jaj); */
  // k direction
  
  /* DMDAVecGetArray(da, user->lKAj, &kaj); */
  for (k=lzs-1; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((nvert[k  ][j][i+1]> solid && nvert[k  ][j][i+1]<innerblank)||
	    (nvert[k+1][j][i+1]> solid && nvert[k+1][j][i+1]<innerblank)) {
	  dudc = (ucat[k+1][j][i  ].x + ucat[k][j][i  ].x -
		  ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.5;
	  dvdc = (ucat[k+1][j][i  ].y + ucat[k][j][i  ].y -
		  ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.5;
	  dwdc = (ucat[k+1][j][i  ].z + ucat[k][j][i  ].z -
		  ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.5;
	  dpdc = (p[k+1][j][i  ]/q[k+1][j][i  ].rho + p[k][j][i  ]/q[k][j][i  ].rho -
		  p[k+1][j][i-1]/q[k+1][j][i-1].rho - p[k][j][i-1]/q[k][j][i-1].rho) * 0.5;
	}
	else if ((nvert[k  ][j][i-1]> solid && nvert[k  ][j][i-1]<innerblank) ||
		 (nvert[k+1][j][i-1]> solid && nvert[k+1][j][i-1]<innerblank)) {
	  dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x -
		  ucat[k+1][j][i  ].x - ucat[k][j][i  ].x) * 0.5;
	  dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y -
		  ucat[k+1][j][i  ].y - ucat[k][j][i  ].y) * 0.5;
	  dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z -
		  ucat[k+1][j][i  ].z - ucat[k][j][i  ].z) * 0.5;
	  dpdc = (p[k+1][j][i+1]/q[k+1][j][i+1].rho + p[k][j][i+1]/q[k][j][i+1].rho -
		  p[k+1][j][i  ]/q[k+1][j][i  ].rho - p[k][j][i  ]/q[k][j][i  ].rho) * 0.5;
	}
	else {
	  dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x -
		  ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.25;
	  dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y -
		  ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.25;
	  dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z -
		  ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.25;
	  dpdc = (p[k+1][j][i+1]/q[k+1][j][i+1].rho + p[k][j][i+1]/q[k][j][i+1].rho -
		  p[k+1][j][i-1]/q[k+1][j][i-1].rho - p[k][j][i-1]/q[k][j][i-1].rho) * 0.25;
	}
	
	if ((nvert[k  ][j+1][i]> solid && nvert[k  ][j+1][i]<innerblank)||
	    (nvert[k+1][j+1][i]> solid && nvert[k+1][j+1][i]<innerblank)) {
	  dude = (ucat[k+1][j  ][i].x + ucat[k][j  ][i].x -
		  ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.5;
	  dvde = (ucat[k+1][j  ][i].y + ucat[k][j  ][i].y -
		  ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.5;
	  dwde = (ucat[k+1][j  ][i].z + ucat[k][j  ][i].z -
		  ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.5;
	  dpde = (p[k+1][j  ][i]/q[k+1][j  ][i].rho + p[k][j  ][i]/q[k][j  ][i].rho -
		  p[k+1][j-1][i]/q[k+1][j-1][i].rho - p[k][j-1][i]/q[k][j-1][i].rho) * 0.5;
	}
	else if ((nvert[k  ][j-1][i]> solid && nvert[k  ][j-1][i]<innerblank) ||
		 (nvert[k+1][j-1][i]> solid && nvert[k+1][j-1][i]<innerblank)){
	  dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x -
		  ucat[k+1][j  ][i].x - ucat[k][j  ][i].x) * 0.5;
	  dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y -
		  ucat[k+1][j  ][i].y - ucat[k][j  ][i].y) * 0.5;
	  dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z -
		  ucat[k+1][j  ][i].z - ucat[k][j  ][i].z) * 0.5;
	  dpde = (p[k+1][j+1][i]/q[k+1][j+1][i].rho + p[k][j+1][i]/q[k][j+1][i].rho -
		  p[k+1][j  ][i]/q[k+1][j  ][i].rho - p[k][j  ][i]/q[k][j  ][i].rho) * 0.5;
	}
	else {
	  dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x -
		  ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.25;
	  dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y -
		  ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.25;
	  dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z -
		  ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.25;
	  dpde = (p[k+1][j+1][i]/q[k+1][j+1][i].rho + p[k][j+1][i]/q[k][j+1][i].rho -
		  p[k+1][j-1][i]/q[k+1][j-1][i].rho - p[k][j-1][i]/q[k][j-1][i].rho) * 0.25;
	  
	}
	
	dudz = ucat[k+1][j][i].x - ucat[k][j][i].x;
	dvdz = ucat[k+1][j][i].y - ucat[k][j][i].y;
	dwdz = ucat[k+1][j][i].z - ucat[k][j][i].z;
	dpdz = p[k+1][j][i]/q[k+1][j][i].rho - p[k][j][i]/q[k][j][i].rho;
	
	
	csi0 = 0.5*(csi[k][j][i].x+csi[k+1][j][i].x); 
	csi1 = 0.5*(csi[k][j][i].y+csi[k+1][j][i].y);
	csi2 = 0.5*(csi[k][j][i].z+csi[k+1][j][i].z);
	
	eta0 = 0.5*(eta[k][j][i].x+eta[k+1][j][i].x); 
	eta1 = 0.5*(eta[k][j][i].y+eta[k+1][j][i].y);
	eta2 = 0.5*(eta[k][j][i].z+eta[k+1][j][i].z);
	
	zet0 = 0.5*(zet[k][j][i].x+zet[k+1][j][i].x); 
	zet1 = 0.5*(zet[k][j][i].y+zet[k+1][j][i].y);
	zet2 = 0.5*(zet[k][j][i].z+zet[k+1][j][i].z);
	
	
	g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
	g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
	g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;
	
	r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
	r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
	r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;
	
	r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
	r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
	r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;
	
	r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
	r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
	r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;
	
	ajc = 0.5*(aj[k][j][i]+aj[k+1][j][i]);

	//double nu = 1./user->ren*0.5*(q[k][j][i].rho+q[k+1][j][i].rho), nu_t =0;

	if(les)nu = 1./user->ren*0.5*(q[k][j][i].rho+q[k][j][i+1].rho);

	if (Sutherland) {
	  //Tr=p/rho should be 1 for free stream to get exactly 1/Re
	  // nu=1/Re * (1+C/Tref)/(T/Tref+C/Tref)(T/Tref)^3/2
	  Tr=user->gamma*user->Ma*user->Ma*
	    (p[k][j][i]+p[k+1][j][i])/(q[k][j][i].rho+q[k+1][j][i].rho); 
	  if (Tr<0) {
	    PetscPrintf(PETSC_COMM_WORLD, "!! k Tr<0 Sutherland %d Tr %f\n", Sutherland, Tr);
	    nu=1./user->ren;
	  } else 
	    nu=1./user->ren*(1+C_sutherland/Tref)/(Tr+C_sutherland/Tref)*pow(Tr, 1.5);
	}

	div = 2./3.*(dudc * csi0 + dude * eta0 + dudz * zet0 +
		     dvdc * csi1 + dvde * eta1 + dvdz * zet1 +
		     dwdc * csi2 + dwde * eta2 + dwdz * zet2 );
	fp3[k][j][i].rhoU += (g11 * dudc + g21 * dude + g31 * dudz + 
			      r11 * zet0 + r21 * zet1 + r31 * zet2) * ajc * (nu);//
	fp3[k][j][i].rhoV += (g11 * dvdc + g21 * dvde + g31 * dvdz + 
			      r12 * zet0 + r22 * zet1 + r32 * zet2) * ajc * (nu);//
	fp3[k][j][i].rhoW += (g11 * dwdc + g21 * dwde + g31 * dwdz + 
			      r13 * zet0 + r23 * zet1 + r33 * zet2) * ajc * (nu);//

	// div= 2/3 div(u) 
	// compressiblity component of tau_ij= -2/3 div(u) delta_ij
	// fpk .i -= div zeta^k_x_i 1/aj (nu+nu_t)
	fp3[k][j][i].rhoU -= div * zet0 * ajc * nu ;
	fp3[k][j][i].rhoV -= div * zet1 * ajc * nu ;
	fp3[k][j][i].rhoW -= div * zet2 * ajc * nu ;

	if (k==0 && user->bctype[4]==1)  // on the boundary;
	  ucat_h=ucat[k][j][i]; 
	else if (k==mz-2 && user->bctype[5]==1)
	  ucat_h=ucat[k+1][j][i];
	else {
	  ucat_h.x=0.5*(ucat[k+1][j][i].x + ucat[k][j][i].x);
	  ucat_h.y=0.5*(ucat[k+1][j][i].y + ucat[k][j][i].y);
	  ucat_h.z=0.5*(ucat[k+1][j][i].z + ucat[k][j][i].z);
	}

	fp3[k][j][i].rhoE  =  (fp3[k][j][i].rhoU * ucat_h.x +
			       fp3[k][j][i].rhoV * ucat_h.y +
			       fp3[k][j][i].rhoW * ucat_h.z)+
	  (g11*dpdc+g21*dpde+g31*dpdz) *ajc *(gamma*nu/pr/(gamma-1.0));

	/* if ((k==1 || k==2) && j==1 && i==1) { */
	/*   PetscReal xx; */
	/*   PetscReal xx2; */
	/*   xx=(gamma*nu/pr/(gamma-1.0)); */
	/*   xx2=(g11*dpdc+g21*dpde+g31*dpdz) *ajc *(gamma*nu/pr/(gamma-1.0)); */
	/*   PetscPrintf(PETSC_COMM_WORLD, "!!Visc TGV rhoV %le rhoW %le gamma %le pr %le calc E %le %le aj %le\n",fp3[k][j][i].rhoV,fp3[k][j][i].rhoW, gamma, pr,xx, xx2, ajc); */
	/* } */

	if( les || rans) { //(rans && ti>0) ) {
	  //nu_t = pow( 0.5 * ( sqrt(lnu_t[k][j][i]) + sqrt(lnu_t[k+1][j][i]) ), 2.0) * Sabs;
	  //nu_t = 0.5 * (lnu_t[k][j][i]*q[k][j][i].rho + lnu_t[k+1][j][i]*q[k+1][j][i].rho);
	  nu_t = 0.5 * (lnu_t[k][j][i]*lsigma[k][j][i].z + lnu_t[k+1][j][i]*lsigma[k+1][j][i].z);
	  rho_c = 0.5* (q[k][j][i].rho+q[k+1][j][i].rho); 
	  sigma_c=0.5* (lsigma[k][j][i].z + lsigma[k+1][j][i].z);
	  if ( (user->bctype[4]==1 && k==0) || (user->bctype[5]==1 && k==mz-2) ) nu_t=0;
	  
	  t.x=  (g11 * dudc + g21 * dude + g31 * dudz + 
		 r11 * zet0 + r21 * zet1 + r31 * zet2 -div *zet0) * ajc * (nu_t*rho_c); 
	  t.y=  (g11 * dvdc + g21 * dvde + g31 * dvdz + 
		 r12 * zet0 + r22 * zet1 + r32 * zet2 -div *zet1) * ajc * (nu_t*rho_c);
	  t.z=  (g11 * dwdc + g21 * dwde + g31 * dwdz + 
		 r13 * zet0 + r23 * zet1 + r33 * zet2 -div *zet2) * ajc * (nu_t*rho_c);
	  S2=.5*((r11-div)*(r11-div)+r12*r12+r13*r13+r21*r21+
		 (r22-div)*(r22-div)+r23*r23+r31*r31+r32*r32+(r33-div)*(r33-div));
	  double filter  = pow( 1./aj[k][j][i],1./3.);
	  
	  Ss.x=-  S2*pow(filter,2.0)*CI*2.0/3.0*zet0*rho_c*sigma_c ;
	  Ss.y=-  S2*pow(filter,2.0)*CI*2.0/3.0*zet1*rho_c*sigma_c ;
	  Ss.z=-  S2*pow(filter,2.0)*CI*2.0/3.0*zet2*rho_c*sigma_c ;	  	  
	  
	  fp3[k][j][i].rhoU += t.x+Ss.x;
	  fp3[k][j][i].rhoV += t.y+Ss.y;
	  fp3[k][j][i].rhoW += t.z+Ss.z;
	  					
	  alpha3[k][j][i].x=(t.x+Ss.x);		
	  alpha3[k][j][i].y=(t.y+Ss.y);		
	  alpha3[k][j][i].z=(t.z+Ss.z);	  	  	  
	}
			
      }
    }
  }
  
  /* DMDAVecRestoreArray(da, user->lKAj, &kaj); */
  if (les){
	DMDAVecGetArray(da, user->lBBeta, &bbeta);
}

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {	
	visc[k][j][i].rho  = 0.;

	// Central Scheme for Viscous Fluxes
	if (twoD!=1)
	visc[k][j][i].rhoU =
	  (fp1[k][j][i].rhoU - fp1[k][j][i-1].rhoU +
	   fp2[k][j][i].rhoU - fp2[k][j-1][i].rhoU +
	   fp3[k][j][i].rhoU - fp3[k-1][j][i].rhoU );
	else visc[k][j][i].rhoU = 0.; 
	visc[k][j][i].rhoV =
	  (fp1[k][j][i].rhoV - fp1[k][j][i-1].rhoV +
	   fp2[k][j][i].rhoV - fp2[k][j-1][i].rhoV +
	   fp3[k][j][i].rhoV - fp3[k-1][j][i].rhoV );
	visc[k][j][i].rhoW =
	  (fp1[k][j][i].rhoW - fp1[k][j][i-1].rhoW +
	   fp2[k][j][i].rhoW - fp2[k][j-1][i].rhoW +
	   fp3[k][j][i].rhoW - fp3[k-1][j][i].rhoW );
	visc[k][j][i].rhoE =
	  (fp1[k][j][i].rhoE - fp1[k][j][i-1].rhoE +
	   fp2[k][j][i].rhoE - fp2[k][j-1][i].rhoE +
	   fp3[k][j][i].rhoE - fp3[k-1][j][i].rhoE );
	if (les && j>2 && j < my-2 && nvert[k][j][i]<0.5){
	  
	//   visc[k][j][i].rhoE += 1/q[k][j][i].rho*( // changed from -= to += ---iman        
	// 			 q[k][j][i].rhoU *
	// 			 (alpha1[k][j][i].x - alpha1[k][j][i-1].x +
	// 			  alpha2[k][j][i].x - alpha2[k][j-1][i].x +
	// 			  alpha3[k][j][i].x - alpha3[k-1][j][i].x ) +
	// 			 q[k][j][i].rhoV *
	// 			 (alpha1[k][j][i].y - alpha1[k][j][i-1].y +
	// 			  alpha2[k][j][i].y - alpha2[k][j-1][i].y +
	// 			  alpha3[k][j][i].y - alpha3[k-1][j][i].y )+
	// 			 q[k][j][i].rhoW *
	// 			 (alpha1[k][j][i].z - alpha1[k][j][i-1].z +
	// 			  alpha2[k][j][i].z - alpha2[k][j-1][i].z +
	// 			  alpha3[k][j][i].z - alpha3[k-1][j][i].z));

	// Advait Added
	/* Calculate Metrics at Edges of All the blocks.
		*Notation: csi_faces-> L,R; eta_faces-> U,D; zeta_faces-> F,B 
		To calculate α = uᵢ ∂	(ρ x ξ𐞥ⱼ x τᵢⱼ)
							∂ξ𐞥
				Simplified to:
					 α = uᵢ ∂ (ρ τᵢⱼ ξ𐞥ⱼ)
							∂ξ𐞥	   J
		The metrics (ξ𐞥ⱼ) are calculated at the face centers for Each of the six faces
			So, csi0,csi1,csi2 is ∂ξ/∂x, ∂ξ/∂y, ∂ξ/∂z the same as csi[k][j][i].x, csi[k][j][i].y, csi[k][j][i].z ;
		*/
		rho.L = (q[k][j][i].rho + q[k][j][i-1].rho)*0.5; rho.R = (q[k][j][i+1].rho + q[k][j][i].rho)*0.5;
		rho.D = (q[k][j][i].rho + q[k][j-1][i].rho)*0.5; rho.U = (q[k][j+1][i].rho + q[k][j][i].rho)*0.5;
		rho.B = (q[k][j][i].rho + q[k-1][j][i].rho)*0.5; rho.F = (q[k+1][j][i].rho + q[k][j][i].rho)*0.5;

		u[0] = q[k][j][i].rhoU/q[k][j][i].rho ; u[1] = q[k][j][i].rhoV/q[k][j][i].rho ; u[2] = q[k][j][i].rhoW/q[k][j][i].rho ;

		visc[k][j][i].rhoE += 	u[0]* (
									rho.R*alpha1[k][j][i].x - rho.L*alpha1[k][j][i-1].x+
									rho.U*alpha2[k][j][i].x - rho.D*alpha2[k][j-1][i].x+
									rho.F* alpha3[k][j][i].x - rho.B*alpha3[k-1][j][i].x) +
							 	u[1]*(		
									rho.R*alpha1[k][j][i].y - rho.L*alpha1[k][j][i-1].y+
									rho.U*alpha2[k][j][i].y - rho.D*alpha2[k][j-1][i].y+
									rho.F* alpha3[k][j][i].y - rho.B*alpha3[k-1][j][i].y) +
								u[2]*(
									rho.R*alpha1[k][j][i].z - rho.L*alpha1[k][j][i-1].z+
									rho.U*alpha2[k][j][i].z - rho.D*alpha2[k][j-1][i].z+
									rho.F* alpha3[k][j][i].z - rho.B*alpha3[k-1][j][i].z);								
		  ////////////////// 1/23/2023
	  //visc[k][j][i].rhoE -=bbeta[k][j][i]/2.0*pow(user->Ma,2.) ;
			}// End if LES
		}
	  }
    }
  
  DMDAVecRestoreArray(cda, user->lQ, &q);
  DMDAVecRestoreArray(cda, Visc,  &visc);
  
  DMDAVecRestoreArray(da, user->lP, &p);
  
  
  DMDAVecRestoreArray(fda, Csi, &csi);
  DMDAVecRestoreArray(fda, Eta, &eta);
  DMDAVecRestoreArray(fda, Zet, &zet);
  
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  
  DMDAVecRestoreArray(cda, Fp1, &fp1);
  DMDAVecRestoreArray(cda, Fp2, &fp2);
  DMDAVecRestoreArray(cda, Fp3, &fp3);
  
  DMDAVecRestoreArray(da, user->lAj, &aj);
  
  DMDAVecRestoreArray(fda, user->lUcat,  &ucat);
  
  
  
  if(les) {
    DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
    
    DMDAVecRestoreArray(da, user->lBBeta , &bbeta);

    DMDAVecRestoreArray(fda, Alpha1, &alpha1);
    DMDAVecRestoreArray(fda, Alpha2, &alpha2);
    DMDAVecRestoreArray(fda, Alpha3, &alpha3);
    DMDAVecRestoreArray(fda, user->lSigma, &lsigma);
  
    VecDestroy(&Alpha1);
    VecDestroy(&Alpha2);
    VecDestroy(&Alpha3);
    
  } else if (rans) {    
    DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
  }
  
  PetscReal norm;

  VecNorm(Visc, NORM_INFINITY, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "!!norm of Visc %le Sutherland %d Tref %f\n", norm, Sutherland, Tref);
  VecDestroy(&Fp1);
  VecDestroy(&Fp2);
  VecDestroy(&Fp3);

  return(0);
}

PetscReal minmod(PetscReal a, PetscReal b)
{
	PetscReal f;
	if (a*b<0) f=0.;
	else 
		f=PetscSign(a)*PetscMin(fabs(a), fabs(b));
	return(f);
}


PetscReal MUSCL_QL(PetscReal Qm, PetscReal Q0, PetscReal Qp)
{
	/* Inputs: Q0 is Q at i
	   Qp is Q at i+1
	   Qm is Q at i-1
MUSCL: creates a linear polynomial in (i-1/2, i+1/2)
Q(csi)=Q0+S(cis-i)
S is the minimum o the slope between (i,i-1) & (i+1,i)
QL=Q0+ S (i+1/2-i)==> QL=Q0 + 1/2 S
Output: returns QL at i+1/2    
*/
	PetscReal S, Splus, Sminus, QL;

	Splus =Qp - Q0;
	Sminus=Q0 - Qm;
	S=minmod(Splus, Sminus);
	QL=Q0 + 0.5*S;

	return(QL);
}

PetscReal MUSCL_QR(PetscReal Qm, PetscReal Q0, PetscReal Qp)
{

	/* Inputs: Q0 is Q at i+1
	   Qp is Q at i+2
	   Qm is Q at i	     
MUSCL: creates a linear polynomial in (i+1/2, i+1+1/2)
Q(csi)=Q0+S(cis-(i+1))
S is the minimum o the slope between (i+1,i) & (i+2,i+1)
QR=Q0+ S (i+1/2-(i+1)) ==> QR=Q0 - 1/2 S
Output: returns QR at i+1/2      
*/
	PetscReal S, Splus, Sminus, QR;

	Splus=Qp-Q0;
	Sminus=Q0-Qm;
	S=minmod(Splus, Sminus);
	QR=Q0 - 0.5*S;

	return(QR);
}

PetscErrorCode MUSCL_QRL(PetscReal Qm, PetscReal Q0, PetscReal Qp, PetscReal QL, PetscReal QR)
{

	/* Inputs: Q0 is Q at i
	   Qp is Q at i+1
	   Qm is Q at i-1	     
MUSCL: creates a linear polynomial in (i-1/2, i+1/2)
Q(csi)=Q0+S(cis - i)
S is the minimum o the slope between (i,i-1) & (i+1,i)
QL=Q0+ S (i+1/2-i) ==> QL=Q0 + 1/2 S
QR=Q0+ S (i-1/2-i) ==> QR=Q0 - 1/2 S
Output: returns QR at i-1/2 
QL at i+1/2    
*/
	PetscReal S, Splus, Sminus;

	Splus=Qp-Q0;
	Sminus=Q0-Qm;
	S=minmod(Splus, Sminus);
	QL=Q0 + 0.5*S;
	QR=Q0 - 0.5*S;

	return(0);
}

PetscReal compute_ri(PetscReal fR, PetscReal f0, PetscReal fL)
{
	PetscReal dfL,dfR,R;
	dfL=fL-f0;
	dfR=f0-fR;
	R=(fabs(2*dfL*dfR)+ks)/(dfL*dfL+dfR*dfR+ks);


	return(R);
}


PetscReal sigma_cent(PetscReal fL, PetscReal f0, PetscReal fR, PetscReal fRR)
{  
	PetscReal rd,r,rL,rR,rmax,sig;
	rL=compute_ri(fL,f0,fR);
	rR=compute_ri(f0,fR,fRR);

	r=PetscMin(rR,rL);
	rd=0.5;
	rmax=PetscMax(rd,r);
	sig=0.5+0.5*tanh(3.0*(r-rd)/rmax)/tanh(3.0);
	return(sig);

}
PetscReal WENO_QL(PetscReal Qm, PetscReal Q0, PetscReal Qp)
{
	/* WENO 3rd order based on Procedure 2.2 of Shu 
	   "ESSENTIALLY NON-OSCILLATORY AND WEIGHTED
	   ESSENTIALLY NON-OSCILLATORY SCHEMES
	   FOR HYPERBOLIC CONSERVATION LAWS"

Inputs: Qm is Q at i-1 
        Q0 is Q at i
        Qp is Q at i+1
WENO: creates a quadratic polynomial in (i-1/2, i+1/2) 
using different stencils of ENO Qr0 and Qr1
combines them to get WENO 
k=2 here, i.e, the accuracy is 2k-1=3

1) Qr(i+1/2)=sum_(j=0 to k-1) [c_rj Q(i-r+j)]
c_rj for k=2          j=0      j=1
r=-1   3/2      -1/2
r= 0   1/2       1/2
r= 1  -1/2       3/2
2) d0=2/3, d1=1/3 such that Q(i+1/2)=sigma_(r=0 to k-1) dr*Qr = Q(x_i+1/2)+O(dx^2k-1)
3) smoothness beta
4) wieghts w_r=alfa_r / sum_s(alfa_s) alfa_r = dr/(eps + beta_r)^2
5) QL=sum_(r=0 to k-1) w_r * Qr
Output: returns QL at i+1/2    
*/

	PetscReal d0=0.6666667, d1=0.3333333;
	PetscReal beta0, beta1;
	PetscReal alfa0, alfa1, eps=1e-6, alfasum;
	PetscReal w0, w1;
	PetscReal QL, Qr0, Qr1;


	Qr0= 0.5*Q0 + 0.5*Qp;
	Qr1=-0.5*Qm + 1.5*Q0;

	beta0=(Qp-Q0)*(Qp-Q0);
	beta1=(Qm-Q0)*(Qm-Q0);

	alfa0=d0/((eps+beta0)*(eps+beta0));
	alfa1=d1/((eps+beta1)*(eps+beta1));
	alfasum=alfa0+alfa1;

	w0=alfa0/alfasum;
	w1=alfa1/alfasum;

	QL=w0*Qr0 + w1*Qr1;

	return(QL);
}

PetscReal WENO_QR(PetscReal Qm, PetscReal Q0, PetscReal Qp)
{
	/* WENO 3rd order based on Procedure 2.2 of Shu 
	   "ESSENTIALLY NON-OSCILLATORY AND WEIGHTED
	   ESSENTIALLY NON-OSCILLATORY SCHEMES
	   FOR HYPERBOLIC CONSERVATION LAWS"

Inputs: Qm is Q at i-1
Q0 is Q at i
Qp is Q at i+1cen
WENO: creates a quadratic polynomial in (i-1/2, i+1/2) 
using different stencils of ENO Qr0 and Qr1
combines them to get WNO
k=2 here, i.e, the accuracy is 2k-1=3

1) Qr(i-1/2)=sum_(j=0 to k-1) [c_(r-1)j Q(i-r+j)]
c_rj for k=2          j=0      j=1
r=-1   3/2      -1/2
r= 0   1/2       1/2
r= 1  -1/2       3/2
2) d0=2/3, d1=1/3 such that Q(i+1/2)=sigma_(r=0 to k-1) dr*Qr = Q(x_i+1/2)+O(dx^2k-1)
3) smoothness beta
4) wieghts w_r=alfa_r / sum_s(alfa_s) alfa_r = dr/(eps + beta_r)^2
5) QL=sum_(r=0 to k-1) w_r * Qr
Output: returns QR at i-1/2    
*/

	PetscReal d0=0.6666667, d1=0.3333333;
	PetscReal beta0, beta1;
	PetscReal alfa0, alfa1, eps=1e-6, alfasum;
	PetscReal w0, w1;
	PetscReal QR, Qr0, Qr1;


	Qr0= 1.5*Q0 - 0.5*Qp;
	Qr1= 0.5*Qm + 0.5*Q0;

	beta0=(Qp-Q0)*(Qp-Q0);
	beta1=(Qm-Q0)*(Qm-Q0);

	alfa0=d0/((eps+beta0)*(eps+beta0));
	alfa1=d1/((eps+beta1)*(eps+beta1));
	alfasum=alfa0+alfa1;

	w0=alfa0/alfasum;
	w1=alfa1/alfasum;

	QR=w0*Qr0 + w1*Qr1;

	return(QR);
}
PetscReal Fp_skew(PetscReal u0, PetscReal u1, PetscReal u2, PetscReal u3, PetscReal v0, PetscReal v1, PetscReal v2, PetscReal v3)
{
  PetscReal uv;
  uv=(u1+u2)*(v1+v2)/3.0 - 1./24.*(u0*v0+u0*v2+u1*v1+u1*v3+u2*v2+u2*v0+u3*v1+u3*v3);
  return(uv);
}
PetscReal Fp_2(PetscReal u[4], PetscReal v[4])
{
  // doi:10.1006/jcph.2000.6492  non-dissipative skew central for convective fluxes
  PetscReal uv;
  uv=(u[1]+u[2])*(v[1]+v[2])/3.0 - 1./24.*(u[0]*v[0]+u[0]*v[2]+u[1]*v[1]+u[1]*v[3]+u[2]*v[2]+u[2]*v[0]+u[3]*v[1]+u[3]*v[3]);
  return(uv);
}

PetscReal Fp_1(PetscReal u[4], PetscReal v[4])
{
  PetscReal uv;
  uv=1.0/12.0*(-u[0]*v[0]+7*u[1]*v[1]+7*u[2]*v[2]-u[3]*v[3]);
  return(uv);
}

PetscReal Fp4(PetscReal u[4])
{
  // Central
  PetscReal uv;
  uv=1.0/12.0*(-u[0]+7*u[1]+7*u[2]-u[3]);
  return(uv);
}

PetscReal Skew_Central_Flux(PetscReal u[4], PetscReal v[4]){
	PetscReal uv;
	// u[0] = u_i-1 , u[1] = u_i , u[2] = u_i+1, u[3] = u_i+2 
	uv = 1.0/12.0 *( -u[3]*v[3] + 7*u[2]*v[2] + 7*u[1]*v[1] -u[0]*v[0] ) + 1./3.*( 0.5* (u[2]*v[2]+ u[1]*v[1]) - (0.25* (u[2]+v[1])*(u[2]+v[1])) );
	return uv;
}

PetscReal ENO_QR(PetscReal Qm, PetscReal Q0, PetscReal Qp, PetscInt k)
{
  /* ENO 3rd order based on Procedure 2.1 (page 17) of Shu 
     "ESSENTIALLY NON-OSCILLATORY AND WEIGHTED
      ESSENTIALLY NON-OSCILLATORY SCHEMES
      FOR HYPERBOLIC CONSERVATION LAWS"

     Inputs: Qm is Q at i-1
             Q0 is Q at i
             Qp is Q at i+1
	     i+2 is wall
	     ENO: creates a polynomial in (i-1/2, i+1/2)
	     k=3 here, i.e, the accuracy is k=3
		   
		   1) Qr(i-1/2)=sum_(j=0 to k-1) [c_(r-1)j Q(i-r+j)]
	    i+1->i or Qr(i+1/2)=sum_(j=0 to k-1) [c_(r-1)j Q(i+1-r+j)] 
		   k=3, need r=2 to get i-1, i, i+1 stencil
		   need c_1j
		   c_rj for k=3          j=0      j=1    j=2
				  r= 1  -1/6      5/6    1/3
		   k=2, need r=1 to get i, i+1 stencil
		   c_rj for k=2          j=0      j=1
		                  r=-1   3/2      -1/2
				  r= 0   1/2       1/2
				  r= 1  -1/2       3/2
     Output: returns QR at i-1/2  (i+1/2) when i+2 is wall

  */

  PetscReal QR;
  if (k==2) 
    QR=0.5*Q0+0.5*Qp; //2nd
  else 
    QR=-Qm*0.1666666666666666 + Q0*0.8333333333333333 + Qp*0.3333333333333333; //3rd
  return(QR);
}

PetscReal ENO_QL(PetscReal Qm, PetscReal Q0, PetscReal Qp, PetscInt k)
{
  /* ENO 3rd order based on Procedure 2.1 (page 17) of Shu 
     "ESSENTIALLY NON-OSCILLATORY AND WEIGHTED
      ESSENTIALLY NON-OSCILLATORY SCHEMES
      FOR HYPERBOLIC CONSERVATION LAWS"

     Inputs: Qm is Q at i
             Q0 is Q at i+1
             Qp is Q at i+2
	     i-1 is wall
	     ENO: k=3 here, i.e, the accuracy is k=3
		   1) Qr(i+1/2)=sum_(j=0 to k-1) [c_rj Q(i-r+j)]
		   k=3, need r=0 to get i, i+1, i+2 stencil
		   need c_1j
		   c_rj for k=3          j=0      j=1    j=2
				  r= 0   1/3      5/6    -1/6
		   k=2, need r=0 to get i, i+1 stencil
		   c_rj for k=2          j=0      j=1
		                  r=-1   3/2      -1/2
				  r= 0   1/2       1/2
				  r= 1  -1/2       3/2
     Output: returns QL at i+1/2 when i-1 is wall

  */

  PetscReal QL;
  if (k==2) 
    QL=0.5*Q0+0.5*Qm; //2nd
  else
    QL=Qm*0.3333333333333333 + Q0*0.8333333333333333 - Qp*0.1666666666666666; //3rd
  
  return(QL);
}

PetscErrorCode formConvection(UserCtx *user, Vec Conv)
{
  // centeral difference partial tranforamtion

  DM		da = user->da, cda = user->cda;
  DMDALocalInfo	info=user->info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k;
  PetscInt      fluxorder=3, order=2; // ENO order==2 or 3rd, fluxorder ==1 Godunov, 2 MUSCL, 3 WENOe
  PetscOptionsGetInt(NULL, NULL, "-ENOorder", &order, NULL);
  PetscOptionsGetInt(NULL, NULL, "-fluxorder", &fluxorder, NULL);

  PetscInt twoD=0;
  PetscOptionsGetInt(NULL, NULL, "-twoD", &twoD, NULL);

  // PetscPrintf(PETSC_COMM_WORLD, "!!Flux order %i ENO order %i TwoD %d\n", fluxorder, order, twoD);

 
  //Vec		Fp1, Fp2, Fp3;
  Vec		lFp1, lFp2, lFp3;
  CompVars	***fp1, ***fp2, ***fp3;
  CompVars	***conv, ***q, QR,QL, qm, qp;
 
  PetscReal	uconR, uconL, pR, pL, ***p, pm, pp;
  PetscReal     gii, CsoundL, CsoundR, saR, saL, SA, gamma=1.4;

  PetscReal	***nvert,***aj, blank=0.1, blankwall=2.9;
  Cmpnts	***csi, ***eta, ***zet, csiH;// csiR, csiL;
  PetscErrorCode ierr;

  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;

  /* VecDuplicate(user->Q, &Fp1); */
  /* VecDuplicate(user->Q, &Fp2); */
  /* VecDuplicate(user->Q, &Fp3); */
   
  /* DMDAVecGetArray(cda, Fp1, &fp1); */
  /* DMDAVecGetArray(cda, Fp2, &fp2); */
  /* DMDAVecGetArray(cda, Fp3, &fp3); */

  VecDuplicate(user->lQ, &lFp1);
  VecDuplicate(user->lQ, &lFp2);
  ierr =  VecDuplicate(user->lQ, &lFp3); CHKERRQ(ierr);
  VecSet(lFp1, 0.);

  DMDAVecGetArray(cda, lFp1, &fp1);
  DMDAVecGetArray(cda, lFp2, &fp2);
  ierr = DMDAVecGetArray(cda, lFp3, &fp3); CHKERRQ(ierr);

  DMDAVecGetArray(cda, user->lQ, &q);

  DMDAVecGetArray(user->fda, user->lCsi, &csi);
  DMDAVecGetArray(user->fda, user->lEta, &eta);
  DMDAVecGetArray(user->fda, user->lZet, &zet);

  DMDAVecGetArray(da, user->lP, &p);
  DMDAVecGetArray(da, user->lNvert, &nvert);
    
  // i+1/2->i
  if (twoD!=1) {
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs-1; i<lxe; i++){
	// calc F1 on i+1/2 nodes... i+1/2 saved in i 
	//-Lax-Freidrich flux: F(QR,QL)=1/2(E(QR)+E(QL))-1/2|S(A)|(QR-QL)
	/* if (i==0) { // on the boundary */
	/*   QL=q[k][j][i]; */
	/*   QR=q[k][j][i]; */

	/*   pL=p[k][j][i]; */
	/*   pR=p[k][j][i]; */
	/* } else if (i==mx-2) { */
	/*   QL=q[k][j][i+1]; */
	/*   QR=q[k][j][i+1]; */

	/*   pL=p[k][j][i+1]; */
	/*   pR=p[k][j][i+1]; */
	/* } else  */
	if (fluxorder==1 || i==0 || i==mx-2 ||
	  ((nvert[k][j][i]+nvert[k][j][i+1])>blankwall) ) { // if i or i+1 are wall nodes
	  QL=q[k][j][i];
	  QR=q[k][j][i+1];

	  pL=p[k][j][i];
	  pR=p[k][j][i+1];
	} else if (fluxorder==2) {
	  /* we know i and i+1 are not wall nodes but i-1 or i+1 might be wall nodes  */
	  if (nvert[k][j][i-1]>blankwall) {
	     /* put Qm equal to Qp to revert WENO to ENO */
	    qm=q[k][j][i+1];
	    pm=p[k][j][i+1];
	  } else {
	    qm=q[k][j][i-1];
	    pm=p[k][j][i-1];
	  }
	  if (nvert[k][j][i+2]>blankwall) {
	     /* put Qp equal to Qm to revert WENO to ENO */
	    qp=q[k][j][i];
	    pp=p[k][j][i];
	  } else {
	    qp=q[k][j][i+2];
	    pp=p[k][j][i+2];
	  }
	  
	  QL.rho =MUSCL_QL(qm.rho , q[k][j][i].rho , q[k][j][i+1].rho);
	  QL.rhoU=MUSCL_QL(qm.rhoU, q[k][j][i].rhoU, q[k][j][i+1].rhoU);
	  QL.rhoV=MUSCL_QL(qm.rhoV, q[k][j][i].rhoV, q[k][j][i+1].rhoV);
	  QL.rhoW=MUSCL_QL(qm.rhoW, q[k][j][i].rhoW, q[k][j][i+1].rhoW);
	  QL.rhoE=MUSCL_QL(qm.rhoE, q[k][j][i].rhoE, q[k][j][i+1].rhoE);
	  pL     =MUSCL_QL(pm     , p[k][j][i],      p[k][j][i+1]);

	  QR.rho =MUSCL_QR(q[k][j][i].rho , q[k][j][i+1].rho , qp.rho);
	  QR.rhoU=MUSCL_QR(q[k][j][i].rhoU, q[k][j][i+1].rhoU, qp.rhoU);
	  QR.rhoV=MUSCL_QR(q[k][j][i].rhoV, q[k][j][i+1].rhoV, qp.rhoV);
	  QR.rhoW=MUSCL_QR(q[k][j][i].rhoW, q[k][j][i+1].rhoW, qp.rhoW);
	  QR.rhoE=MUSCL_QR(q[k][j][i].rhoE, q[k][j][i+1].rhoE, qp.rhoE);
	  pR     =MUSCL_QR(p[k][j][i]     , p[k][j][i+1],      pp);
	} else { // WENO
	   /* we know i and i+1 are not wall nodes but i-1 or i+1 might be wall nodes  */
	  if ((nvert[k][j][i-1]>blankwall)  && (i+2)<mx) {
	     /* put Qm equal to Qp to revert WENO to ENO */	
	  /*   qm=q[k][j][i+1]; */
	  /*   pm=p[k][j][i+1]; */
	  /* QL.rho =MUSCL_QL(qm.rho , q[k][j][i].rho , q[k][j][i+1].rho); */
	  /* QL.rhoU=MUSCL_QL(qm.rhoU, q[k][j][i].rhoU, q[k][j][i+1].rhoU); */
	  /* QL.rhoV=MUSCL_QL(qm.rhoV, q[k][j][i].rhoV, q[k][j][i+1].rhoV); */
	  /* QL.rhoW=MUSCL_QL(qm.rhoW, q[k][j][i].rhoW, q[k][j][i+1].rhoW); */
	  /* QL.rhoE=MUSCL_QL(qm.rhoE, q[k][j][i].rhoE, q[k][j][i+1].rhoE); */
	  /* pL     =MUSCL_QL(pm     , p[k][j][i],      p[k][j][i+1]); */

	    QL.rho =ENO_QL(q[k][j][i].rho , q[k][j][i+1].rho , q[k][j][i+2].rho ,order);
	    QL.rhoU=ENO_QL(q[k][j][i].rhoU, q[k][j][i+1].rhoU, q[k][j][i+2].rhoU,order);
	    QL.rhoV=ENO_QL(q[k][j][i].rhoV, q[k][j][i+1].rhoV, q[k][j][i+2].rhoV,order);
	    QL.rhoW=ENO_QL(q[k][j][i].rhoW, q[k][j][i+1].rhoW, q[k][j][i+2].rhoW,order);
	    QL.rhoE=ENO_QL(q[k][j][i].rhoE, q[k][j][i+1].rhoE, q[k][j][i+2].rhoE,order);
	    pL     =ENO_QL(p[k][j][i]     , p[k][j][i+1],      p[k][j][i+2]     ,order);	 

	  } else {
	    qm=q[k][j][i-1];
	    pm=p[k][j][i-1];       	 
	  
	    QL.rho =WENO_QL(qm.rho , q[k][j][i].rho , q[k][j][i+1].rho);
	    QL.rhoU=WENO_QL(qm.rhoU, q[k][j][i].rhoU, q[k][j][i+1].rhoU);
	    QL.rhoV=WENO_QL(qm.rhoV, q[k][j][i].rhoV, q[k][j][i+1].rhoV);
	    QL.rhoW=WENO_QL(qm.rhoW, q[k][j][i].rhoW, q[k][j][i+1].rhoW);
	    QL.rhoE=WENO_QL(qm.rhoE, q[k][j][i].rhoE, q[k][j][i+1].rhoE);
	    pL     =WENO_QL(pm     , p[k][j][i],      p[k][j][i+1]);
	  }

	  if (nvert[k][j][i+2]>blankwall) {// || (i+2)==mx-1) {
	     /* put Qp equal to Qm to revert WENO to ENO */
	    qp=q[k][j][i+1];
	    pp=p[k][j][i+1];
	    QR.rho =ENO_QR(q[k][j][i-1].rho , q[k][j][i].rho , qp.rho ,order);
	    QR.rhoU=ENO_QR(q[k][j][i-1].rhoU, q[k][j][i].rhoU, qp.rhoU,order);
	    QR.rhoV=ENO_QR(q[k][j][i-1].rhoV, q[k][j][i].rhoV, qp.rhoV,order);
	    QR.rhoW=ENO_QR(q[k][j][i-1].rhoW, q[k][j][i].rhoW, qp.rhoW,order);
	    QR.rhoE=ENO_QR(q[k][j][i-1].rhoE, q[k][j][i].rhoE, qp.rhoE,order);
	    pR     =ENO_QR(p[k][j][i-1]     , p[k][j][i],      pp     ,order);
	  /*   qp=q[k][j][i]; */
	  /*   pp=p[k][j][i]; */
	  /* QR.rho =MUSCL_QR(q[k][j][i].rho , q[k][j][i+1].rho , qp.rho); */
	  /* QR.rhoU=MUSCL_QR(q[k][j][i].rhoU, q[k][j][i+1].rhoU, qp.rhoU); */
	  /* QR.rhoV=MUSCL_QR(q[k][j][i].rhoV, q[k][j][i+1].rhoV, qp.rhoV); */
	  /* QR.rhoW=MUSCL_QR(q[k][j][i].rhoW, q[k][j][i+1].rhoW, qp.rhoW); */
	  /* QR.rhoE=MUSCL_QR(q[k][j][i].rhoE, q[k][j][i+1].rhoE, qp.rhoE); */
	  /* pR     =MUSCL_QR(p[k][j][i]     , p[k][j][i+1],      pp); */
	  } else {
	    qp=q[k][j][i+2];
	    pp=p[k][j][i+2];
	  

	  QR.rho =WENO_QR(q[k][j][i].rho , q[k][j][i+1].rho , qp.rho);
	  QR.rhoU=WENO_QR(q[k][j][i].rhoU, q[k][j][i+1].rhoU, qp.rhoU);
	  QR.rhoV=WENO_QR(q[k][j][i].rhoV, q[k][j][i+1].rhoV, qp.rhoV);
	  QR.rhoW=WENO_QR(q[k][j][i].rhoW, q[k][j][i+1].rhoW, qp.rhoW);
	  QR.rhoE=WENO_QR(q[k][j][i].rhoE, q[k][j][i+1].rhoE, qp.rhoE);
	  pR     =WENO_QR(p[k][j][i]     , p[k][j][i+1],      pp);
	  }
	}

	/* // calc pL & pR from Eq of State */
	/* pL=EqOfStateP(QL); */
	/* pR=EqOfStateP(QR); */

	//csi at half point i+1/2
	//	csiL=csi[k][j][i];
	//	csiR=csi[k][j][i+1];
	csiH.x = 0.5*(csi[k][j][i].x+csi[k][j][i+1].x);
	csiH.y = 0.5*(csi[k][j][i].y+csi[k][j][i+1].y);
	csiH.z = 0.5*(csi[k][j][i].z+csi[k][j][i+1].z);

	if (i==mx-2 && user->bctype[1]==1)
	  uconL=(csiH.x*(QR.rhoU/QR.rho)+
		 csiH.y*(QR.rhoV/QR.rho)+
		 csiH.z*(QR.rhoW/QR.rho));
	else
	  uconL=(csiH.x*(QL.rhoU/QL.rho)+
		 csiH.y*(QL.rhoV/QL.rho)+
		 csiH.z*(QL.rhoW/QL.rho));
	if (i==0 && user->bctype[0]==1) 
	  uconR=uconL;
	else 
	  uconR=(csiH.x*(QR.rhoU/QR.rho)+
		 csiH.y*(QR.rhoV/QR.rho)+
		 csiH.z*(QR.rhoW/QR.rho));

	gii=csiH.x*csiH.x+csiH.y*csiH.y+csiH.z*csiH.z;
	//	giiR=csiR.x*csiR.x+csiR.y*csiR.y+csiR.z*csiR.z;
	CsoundL=gamma*pL/QL.rho;
	CsoundR=gamma*pR/QR.rho;
	saL=fabs(uconL)+sqrt(fabs(CsoundL)*gii);
	saR=fabs(uconR)+sqrt(fabs(CsoundR)*gii);
	SA=PetscMax(saL,saR);

	fp1[k][j][i].rho= 0.5*( QL.rho*uconL  + QR.rho *uconR) - 0.5*SA*(QR.rho-QL.rho);
	fp1[k][j][i].rhoU=0.5*((QL.rhoU*uconL + QR.rhoU*uconR) + 
			       (pL   * csiH.x + pR * csiH.x))  - 0.5*SA*(QR.rhoU-QL.rhoU);
	fp1[k][j][i].rhoV=0.5*((QL.rhoV*uconL + QR.rhoV*uconR) + 
			       (pL   * csiH.y + pR * csiH.y))  - 0.5*SA*(QR.rhoV-QL.rhoV);
	fp1[k][j][i].rhoW=0.5*((QL.rhoW*uconL + QR.rhoW*uconR) + 
			       (pL   * csiH.z + pR * csiH.z))  - 0.5*SA*(QR.rhoW-QL.rhoW);
	fp1[k][j][i].rhoE=0.5*((pL+QL.rhoE)*uconL + 
			       (pR+QR.rhoE)*uconR)             - 0.5*SA*(QR.rhoE-QL.rhoE);
      }
    }
  }
  } // if !twoD
    
  /* j direction */
  for (k=lzs; k<lze; k++){
    for (j=lys-1; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	// calc F2 on j+1/2 nodes... j+1/2 saved in j
	//-Lax-Freidrich flux: F(QR,QL)=1/2(E(QR)+E(QL))-1/2|S(A)|(QR-QL)
	if (j==0 && isothermal && user->bctype[2]==1) { // on the boundary
	  QL=q[k][j][i];
	  QR=q[k][j][i];

	  pL=p[k][j][i];
	  pR=p[k][j][i];
	} else if (j==my-2 && isothermal && user->bctype[3]==1) {
	  QL=q[k][j+1][i];
	  QR=q[k][j+1][i];

	  pL=p[k][j+1][i];
	  pR=p[k][j+1][i];
	} else if (fluxorder==1 ||  j==0 || j==my-2 || //j==1 || j==my-3 ||//
	    ((nvert[k][j][i]+nvert[k][j+1][i])>blankwall) ) { // if j or j+1 are wall nodes
	  QL=q[k][j][i];
	  QR=q[k][j+1][i];
	  pL=p[k][j][i];
	  pR=p[k][j+1][i];
	} else if (fluxorder==2 || (j<30 && k<7)) { // iman 3/6/2024 for leading edge of plate
	  /* we know j and j+1 are not wall nodes but j-1 or j+1 might be wall nodes  */
	  if (nvert[k][j-1][i]>blankwall) {
	    /* put Qm equal to Qp to revert WENO to ENO */
	    qm=q[k][j+1][i];
	    pm=p[k][j+1][i];
	  } else {
	    qm=q[k][j-1][i];
	    pm=p[k][j-1][i];
	  }
	  if (nvert[k][j+2][i]>blankwall) {
	     /* put Qp equal to Qm to revert WENO to ENO */
	    qp=q[k][j][i];
	    pp=p[k][j][i];
	  } else {
	    qp=q[k][j+2][i];
	    pp=p[k][j+2][i];
	  }
	  QL.rho =MUSCL_QL(qm.rho , q[k][j][i].rho , q[k][j+1][i].rho);
	  QL.rhoU=MUSCL_QL(qm.rhoU, q[k][j][i].rhoU, q[k][j+1][i].rhoU);
	  QL.rhoV=MUSCL_QL(qm.rhoV, q[k][j][i].rhoV, q[k][j+1][i].rhoV);
	  QL.rhoW=MUSCL_QL(qm.rhoW, q[k][j][i].rhoW, q[k][j+1][i].rhoW);
	  QL.rhoE=MUSCL_QL(qm.rhoE, q[k][j][i].rhoE, q[k][j+1][i].rhoE);
	  pL     =MUSCL_QL(pm     , p[k][j][i],      p[k][j+1][i]);

	  QR.rho =MUSCL_QR(q[k][j][i].rho , q[k][j+1][i].rho , qp.rho);
	  QR.rhoU=MUSCL_QR(q[k][j][i].rhoU, q[k][j+1][i].rhoU, qp.rhoU);
	  QR.rhoV=MUSCL_QR(q[k][j][i].rhoV, q[k][j+1][i].rhoV, qp.rhoV);
	  QR.rhoW=MUSCL_QR(q[k][j][i].rhoW, q[k][j+1][i].rhoW, qp.rhoW);
	  QR.rhoE=MUSCL_QR(q[k][j][i].rhoE, q[k][j+1][i].rhoE, qp.rhoE);
	  pR     =MUSCL_QR(p[k][j][i]     , p[k][j+1][i],      pp);
	} else { //WENO
	   /* we know j and j+1 are not wall nodes but j-1 or j+2 might be wall nodes  */
	  if (nvert[k][j-1][i]>blankwall) {
	    /* put Qm equal to Qp to revert WENO to ENO */
	    /* qm=q[k][j+1][i]; */
	    /* pm=p[k][j+1][i]; */
	  QL.rho =ENO_QL(q[k][j][i].rho , q[k][j+1][i].rho , q[k][j+2][i].rho ,order);
	  QL.rhoU=ENO_QL(q[k][j][i].rhoU, q[k][j+1][i].rhoU, q[k][j+2][i].rhoU,order);
	  QL.rhoV=ENO_QL(q[k][j][i].rhoV, q[k][j+1][i].rhoV, q[k][j+2][i].rhoV,order);
	  QL.rhoW=ENO_QL(q[k][j][i].rhoW, q[k][j+1][i].rhoW, q[k][j+2][i].rhoW,order);
	  QL.rhoE=ENO_QL(q[k][j][i].rhoE, q[k][j+1][i].rhoE, q[k][j+2][i].rhoE,order);
	  pL     =ENO_QL(p[k][j][i]     , p[k][j+1][i],      p[k][j+2][i]     ,order);
	  } else {
	    qm=q[k][j-1][i];
	    pm=p[k][j-1][i];

	  QL.rhoU=WENO_QL(qm.rhoU, q[k][j][i].rhoU, q[k][j+1][i].rhoU);
	  QL.rhoV=WENO_QL(qm.rhoV, q[k][j][i].rhoV, q[k][j+1][i].rhoV);
	  QL.rhoW=WENO_QL(qm.rhoW, q[k][j][i].rhoW, q[k][j+1][i].rhoW);
	  if (j==1 && isothermal && user->bctype[2]==1) { // do not use values on the boundary
	    QL.rho =ENO_QL(q[k][j][i].rho , q[k][j+1][i].rho , q[k][j+2][i].rho ,order);
	    QL.rhoE=ENO_QL(q[k][j][i].rhoE, q[k][j+1][i].rhoE, q[k][j+2][i].rhoE,order);
	    pL     =ENO_QL(p[k][j][i]     , p[k][j+1][i],      p[k][j+2][i]     ,order);
	  } else {
	  QL.rho =WENO_QL(qm.rho , q[k][j][i].rho , q[k][j+1][i].rho);
	  QL.rhoE=WENO_QL(qm.rhoE, q[k][j][i].rhoE, q[k][j+1][i].rhoE);
	  pL     =WENO_QL(pm     , p[k][j][i],      p[k][j+1][i]);
	  }
	  }
	  if (nvert[k][j+2][i]>blankwall) {
	    /* put Qp equal to Qm to revert WENO to ENO */
	    /* qp=q[k][j][i]; */
	    /* pp=p[k][j][i]; */
	    qp=q[k][j+1][i]; 
	    pp=p[k][j+1][i];
	    QR.rho =ENO_QR(q[k][j-1][i].rho , q[k][j][i].rho , qp.rho ,order);
	    QR.rhoU=ENO_QR(q[k][j-1][i].rhoU, q[k][j][i].rhoU, qp.rhoU,order);
	    QR.rhoV=ENO_QR(q[k][j-1][i].rhoV, q[k][j][i].rhoV, qp.rhoV,order);
	    QR.rhoW=ENO_QR(q[k][j-1][i].rhoW, q[k][j][i].rhoW, qp.rhoW,order);
	    QR.rhoE=ENO_QR(q[k][j-1][i].rhoE, q[k][j][i].rhoE, qp.rhoE,order);
	    pR     =ENO_QR(p[k][j-1][i]     , p[k][j][i],      pp     ,order);
	  } else {
	    qp=q[k][j+2][i];
	    pp=p[k][j+2][i];

	  QR.rho =WENO_QR(q[k][j][i].rho , q[k][j+1][i].rho , qp.rho);
	  QR.rhoU=WENO_QR(q[k][j][i].rhoU, q[k][j+1][i].rhoU, qp.rhoU);
	  QR.rhoV=WENO_QR(q[k][j][i].rhoV, q[k][j+1][i].rhoV, qp.rhoV);
	  QR.rhoW=WENO_QR(q[k][j][i].rhoW, q[k][j+1][i].rhoW, qp.rhoW);
	  QR.rhoE=WENO_QR(q[k][j][i].rhoE, q[k][j+1][i].rhoE, qp.rhoE);
	  pR     =WENO_QR(p[k][j][i]     , p[k][j+1][i],      pp);
	  }
	   
	}

	/* // calc pL & pR from Eq of State */
	/* pL=EqOfStateP(QL); */
	/* pR=EqOfStateP(QR); */

	//csi at half point j+1/2
	//	csiL=eta[k][j][i];
	//	csiR=eta[k][j+1][i];
	csiH.x=0.5*(eta[k][j][i].x+eta[k][j+1][i].x);
	csiH.y=0.5*(eta[k][j][i].y+eta[k][j+1][i].y);
	csiH.z=0.5*(eta[k][j][i].z+eta[k][j+1][i].z);

	if (j==my-2 && user->bctype[3]==1)
	  uconL=(csiH.x*(QR.rhoU/QR.rho)+
		 csiH.y*(QR.rhoV/QR.rho)+
		 csiH.z*(QR.rhoW/QR.rho));
	else
	  uconL=(csiH.x*(QL.rhoU/QL.rho)+
		 csiH.y*(QL.rhoV/QL.rho)+
		 csiH.z*(QL.rhoW/QL.rho));
	if (j==0 && user->bctype[2]==1)
	  uconR=uconL;
	else
	  uconR=(csiH.x*(QR.rhoU/QR.rho)+
		 csiH.y*(QR.rhoV/QR.rho)+
		 csiH.z*(QR.rhoW/QR.rho));

	gii=csiH.x*csiH.x+csiH.y*csiH.y+csiH.z*csiH.z;
	//	giiR=csiR.x*csiR.x+csiR.y*csiR.y+csiR.z*csiR.z;
	CsoundL=gamma*pL/QL.rho;
	CsoundR=gamma*pR/QR.rho;
	saL=fabs(uconL)+sqrt(fabs(CsoundL)*gii);
	saR=fabs(uconR)+sqrt(fabs(CsoundR)*gii);
	SA=PetscMax(saL,saR);
	
	fp2[k][j][i].rho= 0.5*( QL.rho*uconL  + QR.rho *uconR) - 0.5*SA*(QR.rho-QL.rho);
	fp2[k][j][i].rhoU=0.5*((QL.rhoU*uconL + QR.rhoU*uconR) + 
			       (pL   * csiH.x + pR * csiH.x))  - 0.5*SA*(QR.rhoU-QL.rhoU);
	fp2[k][j][i].rhoV=0.5*((QL.rhoV*uconL + QR.rhoV*uconR) + 
			       (pL   * csiH.y + pR * csiH.y))  - 0.5*SA*(QR.rhoV-QL.rhoV);
	fp2[k][j][i].rhoW=0.5*((QL.rhoW*uconL + QR.rhoW*uconR) + 
			       (pL   * csiH.z + pR * csiH.z))  - 0.5*SA*(QR.rhoW-QL.rhoW);
	fp2[k][j][i].rhoE=0.5*((pL+QL.rhoE)*uconL + 
			       (pR+QR.rhoE)*uconR)             - 0.5*SA*(QR.rhoE-QL.rhoE);
      }
    }
  }

  for (k=lzs-1; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	// calc F3 on k+1/2 nodes... k+1/2 saved in k
	//-Lax-Freidrich flux: F(QR,QL)=1/2(E(QR)+E(QL))-1/2|S(A)|(QR-QL)
	/* if (k==0) { // on the boundary */
	/*   QL=q[k][j][i]; */
	/*   QR=q[k][j][i]; */

	/*   pL=p[k][j][i]; */
	/*   pR=p[k][j][i]; */
	/* } else if (k==mz-2) { */
	/*   QL=q[k+1][j][i]; */
	/*   QR=q[k+1][j][i]; */

	/*   pL=p[k+1][j][i]; */
	/*   pR=p[k+1][j][i]; */
	/* } else 	 */
	if (fluxorder==1 ||k==0 || k==mz-2 || // k==1 || k==mz-3 ||// 
	    ((nvert[k][j][i]+nvert[k+1][j][i])>blankwall) ) { // if k or k+1 are wall nodes
	  QL=q[k][j][i];
	  QR=q[k+1][j][i];
	  pL=p[k][j][i];
	  pR=p[k+1][j][i];
	} else if (fluxorder==2 || (j<30 && k<7)) { // iman 3/6/2024 for leading edge of plate
	  /* we know k and k+1 are not wall nodes but k-1 or k+2 might be wall nodes  */
	  if (nvert[k-1][j][i]>blankwall) {
	    /* put Qm equal to Qp to revert WENO to ENO */
	    qm=q[k+1][j][i];
	    pm=p[k+1][j][i];
	  } else {
	    qm=q[k-1][j][i];
	    pm=p[k-1][j][i];
	  }
	  if (nvert[k+2][j][i]>blankwall) {
	    /* put Qp equal to Qm to revert WENO to ENO */
	    qp=q[k][j][i];
	    pp=p[k][j][i];
	  } else {
	    qp=q[k+2][j][i];
	    pp=p[k+2][j][i];
	  }
	  QL.rho =MUSCL_QL(qm.rho , q[k][j][i].rho , q[k+1][j][i].rho);
	  QL.rhoU=MUSCL_QL(qm.rhoU, q[k][j][i].rhoU, q[k+1][j][i].rhoU);
	  QL.rhoV=MUSCL_QL(qm.rhoV, q[k][j][i].rhoV, q[k+1][j][i].rhoV);
	  QL.rhoW=MUSCL_QL(qm.rhoW, q[k][j][i].rhoW, q[k+1][j][i].rhoW);
	  QL.rhoE=MUSCL_QL(qm.rhoE, q[k][j][i].rhoE, q[k+1][j][i].rhoE);
	  pL     =MUSCL_QL(pm     , p[k][j][i],      p[k+1][j][i]);

	  QR.rho =MUSCL_QR(q[k][j][i].rho , q[k+1][j][i].rho , qp.rho);
	  QR.rhoU=MUSCL_QR(q[k][j][i].rhoU, q[k+1][j][i].rhoU, qp.rhoU);
	  QR.rhoV=MUSCL_QR(q[k][j][i].rhoV, q[k+1][j][i].rhoV, qp.rhoV);
	  QR.rhoW=MUSCL_QR(q[k][j][i].rhoW, q[k+1][j][i].rhoW, qp.rhoW);
	  QR.rhoE=MUSCL_QR(q[k][j][i].rhoE, q[k+1][j][i].rhoE, qp.rhoE);
	  pR     =MUSCL_QR(p[k][j][i]     , p[k+1][j][i],      pp);
	} else { // WENO
	  /* we know k and k+1 are not wall nodes but k-1 or k+2 might be wall nodes  */
	  if (nvert[k-1][j][i]>blankwall) {
	    /* put Qm equal to Qp to revert WENO to ENO */
	    /* qm=q[k+1][j][i]; */
	    /* pm=p[k+1][j][i]; */
	  QL.rho =ENO_QL(q[k][j][i].rho , q[k+1][j][i].rho , q[k+2][j][i].rho ,order);
	  QL.rhoU=ENO_QL(q[k][j][i].rhoU, q[k+1][j][i].rhoU, q[k+2][j][i].rhoU,order);
	  QL.rhoV=ENO_QL(q[k][j][i].rhoV, q[k+1][j][i].rhoV, q[k+2][j][i].rhoV,order);
	  QL.rhoW=ENO_QL(q[k][j][i].rhoW, q[k+1][j][i].rhoW, q[k+2][j][i].rhoW,order);
	  QL.rhoE=ENO_QL(q[k][j][i].rhoE, q[k+1][j][i].rhoE, q[k+2][j][i].rhoE,order);
	  pL     =ENO_QL(p[k][j][i]     , p[k+1][j][i],      p[k+2][j][i]     ,order);
	  } else {
	    qm=q[k-1][j][i];
	    pm=p[k-1][j][i];

	  QL.rho =WENO_QL(qm.rho , q[k][j][i].rho , q[k+1][j][i].rho);
	  QL.rhoU=WENO_QL(qm.rhoU, q[k][j][i].rhoU, q[k+1][j][i].rhoU);
	  QL.rhoV=WENO_QL(qm.rhoV, q[k][j][i].rhoV, q[k+1][j][i].rhoV);
	  QL.rhoW=WENO_QL(qm.rhoW, q[k][j][i].rhoW, q[k+1][j][i].rhoW);
	  QL.rhoE=WENO_QL(qm.rhoE, q[k][j][i].rhoE, q[k+1][j][i].rhoE);
	  pL     =WENO_QL(pm     , p[k][j][i],      p[k+1][j][i]);

	  }
	  if (nvert[k+2][j][i]>blankwall) {
	    /* put Qp equal to Qm to revert WENO to ENO */
	    /* qp=q[k][j][i]; */
	    /* pp=p[k][j][i]; */
	    qp=q[k+1][j][i];
	    pp=p[k+1][j][i];
	    QR.rho =ENO_QR(q[k-1][j][i].rho , q[k][j][i].rho , qp.rho ,order);
	    QR.rhoU=ENO_QR(q[k-1][j][i].rhoU, q[k][j][i].rhoU, qp.rhoU,order);
	    QR.rhoV=ENO_QR(q[k-1][j][i].rhoV, q[k][j][i].rhoV, qp.rhoV,order);
	    QR.rhoW=ENO_QR(q[k-1][j][i].rhoW, q[k][j][i].rhoW, qp.rhoW,order);
	    QR.rhoE=ENO_QR(q[k-1][j][i].rhoE, q[k][j][i].rhoE, qp.rhoE,order);
	    pR     =ENO_QR(p[k-1][j][i]     , p[k][j][i],      pp     ,order);

	  } else {
	    qp=q[k+2][j][i];
	    pp=p[k+2][j][i];
	  QR.rho =WENO_QR(q[k][j][i].rho , q[k+1][j][i].rho , qp.rho);
	  QR.rhoU=WENO_QR(q[k][j][i].rhoU, q[k+1][j][i].rhoU, qp.rhoU);
	  QR.rhoV=WENO_QR(q[k][j][i].rhoV, q[k+1][j][i].rhoV, qp.rhoV);
	  QR.rhoW=WENO_QR(q[k][j][i].rhoW, q[k+1][j][i].rhoW, qp.rhoW);
	  QR.rhoE=WENO_QR(q[k][j][i].rhoE, q[k+1][j][i].rhoE, qp.rhoE);
	  pR     =WENO_QR(p[k][j][i]     , p[k+1][j][i],      pp);

	  }

	}

    	/* // calc pL & pR from Eq of State */
	/* pL=EqOfStateP(QL); */
	/* pR=EqOfStateP(QR); */

	// csiH at half node k+1/2
	csiH.x=0.5*(zet[k][j][i].x+zet[k+1][j][i].x);
	csiH.y=0.5*(zet[k][j][i].y+zet[k+1][j][i].y);
	csiH.z=0.5*(zet[k][j][i].z+zet[k+1][j][i].z);

	if (k==mz-2 && user->bctype[5]==1)
	  uconL=(csiH.x*(QR.rhoU/QR.rho)+
		 csiH.y*(QR.rhoV/QR.rho)+
		 csiH.z*(QR.rhoW/QR.rho));
	else
	  uconL=(csiH.x*(QL.rhoU/QL.rho)+
		 csiH.y*(QL.rhoV/QL.rho)+
		 csiH.z*(QL.rhoW/QL.rho));
	if (k==0 && user->bctype[4]==1)
	  uconR=uconL;
	else
	  uconR=(csiH.x*(QR.rhoU/QR.rho)+
		 csiH.y*(QR.rhoV/QR.rho)+
		 csiH.z*(QR.rhoW/QR.rho));
	gii=csiH.x*csiH.x+csiH.y*csiH.y+csiH.z*csiH.z;
	//	giiR=csiR.x*csiR.x+csiR.y*csiR.y+csiR.z*csiR.z;
	CsoundL=gamma*pL/QL.rho;
	CsoundR=gamma*pR/QR.rho;
	saL=fabs(uconL)+sqrt(fabs(CsoundL)*gii);
	saR=fabs(uconR)+sqrt(fabs(CsoundR)*gii);
	SA=PetscMax(saL,saR);

	fp3[k][j][i].rho= 0.5*( QL.rho*uconL  + QR.rho *uconR) - 0.5*SA*(QR.rho -QL.rho );
	fp3[k][j][i].rhoU=0.5*((QL.rhoU*uconL + QR.rhoU*uconR) + 
			       (pL   * csiH.x + pR * csiH.x))  - 0.5*SA*(QR.rhoU-QL.rhoU);
	fp3[k][j][i].rhoV=0.5*((QL.rhoV*uconL + QR.rhoV*uconR) + 
			       (pL   * csiH.y + pR * csiH.y))  - 0.5*SA*(QR.rhoV-QL.rhoV);
	fp3[k][j][i].rhoW=0.5*((QL.rhoW*uconL + QR.rhoW*uconR) + 
			       (pL   * csiH.z + pR * csiH.z))  - 0.5*SA*(QR.rhoW-QL.rhoW);
	fp3[k][j][i].rhoE=0.5*((pL+QL.rhoE)*uconL + 
			       (pR+QR.rhoE)*uconR)             - 0.5*SA*(QR.rhoE-QL.rhoE);
      }
    }
  }

  DMDAVecRestoreArray(cda, user->lQ, &q);

  DMDAVecRestoreArray(user->fda, user->lCsi, &csi);
  DMDAVecRestoreArray(user->fda, user->lEta, &eta);
  DMDAVecRestoreArray(user->fda, user->lZet, &zet);

  DMDAVecRestoreArray(da, user->lP, &p);

  /* DMDAVecRestoreArray(cda, Fp1, &fp1); */
  /* DMDAVecRestoreArray(cda, Fp2, &fp2); */
  /* DMDAVecRestoreArray(cda, Fp3, &fp3); */

  /* VecDuplicate(user->lQ, &lFp1); */
  /* VecDuplicate(user->lQ, &lFp2); */
  /* VecDuplicate(user->lQ, &lFp3); */

  /* DMGlobalToLocalBegin(cda,Fp1, INSERT_VALUES, lFp1); */
  /* DMGlobalToLocalEnd(cda,Fp1, INSERT_VALUES, lFp1); */
  /* DMGlobalToLocalBegin(cda,Fp2, INSERT_VALUES, lFp2); */
  /* DMGlobalToLocalEnd(cda,Fp2, INSERT_VALUES, lFp2); */
  /* DMGlobalToLocalBegin(cda,Fp3, INSERT_VALUES, lFp3); */
  /* DMGlobalToLocalEnd(cda,Fp3, INSERT_VALUES, lFp3); */

  /* VecDestroy(&Fp1); */
  /* VecDestroy(&Fp2); */
  /* VecDestroy(&Fp3); */
  
  /* DMDAVecGetArray(cda, lFp1, &fp1); */
  /* DMDAVecGetArray(cda, lFp2, &fp2); */
  /* DMDAVecGetArray(cda, lFp3, &fp3); */

  DMDAVecGetArray(cda, Conv,  &conv);
  // DMDAVecGetArray(da, user->lAj, &aj);
  
  
  /* Calculate the convective terms under cartesian coordinates */
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i]< blank) {
	conv[k][j][i].rho = (
	  fp1[k][j][i].rho - fp1[k][j][i-1].rho +
	  fp2[k][j][i].rho - fp2[k][j-1][i].rho +
	  fp3[k][j][i].rho - fp3[k-1][j][i].rho);

	if (twoD!=1)
	conv[k][j][i].rhoU = (
	  fp1[k][j][i].rhoU - fp1[k][j][i-1].rhoU +
	  fp2[k][j][i].rhoU - fp2[k][j-1][i].rhoU +
	  fp3[k][j][i].rhoU - fp3[k-1][j][i].rhoU);
	else conv[k][j][i].rhoU = 0.;

	conv[k][j][i].rhoV = (
	  fp1[k][j][i].rhoV - fp1[k][j][i-1].rhoV +
	  fp2[k][j][i].rhoV - fp2[k][j-1][i].rhoV +
	  fp3[k][j][i].rhoV - fp3[k-1][j][i].rhoV);

	conv[k][j][i].rhoW = (
	  fp1[k][j][i].rhoW - fp1[k][j][i-1].rhoW +
	  fp2[k][j][i].rhoW - fp2[k][j-1][i].rhoW +
	  fp3[k][j][i].rhoW - fp3[k-1][j][i].rhoW);

	conv[k][j][i].rhoE = (
	  fp1[k][j][i].rhoE - fp1[k][j][i-1].rhoE +
	  fp2[k][j][i].rhoE - fp2[k][j-1][i].rhoE +
	  fp3[k][j][i].rhoE - fp3[k-1][j][i].rhoE);
	} else {
	  /* if (nvert[k][j][i]<blankwall)  */
	  /* conv[k][j][i].rho = ( */
	  /* fp1[k][j][i].rho - fp1[k][j][i-1].rho + */
	  /* fp2[k][j][i].rho - fp2[k][j-1][i].rho + */
	  /* fp3[k][j][i].rho - fp3[k-1][j][i].rho); */
	  /* else */
	    conv[k][j][i].rho = 0.;
	conv[k][j][i].rhoU = 0.;
	conv[k][j][i].rhoV = 0.;
	conv[k][j][i].rhoW = 0.;
	conv[k][j][i].rhoE = 0.;
	}
	  
      }
    }
  } 

  DMDAVecRestoreArray(cda, Conv,  &conv);
  // DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  PetscReal norm;
  VecNorm(Conv, NORM_INFINITY, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "!!norm of Conv %le Flux order %i ENO order %i TwoD %d\n", norm, fluxorder, order, twoD);
  
  DMDAVecRestoreArray(cda, lFp1, &fp1);
  DMDAVecRestoreArray(cda, lFp2, &fp2);
  DMDAVecRestoreArray(cda, lFp3, &fp3);

  VecDestroy(&lFp1);
  VecDestroy(&lFp2);
  VecDestroy(&lFp3);

  return(0);
}

// Needs to be added
//PetscErrorCode filter_sigma(Vec Sigma)
//{
//}

PetscErrorCode formConvectionLES(UserCtx *user, Vec Conv, PetscInt istage)
{
  // centeral difference partial tranforamtion
  
  DM		da = user->da, cda = user->cda;
  DMDALocalInfo	info=user->info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k;
  PetscInt      fluxorder=3, order=2; // ENO order==2 or 3rd, fluxorder ==1 Godunov, 2 MUSCL, 3 WENOe
  PetscOptionsGetInt(NULL, NULL, "-ENOorder", &order, NULL);
  PetscOptionsGetInt(NULL, NULL, "-fluxorder", &fluxorder, NULL);
  PetscInt twoD=0;
  PetscOptionsGetInt(NULL, NULL, "-twoD", &twoD, NULL);
 
  PetscPrintf(PETSC_COMM_WORLD, "!!Hybrid Flux WENO order %i ENO order %i TwoD %d\n", fluxorder, order, twoD);
   
 
  //Vec		Fp1, Fp2, Fp3;
  Vec		lFp1, lFp2, lFp3, Conv_cent;
  CompVars	***fp1, ***fp2, ***fp3;
  CompVars	***conv, ***q, QR,QL, qm, qp,Fp_cent,fp_WE;
  
  PetscReal	uconR, uconL, pR, pL, ***p, pm, pp;
  PetscReal     gii, CsoundL, CsoundR, saR, saL, SA, gamma=1.4;
  PetscReal     sigma,rho_cent, rhoU_cent, rhoV_cent,rhoW_cent, rhoE_cent, Ucont_cent;
  PetscReal	***nvert,***aj, blank=0.1, blankwall=2.9;
  Cmpnts	***csi, ***eta, ***zet,***sigma1,***lsigma, csiH;// csiR, csiL;
  PetscErrorCode ierr;

  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;

  PetscInt  rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  /* VecDuplicate(user->Q, &Fp1); */
  /* VecDuplicate(user->Q, &Fp2); */
  /* VecDuplicate(user->Q, &Fp3); */
  
  /* DMDAVecGetArray(cda, Fp1, &fp1); */
  /* DMDAVecGetArray(cda, Fp2, &fp2); */
  /* DMDAVecGetArray(cda, Fp3, &fp3); */
  
  VecDuplicate(user->lQ, &lFp1);
  VecDuplicate(user->lQ, &lFp2);
  ierr = VecDuplicate(user->lQ, &lFp3); CHKERRQ(ierr);
  VecSet(lFp1, 0.);

  DMDAVecGetArray(cda, lFp1, &fp1);
  DMDAVecGetArray(cda, lFp2, &fp2);
  
  
  DMDAVecGetArray(cda, user->lQ, &q);  
  
  DMDAVecGetArray(da, user->lP, &p);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  //Amir 
  //if bcype[0]==7 =======> bctype[1]==7. 

  //Amir 
  //if bcype[0]==7 =======> bctype[1]==7. 
  // i+1/2->i

  PetscInt ABC=0; //Absorbing BC
  PetscOptionsGetInt(NULL, NULL, "-ABC", &ABC, NULL);  

  PetscReal s_y=1.0, s_z=1.0, s_x=1.0, s_y_min=.0, s_z_min=.0, s_x_min=.0; 
  //  PetscOptionsInsertFile(NULL, PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  PetscOptionsGetReal(NULL, NULL, "-s_x", &s_x, NULL);
  PetscOptionsGetReal(NULL, NULL, "-s_y", &s_y, NULL);
  PetscOptionsGetReal(NULL, NULL, "-s_z", &s_z, NULL);

  PetscOptionsGetReal(NULL, NULL, "-s_x_min", &s_x_min, NULL);
  PetscOptionsGetReal(NULL, NULL, "-s_y_min", &s_y_min, NULL);
  PetscOptionsGetReal(NULL, NULL, "-s_z_min", &s_z_min, NULL);

  ks=0.00001;
    
  PetscOptionsGetReal(NULL, NULL, "-sigma_thrsh", &ks, NULL);
  
  // Details on istage may be found in the Runge-Kutta subroutine (check-solvers.c). istage=0 is the first stage of RK4 implementation
  if (istage==0){    

    VecStrideSet(user->Sigma, 0, s_x); 
    VecStrideSet(user->Sigma, 1, s_y); 
    VecStrideSet(user->Sigma, 2, s_z); 
    

    if (calc_sigma){
      DMDAVecGetArray(user->fda, user->Sigma, &sigma1);      
      for (k=lzs; k<lze; k++){
  			for (j=lys; j<lye; j++){
  	  			for (i=lxs; i<lxe; i++){
	    
	    			if (ABC!=0 && k>ABC) sigma1[k][j][i].x=0.;
  	    			else if (i<mx-4 && i>2)
  	      			sigma1[k][j][i].x=sigma_cent(p[k][j][i-1],p[k][j][i],p[k][j][i+1],p[k][j][i+2]);
  	    			sigma1[k][j][i].x=PetscMin(s_x,    sigma1[k][j][i].x);
	    			sigma1[k][j][i].x=PetscMax(s_x_min,sigma1[k][j][i].x);	    
	    //if (nvert[k][j][i]>blank) sigma1[k][j][i].x= 0.;	    
  	  			}
  			}
      	}
      //    double M;
      
      for (k=lzs; k<lze; k++){
  			for (j=lys; j<lye; j++){
  	  			for (i=lxs; i<lxe; i++){
	    			if (ABC!=0 && k>ABC) sigma1[k][j][i].z=0.;
  	    			else if (user->bctype[4]==7 || (k<mz-4 && k>4))
  	      			sigma1[k][j][i].z=sigma_cent(p[k-1][j][i],p[k][j][i],p[k+1][j][i],p[k+2][j][i]);
	    //if (nvert[k][j][i]>blank) sigma1[k][j][i].z= 0.;
  	    			sigma1[k][j][i].z=PetscMin(s_z    ,sigma1[k][j][i].z);
	    			sigma1[k][j][i].z=PetscMax(s_z_min,sigma1[k][j][i].z);	    
  	  			}
  			}
      	}
      
      for (k=lzs; k<lze; k++){
  			for (j=lys; j<lye; j++){
  	  			for (i=lxs; i<lxe; i++){
	    			if (ABC!=0 && k>ABC) sigma1[k][j][i].y=0.;
  	    			else if (j<my-4 && j>4 )
  	      			sigma1[k][j][i].y=sigma_cent(p[k][j-1][i],p[k][j][i],p[k][j+1][i],p[k][j+2][i]);
	    //if (nvert[k][j][i]>blank) sigma1[k][j][i].y= 0.;
  	    			sigma1[k][j][i].y=PetscMin(s_y    ,sigma1[k][j][i].y);
  	    			sigma1[k][j][i].y=PetscMax(s_y_min,sigma1[k][j][i].y);	    
  	  			}
  			}
      	}
      
    PetscPrintf(PETSC_COMM_WORLD, "!!Calc Sigma done rank %d!! s_y %le s_y_min %le ABC %d \n", rank, s_y, s_y_min, ABC);
    DMDAVecRestoreArray(user->fda, user->Sigma, &sigma1);
    
    
    DMGlobalToLocalBegin(user->fda, user->Sigma, INSERT_VALUES, user->lSigma);
    DMGlobalToLocalEnd(user->fda, user->Sigma, INSERT_VALUES, user->lSigma);
    
    DMDAVecGetArray(user->fda, user->Sigma, &sigma1);
    DMDAVecGetArray(user->fda, user->lSigma, &lsigma);
    
    //    PetscPrintf(PETSC_COMM_SELF, "!!Start Filetring Sigma rank %d!!\n", rank);
    //filtering shock sensor to smooth shock sensor
    for (k=lzs; k<lze; k++){
      	for (j=lys; j<lye; j++){
    		for (i=lxs; i<lxe; i++){
	  
    	  		sigma1[k][j][i].x *= sqrt(lsigma[k][j-1][i].x*lsigma[k][j][i].x)*
  	    		sqrt(lsigma[k-1][j][i].x*lsigma[k+1][j][i].x)*
				sqrt(lsigma[k][j][i-1].x*lsigma[k][j][i+1].x);
				
				sigma1[k][j][i].y *= sqrt(lsigma[k][j-1][i].y*lsigma[k][j][i].y)*
				sqrt(lsigma[k-1][j][i].y*lsigma[k+1][j][i].y)*
				sqrt(lsigma[k][j][i-1].y*lsigma[k][j][i+1].y);
			
				sigma1[k][j][i].z *=sqrt(lsigma[k][j-1][i].z*lsigma[k][j][i].z)*
				sqrt(lsigma[k-1][j][i].z*lsigma[k+1][j][i].z)*
				sqrt(lsigma[k][j][i-1].z*lsigma[k][j][i+1].z);
    		}
      	}
    }
    //  PetscPrintf(PETSC_COMM_SELF, "!!Done Filetring Sigma rank %d!!\n", rank);
    //filter_sigma(user->lSigma);
    
    DMDAVecRestoreArray(user->fda, user->lSigma, &lsigma);
    DMDAVecRestoreArray(user->fda, user->Sigma, &sigma1);
    } // calc sigma done
    
    DMGlobalToLocalBegin(user->fda, user->Sigma, INSERT_VALUES, user->lSigma);
    DMGlobalToLocalEnd(user->fda, user->Sigma, INSERT_VALUES, user->lSigma);
	
	DMDAVecGetArray(user->fda, user->Sigma, &sigma1);
	DMDAVecGetArray(user->fda, user->lSigma, &lsigma);
	//Advait
	for (k=lzs; k<lze; k++){
		for (j=lys; j<lye; j++){
	  		for (i=lxs; i<lxe; i++){
				if(j<my-4 && j>4 && lsigma[k][j-1][i].y+lsigma[k][j+1][i].y-2.*lsigma[k][j][i].y<0.0){	
					sigma1[k][j][i].y=(lsigma[k][j-1][i].y*lsigma[k][j][i].y*lsigma[k][j][i].y);
				}
				if(k<mz-4 && k>4 && lsigma[k-1][j][i].z+lsigma[k+1][j][i].z-2.*lsigma[k][j][i].z<0.0 ){
					sigma1[k][j][i].z=(lsigma[k-1][j][i].z*lsigma[k+1][j][i].z*lsigma[k][j][i].z);
				}
				sigma1[k][j][i].z=PetscMin(sigma1[k][j][i].z,sigma1[k][j][i].y);
	    		sigma1[k][j][i].y=sigma1[k][j][i].z;
			}
		}
	}

	DMDAVecRestoreArray(user->fda, user->Sigma, &sigma1);
	DMDAVecRestoreArray(user->fda, user->lSigma, &lsigma);

	DMGlobalToLocalBegin(user->fda, user->Sigma, INSERT_VALUES, user->lSigma);
    DMGlobalToLocalEnd(user->fda, user->Sigma, INSERT_VALUES, user->lSigma);
  }//end istage==0

  DMDAVecGetArray(user->fda, user->lSigma, &lsigma); 
  
  
  /* CompVars ***conv_cent; */
  /* VecDuplicate(user->Q, &Conv_cent); */
  /* DMDAVecGetArray(cda, Conv_cent, &conv_cent); */
  //int order=3;
  
  
  
  PetscReal c[5]; c[0]=1.0/12.0; c[1]=-2.0/3.0;c[2]=0.0;c[3]=2.0/3.0;c[4]=-1.0/12.0;
  
  //	double  Ucont[5], Ucontr, Ucontl;
  PetscInt m,a,L;
  
  
  PetscReal p_cent,r[4],u[4],v[4],w[4],e[4],pp_c[4],ep[4],Ucont[4];

 
  DMDAVecGetArray(user->fda, user->lCsi, &csi);
  //	int central=3;
  // i+1/2->i
  if (twoD!=1) {
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs-1; i<lxe; i++){
	
	// 1/J x ∂ξ/∂xi -> calculate metrics at i+1/2
	csiH.x = 0.5*(csi[k][j][i].x+csi[k][j][i+1].x);
	csiH.y = 0.5*(csi[k][j][i].y+csi[k][j][i+1].y);
	csiH.z = 0.5*(csi[k][j][i].z+csi[k][j][i+1].z);
	
	// sigma-> sigma_x for i direction flux calc
	sigma=lsigma[k][j][i].x;
	
	// u0=q[k][j][i-1].rhoU
	// u1=q[k][j][i].rhoU
	// u2=q[k][j][i+1].rhoU
	// u3=q[k][j][i+2].rhoU
	//calculating at i+1/2
	if (i>0 && i<mx-2) {

	if ((nvert[k][j][i]+nvert[k][j][i+1]+nvert[k][j][i+2]+nvert[k][j][i-1])<2.9){	  
	  for (a=-1;a<3;a++){ // {a-> -1,0,1,2}
	    L=i+a; // L -> i-1, i, i+1, i+2
	    m=1+a;
		// Note that the values are central with respect to the i+1/2 node
	    r[m]=q[k][j][L].rho; // r[m]-> 4 values of rho from i-1, i, i+1, i+2 
	    u[m]=q[k][j][L].rhoU;// u[m]-> 4 values of rhoU from i-1, i, i+1, i+2
	    v[m]=q[k][j][L].rhoV;// v[m]-> 4 values of rhoV from i-1, i, i+1, i+2
	    w[m]=q[k][j][L].rhoW;// w[m]-> 4 values of rhoW from i-1, i, i+1, i+2
	    e[m]=q[k][j][L].rhoE;// e[m]-> 4 values of rhoE from i-1, i, i+1, i+2
	    pp_c[m]=p[k][j][L];  // p[m]-> 4 values of p from i-1, i, i+1, i+2
	    ep[m]=pp_c[m]+e[m];	// ep[m]-> 4 values of p+rhoE from i-1, i, i+1, i+2

		// Calculate the same array of contravariant velocity: Uₓᵢ = 1/ρ ( ρuᵢ ξᵢ ) {This loop (i-loop) calculates U_cont in csi direction}
	    Ucont[m]=1./q[k][j][L].rho*(csi[k][j][L].x*q[k][j][L].rhoU+csi[k][j][L].y*q[k][j][L].rhoV+csi[k][j][L].z*q[k][j][L].rhoW);
	  }
	  
	  // Pressure at i+1/2
	  p_cent=0.5*(p[k][j][i]+p[k][j][i+1]);
	  // rho_cent=0.5*(q[k][j][i].rho+q[k][j][i+1].rho);

	  // Calculate the Central Convective Fluxes using different interpolation schemes
	  // Note: this does not seem to follow the same central scheme as in the paper
	  if (central==2){
	    Fp_cent.rho =Fp_2(r,Ucont);
	    Fp_cent.rhoU=Fp_2(u,Ucont)+Fp4(pp_c)*csiH.x;
	    Fp_cent.rhoV=Fp_2(v,Ucont)+Fp4(pp_c)*csiH.y;
	    Fp_cent.rhoW=Fp_2(w,Ucont)+Fp4(pp_c)*csiH.z;
	    Fp_cent.rhoE=Fp_2(ep,Ucont);
	  }
	  else if (central==3){
	    Fp_cent.rho=Fp_1(r,Ucont);
	    Fp_cent.rhoU=Fp_1(u,Ucont)+Fp4(pp_c)*csiH.x;
	    Fp_cent.rhoV=Fp_1(v,Ucont)+Fp4(pp_c)*csiH.y;
	    Fp_cent.rhoW=Fp_1(w,Ucont)+Fp4(pp_c)*csiH.z;
	    Fp_cent.rhoE=Fp_1(ep,Ucont);	      
	  }
	  else {
		// Simple central approximation using averaging!
	    rho_cent=0.5*(q[k][j][i].rho+q[k][j][i+1].rho);
	    rhoU_cent=0.5*(q[k][j][i].rhoU+q[k][j][i+1].rhoU);
	    rhoV_cent=0.5*(q[k][j][i].rhoV+q[k][j][i+1].rhoV);
	    rhoW_cent=0.5*(q[k][j][i].rhoW+q[k][j][i+1].rhoW);
	    rhoE_cent=0.5*(q[k][j][i].rhoE+q[k][j][i+1].rhoE);
	    Ucont_cent=(csiH.x*rhoU_cent+csiH.y*rhoV_cent+csiH.z*rhoW_cent)/rho_cent;
	    
	    Fp_cent.rho=        rho_cent*Ucont_cent;
	    Fp_cent.rhoU =	rhoU_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j][i+1])*csiH.x;
	    Fp_cent.rhoV =	rhoV_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j][i+1])*csiH.y;
	    Fp_cent.rhoW =	rhoW_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j][i+1])*csiH.z;
	    Fp_cent.rhoE =	Ucont_cent*(rhoE_cent+0.5*(p[k][j][i]+p[k][j][i+1]));	    	    
	}
	  
	} else { // if nvert> wall
	  rho_cent =0.5*(q[k][j][i].rho +q[k][j][i+1].rho);
	  rhoU_cent=0.5*(q[k][j][i].rhoU+q[k][j][i+1].rhoU);
	  rhoV_cent=0.5*(q[k][j][i].rhoV+q[k][j][i+1].rhoV);
	  rhoW_cent=0.5*(q[k][j][i].rhoW+q[k][j][i+1].rhoW);
	  rhoE_cent=0.5*(q[k][j][i].rhoE+q[k][j][i+1].rhoE);
	  Ucont_cent=(csiH.x*rhoU_cent+csiH.y*rhoV_cent+csiH.z*rhoW_cent)/rho_cent;
	  
	  Fp_cent.rho  =      rho_cent*Ucont_cent;
	  Fp_cent.rhoU =      rhoU_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j][i+1])*csiH.x;
	  Fp_cent.rhoV =      rhoV_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j][i+1])*csiH.y;
	  Fp_cent.rhoW =      rhoW_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j][i+1])*csiH.z;
	  Fp_cent.rhoE =      Ucont_cent*(rhoE_cent+0.5*(p[k][j][i]+p[k][j][i+1]));	  	  
	}	
	  	  
	}   
	  else { // i=0 or >=mx-2
	  rho_cent =0.5*(q[k][j][i].rho +q[k][j][i+1].rho);
	  rhoU_cent=0.5*(q[k][j][i].rhoU+q[k][j][i+1].rhoU);
	  rhoV_cent=0.5*(q[k][j][i].rhoV+q[k][j][i+1].rhoV);
	  rhoW_cent=0.5*(q[k][j][i].rhoW+q[k][j][i+1].rhoW);
	  rhoE_cent=0.5*(q[k][j][i].rhoE+q[k][j][i+1].rhoE);
	  Ucont_cent=(csiH.x*rhoU_cent+csiH.y*rhoV_cent+csiH.z*rhoW_cent)/rho_cent;
	  
	  Fp_cent.rho  =      rho_cent*Ucont_cent;
	  Fp_cent.rhoU =      rhoU_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j][i+1])*csiH.x;
	  Fp_cent.rhoV =      rhoV_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j][i+1])*csiH.y;
	  Fp_cent.rhoW =      rhoW_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j][i+1])*csiH.z;
	  Fp_cent.rhoE =      Ucont_cent*(rhoE_cent+0.5*(p[k][j][i]+p[k][j][i+1]));	  	  
	}	
	// Done Calculating Central Fluxes!
	
	// calc F1 on i+1/2 nodes... i+1/2 saved in i 
	//-Lax-Freidrich flux: F(QR,QL)=1/2(E(QR)+E(QL))-1/2|S(A)|(QR-QL)
	if (fluxorder==1 || ((i==0 || i==mx-2) && user->bctype[0]!=7) ||
	    ((nvert[k][j][i]+nvert[k][j][i+1])>blankwall) ) { // if i or i+1 are wall nodes
	  QL=q[k][j][i];
	  QR=q[k][j][i+1];
	  pL=p[k][j][i];
	  pR=p[k][j][i+1];
	} else if (fluxorder==2) { // MUSCL
	  /* we know i and i+1 are not wall nodes but i-1 or i+2 might be wall nodes  */
	  if (nvert[k][j][i-1]>blankwall) {
	    /* put Qm equal to Qp to revert WENO to ENO */
	    qm=q[k][j][i+1];
	    pm=p[k][j][i+1];
	  } else {
	    qm=q[k][j][i-1];
	    pm=p[k][j][i-1];
	  }
	  if (nvert[k][j][i+2]>blankwall) {
	    /* put Qp equal to Qm to revert WENO to ENO */
	    qp=q[k][j][i];
	    pp=p[k][j][i];
	  } else {
	    qp=q[k][j][i+2];
	    pp=p[k][j][i+2];
	  }
	  
	  QL.rho =MUSCL_QL(qm.rho , q[k][j][i].rho , q[k][j][i+1].rho);
	  QL.rhoU=MUSCL_QL(qm.rhoU, q[k][j][i].rhoU, q[k][j][i+1].rhoU);
	  QL.rhoV=MUSCL_QL(qm.rhoV, q[k][j][i].rhoV, q[k][j][i+1].rhoV);
	  QL.rhoW=MUSCL_QL(qm.rhoW, q[k][j][i].rhoW, q[k][j][i+1].rhoW);
	  QL.rhoE=MUSCL_QL(qm.rhoE, q[k][j][i].rhoE, q[k][j][i+1].rhoE);
	  pL     =MUSCL_QL(pm     , p[k][j][i],      p[k][j][i+1]);
	  
	  QR.rho =MUSCL_QR(q[k][j][i].rho , q[k][j][i+1].rho , qp.rho);
	  QR.rhoU=MUSCL_QR(q[k][j][i].rhoU, q[k][j][i+1].rhoU, qp.rhoU);
	  QR.rhoV=MUSCL_QR(q[k][j][i].rhoV, q[k][j][i+1].rhoV, qp.rhoV);
	  QR.rhoW=MUSCL_QR(q[k][j][i].rhoW, q[k][j][i+1].rhoW, qp.rhoW);
	  QR.rhoE=MUSCL_QR(q[k][j][i].rhoE, q[k][j][i+1].rhoE, qp.rhoE);
	  pR     =MUSCL_QR(p[k][j][i]     , p[k][j][i+1],      pp);
	} else { // WENO
	  /* we know i and i+1 are not wall nodes but i-1 or i+2 might be wall nodes  */
	  if (nvert[k][j][i-1]>blankwall && (i+2)<mx) {
	    
	    QL.rho =ENO_QL(q[k][j][i].rho , q[k][j][i+1].rho , q[k][j][i+2].rho ,order);
	    QL.rhoU=ENO_QL(q[k][j][i].rhoU, q[k][j][i+1].rhoU, q[k][j][i+2].rhoU,order);
	    QL.rhoV=ENO_QL(q[k][j][i].rhoV, q[k][j][i+1].rhoV, q[k][j][i+2].rhoV,order);
	    QL.rhoW=ENO_QL(q[k][j][i].rhoW, q[k][j][i+1].rhoW, q[k][j][i+2].rhoW,order);
	    QL.rhoE=ENO_QL(q[k][j][i].rhoE, q[k][j][i+1].rhoE, q[k][j][i+2].rhoE,order);
	    pL     =ENO_QL(p[k][j][i]     , p[k][j][i+1],      p[k][j][i+2]     ,order);	 
	    
	  } else {
	    qm=q[k][j][i-1];
	    pm=p[k][j][i-1];       	 
            
	    QL.rho =WENO_QL(qm.rho , q[k][j][i].rho , q[k][j][i+1].rho);
	    QL.rhoU=WENO_QL(qm.rhoU, q[k][j][i].rhoU, q[k][j][i+1].rhoU);
	    QL.rhoV=WENO_QL(qm.rhoV, q[k][j][i].rhoV, q[k][j][i+1].rhoV);
	    QL.rhoW=WENO_QL(qm.rhoW, q[k][j][i].rhoW, q[k][j][i+1].rhoW);
	    QL.rhoE=WENO_QL(qm.rhoE, q[k][j][i].rhoE, q[k][j][i+1].rhoE);
	    pL     =WENO_QL(pm     , p[k][j][i],      p[k][j][i+1]);
	  }
	  
	  if (nvert[k][j][i+2]>blankwall) {
	    /* put Qp equal to Qm to revert WENO to ENO */
	    qp=q[k][j][i+1];
	    pp=p[k][j][i+1];
	    QR.rho =ENO_QR(q[k][j][i-1].rho , q[k][j][i].rho , qp.rho ,order);
	    QR.rhoU=ENO_QR(q[k][j][i-1].rhoU, q[k][j][i].rhoU, qp.rhoU,order);
	    QR.rhoV=ENO_QR(q[k][j][i-1].rhoV, q[k][j][i].rhoV, qp.rhoV,order);
	    QR.rhoW=ENO_QR(q[k][j][i-1].rhoW, q[k][j][i].rhoW, qp.rhoW,order);
	    QR.rhoE=ENO_QR(q[k][j][i-1].rhoE, q[k][j][i].rhoE, qp.rhoE,order);
	    pR     =ENO_QR(p[k][j][i-1]     , p[k][j][i],      pp     ,order);
	    /*   qp=q[k][j][i]; */
	    /*   pp=p[k][j][i]; */
	    /* QR.rho =MUSCL_QR(q[k][j][i].rho , q[k][j][i+1].rho , qp.rho); */
	    /* QR.rhoU=MUSCL_QR(q[k][j][i].rhoU, q[k][j][i+1].rhoU, qp.rhoU); */
	    /* QR.rhoV=MUSCL_QR(q[k][j][i].rhoV, q[k][j][i+1].rhoV, qp.rhoV); */
	    /* QR.rhoW=MUSCL_QR(q[k][j][i].rhoW, q[k][j][i+1].rhoW, qp.rhoW); */
	    /* QR.rhoE=MUSCL_QR(q[k][j][i].rhoE, q[k][j][i+1].rhoE, qp.rhoE); */
	    /* pR     =MUSCL_QR(p[k][j][i]     , p[k][j][i+1],      pp); */
	  } else {
	    qp=q[k][j][i+2];
	    pp=p[k][j][i+2];
            
	    
	    QR.rho =WENO_QR(q[k][j][i].rho , q[k][j][i+1].rho , qp.rho);
	    QR.rhoU=WENO_QR(q[k][j][i].rhoU, q[k][j][i+1].rhoU, qp.rhoU);
	    QR.rhoV=WENO_QR(q[k][j][i].rhoV, q[k][j][i+1].rhoV, qp.rhoV);
	    QR.rhoW=WENO_QR(q[k][j][i].rhoW, q[k][j][i+1].rhoW, qp.rhoW);
	    QR.rhoE=WENO_QR(q[k][j][i].rhoE, q[k][j][i+1].rhoE, qp.rhoE);
	    pR     =WENO_QR(p[k][j][i]     , p[k][j][i+1],      pp);
	  }
	}
	
	
	//csi at half point i+1/2
	//	csiL=csi[k][j][i];
	//	csiR=csi[k][j][i+1];
	/* csiH.x = 0.5*(csi[k][j][i].x+csi[k][j][i+1].x); */
	/* csiH.y = 0.5*(csi[k][j][i].y+csi[k][j][i+1].y); */
	/* csiH.z = 0.5*(csi[k][j][i].z+csi[k][j][i+1].z); */
	
	uconL=(csiH.x*(QL.rhoU/QL.rho)+
	       csiH.y*(QL.rhoV/QL.rho)+
	       csiH.z*(QL.rhoW/QL.rho));
	uconR=(csiH.x*(QR.rhoU/QR.rho)+
	       csiH.y*(QR.rhoV/QR.rho)+
	       csiH.z*(QR.rhoW/QR.rho));
	gii=csiH.x*csiH.x+csiH.y*csiH.y+csiH.z*csiH.z;
	//	giiR=csiR.x*csiR.x+csiR.y*csiR.y+csiR.z*csiR.z;
	CsoundL=gamma*pL/QL.rho;
	CsoundR=gamma*pR/QR.rho;
	saL=fabs(uconL)+sqrt(fabs(CsoundL)*gii);
	saR=fabs(uconR)+sqrt(fabs(CsoundR)*gii);
	SA=PetscMax(saL,saR);		
	
	fp_WE.rho= 0.5*( QL.rho*uconL  + QR.rho *uconR) - 0.5*SA*(QR.rho-QL.rho);
	fp_WE.rhoU=0.5*((QL.rhoU*uconL + QR.rhoU*uconR) + 
			(pL   * csiH.x + pR * csiH.x))  - 0.5*SA*(QR.rhoU-QL.rhoU);
	fp_WE.rhoV=0.5*((QL.rhoV*uconL + QR.rhoV*uconR) + 
			(pL   * csiH.y + pR * csiH.y))  - 0.5*SA*(QR.rhoV-QL.rhoV);
	fp_WE.rhoW=0.5*((QL.rhoW*uconL + QR.rhoW*uconR) + 
			(pL   * csiH.z + pR * csiH.z))  - 0.5*SA*(QR.rhoW-QL.rhoW);
	fp_WE.rhoE=0.5*((pL+QL.rhoE)*uconL + 
			(pR+QR.rhoE)*uconR)             - 0.5*SA*(QR.rhoE-QL.rhoE);

	/* hybrid scheme 2nd order central  */
	fp1[k][j][i].rho =fp_WE.rho *(1.0-sigma)+Fp_cent.rho *sigma;
	fp1[k][j][i].rhoU=fp_WE.rhoU*(1.0-sigma)+Fp_cent.rhoU*sigma;
	fp1[k][j][i].rhoV=fp_WE.rhoV*(1.0-sigma)+Fp_cent.rhoV*sigma;
	fp1[k][j][i].rhoW=fp_WE.rhoW*(1.0-sigma)+Fp_cent.rhoW*sigma;
	fp1[k][j][i].rhoE=fp_WE.rhoE*(1.0-sigma)+Fp_cent.rhoE*sigma;
      }
    }
  }
  } // if !twoD
  DMDAVecRestoreArray(user->fda, user->lCsi, &csi);
  //    PetscPrintf(PETSC_COMM_SELF, "!!Central Done i rank %d!!\n", rank);
  
  
  DMDAVecGetArray(user->fda, user->lEta, &eta);
  /* j direction */
  // j+1/2 -> j
  for (k=lzs; k<lze; k++){
    for (j=lys-1; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	// calc F2 on j+1/2 nodes... j+1/2 saved in j
	//	central=2;
	
	csiH.x = 0.5*(eta[k][j][i].x+eta[k][j+1][i].x);
	csiH.y = 0.5*(eta[k][j][i].y+eta[k][j+1][i].y);
	csiH.z = 0.5*(eta[k][j][i].z+eta[k][j+1][i].z);
	
	sigma=lsigma[k][j][i].y;
	
	// u0=q[k][j][i-1].rhoU
	// u1=q[k][j][i].rhoU
	// u2=q[k][j][i+1].rhoU
	// u3=q[k][j][i+2].rhoU
	//calculating at j+1/2
	if (j>0 && j<my-2){ 
	if ((nvert[k][j-1][i]+nvert[k][j][i]+nvert[k][j+1][i]+nvert[k][j+2][i])<2.9) {
	  
	  for (a=-1;a<3;a++){
	    L=j+a;
	    m=1+a;
	    r[m]=q[k][L][i].rho;
	    u[m]=q[k][L][i].rhoU;
	    v[m]=q[k][L][i].rhoV;
	    w[m]=q[k][L][i].rhoW;
	    e[m]=q[k][L][i].rhoE;
	    pp_c[m]=p[k][L][i];
	    ep[m]=pp_c[m]+e[m];
	    Ucont[m]=1./q[k][L][i].rho*(eta[k][L][i].x*q[k][L][i].rhoU+eta[k][L][i].y*q[k][L][i].rhoV+eta[k][L][i].z*q[k][L][i].rhoW);
	  }
	  p_cent=(p[k][j][i]+p[k][j+1][i])*.5;
	  if (central==2){
	    Fp_cent.rho= Fp_2(r,Ucont);
	    Fp_cent.rhoU=Fp_2(u,Ucont)+Fp4(pp_c)*csiH.x;
	    Fp_cent.rhoV=Fp_2(v,Ucont)+Fp4(pp_c)*csiH.y;
	    Fp_cent.rhoW=Fp_2(w,Ucont)+Fp4(pp_c)*csiH.z;
	    Fp_cent.rhoE=Fp_2(ep,Ucont);
	  }
	  else if (central==3){
	    Fp_cent.rho =Fp_1(r,Ucont);
	    Fp_cent.rhoU=Fp_1(u,Ucont)+Fp4(pp_c)*csiH.x;
	    Fp_cent.rhoV=Fp_1(v,Ucont)+Fp4(pp_c)*csiH.y;
	    Fp_cent.rhoW=Fp_1(w,Ucont)+Fp4(pp_c)*csiH.z;
	    Fp_cent.rhoE=Fp_1(ep,Ucont);
	  }else{
	    rho_cent =0.5*(q[k][j][i].rho +q[k][j+1][i].rho);
	    rhoU_cent=0.5*(q[k][j][i].rhoU+q[k][j+1][i].rhoU);
	    rhoV_cent=0.5*(q[k][j][i].rhoV+q[k][j+1][i].rhoV);
	    rhoW_cent=0.5*(q[k][j][i].rhoW+q[k][j+1][i].rhoW);
	    rhoE_cent=0.5*(q[k][j][i].rhoE+q[k][j+1][i].rhoE);
	    Ucont_cent=(csiH.x*rhoU_cent+csiH.y*rhoV_cent+csiH.z*rhoW_cent)/rho_cent;
	    
	    Fp_cent.rho  =      rho_cent*Ucont_cent;
	    Fp_cent.rhoU =	rhoU_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j+1][i])*csiH.x;
	    Fp_cent.rhoV =	rhoV_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j+1][i])*csiH.y;
	    Fp_cent.rhoW =	rhoW_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j+1][i])*csiH.z;
	    Fp_cent.rhoE =	Ucont_cent*(rhoE_cent+0.5*(p[k][j][i]+p[k][j+1][i]));
	    
	    
	  }
	  
	}
	else{ // if nvert > wall
	  rho_cent =0.5*(q[k][j][i].rho +q[k][j+1][i].rho);
	  rhoU_cent=0.5*(q[k][j][i].rhoU+q[k][j+1][i].rhoU);
	  rhoV_cent=0.5*(q[k][j][i].rhoV+q[k][j+1][i].rhoV);
	  rhoW_cent=0.5*(q[k][j][i].rhoW+q[k][j+1][i].rhoW);
	  rhoE_cent=0.5*(q[k][j][i].rhoE+q[k][j+1][i].rhoE);
	  Ucont_cent=(csiH.x*rhoU_cent+csiH.y*rhoV_cent+csiH.z*rhoW_cent)/rho_cent;
	  
	  Fp_cent.rho  = rho_cent*Ucont_cent;
	  Fp_cent.rhoU = rhoU_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j+1][i])*csiH.x;
	  Fp_cent.rhoV = rhoV_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j+1][i])*csiH.y;
	  Fp_cent.rhoW = rhoW_cent*Ucont_cent+0.5*(p[k][j][i]+p[k][j+1][i])*csiH.z;
	  Fp_cent.rhoE = Ucont_cent*(rhoE_cent+0.5*(p[k][j][i]+p[k][j+1][i]));  
	}

	} else{ // if j=0 or >=my-2
	  rho_cent =0.5*(q[k][j][i].rho +q[k][j+1][i].rho);
	  rhoU_cent=0.5*(q[k][j][i].rhoU+q[k][j+1][i].rhoU);
	  rhoV_cent=0.5*(q[k][j][i].rhoV+q[k][j+1][i].rhoV);
	  rhoW_cent=0.5*(q[k][j][i].rhoW+q[k][j+1][i].rhoW);
	  rhoE_cent=0.5*(q[k][j][i].rhoE+q[k][j+1][i].rhoE);
	  Ucont_cent=(csiH.x*rhoU_cent+csiH.y*rhoV_cent+csiH.z*rhoW_cent)/rho_cent;
	  
	  Fp_cent.rho  = rho_cent *Ucont_cent;
	  Fp_cent.rhoU = rhoU_cent*Ucont_cent +0.5*(p[k][j][i]+p[k][j+1][i])*csiH.x;
	  Fp_cent.rhoV = rhoV_cent*Ucont_cent +0.5*(p[k][j][i]+p[k][j+1][i])*csiH.y;
	  Fp_cent.rhoW = rhoW_cent*Ucont_cent +0.5*(p[k][j][i]+p[k][j+1][i])*csiH.z;
	  Fp_cent.rhoE = Ucont_cent*(rhoE_cent+0.5*(p[k][j][i]+p[k][j+1][i]));	  
	}
	
	//-Lax-Freidrich flux: F(QR,QL)=1/2(E(QR)+E(QL))-1/2|S(A)|(QR-QL)
	if (fluxorder==1 || ((j==0 || j==my-2) && user->bctype[2]!=7) ||
	    ((nvert[k][j][i]+nvert[k][j+1][i])>blankwall) ) { // if j or j+1 are wall nodes
	  QL=q[k][j][i];
	  QR=q[k][j+1][i];
	  pL=p[k][j][i];
	  pR=p[k][j+1][i];
	} else if (fluxorder==2) {
	  /* we know j and j+1 are not wall nodes but j-1 or j+2 might be wall nodes  */
	  if (nvert[k][j-1][i]>blankwall) {
	    /* put Qm equal to Qp to revert WENO to ENO */
	    qm=q[k][j+1][i];
	    pm=p[k][j+1][i];
	  } else {
	    qm=q[k][j-1][i];
	    pm=p[k][j-1][i];
	  }
	  if (nvert[k][j+2][i]>blankwall) {
	    /* put Qp equal to Qm to revert WENO to ENO */
	    qp=q[k][j][i];
	    pp=p[k][j][i];
	  } else {
	    qp=q[k][j+2][i];
	    pp=p[k][j+2][i];
	  }
	  QL.rho =MUSCL_QL(qm.rho , q[k][j][i].rho , q[k][j+1][i].rho);
	  QL.rhoU=MUSCL_QL(qm.rhoU, q[k][j][i].rhoU, q[k][j+1][i].rhoU);
	  QL.rhoV=MUSCL_QL(qm.rhoV, q[k][j][i].rhoV, q[k][j+1][i].rhoV);
	  QL.rhoW=MUSCL_QL(qm.rhoW, q[k][j][i].rhoW, q[k][j+1][i].rhoW);
	  QL.rhoE=MUSCL_QL(qm.rhoE, q[k][j][i].rhoE, q[k][j+1][i].rhoE);
	  pL     =MUSCL_QL(pm     , p[k][j][i],      p[k][j+1][i]);
	  
	  QR.rho =MUSCL_QR(q[k][j][i].rho , q[k][j+1][i].rho , qp.rho);
	  QR.rhoU=MUSCL_QR(q[k][j][i].rhoU, q[k][j+1][i].rhoU, qp.rhoU);
	  QR.rhoV=MUSCL_QR(q[k][j][i].rhoV, q[k][j+1][i].rhoV, qp.rhoV);
	  QR.rhoW=MUSCL_QR(q[k][j][i].rhoW, q[k][j+1][i].rhoW, qp.rhoW);
	  QR.rhoE=MUSCL_QR(q[k][j][i].rhoE, q[k][j+1][i].rhoE, qp.rhoE);
	  pR     =MUSCL_QR(p[k][j][i]     , p[k][j+1][i],      pp);
	  
	} else { //WENO
	  /* we know j and j+1 are not wall nodes but j-1 or j+2 might be wall nodes  */
	  if (nvert[k][j-1][i]>blankwall) {
	    /* put Qm equal to Qp to revert WENO to ENO */
	    /* qm=q[k][j+1][i]; */
	    /* pm=p[k][j+1][i]; */
	    QL.rho =ENO_QL(q[k][j][i].rho , q[k][j+1][i].rho , q[k][j+2][i].rho ,order);
	    QL.rhoU=ENO_QL(q[k][j][i].rhoU, q[k][j+1][i].rhoU, q[k][j+2][i].rhoU,order);
	    QL.rhoV=ENO_QL(q[k][j][i].rhoV, q[k][j+1][i].rhoV, q[k][j+2][i].rhoV,order);
	    QL.rhoW=ENO_QL(q[k][j][i].rhoW, q[k][j+1][i].rhoW, q[k][j+2][i].rhoW,order);
	    QL.rhoE=ENO_QL(q[k][j][i].rhoE, q[k][j+1][i].rhoE, q[k][j+2][i].rhoE,order);
	    pL     =ENO_QL(p[k][j][i]     , p[k][j+1][i],      p[k][j+2][i]     ,order);
	  } else {
	    qm=q[k][j-1][i];
	    pm=p[k][j-1][i];
	    
	    QL.rho =WENO_QL(qm.rho , q[k][j][i].rho , q[k][j+1][i].rho);
	    QL.rhoU=WENO_QL(qm.rhoU, q[k][j][i].rhoU, q[k][j+1][i].rhoU);
	    QL.rhoV=WENO_QL(qm.rhoV, q[k][j][i].rhoV, q[k][j+1][i].rhoV);
	    QL.rhoW=WENO_QL(qm.rhoW, q[k][j][i].rhoW, q[k][j+1][i].rhoW);
	    QL.rhoE=WENO_QL(qm.rhoE, q[k][j][i].rhoE, q[k][j+1][i].rhoE);
	    pL     =WENO_QL(pm     , p[k][j][i],      p[k][j+1][i]);
	  }
	  if (nvert[k][j+2][i]>blankwall) {
	    /* put Qp equal to Qm to revert WENO to ENO */
	    /* qp=q[k][j][i]; */
	    /* pp=p[k][j][i]; */
	    qp=q[k][j+1][i]; 
	    pp=p[k][j+1][i];
	    QR.rho =ENO_QR(q[k][j-1][i].rho , q[k][j][i].rho , qp.rho ,order);
	    QR.rhoU=ENO_QR(q[k][j-1][i].rhoU, q[k][j][i].rhoU, qp.rhoU,order);
	    QR.rhoV=ENO_QR(q[k][j-1][i].rhoV, q[k][j][i].rhoV, qp.rhoV,order);
	    QR.rhoW=ENO_QR(q[k][j-1][i].rhoW, q[k][j][i].rhoW, qp.rhoW,order);
	    QR.rhoE=ENO_QR(q[k][j-1][i].rhoE, q[k][j][i].rhoE, qp.rhoE,order);
	    pR     =ENO_QR(p[k][j-1][i]     , p[k][j][i],      pp     ,order);
	  } else {
	    qp=q[k][j+2][i];
	    pp=p[k][j+2][i];
	    
	    QR.rho =WENO_QR(q[k][j][i].rho , q[k][j+1][i].rho , qp.rho);
	    QR.rhoU=WENO_QR(q[k][j][i].rhoU, q[k][j+1][i].rhoU, qp.rhoU);
	    QR.rhoV=WENO_QR(q[k][j][i].rhoV, q[k][j+1][i].rhoV, qp.rhoV);
	    QR.rhoW=WENO_QR(q[k][j][i].rhoW, q[k][j+1][i].rhoW, qp.rhoW);
	    QR.rhoE=WENO_QR(q[k][j][i].rhoE, q[k][j+1][i].rhoE, qp.rhoE);
	    pR     =WENO_QR(p[k][j][i]     , p[k][j+1][i],      pp);
	  }
          
	}
	
	//csi at half point j+1/2
	//	csiL=eta[k][j][i];
	//	csiR=eta[k][j+1][i];
	/* csiH.x=0.5*(eta[k][j][i].x+eta[k][j+1][i].x); */
	/* csiH.y=0.5*(eta[k][j][i].y+eta[k][j+1][i].y); */
	/* csiH.z=0.5*(eta[k][j][i].z+eta[k][j+1][i].z); */
	
	uconL=(csiH.x*(QL.rhoU/QL.rho)+
	       csiH.y*(QL.rhoV/QL.rho)+
	       csiH.z*(QL.rhoW/QL.rho));
	uconR=(csiH.x*(QR.rhoU/QR.rho)+
	       csiH.y*(QR.rhoV/QR.rho)+
	       csiH.z*(QR.rhoW/QR.rho));
	gii=csiH.x*csiH.x+csiH.y*csiH.y+csiH.z*csiH.z;
	//	giiR=csiR.x*csiR.x+csiR.y*csiR.y+csiR.z*csiR.z;
	CsoundL=gamma*pL/QL.rho;
	CsoundR=gamma*pR/QR.rho;
	saL=fabs(uconL)+sqrt(fabs(CsoundL)*gii);
	saR=fabs(uconR)+sqrt(fabs(CsoundR)*gii);
	SA=PetscMax(saL,saR);
		
	fp_WE.rho= 0.5*( QL.rho*uconL  + QR.rho *uconR) - 0.5*SA*(QR.rho-QL.rho);
	
	fp_WE.rhoU=0.5*((QL.rhoU*uconL + QR.rhoU*uconR) + 
			(pL   * csiH.x + pR * csiH.x))  - 0.5*SA*(QR.rhoU-QL.rhoU);
	fp_WE.rhoV=0.5*((QL.rhoV*uconL + QR.rhoV*uconR) + 
			(pL   * csiH.y + pR * csiH.y))  - 0.5*SA*(QR.rhoV-QL.rhoV);
	fp_WE.rhoW=0.5*((QL.rhoW*uconL + QR.rhoW*uconR) + 
			(pL   * csiH.z + pR * csiH.z))  - 0.5*SA*(QR.rhoW-QL.rhoW);
	fp_WE.rhoE=0.5*((pL+QL.rhoE)*uconL + 
			(pR+QR.rhoE)*uconR)             - 0.5*SA*(QR.rhoE-QL.rhoE);
	
	/* hybrid scheme 4nd order central  */								       
	fp2[k][j][i].rho = (1.-sigma)*fp_WE.rho +sigma*Fp_cent.rho;
	fp2[k][j][i].rhoU= (1.-sigma)*fp_WE.rhoU+sigma*Fp_cent.rhoU;
	fp2[k][j][i].rhoV= (1.-sigma)*fp_WE.rhoV+sigma*Fp_cent.rhoV;
	fp2[k][j][i].rhoW= (1.-sigma)*fp_WE.rhoW+sigma*Fp_cent.rhoW;
	fp2[k][j][i].rhoE= (1.-sigma)*fp_WE.rhoE+sigma*Fp_cent.rhoE;
					
      }
    }
  }

  DMDAVecRestoreArray(user->fda, user->lEta, &eta);
  //  PetscBarrier(NULL);
  //    PetscPrintf(PETSC_COMM_SELF, "!!Central Done j rank %d lzs %d lze %d lys %d lye %d lxs %d lxe %d !!\n", rank, lzs, lze, lys, lye, lxs, lxe);

  ierr = DMDAVecGetArray(user->fda, user->lZet, &zet); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(cda, lFp3, &fp3); CHKERRQ(ierr);
  /* k direction  */
  // k+1/2 -> k
  for (k=lzs-1; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++){
	
	// central=1;
	//	if (k==16) PetscPrintf(PETSC_COMM_SELF, "!! rank %d k %d j %d i %d lye %d lxs %d lxe %d !!\n", rank, k, j, i, lye, lxs, lxe);
	csiH.x = 0.5*(zet[k][j][i].x+zet[k+1][j][i].x);
	csiH.y = 0.5*(zet[k][j][i].y+zet[k+1][j][i].y);
	csiH.z = 0.5*(zet[k][j][i].z+zet[k+1][j][i].z);
	
	sigma=lsigma[k][j][i].z;
	
	// u0=q[k][j][i-1].rhoU
	// u1=q[k][j][i].rhoU
	// u2=q[k][j][i+1].rhoU
	// u3=q[k][j][i+2].rhoU
	//calculating at k+1/2
	if (k>0 && k<mz-2) {
	if ((nvert[k-1][j][i]+nvert[k][j][i]+nvert[k+1][j][i]+nvert[k+2][j][i])<2.9){
	  
	  for (a=-1;a<3;a++){
	    L=k+a;
	    m=1+a;
	    r[m]=q[L][j][i].rho;
	    u[m]=q[L][j][i].rhoU;
	    v[m]=q[L][j][i].rhoV;
	    w[m]=q[L][j][i].rhoW;
	    e[m]=q[L][j][i].rhoE;
	    pp_c[m]=p[L][j][i];
	    ep[m]=e[m]+pp_c[m];
	    Ucont[m]=1./q[L][j][i].rho*(zet[L][j][i].x*q[L][j][i].rhoU+zet[L][j][i].y*q[L][j][i].rhoV+zet[L][j][i].z*q[L][j][i].rhoW);
	  }

	  p_cent=0.5*(p[k][j][i]+p[k+1][j][i]);
	  if (central==2){
	    Fp_cent.rho= Fp_2(r,Ucont);
	    Fp_cent.rhoU=Fp_2(u,Ucont)+Fp4(pp_c)*csiH.x;
	    Fp_cent.rhoV=Fp_2(v,Ucont)+Fp4(pp_c)*csiH.y;
	    Fp_cent.rhoW=Fp_2(w,Ucont)+Fp4(pp_c)*csiH.z;
	    Fp_cent.rhoE=Fp_2(ep,Ucont);
	  } else if (central==3){
	    Fp_cent.rho =Fp_1(r,Ucont);
	    Fp_cent.rhoU=Fp_1(u,Ucont)+Fp4(pp_c)*csiH.x;
	    Fp_cent.rhoV=Fp_1(v,Ucont)+Fp4(pp_c)*csiH.y;
	    Fp_cent.rhoW=Fp_1(w,Ucont)+Fp4(pp_c)*csiH.z;
	    Fp_cent.rhoE=Fp_1(ep,Ucont);
	    /* if (i==1 && j==1 && k==1) */
	    /*   PetscPrintf(PETSC_COMM_SELF, "k  %d  Fp1_cent    %le\n",k,Fp_1(w,w)); */
	  } else {
	    rho_cent =0.5*(q[k][j][i].rho +q[k+1][j][i].rho);
	    rhoU_cent=0.5*(q[k][j][i].rhoU+q[k+1][j][i].rhoU);
	    rhoV_cent=0.5*(q[k][j][i].rhoV+q[k+1][j][i].rhoV);
	    rhoW_cent=0.5*(q[k][j][i].rhoW+q[k+1][j][i].rhoW);
	    rhoE_cent=0.5*(q[k][j][i].rhoE+q[k+1][j][i].rhoE);
	    Ucont_cent=(csiH.x*rhoU_cent+csiH.y*rhoV_cent+csiH.z*rhoW_cent)/rho_cent;
	    
	    Fp_cent.rho=        rho_cent  *Ucont_cent;
	    Fp_cent.rhoU =	rhoU_cent *Ucont_cent+0.5*(p[k][j][i]+p[k+1][j][i])*csiH.x;
	    Fp_cent.rhoV =	rhoV_cent *Ucont_cent+0.5*(p[k][j][i]+p[k+1][j][i])*csiH.y;
	    Fp_cent.rhoW =	rhoW_cent *Ucont_cent+0.5*(p[k][j][i]+p[k+1][j][i])*csiH.z;
	    Fp_cent.rhoE =	Ucont_cent*(rhoE_cent+0.5*(p[k][j][i]+p[k+1][j][i]));
	  }
	  
	}
	else{ //  if nvert> wall
	  rho_cent =0.5*(q[k][j][i].rho +q[k+1][j][i].rho);
	  rhoU_cent=0.5*(q[k][j][i].rhoU+q[k+1][j][i].rhoU);
	  rhoV_cent=0.5*(q[k][j][i].rhoV+q[k+1][j][i].rhoV);
	  rhoW_cent=0.5*(q[k][j][i].rhoW+q[k+1][j][i].rhoW);
	  rhoE_cent=0.5*(q[k][j][i].rhoE+q[k+1][j][i].rhoE);
	  Ucont_cent=(csiH.x*rhoU_cent+csiH.y*rhoV_cent+csiH.z*rhoW_cent)/rho_cent;
	  
	  Fp_cent.rho  =        rho_cent  *Ucont_cent;
	  Fp_cent.rhoU =	rhoU_cent *Ucont_cent+0.5*(p[k][j][i]+p[k+1][j][i])*csiH.x;
	  Fp_cent.rhoV =	rhoV_cent *Ucont_cent+0.5*(p[k][j][i]+p[k+1][j][i])*csiH.y;
	  Fp_cent.rhoW =	rhoW_cent *Ucont_cent+0.5*(p[k][j][i]+p[k+1][j][i])*csiH.z;
	  Fp_cent.rhoE =	Ucont_cent*(rhoE_cent+0.5*(p[k][j][i]+p[k+1][j][i]));
	}

	} else { // if k=0 or = mz-2
	  rho_cent =0.5*(q[k][j][i].rho +q[k+1][j][i].rho);
	  rhoU_cent=0.5*(q[k][j][i].rhoU+q[k+1][j][i].rhoU);
	  rhoV_cent=0.5*(q[k][j][i].rhoV+q[k+1][j][i].rhoV);
	  rhoW_cent=0.5*(q[k][j][i].rhoW+q[k+1][j][i].rhoW);
	  rhoE_cent=0.5*(q[k][j][i].rhoE+q[k+1][j][i].rhoE);
	  Ucont_cent=(csiH.x*rhoU_cent+csiH.y*rhoV_cent+csiH.z*rhoW_cent)/rho_cent;
	  
	  Fp_cent.rho  = rho_cent  *Ucont_cent;
	  Fp_cent.rhoU = rhoU_cent *Ucont_cent+0.5*(p[k][j][i]+p[k+1][j][i])*csiH.x;
	  Fp_cent.rhoV = rhoV_cent *Ucont_cent+0.5*(p[k][j][i]+p[k+1][j][i])*csiH.y;
	  Fp_cent.rhoW = rhoW_cent *Ucont_cent+0.5*(p[k][j][i]+p[k+1][j][i])*csiH.z;
	  Fp_cent.rhoE = Ucont_cent*(rhoE_cent+0.5*(p[k][j][i]+p[k+1][j][i]));
	}
	
	
	
	// calc F3 on k+1/2 nodes... k+1/2 saved in k
	// Lax-Freidrich flux: F(QR,QL)=1/2(E(QR)+E(QL))-1/2|S(A)|(QR-QL)
	if (fluxorder==1 || ((k==0 || k==mz-2) && user->bctype[4]!=7) ||
	    ((nvert[k][j][i]+nvert[k+1][j][1])>blankwall) ) { // if k or k+1 are wall nodes
	  QL=q[k][j][i];
	  QR=q[k+1][j][i];
	  pL=p[k][j][i];
	  pR=p[k+1][j][i];
	} else if (fluxorder==2) {
	  /* we know k and k+1 are not wall nodes but k-1 or k+1 might be wall nodes  */
	  if (nvert[k-1][j][i]>blankwall) {
	    /* put Qm equal to Qp to revert WENO to ENO */
	    qm=q[k+1][j][i];
	    pm=p[k+1][j][i];
	  } else {
	    qm=q[k-1][j][i];
	    pm=p[k-1][j][i];
	  }
	  if (nvert[k+2][j][i]>blankwall) {
	    /* put Qp equal to Qm to revert WENO to ENO */
	    qp=q[k][j][i];
	    pp=p[k][j][i];
	  } else {
	    qp=q[k+2][j][i];
	    pp=p[k+2][j][i];
	  }
	  QL.rho =MUSCL_QL(qm.rho , q[k][j][i].rho , q[k+1][j][i].rho);
	  QL.rhoU=MUSCL_QL(qm.rhoU, q[k][j][i].rhoU, q[k+1][j][i].rhoU);
	  QL.rhoV=MUSCL_QL(qm.rhoV, q[k][j][i].rhoV, q[k+1][j][i].rhoV);
	  QL.rhoW=MUSCL_QL(qm.rhoW, q[k][j][i].rhoW, q[k+1][j][i].rhoW);
	  QL.rhoE=MUSCL_QL(qm.rhoE, q[k][j][i].rhoE, q[k+1][j][i].rhoE);
	  pL     =MUSCL_QL(pm     , p[k][j][i],      p[k+1][j][i]);

	  QR.rho =MUSCL_QR(q[k][j][i].rho , q[k+1][j][i].rho , qp.rho);
	  QR.rhoU=MUSCL_QR(q[k][j][i].rhoU, q[k+1][j][i].rhoU, qp.rhoU);
	  QR.rhoV=MUSCL_QR(q[k][j][i].rhoV, q[k+1][j][i].rhoV, qp.rhoV);
	  QR.rhoW=MUSCL_QR(q[k][j][i].rhoW, q[k+1][j][i].rhoW, qp.rhoW);
	  QR.rhoE=MUSCL_QR(q[k][j][i].rhoE, q[k+1][j][i].rhoE, qp.rhoE);
	  pR     =MUSCL_QR(p[k][j][i]     , p[k+1][j][i],      pp);
	} else { // WENO
	  /* we know k and k+1 are not wall nodes but k-1 or k+1 might be wall nodes  */
	  if (nvert[k-1][j][i]>blankwall) {
	    /* put Qm equal to Qp to revert WENO to ENO */
	    //qm=q[k+1][j][i];
	    //pm=p[k+1][j][i];
	    QL.rho =ENO_QL(q[k][j][i].rho , q[k+1][j][i].rho , q[k+2][j][i].rho ,order);
	    QL.rhoU=ENO_QL(q[k][j][i].rhoU, q[k+1][j][i].rhoU, q[k+2][j][i].rhoU,order);
	    QL.rhoV=ENO_QL(q[k][j][i].rhoV, q[k+1][j][i].rhoV, q[k+2][j][i].rhoV,order);
	    QL.rhoW=ENO_QL(q[k][j][i].rhoW, q[k+1][j][i].rhoW, q[k+2][j][i].rhoW,order);
	    QL.rhoE=ENO_QL(q[k][j][i].rhoE, q[k+1][j][i].rhoE, q[k+2][j][i].rhoE,order);
	    pL     =ENO_QL(p[k][j][i]     , p[k+1][j][i],      p[k+2][j][i]     ,order);
	  } else {
	    qm=q[k-1][j][i];
	    pm=p[k-1][j][i];

	    QL.rho =WENO_QL(qm.rho , q[k][j][i].rho , q[k+1][j][i].rho);
	    QL.rhoU=WENO_QL(qm.rhoU, q[k][j][i].rhoU, q[k+1][j][i].rhoU);
	    QL.rhoV=WENO_QL(qm.rhoV, q[k][j][i].rhoV, q[k+1][j][i].rhoV);
	    QL.rhoW=WENO_QL(qm.rhoW, q[k][j][i].rhoW, q[k+1][j][i].rhoW);
	    QL.rhoE=WENO_QL(qm.rhoE, q[k][j][i].rhoE, q[k+1][j][i].rhoE);
	    pL     =WENO_QL(pm     , p[k][j][i],      p[k+1][j][i]);
	  }
	  if (nvert[k+2][j][i]>blankwall) {
	    /* put Qp equal to Qm to revert WENO to ENO */
	    /* qp=q[k][j][i]; */
	    /* pp=p[k][j][i]; */
	    qp=q[k+1][j][i];
	    pp=p[k+1][j][i];
	    QR.rho =ENO_QR(q[k-1][j][i].rho , q[k][j][i].rho , qp.rho ,order);
	    QR.rhoU=ENO_QR(q[k-1][j][i].rhoU, q[k][j][i].rhoU, qp.rhoU,order);
	    QR.rhoV=ENO_QR(q[k-1][j][i].rhoV, q[k][j][i].rhoV, qp.rhoV,order);
	    QR.rhoW=ENO_QR(q[k-1][j][i].rhoW, q[k][j][i].rhoW, qp.rhoW,order);
	    QR.rhoE=ENO_QR(q[k-1][j][i].rhoE, q[k][j][i].rhoE, qp.rhoE,order);
	    pR     =ENO_QR(p[k-1][j][i]     , p[k][j][i],      pp     ,order);

	  } else {
	    qp=q[k+2][j][i];
	    pp=p[k+2][j][i];
	    QR.rho =WENO_QR(q[k][j][i].rho , q[k+1][j][i].rho , qp.rho);
	    QR.rhoU=WENO_QR(q[k][j][i].rhoU, q[k+1][j][i].rhoU, qp.rhoU);
	    QR.rhoV=WENO_QR(q[k][j][i].rhoV, q[k+1][j][i].rhoV, qp.rhoV);
	    QR.rhoW=WENO_QR(q[k][j][i].rhoW, q[k+1][j][i].rhoW, qp.rhoW);
	    QR.rhoE=WENO_QR(q[k][j][i].rhoE, q[k+1][j][i].rhoE, qp.rhoE);
	    pR     =WENO_QR(p[k][j][i]     , p[k+1][j][i],      pp);
	  }
										
	} // WENO

	// csiH at half node k+1/2
	/* csiH.x=0.5*(zet[k][j][i].x+zet[k+1][j][i].x); */
	/* csiH.y=0.5*(zet[k][j][i].y+zet[k+1][j][i].y); */
	/* csiH.z=0.5*(zet[k][j][i].z+zet[k+1][j][i].z); */

	uconL=(csiH.x*(QL.rhoU/QL.rho)+
	       csiH.y*(QL.rhoV/QL.rho)+
	       csiH.z*(QL.rhoW/QL.rho));
	uconR=(csiH.x*(QR.rhoU/QR.rho)+
	       csiH.y*(QR.rhoV/QR.rho)+
	       csiH.z*(QR.rhoW/QR.rho));
	gii=csiH.x*csiH.x+csiH.y*csiH.y+csiH.z*csiH.z;
	//	giiR=csiR.x*csiR.x+csiR.y*csiR.y+csiR.z*csiR.z;
	CsoundL=gamma*pL/QL.rho;
	CsoundR=gamma*pR/QR.rho;
	saL=fabs(uconL)+sqrt(fabs(CsoundL)*gii);
	saR=fabs(uconR)+sqrt(fabs(CsoundR)*gii);
	SA=PetscMax(saL,saR);
	
	fp_WE.rho =0.5*( QL.rho *uconL + QR.rho *uconR) - 0.5*SA*(QR.rho -QL.rho);
	fp_WE.rhoU=0.5*((QL.rhoU*uconL + QR.rhoU*uconR) +
			(pL   * csiH.x + pR * csiH.x))  - 0.5*SA*(QR.rhoU-QL.rhoU);
	fp_WE.rhoV=0.5*((QL.rhoV*uconL + QR.rhoV*uconR) +
			(pL   * csiH.y + pR * csiH.y))  - 0.5*SA*(QR.rhoV-QL.rhoV);
	fp_WE.rhoW=0.5*((QL.rhoW*uconL + QR.rhoW*uconR) +
			(pL   * csiH.z + pR * csiH.z))  - 0.5*SA*(QR.rhoW-QL.rhoW);
	fp_WE.rhoE=0.5*((pL+QL.rhoE)*uconL +
			(pR+QR.rhoE)*uconR)             - 0.5*SA*(QR.rhoE-QL.rhoE);


	/* hybrid scheme  */
	fp3[k][j][i].rho =(1.-sigma)*fp_WE.rho +sigma*Fp_cent.rho;
	fp3[k][j][i].rhoU=(1.-sigma)*fp_WE.rhoU+sigma*Fp_cent.rhoU;
	fp3[k][j][i].rhoV=(1.-sigma)*fp_WE.rhoV+sigma*Fp_cent.rhoV;
	fp3[k][j][i].rhoW=(1.-sigma)*fp_WE.rhoW+sigma*Fp_cent.rhoW;
	fp3[k][j][i].rhoE=(1.-sigma)*fp_WE.rhoE+sigma*Fp_cent.rhoE;

      }
    
  
  DMDAVecRestoreArray(user->fda, user->lZet, &zet);
  DMDAVecRestoreArray(user->fda, user->lSigma, &lsigma);
  //  PetscBarrier(NULL);
  //  PetscPrintf(PETSC_COMM_SELF, "!!Central Done k rank %d!!\n", rank);

  //  PetscPrintf(PETSC_COMM_WORLD, "!!Central Done K!!\n");		
  /* DMDAVecRestoreArray(cda, Fp1, &fp1); */
  /* DMDAVecRestoreArray(cda, Fp2, &fp2); */
  /* DMDAVecRestoreArray(cda, Fp3, &fp3); */

  /* VecDuplicate(user->lQ, &lFp1); */
  /* VecDuplicate(user->lQ, &lFp2); */
  /* VecDuplicate(user->lQ, &lFp3); */

  /* DMGlobalToLocalBegin(cda,Fp1, INSERT_VALUES, lFp1); */
  /* DMGlobalToLocalEnd(cda,Fp1, INSERT_VALUES, lFp1); */
  /* DMGlobalToLocalBegin(cda,Fp2, INSERT_VALUES, lFp2); */
  /* DMGlobalToLocalEnd(cda,Fp2, INSERT_VALUES, lFp2); */
  /* DMGlobalToLocalBegin(cda,Fp3, INSERT_VALUES, lFp3); */
  /* DMGlobalToLocalEnd(cda,Fp3, INSERT_VALUES, lFp3); */

  /* VecDestroy(&Fp1); */
  /* VecDestroy(&Fp2); */
  /* VecDestroy(&Fp3); */

  /* DMDAVecGetArray(cda, lFp1, &fp1); */
  /* DMDAVecGetArray(cda, lFp2, &fp2); */
  /* DMDAVecGetArray(cda, lFp3, &fp3); */

  DMDAVecGetArray(cda, Conv,  &conv);
  DMDAVecGetArray(da, user->lAj, &aj);


  /* Calculate the convective terms under cartesian coordinates */
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i]< blank) {
				  
	  conv[k][j][i].rho = 
	    fp1[k][j][i].rho - fp1[k][j][i-1].rho +
	    fp2[k][j][i].rho - fp2[k][j-1][i].rho +
	    fp3[k][j][i].rho - fp3[k-1][j][i].rho;

	  if (twoD!=1)
	  conv[k][j][i].rhoU = 
	    fp1[k][j][i].rhoU - fp1[k][j][i-1].rhoU +
	    fp2[k][j][i].rhoU - fp2[k][j-1][i].rhoU +
	    fp3[k][j][i].rhoU - fp3[k-1][j][i].rhoU ;
	  else conv[k][j][i].rhoU = 0.;

	  conv[k][j][i].rhoV = 
	    fp1[k][j][i].rhoV - fp1[k][j][i-1].rhoV +
	    fp2[k][j][i].rhoV - fp2[k][j-1][i].rhoV +
	    fp3[k][j][i].rhoV - fp3[k-1][j][i].rhoV;

	  conv[k][j][i].rhoW = (
				fp1[k][j][i].rhoW - fp1[k][j][i-1].rhoW +
				fp2[k][j][i].rhoW - fp2[k][j-1][i].rhoW +
				fp3[k][j][i].rhoW - fp3[k-1][j][i].rhoW)- user->ratio/aj[k][j][i]*q[k][j][i].rho/user->dt;

	  conv[k][j][i].rhoE = 
	    fp1[k][j][i].rhoE - fp1[k][j][i-1].rhoE +
	    fp2[k][j][i].rhoE - fp2[k][j-1][i].rhoE +
	    fp3[k][j][i].rhoE - fp3[k-1][j][i].rhoE - user->ratio/aj[k][j][i]*q[k][j][i].rhoW/user->dt;
	} else {
	  conv[k][j][i].rho = 0.;
	  conv[k][j][i].rhoU = 0.;
	  conv[k][j][i].rhoV = 0.;
	  conv[k][j][i].rhoW = 0.;
	  conv[k][j][i].rhoE = 0.;
	}

      }
    }
  }
      
      
      //Amir
  //
  //  PetscPrintf(PETSC_COMM_WORLD, "!!CONV Done 1!!\n");	

  /* for (k=lzs; k<lze; k++) { */
  /*   for (j=ys; j<ye; j++) { */
  /*     for (i=xs; i<xe; i++) { */
  /* 	//	if (k==2 && j==15 ) */
  /* 	//	  PetscPrintf(PETSC_COMM_WORLD, "i  %d  rhs!!   %le   %le   %le\n",i,rhs[k][j][i].rhoU,rhs[k][j][i].rhoV,rhs[k][j][i].rhoW); */
  /* 	//	if (j<6 && i==4 ) */

  /* 	//	  PetscPrintf(PETSC_COMM_WORLD, "k %d   j  %d !!   V %le   fp2_W:  %le  ConvW    %le    convV  %le  \n",k,j, q[k][j][i].rhoV,fp2[k][j][i].rhoW,conv[k][j][i].rhoW,conv[k][j][i].rhoV); */


  /*     }	   */
  /*   }   */
  /* } */

  DMDAVecRestoreArray(cda, user->lQ, &q);
  DMDAVecRestoreArray(da, user->lP, &p);

  DMDAVecRestoreArray(cda, Conv,  &conv);
  //  DMDAVecRestoreArray(cda, Conv_cent, &conv_cent);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  PetscReal norm;
  VecNorm(Conv, NORM_INFINITY, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "!!norm of Conv %le\n", norm);

  DMDAVecRestoreArray(cda, lFp1, &fp1);
  DMDAVecRestoreArray(cda, lFp2, &fp2);
  DMDAVecRestoreArray(cda, lFp3, &fp3);
  
 

  VecDestroy(&lFp1);
  VecDestroy(&lFp2);
  VecDestroy(&lFp3);
  //  VecDestroy(&Conv_cent);

  return(0);
}


PetscErrorCode formDissipation(UserCtx *user, Vec Disp)
{

  DM		da = user->da, cda = user->cda;
  DMDALocalInfo	info=user->info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k;

  Vec		Fp1, Fp2, Fp3;
  Vec		lFp1, lFp2, lFp3;
  CompVars	***fp1, ***fp2, ***fp3;
  CompVars	***disp, ***q;
 
  PetscReal	ucon, SA,gii, epsilon4=0.0, epsilon2=0., gamma=1.4, Csound;
  PetscOptionsGetReal(NULL, NULL, "-dissipation4", &(epsilon4), NULL);
  PetscOptionsGetReal(NULL, NULL, "-dissipation2", &(epsilon2), NULL);

  PetscReal	***nvert, ***aj, ***p;
  Cmpnts	***csi, ***eta, ***zet;

  //  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;


  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  VecDuplicate(user->Q, &Fp1);
  VecDuplicate(user->Q, &Fp2);
  VecDuplicate(user->Q, &Fp3);

  VecSet(Fp1, 0.0);
  VecSet(Fp2, 0.0);
  VecSet(Fp3, 0.0);

  DMDAVecGetArray(cda, Disp,  &disp);
 
  DMDAVecGetArray(cda, Fp1, &fp1);
  DMDAVecGetArray(cda, Fp2, &fp2);
  DMDAVecGetArray(cda, Fp3, &fp3);

  DMDAVecGetArray(cda, user->lQ, &q);
  DMDAVecGetArray(da, user->lP, &p);

  DMDAVecGetArray(user->fda, user->lCsi, &csi);
  DMDAVecGetArray(user->fda, user->lEta, &eta);
  DMDAVecGetArray(user->fda, user->lZet, &zet);

 
  if (xe==mx) lxe=xe-2;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;

  // forms 4th order dissipation
  // D(i+1/2)=epsilon*Spectral(A1)( q(i+2)-3q(i+1)+3q(i)-q(i-1) )
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	// calc F1 on i+1/2 nodes... i+1/2 saved in i
	ucon=0.25*((csi[k][j][i].x+csi[k][j][i+1].x)*(q[k][j][i].rhoU/q[k][j][i].rho+q[k][j][i+1].rhoU/q[k][j][i+1].rho)+
		   (csi[k][j][i].y+csi[k][j][i+1].y)*(q[k][j][i].rhoV/q[k][j][i].rho+q[k][j][i+1].rhoV/q[k][j][i+1].rho)+
		   (csi[k][j][i].z+csi[k][j][i+1].z)*(q[k][j][i].rhoW/q[k][j][i].rho+q[k][j][i+1].rhoW/q[k][j][i+1].rho));
	gii=0.25*((csi[k][j][i].x+csi[k][j][i+1].x)*(csi[k][j][i].x+csi[k][j][i+1].x)+
		  (csi[k][j][i].x+csi[k][j][i+1].y)*(csi[k][j][i].y+csi[k][j][i+1].y)+
		  (csi[k][j][i].x+csi[k][j][i+1].z)*(csi[k][j][i].z+csi[k][j][i+1].z));
	Csound=gamma*(p[k][j][i]+p[k][j][i+1])/(q[k][j][i].rho+q[k][j][i+1].rho);
	SA=fabs(ucon)+sqrt(fabs(Csound)*gii);
	
	fp1[k][j][i].rho= -epsilon2*0.5*SA*(q[k][j][i+1].rho-q[k][j][i].rho)+
	  epsilon4*SA*(q[k][j][i+2].rho-3.*q[k][j][i+1].rho+3.*q[k][j][i].rho-q[k][j][i-1].rho);
	fp1[k][j][i].rhoU= -epsilon2*0.5*SA*(q[k][j][i+1].rhoU-q[k][j][i].rhoU)+
	  epsilon4*SA*(q[k][j][i+2].rhoU-3.*q[k][j][i+1].rhoU+3.*q[k][j][i].rhoU-q[k][j][i-1].rhoU);
	fp1[k][j][i].rhoV= -epsilon2*0.5*SA*(q[k][j][i+1].rhoV-q[k][j][i].rhoV)+
	  epsilon4*SA*(q[k][j][i+2].rhoV-3.*q[k][j][i+1].rhoV+3.*q[k][j][i].rhoV-q[k][j][i-1].rhoV);
	fp1[k][j][i].rhoW= -epsilon2*0.5*SA*(q[k][j][i+1].rhoW-q[k][j][i].rhoW)+
	  epsilon4*SA*(q[k][j][i+2].rhoW-3.*q[k][j][i+1].rhoW+3.*q[k][j][i].rhoW-q[k][j][i-1].rhoW);
	fp1[k][j][i].rhoE= -epsilon2*0.5*SA*(q[k][j][i+1].rhoE-q[k][j][i].rhoE)+
	  epsilon4*SA*(q[k][j][i+2].rhoE-3.*q[k][j][i+1].rhoE+3.*q[k][j][i].rhoE-q[k][j][i-1].rhoE);
      }
    }
  }

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      if (xs==0) {
      i=xs;
      fp1[k][j][i].rho=fp1[k][j][i+1].rho;
      fp1[k][j][i].rhoU=fp1[k][j][i+1].rhoU;
      fp1[k][j][i].rhoV=fp1[k][j][i+1].rhoV;
      fp1[k][j][i].rhoW=fp1[k][j][i+1].rhoW;
      fp1[k][j][i].rhoE=fp1[k][j][i+1].rhoE;
      }
      if (mx==xe) {
      i=xe-2;
      fp1[k][j][i].rho=fp1[k][j][i-1].rho;
      fp1[k][j][i].rhoU=fp1[k][j][i-1].rhoU;
      fp1[k][j][i].rhoV=fp1[k][j][i-1].rhoV;
      fp1[k][j][i].rhoW=fp1[k][j][i-1].rhoW;
      fp1[k][j][i].rhoE=fp1[k][j][i-1].rhoE;
      }
    }
  }
  
/*   if (xe==mx) lxe=xe-1; */
/*   if (ye==my) lye=ye-2; */
/*   if (ze==mz) lze=ze-1; */
/*   for (k=lzs; k<lze; k++){ */
/*     for (j=lys; j<lye; j++){ */
/*       for (i=lxs; i<lxe; i++){ */
/* 	// calc F2 on j+1/2 nodes... j+1/2 saved in j */
/* 	ucon=0.25*((eta[k][j][i].x+eta[k][j+1][i].x)*(q[k][j][i].rhoU/q[k][j][i].rho+q[k][j+1][i].rhoU/q[k][j+1][i].rho)+ */
/* 		   (eta[k][j][i].y+eta[k][j+1][i].y)*(q[k][j][i].rhoV/q[k][j][i].rho+q[k][j+1][i].rhoV/q[k][j+1][i].rho)+ */
/* 		   (eta[k][j][i].z+eta[k][j+1][i].z)*(q[k][j][i].rhoW/q[k][j][i].rho+q[k][j+1][i].rhoW/q[k][j+1][i].rho)); */

/* 	gii=0.25*((eta[k][j][i].x+eta[k][j+1][i].x)*(eta[k][j][i].x+eta[k][j+1][i].x)+ */
/* 		  (eta[k][j][i].y+eta[k][j+1][i].y)*(eta[k][j][i].y+eta[k][j+1][i].y)+ */
/* 		  (eta[k][j][i].z+eta[k][j+1][i].z)*(eta[k][j][i].z+eta[k][j+1][i].z)); */
/* 	SA=ucon+sqrt(ucon*ucon+gii); */
	
/* 	fp2[k][j][i].rho =epsilon*SA*(q[k][j+2][i].rho -3.*q[k][j+1][i].rho +3.*q[k][j][i].rho -q[k][j-1][i].rho); */
/* 	fp2[k][j][i].rhoU=epsilon*SA*(q[k][j+2][i].rhoU-3.*q[k][j+1][i].rhoU+3.*q[k][j][i].rhoU-q[k][j-1][i].rhoU); */
/* 	fp2[k][j][i].rhoV=epsilon*SA*(q[k][j+2][i].rhoV-3.*q[k][j+1][i].rhoV+3.*q[k][j][i].rhoV-q[k][j-1][i].rhoV); */
/* 	fp2[k][j][i].rhoW=epsilon*SA*(q[k][j+2][i].rhoW-3.*q[k][j+1][i].rhoW+3.*q[k][j][i].rhoW-q[k][j-1][i].rhoW); */
/* 	fp2[k][j][i].rhoE=epsilon*SA*(q[k][j+2][i].rhoE-3.*q[k][j+1][i].rhoE+3.*q[k][j][i].rhoE-q[k][j-1][i].rhoE); */
/*       } */
/*     } */
/*   } */

/*   for (k=lzs; k<lze; k++){ */
/*     for (i=lxs; i<lxe; i++){ */
/*       j=ys; */
/*       fp2[k][j][i].rho=fp2[k][j+1][i].rho; */
/*       fp2[k][j][i].rhoU=fp2[k][j+1][i].rhoU; */
/*       fp2[k][j][i].rhoV=fp2[k][j+1][i].rhoV; */
/*       fp2[k][j][i].rhoW=fp2[k][j+1][i].rhoW; */
/*       fp2[k][j][i].rhoE=fp2[k][j+1][i].rhoE; */
      
/*       j=ye-2; */
/*       fp2[k][j][i].rho=fp2[k][j-1][i].rho; */
/*       fp2[k][j][i].rhoU=fp2[k][j-1][i].rhoU; */
/*       fp2[k][j][i].rhoV=fp2[k][j-1][i].rhoV; */
/*       fp2[k][j][i].rhoW=fp2[k][j-1][i].rhoW; */
/*       fp2[k][j][i].rhoE=fp2[k][j-1][i].rhoE; */
/*     } */
/*   } */

/*   if (xe==mx) lxe=xe-1; */
/*   if (ye==my) lye=ye-1; */
/*   if (ze==mz) lze=ze-2; */
/*   for (k=lzs; k<lze; k++){ */
/*     for (j=lys; j<lye; j++){ */
/*       for (i=lxs; i<lxe; i++){ */
/* 	// calc F3 on k+1/2 nodes... k+1/2 saved in k */
/* 	ucon=0.25*((zet[k][j][i].x+zet[k+1][j][i].x)*(q[k][j][i].rhoU/q[k][j][i].rho+q[k+1][j][i].rhoU/q[k+1][j][i].rho)+ */
/* 		   (zet[k][j][i].y+zet[k+1][j][i].y)*(q[k][j][i].rhoV/q[k][j][i].rho+q[k+1][j][i].rhoV/q[k+1][j][i].rho)+ */
/* 		   (zet[k][j][i].z+zet[k+1][j][i].z)*(q[k][j][i].rhoW/q[k][j][i].rho+q[k+1][j][i].rhoW/q[k+1][j][i].rho)); */
/* 	gii=0.25*((zet[k][j][i].x+zet[k+1][j][i].x)*(zet[k][j][i].x+zet[k+1][j][i].x)+ */
/* 		  (zet[k][j][i].y+zet[k+1][j][i].y)*(zet[k][j][i].y+zet[k+1][j][i].y)+ */
/* 		  (zet[k][j][i].z+zet[k+1][j][i].z)*(zet[k][j][i].z+zet[k+1][j][i].z)); */
/* 	SA=ucon+sqrt(ucon*ucon+gii); */
	
/* 	fp3[k][j][i].rho =epsilon*SA*(q[k+2][j][i].rho -3.*q[k+1][j][i].rho +3.*q[k][j][i].rho -q[k-1][j][i].rho); */
/* 	fp3[k][j][i].rhoU=epsilon*SA*(q[k+2][j][i].rhoU-3.*q[k+1][j][i].rhoU+3.*q[k][j][i].rhoU-q[k-1][j][i].rhoU); */
/* 	fp3[k][j][i].rhoV=epsilon*SA*(q[k+2][j][i].rhoV-3.*q[k+1][j][i].rhoV+3.*q[k][j][i].rhoV-q[k-1][j][i].rhoV); */
/* 	fp3[k][j][i].rhoW=epsilon*SA*(q[k+2][j][i].rhoW-3.*q[k+1][j][i].rhoW+3.*q[k][j][i].rhoW-q[k-1][j][i].rhoW); */
/* 	fp3[k][j][i].rhoE=epsilon*SA*(q[k+2][j][i].rhoE-3.*q[k+1][j][i].rhoE+3.*q[k][j][i].rhoE-q[k-1][j][i].rhoE); */
	
/*       } */
/*     } */
/*   } */

/*   for (j=lys; j<lye; j++){ */
/*     for (i=lxs; i<lxe; i++){ */
/*       k=zs; */
/*       fp3[k][j][i].rho=fp3[k+1][j][i].rho; */
/*       fp3[k][j][i].rhoU=fp3[k+1][j][i].rhoU; */
/*       fp3[k][j][i].rhoV=fp3[k+1][j][i].rhoV; */
/*       fp3[k][j][i].rhoW=fp3[k+1][j][i].rhoW; */
/*       fp3[k][j][i].rhoE=fp3[k+1][j][i].rhoE; */

/*       k=ze-2; */
/*       fp3[k][j][i].rho=fp3[k-1][j][i].rho; */
/*       fp3[k][j][i].rhoU=fp3[k-1][j][i].rhoU; */
/*       fp3[k][j][i].rhoV=fp3[k-1][j][i].rhoV; */
/*       fp3[k][j][i].rhoW=fp3[k-1][j][i].rhoW; */
/*       fp3[k][j][i].rhoE=fp3[k-1][j][i].rhoE; */
/*     } */
/*   } */

  DMDAVecRestoreArray(cda, user->lQ, &q);
  DMDAVecRestoreArray(da, user->lP, &p);

  DMDAVecRestoreArray(user->fda, user->lCsi, &csi);
  DMDAVecRestoreArray(user->fda, user->lEta, &eta);
  DMDAVecRestoreArray(user->fda, user->lZet, &zet);

 

  DMDAVecRestoreArray(cda, Fp1, &fp1);
  DMDAVecRestoreArray(cda, Fp2, &fp2);
  DMDAVecRestoreArray(cda, Fp3, &fp3);

  VecDuplicate(user->lQ, &lFp1);
  VecDuplicate(user->lQ, &lFp2);
  VecDuplicate(user->lQ, &lFp3);

  DMGlobalToLocalBegin(cda,Fp1, INSERT_VALUES, lFp1);
  DMGlobalToLocalEnd(cda,Fp1, INSERT_VALUES, lFp1);
  DMGlobalToLocalBegin(cda,Fp2, INSERT_VALUES, lFp2);
  DMGlobalToLocalEnd(cda,Fp2, INSERT_VALUES, lFp2);
  DMGlobalToLocalBegin(cda,Fp3, INSERT_VALUES, lFp3);
  DMGlobalToLocalEnd(cda,Fp2, INSERT_VALUES, lFp2);

  DMDAVecGetArray(cda, lFp1, &fp1);
  DMDAVecGetArray(cda, lFp2, &fp2);
  DMDAVecGetArray(cda, lFp3, &fp3);
  DMDAVecGetArray(da, user->lAj, &aj);

  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;
  /* Calculate the dissipation terms under cartesian coordinates */
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	disp[k][j][i].rho = (
	  fp1[k][j][i].rho - fp1[k][j][i-1].rho +
	  fp2[k][j][i].rho - fp2[k][j-1][i].rho +
	  fp3[k][j][i].rho - fp3[k-1][j][i].rho);

	disp[k][j][i].rhoU =  (
	  fp1[k][j][i].rhoU - fp1[k][j][i-1].rhoU +
	  fp2[k][j][i].rhoU - fp2[k][j-1][i].rhoU +
	  fp3[k][j][i].rhoU - fp3[k-1][j][i].rhoU);

	disp[k][j][i].rhoV =  (
	  fp1[k][j][i].rhoV - fp1[k][j][i-1].rhoV +
	  fp2[k][j][i].rhoV - fp2[k][j-1][i].rhoV +
	  fp3[k][j][i].rhoV - fp3[k-1][j][i].rhoV);

	disp[k][j][i].rhoW =  (
	  fp1[k][j][i].rhoW - fp1[k][j][i-1].rhoW +
	  fp2[k][j][i].rhoW - fp2[k][j-1][i].rhoW +
	  fp3[k][j][i].rhoW - fp3[k-1][j][i].rhoW);

	disp[k][j][i].rhoE =  (
	  fp1[k][j][i].rhoE - fp1[k][j][i-1].rhoE +
	  fp2[k][j][i].rhoE - fp2[k][j-1][i].rhoE +
	  fp3[k][j][i].rhoE - fp3[k-1][j][i].rhoE);
      }
    }
  } 

  DMDAVecRestoreArray(cda, Disp,  &disp);
  DMDAVecRestoreArray(da, user->lAj, &aj);

  DMDAVecRestoreArray(cda, lFp1, &fp1);
  DMDAVecRestoreArray(cda, lFp2, &fp2);
  DMDAVecRestoreArray(cda, lFp3, &fp3);

  VecDestroy(&Fp1);
  VecDestroy(&Fp2);
  VecDestroy(&Fp3);

  VecDestroy(&lFp1);
  VecDestroy(&lFp2);
  VecDestroy(&lFp3);

  return(0);
}

PetscErrorCode MultAj(UserCtx *user, Vec RHS, Vec Aj)
{
  DM		da = user->da, cda=user->cda;
 
  CompVars	***rhs, ***lrhs;
  PetscReal	***nvert, ***aj, blank=0.1;
  
  DMDALocalInfo	info=user->info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;


  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;

  
  PetscInt	i, j, k;

  DMDAVecGetArray(da, Aj, &aj);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(cda, RHS,  &rhs);
  
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i]< blank) {
	  rhs[k][j][i].rho  = rhs[k][j][i].rho *aj[k][j][i];
	  rhs[k][j][i].rhoU = rhs[k][j][i].rhoU*aj[k][j][i];
	  rhs[k][j][i].rhoV = rhs[k][j][i].rhoV*aj[k][j][i];
	  rhs[k][j][i].rhoW = rhs[k][j][i].rhoW*aj[k][j][i];
	  rhs[k][j][i].rhoE = rhs[k][j][i].rhoE*aj[k][j][i];	  
	} else {
	  rhs[k][j][i].rho  = 0.;
	  rhs[k][j][i].rhoU = 0.;
	  rhs[k][j][i].rhoV = 0.;
	  rhs[k][j][i].rhoW = 0.;
	  rhs[k][j][i].rhoE = 0.;
	}	  
      }
    }
  }

if (TwoD==3){
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	  rhs[k][j][i].rhoW = 0.0;
	}	  
      }
    }
  }


/*
 for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
       for (i=lxs; i<lxe; i++) {
	if (k==2 && j==15 )
	  PetscPrintf(PETSC_COMM_WORLD, "i  %d  rhs!!   %le   %le   %le\n",i,rhs[k][j][i].rhoU,rhs[k][j][i].rhoV,rhs[k][j][i].rhoW);
 	if (k==2 && i==15 )
	
	  PetscPrintf(PETSC_COMM_WORLD, "j  %d  rhs!!   %le   %le   %le\n",j,rhs[k][j][i].rhoU,rhs[k][j][i].rhoV,rhs[k][j][i].rhoW);


	}	  
      }  
  }
*/
  DMDAVecRestoreArray(da, Aj, &aj);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(cda, RHS,  &rhs);
/*
 if (bctype[0]==7){ 
 
  DMDAVecGetArray(cda, RHS,  &rhs);
  DMDAVecGetArray(cda, lRHS,  &lrhs);
 if (xs==0){
  i=0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
	rhs[k][j][i] = lrhs[k][j][i-2];
	}	  
      }
    }

  DMDAVecRestoreArray(cda, RHS,  &rhs);
  DMDAVecRestoreArray(cda, lRHS,  &lrhs);
  DMGlobalToLocalBegin(cda, user->RHS, INSERT_VALUES, user->lRHS);
  DMGlobalToLocalEnd(cda, user->RHS, INSERT_VALUES, user->lRHS);
 

  DMDAVecGetArray(cda, RHS,  &rhs);
  DMDAVecGetArray(cda, lRHS,  &lrhs);
  
  if (xe==xm){
  i=xm-1;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
	rhs[k][j][i] = lrhs[k][j][i+2];
	}	  
      }
    }

  DMDAVecRestoreArray(cda, RHS,  &rhs);
  DMDAVecRestoreArray(cda, lRHS,  &lrhs);
  DMGlobalToLocalBegin(da, user->RHS, INSERT_VALUES, user->lRHS);
  DMGlobalToLocalEnd(da, user->RHS, INSERT_VALUES, user->lRHS);
 }

*/
  
  return(0);
}

PetscErrorCode formRHS(UserCtx *user, Vec Rhs, PetscInt istage)
{
  Vec		Conv, Visc, Disp;
  VecDuplicate(user->Q, &Conv);
  VecDuplicate(user->Q, &Visc);
  //  VecDuplicate(user->Q, &Disp);
  VecSet(Conv, 0.0);
  VecSet(Visc, 0.0);
  //  VecSet(Disp, 0.0);

  if (les) 
    formConvectionLES(user, Conv, istage);
  else 
    formConvection(user, Conv);
  
  if (!invicid) {
    CalcUCat(user);
 //   PetscPrintf(PETSC_COMM_WORLD, "Uca
    formViscous(user, Visc);
  }
//   formDissipation(user, Disp);
  
  VecWAXPY(Rhs,-1.0, Conv, Visc); //Rhs=-conv+visc
  //  VecAXPY(Rhs, -1.0, Disp);

  //mulyiply Rhs by aj (jaocibian)
  MultAj(user, Rhs, user->lAj);

  VecDestroy(&Conv);
  VecDestroy(&Visc);
  // VecDestroy(&Disp);

  return(0);
}


PetscErrorCode Spectral(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;

  PetscReal cfl = 0.5;
  PetscReal vnn = 0.5;

  PetscOptionsGetReal(NULL, NULL, "-cfl", &cfl, NULL);
  PetscInt i, j, k;
  PetscReal abu, abv, abw, ren = user->ren;
  Cmpnts ***ucont, ***csi, ***eta, ***zet;
  PetscReal ***dt, ***aj;
  PetscReal g1, g2, g3, eigen1, eigen2, eigen3, temp, dtii, dtij, dtik;
  PetscReal dtvi, dtvj, dtvk;
  DMDAVecGetArray(fda, user->lUcont, &ucont);
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lAj, &aj);
  DMDAVecGetArray(da, user->Dt, &dt);
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	g1 = (csi[k][j][i].x * csi[k][j][i].x +
	      csi[k][j][i].y * csi[k][j][i].y +
	      csi[k][j][i].z * csi[k][j][i].z) * aj[k][j][i] * aj[k][j][i];

	g2 = (eta[k][j][i].x * eta[k][j][i].x +
	      eta[k][j][i].y * eta[k][j][i].y +
	      eta[k][j][i].z * eta[k][j][i].z) * aj[k][j][i] * aj[k][j][i];

	g3 = (zet[k][j][i].x * zet[k][j][i].x +
	      zet[k][j][i].y * zet[k][j][i].y +
	      zet[k][j][i].z * zet[k][j][i].z) * aj[k][j][i] * aj[k][j][i];

	abu = fabs(ucont[k][j][i].x + ucont[k][j][i-1].x) * 0.5 * aj[k][j][i];
	abv = fabs(ucont[k][j][i].y + ucont[k][j][i-1].y) * 0.5 * aj[k][j][i];
	abw = fabs(ucont[k][j][i].z + ucont[k][j][i-1].z) * 0.5 * aj[k][j][i];

	if (g1<1.e-10) g1 = 1.;
	if (g2<1.e-10) g2 = 1.;
	if (g3<1.e-10) g3 = 1.;
	eigen1 = abu + sqrt(abu*abu + g1);
	eigen2 = abv + sqrt(abv*abv + g2);
	eigen3 = abw + sqrt(abw*abw + g3);

	dtii = cfl / eigen1;
	dtij = cfl / eigen2;
	dtik = cfl / eigen3;

	temp = vnn * ren;
	dtvi = temp / g1;
	dtvj = temp / g2;
	dtvk = temp / g3;

	temp = PetscMin(dtii, dtij);
	temp = PetscMin(temp, dtik);
	temp = PetscMin(temp, dtvi);
	temp = PetscMin(temp, dtvj);
	dt[k][j][i] = PetscMin(temp, dtvk);
      }
    }
  }
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(da, user->Dt, &dt);
  return 0;
}

/* ==================================================================================             */
/*   Equation of state! */
/* ==================================================================================             */

PetscErrorCode EqOfState(UserCtx *user)
{

  DM		da = user->da, cda = user->cda;
  DMDALocalInfo	info=user->info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k;
  
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;

  PetscReal	***p,***lp, gamma=user->gamma;
  CompVars	***q;

  DMDAVecGetArray(cda, user->lQ, &q);
  DMDAVecGetArray(da, user->P, &p);
  
  for (k=zs; k<ze; k++){
    for (j=ys; j<ye; j++){
      for (i=xs; i<xe; i++){
	 // ideal gas law: rhoE=1/2 rho|u|^2+p/(gamma-1)
	p[k][j][i]=(gamma-1)*(q[k][j][i].rhoE-0.5*(q[k][j][i].rhoU*q[k][j][i].rhoU+
						   q[k][j][i].rhoV*q[k][j][i].rhoV+
						   q[k][j][i].rhoW*q[k][j][i].rhoW)/q[k][j][i].rho);

	if (p[k][j][i]<0) PetscPrintf(PETSC_COMM_SELF, "!!! p<0 at i,j,k %d %d %d\n", i,j,k);
      }
    }
  }

  DMDAVecRestoreArray(cda, user->lQ, &q);
  DMDAVecRestoreArray(da, user->P, &p); 
  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
 
  if (user->bctype[0]==7 || user->bctype[1]==7 ||
      user->bctype[2]==7 || user->bctype[3]==7 || 
      user->bctype[4]==7 || user->bctype[5]==7){
    
    
    
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->P, &p);
    
    //Loop Executed for i_periodic
    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xs==0){
	i=xs;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    if(k>0 && k<user->KM && j>0 && j<user->JM){
	      p[k][j][i]=lp[k][j][i-2];
	    }
	  }
	}
      }
      
      if (xe==mx){
	i=mx-1;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    if(k>0 && k<user->KM && j>0 && j<user->JM){
	      p[k][j][i]=lp[k][j][i+2];
	      
	    }
	  }
	}
      }
    }
 
     /* really needed below?  
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
  
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->lP, &lp);
     */
	// 
	if (user->bctype[2]==7||user->bctype[3]==7){		// Advait: Brought back Amir's implementation of conditional
    	if (ys==0){
      		j=ys;
      		for (k=zs; k<ze; k++) {
				for (i=xs; i<xe; i++) {
	  				if(k>0 && k<user->KM){
	    				p[k][j][i]=lp[k][j-2][i];           
	  				}
				}
      		}
    	}
    
    	if (ye==my){
      		j=my-1;
      		for (k=zs; k<ze; k++) {
				for (i=xs; i<xe; i++) {
	  				if(k>0 && k<user->KM ){
	    				p[k][j][i]=lp[k][j+2][i];
	    			}
				}
      		}
    	}
	}
     /* Are these below needed? 
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->lP, &lp);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->P, &p);
     */

    if (user->bctype[4]==7 || user->bctype[5]==7){
      if (zs==0){
	k=zs;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && i<user->IM && j>0 && j<user->JM){
	      p[k][j][i]=lp[k-2][j][i];
	      
	    }
	  }
	}
      }
    
      if (ze==mz){
	k=mz-1;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && i<user->IM && j>0 && j<user->JM){
	      p[k][j][i]=lp[k+2][j][i];
	    }
	  }
	}
      }
    }
  
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);	
    } // if periodic 7
  
  return(0);

}

PetscReal EqOfStateP(CompVars Qp)
{
  PetscReal p, gamma=1.4;
  p=(gamma-1)*(Qp.rhoE-0.5*(Qp.rhoU*Qp.rhoU+
			    Qp.rhoV*Qp.rhoV+
			    Qp.rhoW*Qp.rhoW)/Qp.rho);
  return(p);

}


