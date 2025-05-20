#include "variables.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//extern PetscInt ti,tistart,rstart_flg,fish;
extern PetscInt tistart;
extern PetscBool rstart_flg;
double eps=0.000001;	

void Compute_dp_dxyz (	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,double dpdc, double dpde, double dpdz,double *dp_dx, double *dp_dy, double *dp_dz )
{
  *dp_dx = (dpdc * csi0 + dpde * eta0 + dpdz * zet0) * ajc;
  *dp_dy = (dpdc * csi1 + dpde * eta1 + dpdz * zet1) * ajc;
  *dp_dz = (dpdc * csi2 + dpde * eta2 + dpdz * zet2) * ajc;
 };

void Compute_dp_center (int i, int j, int k,  int mx, int my, int mz, PetscReal ***lt, PetscReal ***nvert, 
			double *dpdc, double *dpde, double *dpdz)

//issues with indices ????

{
       *dpdc = ( lt[k][j][i+1] - lt[k][j][i-1] ) * 0.5;
       *dpde = ( lt[k][j+1][i] - lt[k][j-1][i] ) * 0.5;
       *dpdz = ( lt[k+1][j][i] - lt[k-1][j][i] ) * 0.5;
 
}

 	

void Compute_du_dxyz (	double csi0, double csi1, double csi2, double eta0, double eta1, double eta2, double zet0, double zet1, double zet2, double ajc,
			double dudc, double dvdc, double dwdc, double dude, double dvde, double dwde, double dudz, double dvdz, double dwdz,
			double *du_dx, double *dv_dx, double *dw_dx, double *du_dy, double *dv_dy, double *dw_dy, double *du_dz, double *dv_dz, double *dw_dz )
{
  *du_dx = (dudc * csi0 + dude * eta0 + dudz * zet0) * ajc;
  *du_dy = (dudc * csi1 + dude * eta1 + dudz * zet1) * ajc;
  *du_dz = (dudc * csi2 + dude * eta2 + dudz * zet2) * ajc;
  *dv_dx = (dvdc * csi0 + dvde * eta0 + dvdz * zet0) * ajc;
  *dv_dy = (dvdc * csi1 + dvde * eta1 + dvdz * zet1) * ajc;
  *dv_dz = (dvdc * csi2 + dvde * eta2 + dvdz * zet2) * ajc;
  *dw_dx = (dwdc * csi0 + dwde * eta0 + dwdz * zet0) * ajc;
  *dw_dy = (dwdc * csi1 + dwde * eta1 + dwdz * zet1) * ajc;	
  *dw_dz = (dwdc * csi2 + dwde * eta2 + dwdz * zet2) * ajc;
};

void Compute_du_center (int i, int j, int k,  int mx, int my, int mz, Cmpnts ***ucat, PetscReal ***nvert, 
			double *dudc, double *dvdc, double *dwdc, 
			double *dude, double *dvde, double *dwde,
			double *dudz, double *dvdz, double *dwdz)

//issues with indices ????

{
  *dudc = ( ucat[k][j][i+1].x - ucat[k][j][i-1].x ) * 0.5;
  *dvdc = ( ucat[k][j][i+1].y - ucat[k][j][i-1].y ) * 0.5;
  *dwdc = ( ucat[k][j][i+1].z - ucat[k][j][i-1].z ) * 0.5;
  *dudz = ( ucat[k+1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
       *dvdz = ( ucat[k+1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
       *dwdz = ( ucat[k+1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
       *dude = ( ucat[k][j+1][i].x - ucat[k][j-1][i].x ) * 0.5;
       *dvde = ( ucat[k][j+1][i].y - ucat[k][j-1][i].y ) * 0.5;
       *dwde = ( ucat[k][j+1][i].z - ucat[k][j-1][i].z ) * 0.5;
 
 
}
double integrate_t(double val[3][3][3])
{ 
  int k,j,i;
  double sum=0.0;
  for (k=0;k<3;k++){
    for (j=0;j<3;j++){
      for (i=0;i<3;i++){
	
	sum+=val[k][j][i]/27.0;
      }
    }
  }
  return(sum);
}
double integrate_testfilter(double val[3][3][3], double w[3][3][3])
{
  /*	double v1, v2, v3, v4, v5, v6, v7, v8;
	double w1, w2, w3, w4, w5, w6, w7, w8;
	if (channelz & !fish) testfilter_ik=1;	
	if( testfilter_ik ) {
		return integrate_testfilter_ik(val, w);
        	}
       

	v1 = ( val[0][0][0]*w[0][0][0] + val[1][0][0]*w[1][0][0] + val[0][0][1]*w[0][0][1] + val[1][0][1]*w[1][0][1] + val[0][1][0]*w[0][1][0] + val[1][1][0]*w[1][1][0] + val[0][1][1]*w[0][1][1] + val[1][1][1]*w[1][1][1] );
	v2 = ( val[1][0][0]*w[1][0][0] + val[2][0][0]*w[2][0][0] + val[1][0][1]*w[1][0][1] + val[2][0][1]*w[2][0][1] + val[1][1][0]*w[1][1][0] + val[2][1][0]*w[2][1][0] + val[1][1][1]*w[1][1][1] + val[2][1][1]*w[2][1][1] );
	v3 = ( val[0][0][1]*w[0][0][1] + val[1][0][1]*w[1][0][1] + val[0][0][2]*w[0][0][2] + val[1][0][2]*w[1][0][2] + val[0][1][1]*w[0][1][1] + val[1][1][1]*w[1][1][1] + val[0][1][2]*w[0][1][2] + val[1][1][2]*w[1][1][2] );
	v4 = ( val[1][0][1]*w[1][0][1] + val[2][0][1]*w[2][0][1] + val[1][0][2]*w[1][0][2] + val[2][0][2]*w[2][0][2] + val[1][1][1]*w[1][1][1] + val[2][1][1]*w[2][1][1] + val[1][1][2]*w[1][1][2] + val[2][1][2]*w[2][1][2] );
	
	//top : j=1, 2
	v5 = ( val[0][1][0]*w[0][1][0] + val[1][1][0]*w[1][1][0] + val[0][1][1]*w[0][1][1] + val[1][1][1]*w[1][1][1] + val[0][2][0]*w[0][2][0] + val[1][2][0]*w[1][2][0] + val[0][2][1]*w[0][2][1] + val[1][2][1]*w[1][2][1] );
	v6 = ( val[1][1][0]*w[1][1][0] + val[2][1][0]*w[2][1][0] + val[1][1][1]*w[1][1][1] + val[2][1][1]*w[2][1][1] + val[1][2][0]*w[1][2][0] + val[2][2][0]*w[2][2][0] + val[1][2][1]*w[1][2][1] + val[2][2][1]*w[2][2][1] );
	v7 = ( val[0][1][1]*w[0][1][1] + val[1][1][1]*w[1][1][1] + val[0][1][2]*w[0][1][2] + val[1][1][2]*w[1][1][2] + val[0][2][1]*w[0][2][1] + val[1][2][1]*w[1][2][1] + val[0][2][2]*w[0][2][2] + val[1][2][2]*w[1][2][2] );
	v8 = ( val[1][1][1]*w[1][1][1] + val[2][1][1]*w[2][1][1] + val[1][1][2]*w[1][1][2] + val[2][1][2]*w[2][1][2] + val[1][2][1]*w[1][2][1] + val[2][2][1]*w[2][2][1] + val[1][2][2]*w[1][2][2] + val[2][2][2]*w[2][2][2] );
	
	
	// bottom
	w1 = ( w[0][0][0] + w[1][0][0] + w[0][0][1] + w[1][0][1] + w[0][1][0] + w[1][1][0] + w[0][1][1] + w[1][1][1] );
	w2 = ( w[1][0][0] + w[2][0][0] + w[1][0][1] + w[2][0][1] + w[1][1][0] + w[2][1][0] + w[1][1][1] + w[2][1][1] );
	w3 = ( w[0][0][1] + w[1][0][1] + w[0][0][2] + w[1][0][2] + w[0][1][1] + w[1][1][1] + w[0][1][2] + w[1][1][2] );
	w4 = ( w[1][0][1] + w[2][0][1] + w[1][0][2] + w[2][0][2] + w[1][1][1] + w[2][1][1] + w[1][1][2] + w[2][1][2] );
	
	//top : j=1, 2
	w5 = ( w[0][1][0] + w[1][1][0] + w[0][1][1] + w[1][1][1] + w[0][2][0] + w[1][2][0] + w[0][2][1] + w[1][2][1] );
	w6 = ( w[1][1][0] + w[2][1][0] + w[1][1][1] + w[2][1][1] + w[1][2][0] + w[2][2][0] + w[1][2][1] + w[2][2][1] );
	w7 = ( w[0][1][1] + w[1][1][1] + w[0][1][2] + w[1][1][2] + w[0][2][1] + w[1][2][1] + w[0][2][2] + w[1][2][2] );
	w8 = ( w[1][1][1] + w[2][1][1] + w[1][1][2] + w[2][1][2] + w[1][2][1] + w[2][2][1] + w[1][2][2] + w[2][2][2] );
	return (v1+v2+v3+v4+v5+v6+v7+v8)/(w1+w2+w3+w4+w5+w6+w7+w8);
*/
  int k,j,i;
  double vsum=0.0,wsum=0.0;
	for (k=0;k<3;k++){
	  for (j=0;j<3;j++){
	    for (i=0;i<3;i++){
	      
	      vsum+=w[k][j][i]*val[k][j][i];
	      wsum+=w[k][j][i];
	    }
	  }
	}
	return (vsum/wsum);
};







void Compute_Smagorinsky_Constant_1(UserCtx *user)
{
	PetscPrintf(PETSC_COMM_WORLD,"!!!Compute Smagorinsky!!!\n");
  //if((user->ti<2 && user->tistart==0 && !user->rstart_flg) || les==10) { //New Comment // && tistart==0 && !rstart_flg){
    if((user->ti<2 && tistart==0 && !rstart_flg)|| les==10){
	VecSet(user->CS, 0.0);
    return;
  }
  
  if(les==1) {
	// 0.2
	PetscReal C_s = 0.03;
    VecSet(user->CS, C_s);
	PetscPrintf(PETSC_COMM_WORLD,"Constant Smagorinsky Cs: %f\n",C_s);
    return;
  }
	double max_cs=0.05;
	//	PetscOptionsInsertFile(NULL, PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  	PetscOptionsGetReal(NULL, PETSC_NULL, "-max_cs", &max_cs, PETSC_NULL);
 
	Vec Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;
	Cmpnts	 ***ucat;
	Cmpnts	***csi, ***eta, ***zet;
	CompVars	***qq;

	PetscReal	***nvert,***Cs, ***lp;//, ***lnu_t;

	DM		da = user->da, fda = user->fda, cda=user->cda;
	DMDALocalInfo	info;
	PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
	PetscInt	mx, my, mz; // Dimensions in three directions
	PetscInt	i, j, k;
  
	PetscReal	***aj;
	Cmpnts ***Ax, ***Ay, ***Az, ***cent, ***beta;//, ;
	double ***LM, ***MM;//, ***NM;
	PetscReal ***Sabs;//, ***Cs, ***lCs_o;//, ***nu_t;

	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	PetscReal	ajc;

	PetscReal	dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;
	
	DMDAGetLocalInfo(da, &info);
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

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

	Vec lSx;		// Sx :	 Sxx, Sxy, Sxz
	Vec lSy;		// Sy :	 Syx, Syy, Syz
	Vec lSz;		// Sz :	 Szx, Szy, Szz
	
	Vec lS;
	//DMDAVecGetArray(fda, user->lUcont, &ucont);
	DMDAVecGetArray(cda, user->lQ,  &qq);
	DMDAVecGetArray(fda, Csi, &csi);
	DMDAVecGetArray(fda, Eta, &eta);
	DMDAVecGetArray(fda, Zet, &zet);
	DMDAVecGetArray(da, user->lNvert, &nvert);
	DMDAVecGetArray(da, user->lAj, &aj); 

	DMDAVecGetArray(fda, user->lBeta, &beta); 
	

	
	DMDAVecGetArray(da, user->CS, &Cs);
	DMDAVecGetArray(da, user->lP, &lp);
	
	DMDAVecGetArray(fda, user->lUcat, &ucat);


  
	VecDuplicate(user->lNvert, &user->lLM);
	VecDuplicate(user->lNvert, &user->lMM);
	//VecDuplicate(user->lP, &user->lNM);
	
	VecSet(user->lLM, 0);
	VecSet(user->lMM, 0);
	//VecSet(user->lNM, 0);
	
       	VecDuplicate(user->lUcat, &lSx);
	VecDuplicate(user->lUcat, &lSy);
	VecDuplicate(user->lUcat, &lSz);
	
	VecDuplicate(user->lNvert, &lS);

	VecSet(lSx, 0);  
	VecSet(lSy, 0);  
	VecSet(lSz, 0);	
	VecSet(lS, 0);
	
	DMDAVecGetArray(da, user->lLM, &LM);//
	DMDAVecGetArray(da, user->lMM, &MM);//
	//DMDAVecGetArray(da, user->lNM, &NM);//
  	
	DMDAVecGetArray(fda, lSx, &Ax);//
	DMDAVecGetArray(fda, lSy, &Ay);//
	DMDAVecGetArray(fda, lSz, &Az);//
	DMDAVecGetArray(da, lS, &Sabs);//
  	
	
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if( nvert[k][j][i]>1.1) continue;

		ucat[k][j][i].x=qq[k][j][i].rhoU/qq[k][j][i].rho;
		ucat[k][j][i].y=qq[k][j][i].rhoV/qq[k][j][i].rho;
		ucat[k][j][i].z=qq[k][j][i].rhoW/qq[k][j][i].rho;
		
		ajc = aj[k][j][i];
		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

		double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
		
		Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
		Compute_du_dxyz (	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
							&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );
		double Sdd=2./3.*(du_dx+dv_dy+dw_dz);	
		double Sxx = 0.5*( du_dx + du_dx )-Sdd , Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
		double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy)-Sdd,	Syz = 0.5*(dv_dz + dw_dy);
		double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz)-Sdd;
	
		Sabs[k][j][i] = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz ));
		
		//debug
		/*
		double dx, dy, dz;
		Calculate_dxdydz(ajc, csi[k][j][i], eta[k][j][i], zet[k][j][i], &dx, &dy, &dz);
		printf("%f %f %f\n", dx, dy, dz);
		*/
		//	printf("%e %e\n", du_dx * csi0, Sx[k][j][i].x);	
		/*
			du_dx * csi0 ~ Sx[k][j][i].x
			
			Uu ~ u^2 * area
			Sij ~ area * Sxx ~ csi * Sxx
			|S| Sij ~ du/dx * ( du/dx ) * (volume * dc/dx )
		*/
		

		/* 06.08.2009 ; has significant effects on subrid model performance */
		/*
		double trace =  Sxx + Syy + Szz;
		Sxx -= trace / 3.;
		Syy -= trace / 3.;
		Szz -= trace / 3.;
		absorbed in the pressure gradient
		*/
		/*
		Sx[k][j][i].x = Sxx;	Sx[k][j][i].y = Sxy;	Sx[k][j][i].z = Sxz;	// in physical space
		Sy[k][j][i].x = Syx;	Sy[k][j][i].y = Syy;	Sy[k][j][i].z = Syz;
		Sz[k][j][i].x = Szx ;	Sz[k][j][i].y = Szy ;	Sz[k][j][i].z = Szz;
		*/
		
		Ax[k][j][i].x = du_dx;	Ax[k][j][i].y = du_dy;	Ax[k][j][i].z = du_dz;
		Ay[k][j][i].x = dv_dx;	Ay[k][j][i].y = dv_dy;	Ay[k][j][i].z = dv_dz;
		Az[k][j][i].x = dw_dx;	Az[k][j][i].y = dw_dy;	Az[k][j][i].z = dw_dz;
	}
	
	DMDAVecRestoreArray(fda, lSx, &Ax);//
	DMDAVecRestoreArray(fda, lSy, &Ay);//
	DMDAVecRestoreArray(fda, lSz, &Az);//
	DMDAVecRestoreArray(da, lS, &Sabs);//
		
	DMLocalToLocalBegin(fda, lSx, INSERT_VALUES, lSx);
	DMLocalToLocalEnd(fda, lSx, INSERT_VALUES, lSx);
 
	DMLocalToLocalBegin(fda, lSy, INSERT_VALUES, lSy);
	DMLocalToLocalEnd(fda, lSy, INSERT_VALUES, lSy);
 
	DMLocalToLocalBegin(fda, lSz, INSERT_VALUES, lSz);
	DMLocalToLocalEnd(fda, lSz, INSERT_VALUES, lSz);
 
	DMLocalToLocalBegin(da, lS, INSERT_VALUES, lS);
	DMLocalToLocalEnd(da, lS, INSERT_VALUES, lS);
	
	DMDAVecGetArray(fda, lSx, &Ax);//
	DMDAVecGetArray(fda, lSy, &Ay);//
	DMDAVecGetArray(fda, lSz, &Az);//
	DMDAVecGetArray(da, lS, &Sabs);//
	
//amir

	if (user->bctype[0]==7 || user->bctype[1]==7){
	  if (xs==0){
	    i=xs;
	    for (k=lzs; k<lze; k++) {
	      for (j=lys; j<lye; j++) {
		Ax[k][j][i]=Ax[k][j][i-2];
		Ay[k][j][i]=Ay[k][j][i-2];	
		Az[k][j][i]=Az[k][j][i-2];
	      }
	    }
	  }
	  if (xe==mx){
	    i=mx-1;
	    for (k=lzs; k<lze; k++) {
	      for (j=lys; j<lye; j++) {
		Ax[k][j][i]=Ax[k][j][i+2];
		Ay[k][j][i]=Ay[k][j][i+2];	
		Az[k][j][i]=Az[k][j][i+2];
	      }
	    }
	  }
	}
	if (user->bctype[2]==7 || user->bctype[3]==7){
	  if (ys==0){
	    j=ys;
	    for (k=lzs; k<lze; k++) {
	      for (i=lxs; i<lxe; i++) {
		Ax[k][j][i]=Ax[k][j-2][i];
		Ay[k][j][i]=Ay[k][j-2][i];	
		Az[k][j][i]=Az[k][j-2][i];
	      }
	    }
	  }
	  if (ye==my){
	    j=my-1;
	    for (k=lzs; k<lze; k++) {
	      for (i=lxs; i<lxe; i++) {
		Ax[k][j][i]=Ax[k][j+2][i];
		Ay[k][j][i]=Ay[k][j+2][i];	
		Az[k][j][i]=Az[k][j+2][i];
	      }
	    }
	  }
	}
	if (user->bctype[4]==7 || user->bctype[5]==7){
	  if (zs==0){
	    k=zs;
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		Ax[k][j][i]=Ax[k-2][j][i];
		Ay[k][j][i]=Ay[k-2][j][i];	
		Az[k][j][i]=Az[k-2][j][i];
	      }
	    }
	  }
	  if (ze==mz){
	    k=mz-1;
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		Ax[k][j][i]=Ax[k+2][j][i];
		Ay[k][j][i]=Ay[k+2][j][i];	
		Az[k][j][i]=Az[k+2][j][i];
	      }
	    }
	  }
	}
	DMDAVecRestoreArray(fda, Csi, &csi);
	DMDAVecRestoreArray(fda, Eta, &eta);
	DMDAVecRestoreArray(fda, Zet, &zet);
	


	
 	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]>1.1) {
			LM[k][j][i]=MM[k][j][i]=0;//NM[k][j][i]=0;
			continue;
		}
		ajc = aj[k][j][i];
	
		double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
		double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
		double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
	
		int a, b;
		double Lij[3][3], Sij_hat[3][3], SSij_hat[3][3], Mij[3][3], Nij[3][3], Nij_cat[3][3];
		double Mij_cat[3][3],Sxyz[3][3];
		
		double filter, test_filter;
		double S[3][3][3], S_hat;
		double rhou[3][3][3], rhov[3][3][3], rhow[3][3][3], rho[3][3][3];
		double u[3][3][3], v[3][3][3], w[3][3][3];
		//	double U[3][3][3], V[3][3][3], W[3][3][3];
		//	double Uu[3][3][3], Uv[3][3][3], Uw[3][3][3];
		//	double Vu[3][3][3], Vv[3][3][3], Vw[3][3][3];
		//	double Wu[3][3][3], Wv[3][3][3], Ww[3][3][3];
		
		double rhouu[3][3][3], rhouv[3][3][3], rhouw[3][3][3];
		double rhovv[3][3][3], rhovw[3][3][3], rhoww[3][3][3];
		double S11[3][3][3], S12[3][3][3], S13[3][3][3], S21[3][3][3], S22[3][3][3], S23[3][3][3], S31[3][3][3], S32[3][3][3], S33[3][3][3];
		double SS11[3][3][3], SS12[3][3][3], SS13[3][3][3], SS21[3][3][3], SS22[3][3][3], SS23[3][3][3], SS31[3][3][3], SS32[3][3][3], SS33[3][3][3];
		
		//		double dU_dx[3][3][3], dU_dy[3][3][3], dU_dz[3][3][3];
		//	double dV_dx[3][3][3], dV_dy[3][3][3], dV_dz[3][3][3];
		//	double dW_dx[3][3][3], dW_dy[3][3][3], dW_dz[3][3][3];
		//	double T11[3][3][3], T12[3][3][3], T13[3][3][3], T21[3][3][3], T22[3][3][3], T23[3][3][3], T31[3][3][3], T32[3][3][3], T33[3][3][3];	// Clark model term
		
		double weight[3][3][3];
		double pu[3][3][3];
		double pv[3][3][3];
		double pw[3][3][3];
		double pp[3][3][3];


		double sum_weight=0;
		double coef[3][3][3]={
			0.125, 0.250, 0.125, 
			0.250, 0.500, 0.250, 
			0.125, 0.250, 0.125, 
				
			0.250, 0.500, 0.250,
			0.500, 1.000, 0.500,
			0.250, 0.500, 0.250,
				
			0.125, 0.250, 0.125, 
			0.250, 0.500, 0.250,
			0.125, 0.250, 0.125
		};
//test	
/*	double coef[3][3][3]={
			0.0, 0.00, 0.0, 
			0.00, 0.0, 0.0, 
			0.00, 0.00, 0.0, 
				
			0.000, 0.000, 0.000,
			0.0, 1.000, 0.000,
			0.000, 0.00, 0.000,
				
			0.00, 0.00, 0.00, 
			0.00, 0.00, 0.000,
			0.0, 0.0, 0.0
		};
*/


		double dx,dy,dz;
		//	Calculate_dxdydz(aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i], &dx, &dy, &dz);
		//	double dx2=dx*dx, dy2=dy*dy, dz2=dz*dz;
		

		int p,q,r;
		for(p=-1; p<=1; p++)
		for(q=-1; q<=1; q++)
		for(r=-1; r<=1; r++) {
			int R=r+1, Q=q+1, P=p+1;
			int KK=k+r, JJ=j+q, II=i+p;

			u[R][Q][P] = v[R][Q][P] = w[R][Q][P] = 0;
			rhou[R][Q][P] = rhov[R][Q][P] = rhow[R][Q][P] = 0;
			rhouu[R][Q][P] = rhouv[R][Q][P] = rhouw[R][Q][P] = 0;
			rhovv[R][Q][P] = rhovw[R][Q][P] = rhoww[R][Q][P] = 0;
			
			S11[R][Q][P] = S12[R][Q][P] = S13[R][Q][P] = 0;
			S21[R][Q][P] = S22[R][Q][P] = S23[R][Q][P] = 0;
			S31[R][Q][P] = S32[R][Q][P] = S33[R][Q][P] = 0;
						
			S[R][Q][P] = 0;
				
			SS11[R][Q][P] = SS12[R][Q][P] = SS13[R][Q][P] = 0;
			SS21[R][Q][P] = SS22[R][Q][P] = SS23[R][Q][P] = 0;
			SS31[R][Q][P] = SS32[R][Q][P] = SS33[R][Q][P] = 0;
			
	//		weight[R][Q][P] = 1./aj[KK][JJ][II]*qq[KK][JJ][II].rho*coef[R][Q][P];
			weight[R][Q][P] =1./(fabs(R-1.0)+1)*1./(fabs(Q-1.0)+1)*1./(fabs(P-1.0)+1);
		
			//rhou[R][Q][P]
			u[R][Q][P] = ucat[KK][JJ][II].x;
                        v[R][Q][P] = ucat[KK][JJ][II].y;
                        w[R][Q][P] = ucat[KK][JJ][II].z;

			  
			rhou[R][Q][P] = qq[KK][JJ][II].rhoU;
			rhov[R][Q][P] = qq[KK][JJ][II].rhoV;
			rhow[R][Q][P] = qq[KK][JJ][II].rhoW;
			 rho[R][Q][P] = qq[KK][JJ][II].rho;
			// metric tensors are also test-filtered : Big difference
			
					
			/*
			Uu[R][Q][P] = U[R][Q][P]*u[R][Q][P];
			Uv[R][Q][P] = U[R][Q][P]*v[R][Q][P];
			Uw[R][Q][P] = U[R][Q][P]*w[R][Q][P];
			
			Vu[R][Q][P] = V[R][Q][P]*u[R][Q][P];
			Vv[R][Q][P] = V[R][Q][P]*v[R][Q][P];
			Vw[R][Q][P] = V[R][Q][P]*w[R][Q][P];
			
			Wu[R][Q][P] = W[R][Q][P]*u[R][Q][P];
			Wv[R][Q][P] = W[R][Q][P]*v[R][Q][P];
			Ww[R][Q][P] = W[R][Q][P]*w[R][Q][P];
			*/ 

			// rho u_i u_j is computed.
			rhouu[R][Q][P] = rhou[R][Q][P]* u[R][Q][P];	
			rhouv[R][Q][P] = rhou[R][Q][P]*v[R][Q][P];	
			rhouw[R][Q][P] = rhou[R][Q][P]*w[R][Q][P];
			rhovv[R][Q][P] = rhov[R][Q][P]*v[R][Q][P];	
			rhovw[R][Q][P] = rhov[R][Q][P]*w[R][Q][P];	
			rhoww[R][Q][P] = rhow[R][Q][P]*w[R][Q][P];
			
			const double du_dx = Ax[KK][JJ][II].x, du_dy = Ax[KK][JJ][II].y, du_dz = Ax[KK][JJ][II].z;
			const double dv_dx = Ay[KK][JJ][II].x, dv_dy = Ay[KK][JJ][II].y, dv_dz = Ay[KK][JJ][II].z;
			const double dw_dx = Az[KK][JJ][II].x, dw_dy = Az[KK][JJ][II].y, dw_dz = Az[KK][JJ][II].z;
		
			const double Sxx = 0.5*( du_dx + du_dx ), Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
			const double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy),	Syz = 0.5*(dv_dz + dw_dy);
			const double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz);
			
			
	/*		dU_dx[R][Q][P] = du_dx, dU_dy[R][Q][P] = du_dy, dU_dz[R][Q][P] = du_dz;
			dV_dx[R][Q][P] = dv_dx, dV_dy[R][Q][P] = dv_dy, dV_dz[R][Q][P] = dv_dz;
			dW_dx[R][Q][P] = dw_dx, dW_dy[R][Q][P] = dw_dy, dW_dz[R][Q][P] = dw_dz;
			
			T11[R][Q][P] = ( du_dx * du_dx * dx2 + du_dy * du_dy * dy2 + du_dz * du_dz * dz2 ) / 12.;
			T12[R][Q][P] = ( du_dx * dv_dx * dx2 + du_dy * dv_dy * dy2 + du_dz * dv_dz * dz2 ) / 12.;
			T13[R][Q][P] = ( du_dx * dw_dx * dx2 + du_dy * dw_dy * dy2 + du_dz * dw_dz * dz2 ) / 12.;
			T21[R][Q][P] = T12[R][Q][P];
			T22[R][Q][P] = ( dv_dx * dv_dx * dx2 + dv_dy * dv_dy * dy2 + dv_dz * dv_dz * dz2 ) / 12.;
			T23[R][Q][P] = ( dv_dx * dw_dx * dx2 + dv_dy * dw_dy * dy2 + dv_dz * dw_dz * dz2 ) / 12.;
			T31[R][Q][P] = T13[R][Q][P];
			T32[R][Q][P] = T23[R][Q][P];
			T33[R][Q][P] = ( dw_dx * dw_dx * dx2 + dw_dy * dw_dy * dy2 + dw_dz * dw_dz * dz2 ) / 12.;
			*/
			S11[R][Q][P] = Sxx, S12[R][Q][P] = Sxy, S13[R][Q][P] = Sxz;
			S21[R][Q][P] = Syx, S22[R][Q][P] = Syy, S23[R][Q][P] = Syz;
			S31[R][Q][P] = Szx, S32[R][Q][P] = Szy, S33[R][Q][P] = Szz;
			
			S[R][Q][P] = Sabs[KK][JJ][II];
			
			SS11[R][Q][P] = S11[R][Q][P]*S[R][Q][P]*rho[R][Q][P], SS12[R][Q][P] = S12[R][Q][P]*S[R][Q][P]*rho[R][Q][P], SS13[R][Q][P] = S13[R][Q][P]*S[R][Q][P]*rho[R][Q][P];
			SS21[R][Q][P] = S21[R][Q][P]*S[R][Q][P]*rho[R][Q][P], SS22[R][Q][P] = S22[R][Q][P]*S[R][Q][P]*rho[R][Q][P], SS23[R][Q][P] = S23[R][Q][P]*S[R][Q][P]*rho[R][Q][P];
			SS31[R][Q][P] = S31[R][Q][P]*S[R][Q][P]*rho[R][Q][P], SS32[R][Q][P] = S32[R][Q][P]*S[R][Q][P]*rho[R][Q][P], SS33[R][Q][P] = S33[R][Q][P]*S[R][Q][P]*rho[R][Q][P];
///////////////////for energy subgrid//////
			  pu[R][Q][P]=ucat[KK][JJ][II].x*lp[KK][JJ][II];
			  pv[R][Q][P]=ucat[KK][JJ][II].y*lp[KK][JJ][II]; 
			  pw[R][Q][P]=ucat[KK][JJ][II].z*lp[KK][JJ][II]; 
			  pp[R][Q][P]=lp[KK][JJ][II]; 
			



		}
		
		//Amir
		
		for(p=-1; p<=1; p++)
		for(q=-1; q<=1; q++)
		for(r=-1; r<=1; r++) sum_weight += weight[r+1][q+1][p+1] * coef[r+1][q+1][p+1];
		
		filter = pow( 1./aj[k][j][i], 1./3. );
		
		if(testfilter_ik) test_filter = pow(5.0, 1./3.) * filter;
		else {
			test_filter = 2.0 * filter;
		//	test_filter = pow( sum_weight, 1./3. );
		}
		
		//		double _U=integrate_testfilter(U, weight);
		//		double _V=integrate_testfilter(V, weight);
		//		double _W=integrate_testfilter(W, weight);
		double _u=integrate_testfilter(u, weight);
		double _v=integrate_testfilter(v, weight);
		double _w=integrate_testfilter(w, weight);	
		double _rho=integrate_testfilter(rho, weight);
		double _rhou=integrate_testfilter(rhou, weight);
		double _rhov=integrate_testfilter(rhov, weight);
		double _rhow=integrate_testfilter(rhow, weight);
// there is issue with the model
		double _p=integrate_testfilter(pp, weight);
		double _pu=integrate_testfilter(pu, weight);
		double _pv=integrate_testfilter(pv, weight);
		double _pw=integrate_testfilter(pw, weight);
	 	beta[k][j][i].x = _pu - _p*_u;	
		beta[k][j][i].y = _pv - _p*_v;//_pv -_v*_p;
		beta[k][j][i].z = _pw - _p*_w;//_pw -_w*_p;
		/*		
		double _T11=integrate_testfilter(T11, weight);
		double _T12=integrate_testfilter(T12, weight);
		double _T13=integrate_testfilter(T13, weight);
		double _T21=integrate_testfilter(T21, weight);
		double _T22=integrate_testfilter(T22, weight);
		double _T23=integrate_testfilter(T23, weight);
		double _T31=integrate_testfilter(T31, weight);
		double _T32=integrate_testfilter(T32, weight);
		double _T33=integrate_testfilter(T33, weight);
		
		double _dU_dx=integrate_testfilter(dU_dx, weight);
		double _dV_dx=integrate_testfilter(dV_dx, weight);
		double _dW_dx=integrate_testfilter(dW_dx, weight);
		double _dU_dy=integrate_testfilter(dU_dy, weight);
		double _dV_dy=integrate_testfilter(dV_dy, weight);
		double _dW_dy=integrate_testfilter(dW_dy, weight);
		double _dU_dz=integrate_testfilter(dU_dz, weight);
		double _dV_dz=integrate_testfilter(dV_dz, weight);
		double _dW_dz=integrate_testfilter(dW_dz, weight);
		
		double _R11 = ( _dU_dx * _dU_dx * dx2 + _dU_dy * _dU_dy * dy2 + _dU_dz * _dU_dz * dz2 ) / 12.;
		double _R12 = ( _dU_dx * _dV_dx * dx2 + _dU_dy * _dV_dy * dy2 + _dU_dz * _dV_dz * dz2 ) / 12.;
		double _R13 = ( _dU_dx * _dW_dx * dx2 + _dU_dy * _dW_dy * dy2 + _dU_dz * _dW_dz * dz2 ) / 12.;
		double _R21 = _R12;
		double _R22 = ( _dV_dx * _dV_dx * dx2 + _dV_dy * _dV_dy * dy2 + _dV_dz * _dV_dz * dz2 ) / 12.;
		double _R23 = ( _dV_dx * _dW_dx * dx2 + _dV_dy * _dW_dy * dy2 + _dV_dz * _dW_dz * dz2 ) / 12.;
		double _R31 = _R13;
		double _R32 = _R23;
		double _R33 = ( _dW_dx * _dW_dx * dx2 + _dW_dy * _dW_dy * dy2 + _dW_dz * _dW_dz * dz2 ) / 12.;
			
		Lij[0][0] = integrate_testfilter(Uu, weight) - _U*_u;
		Lij[0][1] = integrate_testfilter(Uv, weight) - _U*_v;
		Lij[0][2] = integrate_testfilter(Uw, weight) - _U*_w;
		Lij[1][0] = integrate_testfilter(Vu, weight) - _V*_u;
		Li
j[1][1] = integrate_testfilter(Vv, weight) - _V*_v;
		Lij[1][2] = integrate_testfilter(Vw, weight) - _V*_w;
		Lij[2][0] = integrate_testfilter(Wu, weight) - _W*_u;
		Lij[2][1] = integrate_testfilter(Wv, weight) - _W*_v;
		Lij[2][2] = integrate_testfilter(Ww, weight) - _W*_w;
		*/

	//amir 
		Lij[0][0] = integrate_testfilter(rhouu, weight) - _rhou*_rhou/_rho;
		Lij[0][1] = integrate_testfilter(rhouv, weight) - _rhou*_rhov/_rho;
		Lij[0][2] = integrate_testfilter(rhouw, weight) - _rhou*_rhow/_rho;
		Lij[1][0] = integrate_testfilter(rhouv, weight) - _rhov*_rhou/_rho;
		Lij[1][1] = integrate_testfilter(rhovv, weight) - _rhov*_rhov/_rho;
		Lij[1][2] = integrate_testfilter(rhovw, weight) - _rhov*_rhow/_rho;
		Lij[2][0] = integrate_testfilter(rhouw, weight) - _rhow*_rhou/_rho;
		Lij[2][1] = integrate_testfilter(rhovw, weight) - _rhow*_rhov/_rho;
		Lij[2][2] = integrate_testfilter(rhoww, weight) - _rhow*_rhow/_rho;
	
		Sxyz[0][0] = S11[1][1][1];
		Sxyz[0][1] = S12[1][1][1];
		Sxyz[0][2] = S13[1][1][1];
		Sxyz[1][0] = S21[1][1][1];
		Sxyz[1][1] = S22[1][1][1];
		Sxyz[1][2] = S23[1][1][1];
		Sxyz[2][0] = S31[1][1][1];
		Sxyz[2][1] = S32[1][1][1];
		Sxyz[2][2] = S33[1][1][1];
		/*			
		Nij_cat[0][0] = 4.0 * _R11  - _T11;
		Nij_cat[0][1] = 4.0 * _R12  - _T12;
		Nij_cat[0][2] = 4.0 * _R13  - _T13;
		Nij_cat[1][0] = 4.0 * _R21  - _T21;
		Nij_cat[1][1] = 4.0 * _R22  - _T22;
		Nij_cat[1][2] = 4.0 * _R23  - _T23;
		Nij_cat[2][0] = 4.0 * _R31  - _T31;
		Nij_cat[2][1] = 4.0 * _R32  - _T32;
		Nij_cat[2][2] = 4.0 * _R33  - _T33;
		
		Nij[0][0] = Nij_cat[0][0] * csi0 + Nij_cat[0][1] * csi1 + Nij_cat[0][2] * csi2;
		Nij[0][1] = Nij_cat[0][0] * eta0 + Nij_cat[0][1] * eta1 + Nij_cat[0][2] * eta2;
		Nij[0][2] = Nij_cat[0][0] * zet0 + Nij_cat[0][1] * zet1 + Nij_cat[0][2] * zet2;
		Nij[1][0] = Nij_cat[1][0] * csi0 + Nij_cat[1][1] * csi1 + Nij_cat[1][2] * csi2;
		Nij[1][1] = Nij_cat[1][0] * eta0 + Nij_cat[1][1] * eta1 + Nij_cat[1][2] * eta2;
		Nij[1][2] = Nij_cat[1][0] * zet0 + Nij_cat[1][1] * zet1 + Nij_cat[1][2] * zet2;
		Nij[2][0] = Nij_cat[2][0] * csi0 + Nij_cat[2][1] * csi1 + Nij_cat[2][2] * csi2;
		Nij[2][1] = Nij_cat[2][0] * eta0 + Nij_cat[2][1] * eta1 + Nij_cat[2][2] * eta2;
		Nij[2][2] = Nij_cat[2][0] * zet0 + Nij_cat[2][1] * zet1 + Nij_cat[2][2] * zet2;
		*/

		Sij_hat[0][0] = integrate_testfilter(S11, weight);	
		Sij_hat[0][1] = integrate_testfilter(S12, weight);	
		Sij_hat[0][2] = integrate_testfilter(S13, weight);
		Sij_hat[1][0] = integrate_testfilter(S21, weight);	
		Sij_hat[1][1] = integrate_testfilter(S22, weight);	
		Sij_hat[1][2] = integrate_testfilter(S23, weight);
		Sij_hat[2][0] = integrate_testfilter(S31, weight);	
		Sij_hat[2][1] = integrate_testfilter(S32, weight);	
		Sij_hat[2][2] = integrate_testfilter(S33, weight);
		SSij_hat[0][0] = integrate_testfilter(SS11, weight);	
		SSij_hat[0][1] = integrate_testfilter(SS12, weight);	
		SSij_hat[0][2] = integrate_testfilter(SS13, weight);
		SSij_hat[1][0] = integrate_testfilter(SS21, weight);
		SSij_hat[1][1] = integrate_testfilter(SS22, weight);	
		SSij_hat[1][2] = integrate_testfilter(SS23, weight);
		SSij_hat[2][0] = integrate_testfilter(SS31, weight);	
		SSij_hat[2][1] = integrate_testfilter(SS32, weight);	
		SSij_hat[2][2] = integrate_testfilter(SS33, weight);
		
		
		S_hat=0;
		for(a=0; a<3; a++)
		for(b=0; b<3; b++) {
			S_hat += pow( Sij_hat[a][b], 2. );
		}
		S_hat = sqrt ( 2 * S_hat );
		
		//S_hat = integrate_testfilter(S, weight);
	/*	
		
		double gg[3][3], ggc[3][3], G[3][3];
		double xcsi, xeta, xzet, ycsi, yeta, yzet, zcsi, zeta, zzet;

		gg[0][0]=csi0, gg[0][1]=csi1, gg[0][2]=csi2;
		gg[1][0]=eta0, gg[1][1]=eta1, gg[1][2]=eta2;
		gg[2][0]=zet0, gg[2][1]=zet1, gg[2][2]=zet2;
		Calculate_Covariant_metrics(gg, ggc);
		xcsi=ggc[0][0], xeta=ggc[0][1], xzet=ggc[0][2];
		ycsi=ggc[1][0], yeta=ggc[1][1], yzet=ggc[1][2];
		zcsi=ggc[2][0], zeta=ggc[2][1], zzet=ggc[2][2];
		G[0][0] = xcsi * xcsi + ycsi * ycsi + zcsi * zcsi;
		G[1][1] = xeta * xeta + yeta * yeta + zeta * zeta;
		G[2][2] = xzet * xzet + yzet * yzet + zzet * zzet;
		G[0][1] = G[1][0] = xeta * xcsi + yeta * ycsi + zeta * zcsi;
		G[0][2] = G[2][0] = xzet * xcsi + yzet * ycsi + zzet * zcsi;
		G[1][2] = G[2][1] = xeta * xzet + yeta * yzet + zeta * zzet;*/
		double L_cat=0.0, M_cat=0.0;
		for(a=0; a<3; a++)
		for(b=0; b<3; b++) {
		  M_cat += -2.0*_rho*pow( test_filter, 2. ) * S_hat * (Sij_hat[a][b] * Sxyz[a][b] -0.33* Sij_hat[a][a]*Sxyz[b][b])
		    + 2.0*pow( filter, 2. ) * (SSij_hat[a][b]*Sxyz[a][b]-0.33*SSij_hat[a][a]*Sxyz[b][b]);
		  L_cat += Lij[a][b]*Sxyz[a][b]-0.33*Sxyz[a][a]*Lij[b][b];
			/*
			Mij_cat_i[a][b] = -4. * S_hat_i * Sij_hat_i[a][b] + SSij_hat_i[a][b];
			Mij_cat_j[a][b] = -4. * S_hat_j * Sij_hat_j[a][b] + SSij_hat_j[a][b];
			Mij_cat_k[a][b] = -4. * S_hat_k * Sij_hat_k[a][b] + SSij_hat_k[a][b];*/
		}
		    /*
		Mij[0][0] = Mij_cat[0][0] * csi0 + Mij_cat[0][1] * csi1 + Mij_cat[0][2] * csi2;
		Mij[0][1] = Mij_cat[0][0] * eta0 + Mij_cat[0][1] * eta1 + Mij_cat[0][2] * eta2;
		Mij[0][2] = Mij_cat[0][0] * zet0 + Mij_cat[0][1] * zet1 + Mij_cat[0][2] * zet2;
		Mij[1][0] = Mij_cat[1][0] * csi0 + Mij_cat[1][1] * csi1 + Mij_cat[1][2] * csi2;
		Mij[1][1] = Mij_cat[1][0] * eta0 + Mij_cat[1][1] * eta1 + Mij_cat[1][2] * eta2;
		Mij[1][2] = Mij_cat[1][0] * zet0 + Mij_cat[1][1] * zet1 + Mij_cat[1][2] * zet2;
		Mij[2][0] = Mij_cat[2][0] * csi0 + Mij_cat[2][1] * csi1 + Mij_cat[2][2] * csi2;
		Mij[2][1] = Mij_cat[2][0] * eta0 + Mij_cat[2][1] * eta1 + Mij_cat[2][2] * eta2;
		Mij[2][2] = Mij_cat[2][0] * zet0 + Mij_cat[2][1] * zet1 + Mij_cat[2][2] * zet2;
		    		
		double num=0, denom=0;
		int m, n, l;
		
		
			g11 ~ csi*csi ~ dx^4
			G11 ~ dx^-4
		
	
	
		for(q=0; q<3; q++)
		for(a=0; a<3; a++)
		for(b=0; b<3; b++) {
			num += Lij[b][a] * Mij[a][q] * G[b][q];
			//num += Lij[b][a] * Mij_cat[a][b];
			//if(clark) num -= Nij[a][b] * Mij[a][q] * G[b][q];
		}
		
         	for(m=0; m<3; m++)
		for(n=0; n<3; n++)
		for(l=0; l<3; l++) {
			denom += Mij[n][m] * Mij[n][l] * G[m][l];
			//denom += Mij[l][n] * Mij[l][n];
		}
	
		//printf("%f %f\n", num, denom);
		*/
		LM[k][j][i] = L_cat;
		MM[k][j][i] = M_cat;
    	}
	
	DMDAVecRestoreArray(da, user->lP, &lp);

	DMDAVecRestoreArray(da, user->lLM, &LM);//
	DMDAVecRestoreArray(da, user->lMM, &MM);//
	//DMDAVecRestoreArray(da, user->lNM, &NM);//
	
	DMLocalToLocalBegin(da, user->lLM, INSERT_VALUES, user->lLM);
	DMLocalToLocalEnd(da, user->lLM, INSERT_VALUES, user->lLM);
	DMLocalToLocalBegin(da, user->lMM, INSERT_VALUES, user->lMM);
	DMLocalToLocalEnd(da, user->lMM, INSERT_VALUES, user->lMM);
	//DALocalToLocalBegin(da, user->lNM, INSERT_VALUES, user->lNM);
	//DALocalToLocalEnd(da, user->lNM, INSERT_VALUES, user->lNM);
	  
	DMDAVecGetArray(da, user->lLM, &LM);//
	DMDAVecGetArray(da, user->lMM, &MM);//
	//DMDAVecGetArray(da, user->lNM, &NM);//
	
/*	
	if (user->bctype[0]==7 || user->bctype[1]==7){
	  if (xs==0){
	    i=xs;
	    for (k=lzs; k<lze; k++) {
	      for (j=lys; j<lye; j++) {
		LM[k][j][i]=LM[k][j][i-2];
		MM[k][j][i]=MM[k][j][i-2];	
	
	      }
	    }
	  }
	  if (xe==mx){
	    i=mx-1;
	    for (k=lzs; k<lze; k++) {
	      for (j=lys; j<lye; j++) {
		LM[k][j][i]=LM[k][j][i+2];
		MM[k][j][i]=MM[k][j][i+2];	

	      }
	    }
	  }
	}
	if (user->bctype[2]==7 || user->bctype[3]==7){
	  if (ys==0){
	    j=ys;
	    for (k=lzs; k<lze; k++) {
	      for (i=lxs; i<lxe; i++) {
		LM[k][j][i]=LM[k][j-2][i];
		MM[k][j][i]=MM[k][j-2][i];	

	      }
	    }
	  }
	  if (ye==my){
	    j=my-1;
	    for (k=lzs; k<lze; k++) {
	      for (i=lxs; i<lxe; i++) {
		LM[k][j][i]=LM[k][j+2][i];
		MM[k][j][i]=MM[k][j+2][i];	

	      }
	    }
	  }
	}
	if (user->bctype[4]==7 || user->bctype[5]==7){
	  if (zs==0){
	    k=zs;
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		LM[k][j][i]=LM[k-2][j][i];
		MM[k][j][i]=MM[k-2][j][i];	
	      }
	    }
	  }
	  if (ze==mz){
	    k=mz-1;
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		LM[k][j][i]=LM[k+2][j][i];
		MM[k][j][i]=MM[k+2][j][i];	
	      }
	    }
	  }
	}

*/
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i]>1.1) {
			Cs[k][j][i] = 0;
			continue;
		}
				
		double weight[3][3][3];
		double LM0[3][3][3], MM0[3][3][3];
		int a, b, c, a0=1;
			
		for(a=-a0; a<=a0; a++)
		for(b=-a0; b<=a0; b++)
		for(c=-a0; c<=a0; c++) {
			int R=c+a0, Q=b+a0, P=a+a0;
			int KK=k+c, JJ=j+b, II=i+a;
			
		
			LM0[R][Q][P] = LM[KK][JJ][II];
			MM0[R][Q][P] = MM[KK][JJ][II];
		}
			
		double C=0;
			
		double LM_avg, MM_avg;//, NM_avg;

		if ( i_homo_filter || j_homo_filter || k_homo_filter || les==3 ) {
			LM_avg = LM[k][j][i];
			MM_avg = MM[k][j][i];
		}
		else {
			LM_avg = integrate_t(LM0);
			MM_avg = integrate_t(MM0);
		}
	// ............................................ positive 	
		C =  LM_avg / (MM_avg + eps );
	//.............................................	
	
		if(les==1) Cs[k][j][i] = 0.01;
		else Cs[k][j][i] = PetscMax(C, 0.0001);
	}
	

	
 	if ( i_homo_filter && k_homo_filter ){
	
	double *LM_tmp,*MM_tmp,*JJ_LM,*JJ_MM;
	int *count,*total_count;
	LM_tmp = malloc (my*sizeof(double));
	MM_tmp = malloc (my*sizeof(double));
	count = malloc (my*sizeof(int));
	total_count = malloc (my*sizeof(int));
	JJ_LM=malloc (my*sizeof(double));
	JJ_MM=malloc (my*sizeof(double));
 		for(j=0; j<my; j++) { 
 			LM_tmp[j] = 0.; 
 			MM_tmp[j] = 0.; 
 			count[j] = 0; 
 		} 
		
		for (k=lzs; k<lze; k++) 
 		for (j=lys; j<lye; j++) 
 		for (i=lxs; i<lxe; i++) { 
 			if( nvert[k][j][i]<0.1 ) { 
 				LM_tmp[j] += LM[k][j][i]; 
 				MM_tmp[j] += MM[k][j][i]; 
 				count[j] ++; 
 			} 
 		} 
		
 		MPI_Allreduce( &LM_tmp[0], &JJ_LM[0], my, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
 		MPI_Allreduce( &MM_tmp[0], &JJ_MM[0], my, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
 		MPI_Allreduce( &count[0], &total_count[0], my, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);	 
		
 		for(j=0; j<my; j++) { 
 			if( total_count[j]>0) { 
 				JJ_LM[j] /= (double) (total_count[j]); 
 				JJ_MM[j] /= (double) (total_count[j]); 

//		PetscPrintf(PETSC_COMM_WORLD, "count = %d ze %d lze %d\n", total_count[j],ze,lze);
 			} 
 		} 
	for (j=lys; j<lye; j++) 
 		for (k=lzs; k<lze; k++) 
 		for (i=lxs; i<lxe; i++) { 
 			Cs[k][j][i] = 0.5 * JJ_LM[j] / ( JJ_MM[j]+eps); 
		} 
	
	free(LM_tmp);free(count);
	free(MM_tmp);free(JJ_MM);
	free(JJ_LM);free(total_count);
 	} //if i_homo

		
 		if(i_homo_filter && !k_homo_filter){
		int plane_size = my*mz;
  		 
		double *LM_tmp,*MM_tmp,*JJ_LM,*JJ_MM;
		PetscInt *count,*total_count;
		LM_tmp = malloc (plane_size*sizeof(double));
		MM_tmp = malloc (plane_size*sizeof(double));
		count = malloc (plane_size*sizeof(PetscInt));
		total_count = malloc (plane_size*sizeof(PetscInt));
		JJ_LM= malloc (plane_size*sizeof(double));
		JJ_MM= malloc (plane_size*sizeof(double));
		
 		int pos; 
 		pos=0; 
 		for(pos=0; pos<plane_size; pos++) { 
 			LM_tmp[pos] = 0.0; 
 			MM_tmp[pos] = 0.0;
			JJ_LM[pos]=0.0; 
			JJ_MM[pos]=0.0; 
			count[pos]= 0;
			total_count[pos] = 0; 
		} 
		
		pos=0;

 		for(k=0; k<mz; k++) 
 		for(j=0; j<my; j++) { 
			if( k>=lzs && k<lze && j>=lys && j<lye) { 
				for (i=lxs; i<lxe; i++) { 
 					if( nvert[k][j][i]<0.1 ) { 
 						LM_tmp[pos] += LM[k][j][i]; 
 						MM_tmp[pos] += MM[k][j][i]; 
 						count[pos] ++;
		
 					} 
 				} 
			} 
 		 	 	 pos++; 
 		} 
 	 
/* 		else if(j_homo_filter)  { */
/* 			for(k=0; k<mz; k++) */
/* 			for(i=0; i<mx; i++) { */
/* 				if( i>=lxs && i<lxe && k>=lzs && k<lze) { */
/* 					for (j=lys; j<lye; j++) { */
/* 						if( nvert[k][j][i]<0.1 ) { */
/* 							LM_tmp[pos] += LM[k][j][i]; */
/* 							MM_tmp[pos] += MM[k][j][i]; */
/* 							count[pos] ++; */
/* 						} */
/* 					} */
/* 				} */
/* 				pos++; */
/* 			} */
/* 		} */
/* 		else if(k_homo_filter)  { */
/* 			for(j=0; j<my; j++) */
/* 			for(i=0; i<mx; i++) { */
/* 				if( i>=lxs && i<lxe && j>=lys && j<lye) { */
/* 					for (k=lzs; k<lze; k++) { */
/* 						if( nvert[k][j][i]<0.1 ) { */
/* 							LM_tmp[pos] += LM[k][j][i]; */
/* 							MM_tmp[pos] += MM[k][j][i]; */
/* 							count[pos] ++; */
/* 						} */
/* 					} */
/* 				} */
/* 				pos++; */
/* 			} */
/* 		} */
		
 		MPI_Allreduce( &LM_tmp[0], &JJ_LM[0], plane_size, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
 		MPI_Allreduce( &MM_tmp[0], &JJ_MM[0], plane_size, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
 		MPI_Allreduce( &count[0], &total_count[0], plane_size, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);	 
		
		pos=0;
		for(pos=0; pos<plane_size; pos++) { 
 			if( total_count[pos]>0) { 
 				double N = (double) (total_count[pos]); 
 				JJ_LM[pos] /= N; 
 				JJ_MM[pos] /= N;

 			} 
 		} 
		
		pos=0; 
 			for(k=zs; k<ze; k++)  
			for(j=ys; j<ye; j++) { 
 				if( k>=lzs && k<lze && j>=lys && j<lye) { 
 					for (i=lxs; i<lxe; i++) { 
 						if( nvert[k][j][i]<0.1 ) { 
 							Cs[k][j][i] = 0.5 * JJ_LM[pos] / ( JJ_MM[pos] + eps ); 
 						} 
 					} 
 				} 
 				
			pos++;
 			} 
	free(LM_tmp);free(count);
	free(MM_tmp);free(JJ_MM);
	free(JJ_LM);free(total_count);
 	}
//	i_homo_filter=j_homo_filter=k_homo_filter=1;
 	if ( i_homo_filter && j_homo_filter && k_homo_filter ){
	
	 double MM_temp=0.0,LM_temp=0.0,MM_temp_sum=0.0,LM_temp_sum=0.0;
	 
		for (k=lzs; k<lze; k++){
		 for (j=lys; j<lye; j++){
		  for (i=lxs; i<lxe; i++){
		    //if (i==1)
		    //	printf("LM %f MM %f  \n",MM[k][j][i],LM[k][j][i]);
		   MM_temp+=MM[k][j][i]; 
		   LM_temp+=LM[k][j][i];
		 }
		}
	       }

 		MPI_Allreduce( &LM_temp, &LM_temp_sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
 		MPI_Allreduce( &MM_temp, &MM_temp_sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); 
		printf("LM %f MM %f  Cs %f\n",LM_temp_sum,MM_temp_sum,LM_temp_sum/MM_temp_sum); 
		for (k=zs; k<ze; k++){
		 for (j=ys; j<ye; j++){
		  for (i=xs; i<xe; i++){
		   
		   Cs[k][j][i]=LM_temp_sum/MM_temp_sum; 
		   
		 }
		}
	       }


 }
/* 		else if(j_homo_filter)  { */
/* 			for(k=0; k<mz; k++) */
/* 			for(i=0; i<mx; i++) { */
/* 				if( i>=lxs && i<lxe && k>=lzs && k<lze) { */
/* 					for (j=lys; j<lye; j++) { */
/* 						if( nvert[k][j][i]<0.1 ) { */
/* 							Cs[k][j][i] = 0.5 * JJ_LM[pos] / ( JJ_MM[pos] + eps ); */
/* 						} */
/* 					} */
/* 				} */
/* 				pos++; */
/* 			} */
/* 		} */
/* 		else if(k_homo_filter)  { */
/* 			for(j=0; j<my; j++) */
/* 			for(i=0; i<mx; i++) { */
/* 				if( i>=lxs && i<lxe && j>=lys && j<lye) { */
/* 					for (k=lzs; k<lze; k++) { */
/* 						if( nvert[k][j][i]<0.1 ) { */
/* 							Cs[k][j][i] = 0.5 * JJ_LM[pos] / ( JJ_MM[pos] + eps ); */
/* 						} */
/* 					} */
/* 				} */
/* 				pos++; */
/* 			} */
/* 		} */
/* 	} */
	

	for (k=zs; k<ze; k++)
	for (j=ys; j<ye; j++)
	for (i=xs; i<xe; i++) {	  
	  if((nvert[k][j][i]>1.1 && nvert[k][j][i]<7.) || (k==0 && !user->bctype[4]==7) || (k==mz-1 && !user->bctype[5]==7) || (j==0 && !user->bctype[2]==7) || (j==my-1 && !user->bctype[3]==7)  || (i==0 && !user->bctype[0]==7) || (i==mx-1 && !user->bctype[1]==7)) {
	    Cs[k][j][i] = 0.0;
	  }
		

	  else {
	    if(nvert[k][j][i]>0.1 && nvert[k][j][i]<1.1) {
	      Cs[k][j][i] = PetscMax(0.01, Cs[k][j][i]);	// stabilize at high Re
	    }

		Cs[k][j][i] = PetscMin( PetscMax ( Cs[k][j][i], 0.000 ), max_cs);
	  }
	}
	

	DMDAVecRestoreArray(fda, lSx, &Ax);//
	DMDAVecRestoreArray(fda, lSy, &Ay);//
	DMDAVecRestoreArray(fda, lSz, &Az);//
	DMDAVecRestoreArray(da, lS, &Sabs);//
	DMDAVecRestoreArray(da, user->lLM, &LM);//
	DMDAVecRestoreArray(da, user->lMM, &MM);//
	

//	DMDAVecRestoreArray(fda, user->lUcont, &ucont);/	
	DMDAVecRestoreArray(cda,user->lQ,  &qq);
	DMDAVecRestoreArray(fda, user->lUcat,  &ucat);
//	DMDAVecRestoreArray(fda, Csi, &csi);
//	DMDAVecRestoreArray(fda, Eta, &eta);
//	DMDAVecRestoreArray(fda, Zet, &zet);
	DMDAVecRestoreArray(da, user->lNvert, &nvert);
	DMDAVecRestoreArray(da, user->lAj, &aj); 
	DMDAVecRestoreArray(da, user->CS, &Cs);
	

	//	DMDAVecRestoreArray(fda, user->lCent, &cent);
	VecDestroy(&user->lLM);//
	VecDestroy(&user->lMM);//
  
	VecDestroy(&lSx);//
	VecDestroy(&lSy);//
	VecDestroy(&lSz);//
	VecDestroy(&lS);//
	/*
	DMDAVecRestoreArray(da, user->lIAj, &iaj);  
	DMDAVecRestoreArray(da, user->lJAj, &jaj);  
	DMDAVecRestoreArray(da, user->lKAj, &kaj);  
	
	DMDAVecRestoreArray(fda, user->lICsi, &icsi);
	DMDAVecRestoreArray(fda, user->lIEta, &ieta);
	DMDAVecRestoreArray(fda, user->lIZet, &izet);

	DMDAVecRestoreArray(fda, user->lJCsi, &jcsi);
	DMDAVecRestoreArray(fda, user->lJEta, &jeta);
	DMDAVecRestoreArray(fda, user->lJZet, &jzet);

	DMDAVecRestoreArray(fda, user->lKCsi, &kcsi);
	DMDAVecRestoreArray(fda, user->lKEta, &keta);
	DMDAVecRestoreArray(fda, user->lKZet, &kzet);
	*/
	DMGlobalToLocalBegin(da, user->CS, INSERT_VALUES, user->lCS);
	DMGlobalToLocalEnd(da, user->CS, INSERT_VALUES, user->lCS);
	

	DMDAVecGetArray(da, user->lCS, &Cs);
/* 	for(k=zs; k<ze; k++) */
/* 	for(j=ys; j<ye; j++) */
/* 	for(i=xs; i<xe; i++) { */
/* 		int flag=0, a=i, b=j, c=k; */
		
/* 		if(i_periodic && i==0) a=mx-2, flag=1; */
/* 		else if(i_periodic && i==mx-1) a=1, flag=1; */
		
/* 		if(j_periodic && j==0) b=my-2, flag=1; */
/* 		else if(j_periodic && j==my-1) b=1, flag=1; */
		
/* 		if(k_periodic && k==0) c=mz-2, flag=1; */
/* 		else if(k_periodic && k==mz-1) c=1, flag=1; */
		
		
/* 		/\* if(ii_periodic && i==0) a=-2, flag=1; *\/ */
/* /\* 		else if(ii_periodic && i==mx-1) a=mx+1, flag=1; *\/ */
		
/* /\* 		if(jj_periodic && j==0) b=-2, flag=1; *\/ */
/* /\* 		else if(jj_periodic && j==my-1) b=my+1, flag=1; *\/ */
		
/* /\* 		if(kk_periodic && k==0) c=-2, flag=1; *\/ */
/* /\* 		else if(kk_periodic && k==mz-1) c=mz+1, flag=1; *\/ */


/* 		if(flag) Cs[k][j][i] = Cs[c][b][a]; */
/* 	} */

	if (user->bctype[0]==7 || user->bctype[1]==7){
	  if (xs==0){
	    i=xs;
	    for (k=lzs; k<lze; k++) {
	      for (j=lys; j<lye; j++) {
		Cs[k][j][i]=Cs[k][j][i-2];
	
	      }
	    }
	  }
	  if (xe==mx){
	    i=mx-1;
	    for (k=lzs; k<lze; k++) {
	      for (j=lys; j<lye; j++) {
		Cs[k][j][i]=Cs[k][j][i+2];

	      }
	    }
	  }
	}
	if (user->bctype[2]==7 || user->bctype[3]==7){
	  if (ys==0){
	    j=ys;
	    for (k=lzs; k<lze; k++) {
	      for (i=lxs; i<lxe; i++) {
		Cs[k][j][i]=Cs[k][j-2][i];

	      }
	    }
	  }
	  if (ye==my){
	    j=my-1;
	    for (k=lzs; k<lze; k++) {
	      for (i=lxs; i<lxe; i++) {
		Cs[k][j][i]=Cs[k][j+2][i];

	      }
	    }
	  }
	}
	if (user->bctype[4]==7 || user->bctype[5]==7){
	  if (zs==0){
	    k=zs;
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		Cs[k][j][i]=Cs[k-2][j][i];
	      }
	    }
	  }
	  if (ze==mz){
	    k=mz-1;
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		Cs[k][j][i]=Cs[k+2][j][i];
	      }
	    }
	  }
	}

	DMDAVecRestoreArray(da, user->lCS, &Cs);
	
	if ( ! i_homo_filter && ! j_homo_filter && ! k_homo_filter && les!=3 ) {
		Vec lCs_tmp;
		PetscReal ***Cs_tmp;
		
		VecDuplicate(user->lCS, &lCs_tmp);
		VecCopy(user->lCS, lCs_tmp);
		
		DMDAVecGetArray(da, user->lAj, &aj);
		DMDAVecGetArray(da, user->lNvert, &nvert);
		DMDAVecGetArray(da, user->CS, &Cs);
		DMDAVecGetArray(da, lCs_tmp, &Cs_tmp);
		for(k=lzs; k<lze; k++)
		for(j=lys; j<lye; j++)
		for(i=lxs; i<lxe; i++) {
			if(nvert[k][j][i]<0.1) {
				double weight[3][3][3], C[3][3][3];
				int a, b, c;

				for(a=-1; a<=1; a++)
				for(b=-1; b<=1; b++)
				for(c=-1; c<=1; c++) {
					int R=c+1, Q=b+1, P=a+1;
					int KK=k+c, JJ=j+b, II=i+a;
					
					weight[R][Q][P] = 1./aj[KK][JJ][II];
					
/* 					if( nvert[KK][JJ][II]>1.1 ) weight[R][Q][P]=0; */
					
/* 					if( i_periodic ) { */
/* 						if( II==0 ) II=mx-2; */
/* 						else if( II==mx-1 ) II=1; */
/* 					} */
/* 					/\* else if( ii_periodic ) { *\/ */
/* /\* 						if( II==0 ) II=-2; *\/ */
/* /\* 						else if( II==mx-1 ) II=mx+1; *\/ */
/* /\* 					} *\/ */
/* 					else if( II==0 || II==mx-1) weight[R][Q][P]=0; */
					
/* 					if( j_periodic ) { */
/* 						if( JJ==0 ) JJ=my-2; */
/* 						else if( JJ==my-1 ) JJ=1; */
/* 					} */
/* 					/\* else if( jj_periodic ) { *\/ */
/* /\* 						if( JJ==0 ) JJ=-2; *\/ */
/* /\* 						else if( JJ==my-1 ) JJ=my+1; *\/ */
/* /\* 					} *\/ */
/* 					else if( JJ==0 || j==my-1) weight[R][Q][P]=0; */
					
/* 					if( k_periodic ) { */
/* 						if( KK==0 ) KK=mz-2; */
/* 						else if( KK==mz-1 ) KK=1; */
/* 					} */
/* 				/\* 	else if( kk_periodic ) { *\/ */
/* /\* 						if( KK==0 ) KK=-2; *\/ */
/* /\* 						else if( KK==mz-1 ) KK=mz+1; *\/ */
/* /\* 					} *\/ */
/* 					else if( KK==0 || KK==mz-1) weight[R][Q][P]=0; */
					
					C[R][Q][P] = Cs_tmp[KK][JJ][II];
				}
				Cs[k][j][i] = integrate_testfilter(C,weight);
			}
		}
		DMDAVecRestoreArray(da, user->lAj, &aj);
		DMDAVecRestoreArray(da, user->lNvert, &nvert);
		DMDAVecRestoreArray(da, user->CS, &Cs);
		DMDAVecRestoreArray(da, lCs_tmp, &Cs_tmp);
		
		VecDestroy(&lCs_tmp);
	 }   
         	
	DMGlobalToLocalBegin(da, user->CS, INSERT_VALUES, user->lCS);
	DMGlobalToLocalEnd(da, user->CS, INSERT_VALUES, user->lCS);

	
	DMDAVecGetArray(da, user->lCS, &Cs);
		if (user->bctype[0]==7 || user->bctype[1]==7){
	  if (xs==0){
	    i=xs;
	    for (k=lzs; k<lze; k++) {
	      for (j=lys; j<lye; j++) {
		Cs[k][j][i]=Cs[k][j][i-2];
		beta[k][j][i]=beta[k][j][i-2];
		
	      }
	    }
	  }
	  if (xe==mx){
	    i=mx-1;
	    for (k=lzs; k<lze; k++) {
	      for (j=lys; j<lye; j++) {
		Cs[k][j][i]=Cs[k][j][i+2];
		beta[k][j][i]=beta[k][j][i+2];


	      }
	    }
	  }
	}
	if (user->bctype[2]==7 || user->bctype[3]==7){
	  if (ys==0){
	    j=ys;
	    for (k=lzs; k<lze; k++) {
	      for (i=lxs; i<lxe; i++) {
		Cs[k][j][i]=Cs[k][j-2][i];
		beta[k][j][i]=beta[k][j-2][i];


	      }
	    }
	  }
	  if (ye==my){
	    j=my-1;
	    for (k=lzs; k<lze; k++) {
	      for (i=lxs; i<lxe; i++) {
		Cs[k][j][i]=Cs[k][j+2][i];
		beta[k][j][i]=beta[k][j+2][i];


	      }
	    }
	  }
	}
	if (user->bctype[4]==7 || user->bctype[5]==7){
	  if (zs==0){
	    k=zs;
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		Cs[k][j][i]=Cs[k-2][j][i];
		beta[k][j][i]=beta[k-2][j][i];
	
	      }
	    }
	  }
	  if (ze==mz){
	    k=mz-1;
	    for (j=lys; j<lye; j++) {
	      for (i=lxs; i<lxe; i++) {
		Cs[k][j][i]=Cs[k+2][j][i];
		beta[k][j][i]=beta[k+2][j][i];
	
	      }
	    }
	  }
	}
	
	
	DMDAVecRestoreArray(fda, user->lBeta, &beta); 
	
	DMDAVecRestoreArray(da, user->lCS, &Cs);
	//	DMDAVecRestoreArray(da, user->lP, &lp);
	PetscReal lmax_norm=0, max_norm;
	int p;

	if(i_homo_filter || j_homo_filter || k_homo_filter) PetscPrintf(PETSC_COMM_WORLD, "\nFilter type : Box filter homogeneous\n");
	else PetscPrintf(PETSC_COMM_WORLD, "Filter type : Box filter 3D\n");
	if(clark) PetscPrintf(PETSC_COMM_WORLD, "Clark model\n");
	
	VecMax(user->CS, &p, &lmax_norm);
	MPI_Allreduce(&lmax_norm,&max_norm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
	//	PetscGlobalMax(&lmax_norm, &max_norm, PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "Max Cs  = %e \n", max_norm);
	
	VecMax(user->lBeta, &p, &lmax_norm);
	MPI_Allreduce(&lmax_norm,&max_norm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
	//	PetscGlobalMax(&lmax_norm, &max_norm, PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "Max beta  = %e \n", max_norm);
	


}




void Compute_eddy_viscosity_LES(UserCtx *user)
{	
  DM		da = user->da, fda = user->fda, cda = user->cda;
  DMDALocalInfo	info=user->info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k;
	
  PetscReal    ***nvert, ***aj,***nu_t,***Cs, ***lt,***lp,***bbeta;//***lnu_t,
  Cmpnts      ***csi, ***eta, ***zet, ***ucat, ***beta;
  CompVars    ***q;
  double      dT_dx,dT_dy,dT_dz,dT_x,dT_y,dT_z;	
  //	DMDAGetLocalInfo(da, &info);
  mx = info.mx, my = info.my, mz = info.mz;
  xs = info.xs, xe = xs + info.xm;
  ys = info.ys, ye = ys + info.ym;
  zs = info.zs, ze = zs + info.zm;

  int lxs = xs, lxe = xe;
  int lys = ys, lye = ye;
  int lzs = zs, lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

	
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
 
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->lAj, &aj);
  //	DMDAVecGetArray(da, user->lNu_t, &lnu_t);	
  DMDAVecGetArray(da, user->Nu_t, &nu_t);

  DMDAVecGetArray(da, user->CS, &Cs);

  DMDAVecGetArray(fda, user->lUcat,  &ucat);
  DMDAVecGetArray(cda, user->lQ, &q);
	
  //

  Vec lT;
  VecDuplicate(user->lP,&lT);
  DMDAVecGetArray(da, lT, &lt);

  int tt;
  tt = rand()%100;	

  DMDAVecGetArray(da, user->lP, &lp);

  double M=user->Ma;
  double g=user->gamma;
	
  for (k=zs; k<ze; k++)
    for (j=ys; j<ye; j++)
      for (i=xs; i<xe; i++) {
	
	lt[k][j][i]= g*M*M*lp[k][j][i]/q[k][j][i].rho;

      }

  DMDAVecRestoreArray(da, user->lP, &lp);

  DMDAVecGetArray(fda, user->lBeta, &beta); 

  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {
	if(nvert[k][j][i]>0.1) {
	  nu_t[k][j][i]=0;
	  continue;
	}
	double ajc = aj[k][j][i];
	double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
	double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
	double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
	double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,dpdc,dpde,dpdz;
	double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz,dp_dc,dp_de,dp_dz;
	Compute_du_center (i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, 
			   &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);
	Compute_du_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, 
			 ajc, dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, 
			 dwdz,&du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, 
			 &du_dz, &dv_dz, &dw_dz );
	double Sdd=2./3.*(du_dx+dv_dy+dw_dz);	
	double Sxx = 0.5*( du_dx + du_dx )-Sdd, Sxy = 0.5*(du_dy + dv_dx), Sxz = 0.5*(du_dz + dw_dx);
	double Syx = Sxy, Syy = 0.5*(dv_dy + dv_dy)-Sdd,	Syz = 0.5*(dv_dz + dw_dy);
	double Szx = Sxz, Szy=Syz, Szz = 0.5*(dw_dz + dw_dz)-Sdd;
	
	double Sabs = sqrt( 2.0*( Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Syx*Syx + Syy*Syy + Syz*Syz + Szx*Szx + Szy*Szy + Szz*Szz ) );

	double filter  = pow( 1./aj[k][j][i],1./3.);
	nu_t[k][j][i] = Cs[k][j][i] * pow ( filter, 2.0 ) * Sabs;


	//////////////////////// Energy beta Vreman 1995

	Compute_dp_center (i, j, k, mx, my, mz, lt, nvert, &dpdc, &dpde, &dpdz);
	Compute_dp_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, 
			 ajc, dpdc,dpde,dpdz, &dT_dx, &dT_dy, &dT_dz);
	beta[k][j][i].x = -1.0 * dT_x * q[k][j][i].rho*nu_t[k][j][i]/((g-1.0)*M*M*0.71);	
	beta[k][j][i].y =- 1.0 * dT_dy *q[k][j][i].rho * nu_t[k][j][i]/((g-1.0)*M*M*0.71);		
	beta[k][j][i].z = -1.0 * dT_dz *q[k][j][i].rho * nu_t[k][j][i]*((g-1.0)*M*M*0.71);
               		

      }
  //	DMDAVecRestoreArray(fda, user->lBeta, &beta); 


  //	DMDAVecGetArray(fda, user->lBeta, &beta);
  DMDAVecRestoreArray(fda, user->lUcat, &ucat); 
  DMDAVecRestoreArray(cda, user->lQ, &q);

  VecSet(user->lBBeta, 0.0);
  DMDAVecGetArray(da, user->lBBeta, &bbeta);
		
  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {
	if(nvert[k][j][i]>0.1) {
	  nu_t[k][j][i]=0;
	  bbeta[k][j][i]=0.0;
	  continue;
	}
	double ajc = aj[k][j][i];
	double csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
	double eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
	double zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
	double dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
	double du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
	Compute_du_center (i, j, k, mx, my, mz, beta, nvert, &dudc, &dvdc, 
			   &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz); //should this be beta or bbeta?
	Compute_du_dxyz (csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc, 
			 dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz,
			 &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz );

	bbeta[k][j][i]=0.0;
	if ((nvert[k][j+1][i]+nvert[k][j-1][i]+nvert[k+1][j][i]+nvert[k-1][j][i]+ 
	     nvert[k][j][i+1]+nvert[k][j-1][i-1])<0.1) {
	  bbeta[k][j][i]=.5*(du_dx+dv_dy+dw_dz)/ajc;
	  //	bbeta[k][j][i]=0.0;
	  //PetscPrintf(PETSC_COMM_WORLD, " i  %d    j  %d     k  %d      Max beta = %le   beta  %le     dw_dz  %le \n",i,j,k, bbeta[k][j][i],beta[k][j][i].z,dw_dz);

	}
      }

  // 		  PetscPrintf(PETSC_COMM_SELF," LES\n");
  DMDAVecRestoreArray(da, lT, &lt);
  VecDestroy(&lT);//

  DMDAVecRestoreArray(fda, user->lBeta, &beta); 
  DMDAVecRestoreArray(da, user->lBBeta, &bbeta);
 	
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
		
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(da, user->Nu_t, &nu_t);
  //	DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
	
  DMDAVecRestoreArray(da, user->CS, &Cs);
	
  DMGlobalToLocalBegin(da, user->Nu_t, INSERT_VALUES, user->lNu_t);
  DMGlobalToLocalEnd(da, user->Nu_t, INSERT_VALUES, user->lNu_t);

  PetscReal lmax_norm=0.0, lbetamax=0.0, betamax, max_norm;
  int p,p1;
	
  VecMax(user->lNu_t, &p, &lmax_norm);
  VecMax(user->lBBeta, &p1, &lbetamax);
	
  MPI_Allreduce(&lbetamax,&betamax,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
	
  MPI_Allreduce(&lmax_norm,&max_norm,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
  //	PetscGlobalMax(&lmax_norm, &max_norm, PETSC_COMM_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "Max Nu_t = %e    betamax =%e\n", max_norm,betamax);
}


PetscErrorCode Do_averaging(UserCtx *user)
{
	DMDALocalInfo info = user->info;
	PetscInt	xs = info.xs, xe = info.xs + info.xm;
	PetscInt  	ys = info.ys, ye = info.ys + info.ym;
	PetscInt	zs = info.zs, ze = info.zs + info.zm;
	PetscInt	mx = info.mx, my = info.my, mz = info.mz;

	PetscInt i, j, k;
  
	PetscInt	lxs, lxe, lys, lye, lzs, lze;

	lxs = xs; lxe = xe;
	lys = ys; lye = ye;
	lzs = zs; lze = ze;

	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;

	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;
  
	Cmpnts ***ucat, ***u_cross_sum;
	Cmpnts ***u2sum;
	CompVars ***q;
	PetscInt period=1,travel=0, fish=0;
	  
	//	PetscOptionsInsertFile(NULL, PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
	PetscOptionsGetInt(NULL, PETSC_NULL, "-travel", &travel, PETSC_NULL);
	
	PetscOptionsGetInt(NULL, PETSC_NULL, "-fish", &fish, PETSC_NULL);
 
	if (fish)	  
	PetscOptionsGetInt(NULL, PETSC_NULL, "-N_period_fish", &period, PETSC_NULL);
 
	
	if ( user->ti%period==0)
	  {
	    VecAXPY(user->Q_sum, 1., user->Q);
	    //VecAXPY(user->lUstar_sum, 1., user->lUstar);
	    
	    VecAXPY(user->P_sum, 1., user->P);
	    double max_norm;
	    VecMax(user->Q_sum, &i, &max_norm);
	    PetscPrintf(PETSC_COMM_WORLD, "\n*** Max Qcat_avg = %e \n", max_norm / user->ti);
	}	    
	    //user->utawsum+=user->utaw;
	    //user->utaw1sum+=user->utaw1;
	    //PetscPrintf(PETSC_COMM_WORLD, "utawsum = %le  utaw1sum %le  %d\n",user->utawsum,user->utaw1sum,ti );

	    DMDAVecGetArray(user->fda, user->Ucat, &ucat);
	    //DMDAVecGetArray(user->fda, user->Ucat_sum, &ucat_sum);
	    DMDAVecGetArray(user->fda, user->Ucat_square_sum, &u2sum);
	    DMDAVecGetArray(user->fda, user->Ucat_cross_sum, &u_cross_sum);
	    //DMDAVecGetArray(user->da, user->P_sum, &p_sum);
	    //DMDAVecGetArray(user->da, user->P, &p);
	
	    for (k=zs; k<ze; k++)
	      for (j=ys; j<ye; j++)
		for (i=xs; i<xe; i++) {
		  
//		    ucat_sum[k][j][i].x += ucat[k][j][i].x;
//		    ucat_sum[k][j][i].y += ucat[k][j][i].y;
//		    ucat_sum[k][j][i].z += ucat[k][j][i].z;
//		    p_sum[k][j][i] += p[k][j][i];
		
		  u2sum[k][j][i].x += ucat[k][j][i].x * ucat[k][j][i].x;
		  u2sum[k][j][i].y += ucat[k][j][i].y * ucat[k][j][i].y;
		  u2sum[k][j][i].z += ucat[k][j][i].z * ucat[k][j][i].z;
		  
			u_cross_sum[k][j][i].x += ucat[k][j][i].x * ucat[k][j][i].y;	// uv
			u_cross_sum[k][j][i].y += ucat[k][j][i].y * ucat[k][j][i].z;	// vw
			u_cross_sum[k][j][i].z += ucat[k][j][i].z * ucat[k][j][i].x;	// wu
		}
	
	    DMDAVecRestoreArray(user->fda, user->Ucat, &ucat);
	//DMDAVecRestoreArray(user->fda, user->Ucat_sum, &ucat_sum);
	    DMDAVecRestoreArray(user->fda, user->Ucat_square_sum, &u2sum);
	    DMDAVecRestoreArray(user->fda, user->Ucat_cross_sum, &u_cross_sum);
  
	return 0;
}
