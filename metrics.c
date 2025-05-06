#include "variables.h"

PetscErrorCode FormMetrics(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  Cmpnts	***csi, ***eta, ***zet;
  PetscScalar	***aj;
  Vec		coords;
  Cmpnts	***coor;

  PetscInt	xs, ys, zs, xe, ye, ze,rank;
  DMDALocalInfo	info;
  PetscInt	mx, my, mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
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

  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;

  /* PetscInt	gxs, gxe, gys, gye, gzs, gze; */
  /* gxs = info.gxs; gxe = gxs + info.gxm; */
  /* gys = info.gys; gye = gys + info.gym; */
  /* gzs = info.gzs; gze = gzs + info.gzm; */

  PetscScalar	dxdc, dydc, dzdc, dxde, dyde, dzde, dxdz, dydz, dzdz;
  PetscInt	i, j, k;

  DMGetCoordinatesLocal(da, &coords);
  DMDAVecGetArray(fda, coords, &coor);

  VecSet(user->Csi,0.0);
  VecSet(user->Eta,0.0);
  VecSet(user->Zet,0.0);
  DMDAVecGetArray(fda, user->Csi, &csi);
  DMDAVecGetArray(fda, user->Eta, &eta);
  DMDAVecGetArray(fda, user->Zet, &zet);
  DMDAVecGetArray(da, user->Aj,  &aj);

 /* Calculating transformation metrics in i direction at cell center */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=xs; i<xe; i++){
	/* csi = de X dz */
	dxde = 0.5 * (coor[k  ][j+1][i  ].x - coor[k  ][j-1][i  ].x);
	dyde = 0.5 * (coor[k  ][j+1][i  ].y - coor[k  ][j-1][i  ].y); 
	dzde = 0.5 * (coor[k  ][j+1][i  ].z - coor[k  ][j-1][i  ].z);
				       		    	    	 
	dxdz = 0.5 * (coor[k+1][j  ][i  ].x - coor[k-1][j  ][i  ].x);
	dydz = 0.5 * (coor[k+1][j  ][i  ].y - coor[k-1][j  ][i  ].y);
	dzdz = 0.5 * (coor[k+1][j  ][i  ].z - coor[k-1][j  ][i  ].z);
	  
	csi[k][j][i].x = dyde * dzdz - dzde * dydz;
	csi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
	csi[k][j][i].z = dxde * dydz - dyde * dxdz;		
      }
    }
  }

  /* calculating j direction metrics at cell center */
  for (k=lzs; k<lze; k++){
    for (j=ys; j< ye; j++){
      for (i=lxs; i<lxe; i++){

	/* eta = dz X de */
	dxdc = 0.5 * (coor[k  ][j  ][i+1].x - coor[k  ][j  ][i-1].x);
	dydc = 0.5 * (coor[k  ][j  ][i+1].y - coor[k  ][j  ][i-1].y);
	dzdc = 0.5 * (coor[k  ][j  ][i+1].z - coor[k  ][j  ][i-1].z);
	
	dxdz = 0.5 * (coor[k+1][j  ][i  ].x - coor[k-1][j  ][i  ].x);
	dydz = 0.5 * (coor[k+1][j  ][i  ].y - coor[k-1][j  ][i  ].y);
	dzdz = 0.5 * (coor[k+1][j  ][i  ].z - coor[k-1][j  ][i  ].z);    		         	 		   	 
	  
	eta[k][j][i].x = dydz * dzdc - dzdz * dydc;
	eta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
	eta[k][j][i].z = dxdz * dydc - dydz * dxdc;

      }
    }
  }
  
  /* calculating k direction metrics at cell center */
  for (k=zs; k<ze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	dxdc = 0.5 * (coor[k  ][j  ][i+1].x - coor[k  ][j  ][i-1].x);
	dydc = 0.5 * (coor[k  ][j  ][i+1].y - coor[k  ][j  ][i-1].y);
	dzdc = 0.5 * (coor[k  ][j  ][i+1].z - coor[k  ][j  ][i-1].z);
			    		    	     	     	 
	dxde = 0.5 * (coor[k  ][j+1][i  ].x - coor[k  ][j-1][i  ].x);
	dyde = 0.5 * (coor[k  ][j+1][i  ].y - coor[k  ][j-1][i  ].y); 
	dzde = 0.5 * (coor[k  ][j+1][i  ].z - coor[k  ][j-1][i  ].z);
	  
	zet[k][j][i].x = dydc * dzde - dzdc * dyde;
	zet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
	zet[k][j][i].z = dxdc * dyde - dydc * dxde;
	
      }
    }
  }
 
  /* calculating Jacobian of the transformation at cell center */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	dxdc = 0.5 * (coor[k  ][j  ][i+1].x - coor[k  ][j  ][i-1].x);
	dydc = 0.5 * (coor[k  ][j  ][i+1].y - coor[k  ][j  ][i-1].y);
	dzdc = 0.5 * (coor[k  ][j  ][i+1].z - coor[k  ][j  ][i-1].z);
			    		    	     	     	 
	dxde = 0.5 * (coor[k  ][j+1][i  ].x - coor[k  ][j-1][i  ].x);
	dyde = 0.5 * (coor[k  ][j+1][i  ].y - coor[k  ][j-1][i  ].y); 
	dzde = 0.5 * (coor[k  ][j+1][i  ].z - coor[k  ][j-1][i  ].z);

	dxdz = 0.5 * (coor[k+1][j  ][i  ].x - coor[k-1][j  ][i  ].x);
	dydz = 0.5 * (coor[k+1][j  ][i  ].y - coor[k-1][j  ][i  ].y);
	dzdz = 0.5 * (coor[k+1][j  ][i  ].z - coor[k-1][j  ][i  ].z);  
	  
	aj[k][j][i] = dxdc * (dyde * dzdz - dzde * dydz) -
	  dydc * (dxde * dzdz - dzde * dxdz) +
	  dzdc * (dxde * dydz - dyde * dxdz);
	aj[k][j][i] = 1./aj[k][j][i];
	
      }
    }
  }

  if (xs==0) {
    i = xs;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	eta[k][j][i].x = eta[k][j][i+1].x;
	eta[k][j][i].y = eta[k][j][i+1].y;
	eta[k][j][i].z = eta[k][j][i+1].z;
	
	zet[k][j][i].x = zet[k][j][i+1].x;
	zet[k][j][i].y = zet[k][j][i+1].y;
	zet[k][j][i].z = zet[k][j][i+1].z;
	
	aj[k][j][i] = aj[k][j][i+1];
	
      }
    }
  }
  
  
  if (xe==mx) {
    i = xe-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	
	eta[k][j][i].x = eta[k][j][i-1].x;
	eta[k][j][i].y = eta[k][j][i-1].y;
	eta[k][j][i].z = eta[k][j][i-1].z;
	  
	zet[k][j][i].x = zet[k][j][i-1].x;
	zet[k][j][i].y = zet[k][j][i-1].y;
	zet[k][j][i].z = zet[k][j][i-1].z;
	  
	aj[k][j][i] = aj[k][j][i-1];
	
      }
    }
  }
 

  if (ys==0) {
    j = ys;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	csi[k][j][i].x = csi[k][j+1][i].x;
	csi[k][j][i].y = csi[k][j+1][i].y;
	csi[k][j][i].z = csi[k][j+1][i].z;
	
	zet[k][j][i].x = zet[k][j+1][i].x;
	zet[k][j][i].y = zet[k][j+1][i].y;
	zet[k][j][i].z = zet[k][j+1][i].z;
	
	aj[k][j][i] = aj[k][j+1][i];
      }
    }
  }
  
  
  if (ye==my) {
    j = ye-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	csi[k][j][i].x = csi[k][j-1][i].x;
	csi[k][j][i].y = csi[k][j-1][i].y;
	csi[k][j][i].z = csi[k][j-1][i].z;
	
	zet[k][j][i].x = zet[k][j-1][i].x;
	zet[k][j][i].y = zet[k][j-1][i].y;
	zet[k][j][i].z = zet[k][j-1][i].z;
	
	aj[k][j][i] = aj[k][j-1][i];
      }
    }
  }
   

  if (zs==0) {
    k = zs;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	eta[k][j][i].x = eta[k+1][j][i].x;
	eta[k][j][i].y = eta[k+1][j][i].y;
	eta[k][j][i].z = eta[k+1][j][i].z;
	
	csi[k][j][i].x = csi[k+1][j][i].x;
	csi[k][j][i].y = csi[k+1][j][i].y;
	csi[k][j][i].z = csi[k+1][j][i].z;
	
	aj[k][j][i] = aj[k+1][j][i];
      }
    }
  }
  
 
  if (ze==mz){
    k = ze-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	eta[k][j][i].x = eta[k-1][j][i].x;
	eta[k][j][i].y = eta[k-1][j][i].y;
	eta[k][j][i].z = eta[k-1][j][i].z;
	
	csi[k][j][i].x = csi[k-1][j][i].x;
	csi[k][j][i].y = csi[k-1][j][i].y;
	csi[k][j][i].z = csi[k-1][j][i].z;
	
	aj[k][j][i] = aj[k-1][j][i];
      }
    }
  }

  /* PetscPrintf(PETSC_COMM_WORLD, "metrics csi...aj done!\n"); */
  /* /\* calculating the metrics on other faces/locations *\/ */
  /* Vec		ICsi = user->ICsi, IEta = user->IEta, IZet = user->IZet; */
  /* Vec		JCsi = user->JCsi, JEta = user->JEta, JZet = user->JZet; */
  /* Vec		KCsi = user->KCsi, KEta = user->KEta, KZet = user->KZet; */
  /* Vec		IAj = user->IAj, JAj = user->JAj, KAj = user->KAj; */
  /* Vec           Centx, Centy, Centz; */
  
  /* Cmpnts	***icsi, ***ieta, ***izet; */
  /* Cmpnts	***jcsi, ***jeta, ***jzet; */
  /* Cmpnts	***kcsi, ***keta, ***kzet; */
  /* Cmpnts	***gs, ***cent; */
  /* Cmpnts	***centx, ***centy, ***centz; */
  /* PetscReal	***iaj, ***jaj, ***kaj; */

  /* PetscReal	xcp, ycp, zcp, xcm, ycm, zcm; */
  
  /* DMDAVecGetArray(fda, ICsi, &icsi); */
  /* DMDAVecGetArray(fda, IEta, &ieta); */
  /* DMDAVecGetArray(fda, IZet, &izet); */
  /* DMDAVecGetArray(da, IAj,  &iaj); */

  /* DMDAVecGetArray(fda, user->Cent, &cent); */
  /* DMDAVecGetArray(fda, user->GridSpace, &gs); */

  /* for (k=lzs; k<lze; k++) { */
  /*   for (j=lys; j<lye; j++) { */
  /*     for (i=lxs; i<lxe; i++) { */
  /* 	cent[k][j][i].x = 0.125 * */
  /* 	  (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x + */
  /* 	   coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x + */
  /* 	   coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x + */
  /* 	   coor[k-1][j  ][i-1].x + coor[k-1][j-1][i-1].x); */
  /* 	cent[k][j][i].y = 0.125 * */
  /* 	  (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y + */
  /* 	   coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y + */
  /* 	   coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y + */
  /* 	   coor[k-1][j  ][i-1].y + coor[k-1][j-1][i-1].y); */
  /* 	cent[k][j][i].z = 0.125 * */
  /* 	  (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z + */
  /* 	   coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z + */
  /* 	   coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z + */
  /* 	   coor[k-1][j  ][i-1].z + coor[k-1][j-1][i-1].z); */

  /* 	xcp = 0.25 * */
  /* 	  (coor[k][j][i].x + coor[k][j-1][i].x + */
  /* 	   coor[k-1][j-1][i].x + coor[k-1][j][i].x); */
  /* 	ycp = 0.25 * */
  /* 	  (coor[k][j][i].y + coor[k][j-1][i].y + */
  /* 	   coor[k-1][j-1][i].y + coor[k-1][j][i].y); */
  /* 	zcp = 0.25 * */
  /* 	  (coor[k][j][i].z + coor[k][j-1][i].z + */
  /* 	   coor[k-1][j-1][i].z + coor[k-1][j][i].z); */

  /* 	xcm = 0.25 * */
  /* 	  (coor[k][j][i-1].x + coor[k][j-1][i-1].x + */
  /* 	   coor[k-1][j-1][i-1].x + coor[k-1][j][i-1].x); */
  /* 	ycm = 0.25 * */
  /* 	  (coor[k][j][i-1].y + coor[k][j-1][i-1].y + */
  /* 	   coor[k-1][j-1][i-1].y + coor[k-1][j][i-1].y); */
  /* 	zcm = 0.25 * */
  /* 	  (coor[k][j][i-1].z + coor[k][j-1][i-1].z + */
  /* 	   coor[k-1][j-1][i-1].z + coor[k-1][j][i-1].z); */

  /* 	gs[k][j][i].x = sqrt((xcp-xcm) * (xcp-xcm) + */
  /* 			     (ycp-ycm) * (ycp-ycm) + */
  /* 			     (zcp-zcm) * (zcp-zcm)); */

  /* 	xcp = 0.25 * */
  /* 	  (coor[k][j][i].x + coor[k][j][i-1].x + */
  /* 	   coor[k-1][j][i].x + coor[k-1][j][i-1].x); */
  /* 	ycp = 0.25 * */
  /* 	  (coor[k][j][i].y + coor[k][j][i-1].y + */
  /* 	   coor[k-1][j][i].y + coor[k-1][j][i-1].y); */
  /* 	zcp = 0.25 * */
  /* 	  (coor[k][j][i].z + coor[k][j][i-1].z + */
  /* 	   coor[k-1][j][i].z + coor[k-1][j][i-1].z); */

  /* 	xcm = 0.25 * */
  /* 	  (coor[k][j-1][i].x + coor[k][j-1][i-1].x + */
  /* 	   coor[k-1][j-1][i].x + coor[k-1][j-1][i-1].x); */
  /* 	ycm = 0.25 * */
  /* 	  (coor[k][j-1][i].y + coor[k][j-1][i-1].y + */
  /* 	   coor[k-1][j-1][i].y + coor[k-1][j-1][i-1].y); */
  /* 	zcm = 0.25 * */
  /* 	  (coor[k][j-1][i].z + coor[k][j-1][i-1].z + */
  /* 	   coor[k-1][j-1][i].z + coor[k-1][j-1][i-1].z); */

  /* 	gs[k][j][i].y = sqrt((xcp-xcm) * (xcp-xcm) + */
  /* 			     (ycp-ycm) * (ycp-ycm) + */
  /* 			     (zcp-zcm) * (zcp-zcm)); */

  /* 	xcp = 0.25 * */
  /* 	  (coor[k][j][i].x + coor[k][j][i-1].x + */
  /* 	   coor[k][j-1][i].x + coor[k][j-1][i-1].x); */
  /* 	ycp = 0.25 * */
  /* 	  (coor[k][j][i].y + coor[k][j][i-1].y + */
  /* 	   coor[k][j-1][i].y + coor[k][j-1][i-1].y); */
  /* 	zcp = 0.25 * */
  /* 	  (coor[k][j][i].z + coor[k][j][i-1].z + */
  /* 	   coor[k][j-1][i].z + coor[k][j-1][i-1].z); */

  /* 	xcm = 0.25 * */
  /* 	  (coor[k-1][j][i].x + coor[k-1][j][i-1].x + */
  /* 	   coor[k-1][j-1][i].x + coor[k-1][j-1][i-1].x); */
  /* 	ycm = 0.25 * */
  /* 	  (coor[k-1][j][i].y + coor[k-1][j][i-1].y + */
  /* 	   coor[k-1][j-1][i].y + coor[k-1][j-1][i-1].y); */
  /* 	zcm = 0.25 * */
  /* 	  (coor[k-1][j][i].z + coor[k-1][j][i-1].z + */
  /* 	   coor[k-1][j-1][i].z + coor[k-1][j-1][i-1].z); */

  /* 	gs[k][j][i].z = sqrt((xcp-xcm) * (xcp-xcm) + */
  /* 			     (ycp-ycm) * (ycp-ycm) + */
  /* 			     (zcp-zcm) * (zcp-zcm)); */

  /*     } */
  /*   } */
  /* } */

  /* DMDAVecRestoreArray(fda, user->GridSpace, &gs); */

  /* VecDuplicate(user->lCsi, &Centx); */
  /* DMDAVecGetArray(fda, Centx, &centx); */
  /* for(k=gzs+1; k<gze; k++) { */
  /*   for (j=gys+1; j<gye; j++) { */
  /*     for (i=gxs; i<gxe; i++) { */
  /* 	centx[k][j][i].x = (coor[k  ][j  ][i].x + coor[k-1][j  ][i].x + */
  /* 			    coor[k  ][j-1][i].x + coor[k-1][j-1][i].x) * 0.25; */
  /* 	centx[k][j][i].y = (coor[k  ][j  ][i].y + coor[k-1][j  ][i].y + */
  /* 			    coor[k  ][j-1][i].y + coor[k-1][j-1][i].y) * 0.25; */
  /* 	centx[k][j][i].z = (coor[k  ][j  ][i].z + coor[k-1][j  ][i].z + */
  /* 			    coor[k  ][j-1][i].z + coor[k-1][j-1][i].z) * 0.25; */

  /*     } */
  /*   } */
  /* } */

  /* for (k=lzs; k<lze; k++) { */
  /*   for (j=lys; j<lye; j++) { */
  /*     for (i=xs; i<lxe; i++) {   */
  /* 	if (i>0 && i<mx-2) { */
  /* 	  dxdc = (centx[k][j][i+1].x - centx[k][j][i-1].x) * 0.5; */
  /* 	  dydc = (centx[k][j][i+1].y - centx[k][j][i-1].y) * 0.5; */
  /* 	  dzdc = (centx[k][j][i+1].z - centx[k][j][i-1].z) * 0.5; */
  /* 	} */
  /* 	else if (i==0) { */
  /* 	  if (user->bctype[0]==7) { */
  /* 	  dxdc = centx[k][j][i+1].x - centx[k][j][mx-2].x; */
  /* 	  dydc = centx[k][j][i+1].y - centx[k][j][mx-2].y; */
  /* 	  dzdc = centx[k][j][i+1].z - centx[k][j][mx-2].z; */
  /* 	  } else { */
  /* 	  dxdc = centx[k][j][i+1].x - centx[k][j][i].x; */
  /* 	  dydc = centx[k][j][i+1].y - centx[k][j][i].y; */
  /* 	  dzdc = centx[k][j][i+1].z - centx[k][j][i].z; */
  /* 	  } */
  /* 	} */
  /* 	else if (i==mx-2) { */
  /* 	  if (user->bctype[0]==7) { */
  /* 	  dxdc = centx[k][j][i].x - centx[k][j][1].x; */
  /* 	  dydc = centx[k][j][i].y - centx[k][j][1].y; */
  /* 	  dzdc = centx[k][j][i].z - centx[k][j][1].z; */
  /* 	  } else { */
  /* 	  dxdc = centx[k][j][i].x - centx[k][j][i-1].x; */
  /* 	  dydc = centx[k][j][i].y - centx[k][j][i-1].y; */
  /* 	  dzdc = centx[k][j][i].z - centx[k][j][i-1].z; */
  /* 	  } */
  /* 	} */

  /* 	if (j>1 && j<my-2) { */
  /* 	  dxde = (centx[k][j+1][i].x - centx[k][j-1][i].x) * 0.5; */
  /* 	  dyde = (centx[k][j+1][i].y - centx[k][j-1][i].y) * 0.5; */
  /* 	  dzde = (centx[k][j+1][i].z - centx[k][j-1][i].z) * 0.5; */
  /* 	} */
  /* 	else if (j==1) { */
  /* 	  dxde = centx[k][j+1][i].x - centx[k][j][i].x; */
  /* 	  dyde = centx[k][j+1][i].y - centx[k][j][i].y; */
  /* 	  dzde = centx[k][j+1][i].z - centx[k][j][i].z; */
  /* 	} */
  /* 	else if (j==my-2) { */
  /* 	  dxde = centx[k][j][i].x - centx[k][j-1][i].x; */
  /* 	  dyde = centx[k][j][i].y - centx[k][j-1][i].y; */
  /* 	  dzde = centx[k][j][i].z - centx[k][j-1][i].z; */
  /* 	} */

  /* 	if (k>1 && k<mz-2) { */
  /* 	  dxdz = (centx[k+1][j][i].x - centx[k-1][j][i].x) * 0.5; */
  /* 	  dydz = (centx[k+1][j][i].y - centx[k-1][j][i].y) * 0.5; */
  /* 	  dzdz = (centx[k+1][j][i].z - centx[k-1][j][i].z) * 0.5; */
  /* 	} */
  /* 	else if (k==1) { */
  /* 	  dxdz = (centx[k+1][j][i].x - centx[k][j][i].x); */
  /* 	  dydz = (centx[k+1][j][i].y - centx[k][j][i].y); */
  /* 	  dzdz = (centx[k+1][j][i].z - centx[k][j][i].z); */
  /* 	} */
  /* 	else if (k==mz-2) { */
  /* 	  dxdz = (centx[k][j][i].x - centx[k-1][j][i].x); */
  /* 	  dydz = (centx[k][j][i].y - centx[k-1][j][i].y); */
  /* 	  dzdz = (centx[k][j][i].z - centx[k-1][j][i].z); */
  /* 	} */

  /* 	icsi[k][j][i].x = dyde * dzdz - dzde * dydz; */
  /* 	icsi[k][j][i].y =-dxde * dzdz + dzde * dxdz; */
  /* 	icsi[k][j][i].z = dxde * dydz - dyde * dxdz; */

  /* 	ieta[k][j][i].x = dydz * dzdc - dzdz * dydc; */
  /* 	ieta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc; */
  /* 	ieta[k][j][i].z = dxdz * dydc - dydz * dxdc; */

  /* 	izet[k][j][i].x = dydc * dzde - dzdc * dyde; */
  /* 	izet[k][j][i].y =-dxdc * dzde + dzdc * dxde; */
  /* 	izet[k][j][i].z = dxdc * dyde - dydc * dxde; */

  /* 	iaj[k][j][i] = dxdc * (dyde * dzdz - dzde * dydz) - */
  /* 	  dydc * (dxde * dzdz - dzde * dxdz) + */
  /* 	  dzdc * (dxde * dydz - dyde * dxdz); */
  /* 	iaj[k][j][i] = 1./iaj[k][j][i]; */


  /*     } */
  /*   } */
  /* } */

  /* DMDAVecRestoreArray(fda, Centx, &centx); */
  /* DMDAVecRestoreArray(fda, ICsi, &icsi); */
  /* DMDAVecRestoreArray(fda, IEta, &ieta); */
  /* DMDAVecRestoreArray(fda, IZet, &izet); */
  /* DMDAVecRestoreArray(da, IAj,  &iaj); */
  /* VecDestroy(&Centx); */
  
  /* PetscPrintf(PETSC_COMM_WORLD, "metrics i faces done!\n"); */

  /*   // j direction */
  /* DMDAVecGetArray(fda, JCsi, &jcsi); */
  /* DMDAVecGetArray(fda, JEta, &jeta); */
  /* DMDAVecGetArray(fda, JZet, &jzet); */
  /* DMDAVecGetArray(da, JAj,  &jaj); */

  /* VecDuplicate(user->lCsi, &Centy); */
  /* DMDAVecGetArray(fda, Centy, &centy); */
  /* for(k=gzs+1; k<gze; k++) { */
  /*   for (j=gys; j<gye; j++) { */
  /*     for (i=gxs+1; i<gxe; i++) { */
  /* 	centy[k][j][i].x = (coor[k  ][j][i  ].x + coor[k-1][j][i  ].x + */
  /* 			    coor[k  ][j][i-1].x + coor[k-1][j][i-1].x) * 0.25; */
  /* 	centy[k][j][i].y = (coor[k  ][j][i  ].y + coor[k-1][j][i  ].y + */
  /* 			    coor[k  ][j][i-1].y + coor[k-1][j][i-1].y) * 0.25; */
  /* 	centy[k][j][i].z = (coor[k  ][j][i  ].z + coor[k-1][j][i  ].z + */
  /* 			    coor[k  ][j][i-1].z + coor[k-1][j][i-1].z) * 0.25; */
  /*     } */
  /*   } */
  /* } */

  /* for (k=lzs; k<lze; k++) { */
  /*   for (j=ys; j<lye; j++) { */
  /*     for (i=lxs; i<lxe; i++) { */
  /* 	if (i>1 && i<mx-2) { */
  /* 	  dxdc = (centy[k][j][i+1].x - centy[k][j][i-1].x) * 0.5; */
  /* 	  dydc = (centy[k][j][i+1].y - centy[k][j][i-1].y) * 0.5; */
  /* 	  dzdc = (centy[k][j][i+1].z - centy[k][j][i-1].z) * 0.5; */
  /* 	} */
  /* 	else if (i==1) { */
  /* 	  dxdc = centy[k][j][i+1].x - centy[k][j][i].x; */
  /* 	  dydc = centy[k][j][i+1].y - centy[k][j][i].y; */
  /* 	  dzdc = centy[k][j][i+1].z - centy[k][j][i].z; */
  /* 	} */
  /* 	else if (i==mx-2) { */
  /* 	  dxdc = centy[k][j][i].x - centy[k][j][i-1].x; */
  /* 	  dydc = centy[k][j][i].y - centy[k][j][i-1].y; */
  /* 	  dzdc = centy[k][j][i].z - centy[k][j][i-1].z; */
  /* 	} */

  /* 	if (j>0 && j<my-2) { */
  /* 	  dxde = (centy[k][j+1][i].x - centy[k][j-1][i].x) * 0.5; */
  /* 	  dyde = (centy[k][j+1][i].y - centy[k][j-1][i].y) * 0.5; */
  /* 	  dzde = (centy[k][j+1][i].z - centy[k][j-1][i].z) * 0.5; */
  /* 	} */
  /* 	else if (j==0) { */
  /* 	  dxde = centy[k][j+1][i].x - centy[k][j][i].x; */
  /* 	  dyde = centy[k][j+1][i].y - centy[k][j][i].y; */
  /* 	  dzde = centy[k][j+1][i].z - centy[k][j][i].z; */
  /* 	} */
  /* 	else if (j==my-2) { */
  /* 	  dxde = centy[k][j][i].x - centy[k][j-1][i].x; */
  /* 	  dyde = centy[k][j][i].y - centy[k][j-1][i].y; */
  /* 	  dzde = centy[k][j][i].z - centy[k][j-1][i].z; */
  /* 	} */

  /* 	if (k>1 && k<mz-2) { */
  /* 	  dxdz = (centy[k+1][j][i].x - centy[k-1][j][i].x) * 0.5; */
  /* 	  dydz = (centy[k+1][j][i].y - centy[k-1][j][i].y) * 0.5; */
  /* 	  dzdz = (centy[k+1][j][i].z - centy[k-1][j][i].z) * 0.5; */
  /* 	} */
  /* 	else if (k==1) { */
  /* 	  dxdz = (centy[k+1][j][i].x - centy[k][j][i].x); */
  /* 	  dydz = (centy[k+1][j][i].y - centy[k][j][i].y); */
  /* 	  dzdz = (centy[k+1][j][i].z - centy[k][j][i].z); */
  /* 	} */
  /* 	else if (k==mz-2) { */
  /* 	  dxdz = (centy[k][j][i].x - centy[k-1][j][i].x); */
  /* 	  dydz = (centy[k][j][i].y - centy[k-1][j][i].y); */
  /* 	  dzdz = (centy[k][j][i].z - centy[k-1][j][i].z); */
  /* 	} */

  /* 	jcsi[k][j][i].x = dyde * dzdz - dzde * dydz; */
  /* 	jcsi[k][j][i].y =-dxde * dzdz + dzde * dxdz; */
  /* 	jcsi[k][j][i].z = dxde * dydz - dyde * dxdz; */
	
  /* 	jeta[k][j][i].x = dydz * dzdc - dzdz * dydc; */
  /* 	jeta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc; */
  /* 	jeta[k][j][i].z = dxdz * dydc - dydz * dxdc; */

  /* 	jzet[k][j][i].x = dydc * dzde - dzdc * dyde; */
  /* 	jzet[k][j][i].y =-dxdc * dzde + dzdc * dxde; */
  /* 	jzet[k][j][i].z = dxdc * dyde - dydc * dxde; */

	
  /* 	jaj[k][j][i] = dxdc * (dyde * dzdz - dzde * dydz) - */
  /* 	  dydc * (dxde * dzdz - dzde * dxdz) + */
  /* 	  dzdc * (dxde * dydz - dyde * dxdz); */
  /* 	jaj[k][j][i] = 1./jaj[k][j][i]; */
	
  /*     } */
  /*   } */
  /* } */


  /* DMDAVecRestoreArray(fda, JCsi, &jcsi); */
  /* DMDAVecRestoreArray(fda, JEta, &jeta); */
  /* DMDAVecRestoreArray(fda, JZet, &jzet); */
  /* DMDAVecRestoreArray(da, JAj,  &jaj); */

  /* DMDAVecRestoreArray(fda, Centy, &centy); */
  /* VecDestroy(&Centy); */

  /* PetscPrintf(PETSC_COMM_WORLD, "metrics j direction done!"); */

  /* // k direction */
  /* DMDAVecGetArray(fda, KCsi, &kcsi); */
  /* DMDAVecGetArray(fda, KEta, &keta); */
  /* DMDAVecGetArray(fda, KZet, &kzet); */
  /* DMDAVecGetArray(da, KAj,  &kaj); */

  /* VecDuplicate(user->lCsi, &Centz); */
  /* DMDAVecGetArray(fda, Centz, &centz); */
  /* for(k=gzs; k<gze; k++) { */
  /*   for (j=gys+1; j<gye; j++) { */
  /*     for (i=gxs+1; i<gxe; i++) { */
  /* 	centz[k][j][i].x = (coor[k  ][j][i  ].x + coor[k][j-1][i  ].x + */
  /* 			    coor[k  ][j][i-1].x + coor[k][j-1][i-1].x) * 0.25; */
  /* 	centz[k][j][i].y = (coor[k  ][j][i  ].y + coor[k][j-1][i  ].y + */
  /* 			    coor[k  ][j][i-1].y + coor[k][j-1][i-1].y) * 0.25; */
  /* 	centz[k][j][i].z = (coor[k  ][j][i  ].z + coor[k][j-1][i  ].z + */
  /* 			    coor[k  ][j][i-1].z + coor[k][j-1][i-1].z) * 0.25; */
  /*     } */
  /*   } */
  /* } */

  /* for (k=zs; k<lze; k++) { */
  /*   for (j=lys; j<lye; j++) { */
  /*     for (i=lxs; i<lxe; i++) {   */
  /* 	if (i>1 && i<mx-2) { */
  /* 	  dxdc = (centz[k][j][i+1].x - centz[k][j][i-1].x) * 0.5; */
  /* 	  dydc = (centz[k][j][i+1].y - centz[k][j][i-1].y) * 0.5; */
  /* 	  dzdc = (centz[k][j][i+1].z - centz[k][j][i-1].z) * 0.5; */
  /* 	} */
  /* 	else if (i==1) { */
  /* 	  dxdc = centz[k][j][i+1].x - centz[k][j][i].x; */
  /* 	  dydc = centz[k][j][i+1].y - centz[k][j][i].y; */
  /* 	  dzdc = centz[k][j][i+1].z - centz[k][j][i].z; */
  /* 	} */
  /* 	else if (i==mx-2) { */
  /* 	  dxdc = centz[k][j][i].x - centz[k][j][i-1].x; */
  /* 	  dydc = centz[k][j][i].y - centz[k][j][i-1].y; */
  /* 	  dzdc = centz[k][j][i].z - centz[k][j][i-1].z; */
  /* 	} */

  /* 	if (j>1 && j<my-2) { */
  /* 	  dxde = (centz[k][j+1][i].x - centz[k][j-1][i].x) * 0.5; */
  /* 	  dyde = (centz[k][j+1][i].y - centz[k][j-1][i].y) * 0.5; */
  /* 	  dzde = (centz[k][j+1][i].z - centz[k][j-1][i].z) * 0.5; */
  /* 	} */
  /* 	else if (j==1) { */
  /* 	  dxde = centz[k][j+1][i].x - centz[k][j][i].x; */
  /* 	  dyde = centz[k][j+1][i].y - centz[k][j][i].y; */
  /* 	  dzde = centz[k][j+1][i].z - centz[k][j][i].z; */
  /* 	} */
  /* 	else if (j==my-2) { */
  /* 	  dxde = centz[k][j][i].x - centz[k][j-1][i].x; */
  /* 	  dyde = centz[k][j][i].y - centz[k][j-1][i].y; */
  /* 	  dzde = centz[k][j][i].z - centz[k][j-1][i].z; */
  /* 	} */

  /* 	if (k>0 && k<mz-2) { */
  /* 	  dxdz = (centz[k+1][j][i].x - centz[k-1][j][i].x) * 0.5; */
  /* 	  dydz = (centz[k+1][j][i].y - centz[k-1][j][i].y) * 0.5; */
  /* 	  dzdz = (centz[k+1][j][i].z - centz[k-1][j][i].z) * 0.5; */
  /* 	} */
  /* 	else if (k==0) { */
  /* 	  dxdz = (centz[k+1][j][i].x - centz[k][j][i].x); */
  /* 	  dydz = (centz[k+1][j][i].y - centz[k][j][i].y); */
  /* 	  dzdz = (centz[k+1][j][i].z - centz[k][j][i].z); */
  /* 	} */
  /* 	else if (k==mz-2) { */
  /* 	  dxdz = (centz[k][j][i].x - centz[k-1][j][i].x); */
  /* 	  dydz = (centz[k][j][i].y - centz[k-1][j][i].y); */
  /* 	  dzdz = (centz[k][j][i].z - centz[k-1][j][i].z); */
  /* 	} */

  /* 	kcsi[k][j][i].x = dyde * dzdz - dzde * dydz; */
  /* 	kcsi[k][j][i].y =-dxde * dzdz + dzde * dxdz; */
  /* 	kcsi[k][j][i].z = dxde * dydz - dyde * dxdz; */

  /* 	keta[k][j][i].x = dydz * dzdc - dzdz * dydc; */
  /* 	keta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc; */
  /* 	keta[k][j][i].z = dxdz * dydc - dydz * dxdc; */

  /* 	kzet[k][j][i].x = dydc * dzde - dzdc * dyde; */
  /* 	kzet[k][j][i].y =-dxdc * dzde + dzdc * dxde; */
  /* 	kzet[k][j][i].z = dxdc * dyde - dydc * dxde; */


  /* 	kaj[k][j][i] = dxdc * (dyde * dzdz - dzde * dydz) - */
  /* 	  dydc * (dxde * dzdz - dzde * dxdz) + */
  /* 	  dzdc * (dxde * dydz - dyde * dxdz); */
  /* 	kaj[k][j][i] = 1./kaj[k][j][i]; */

  /* 	/\*	if ((i==1 || i==2) && j==1 && k==5) { */
  /* 	  PetscPrintf(PETSC_COMM_WORLD, "\nDI %e %e %e \n %e %e %e \n%e", dxdc, dydc, dzdc, dxde, dyde, dzde, kzet[k][j][i].z); */
  /* 	  } *\/ */
  /*     } */
  /*   } */
  /* } */
  /* DMDAVecRestoreArray(fda, Centz, &centz); */
  /* VecDestroy(&Centz); */
  
  /* DMDAVecRestoreArray(fda, KCsi, &kcsi); */
  /* DMDAVecRestoreArray(fda, KEta, &keta); */
  /* DMDAVecRestoreArray(fda, KZet, &kzet); */
  /* DMDAVecRestoreArray(da, KAj,  &kaj); */
  /* DMDAVecRestoreArray(fda, user->Cent, &cent); */

  DMDAVecRestoreArray(fda, coords, &coor);

  DMDAVecRestoreArray(fda, user->Csi, &csi);
  DMDAVecRestoreArray(fda, user->Eta, &eta);
  DMDAVecRestoreArray(fda, user->Zet, &zet);
  DMDAVecRestoreArray(da, user->Aj,  &aj);

  DMGlobalToLocalBegin(fda, user->Csi, INSERT_VALUES, user->lCsi);
  DMGlobalToLocalEnd(fda, user->Csi, INSERT_VALUES, user->lCsi);

  DMGlobalToLocalBegin(fda, user->Eta, INSERT_VALUES, user->lEta);
  DMGlobalToLocalEnd(fda, user->Eta, INSERT_VALUES, user->lEta);

  DMGlobalToLocalBegin(fda, user->Zet, INSERT_VALUES, user->lZet);
  DMGlobalToLocalEnd(fda, user->Zet, INSERT_VALUES, user->lZet);

  DMGlobalToLocalBegin(da, user->Aj, INSERT_VALUES, user->lAj);
  DMGlobalToLocalEnd(da, user->Aj, INSERT_VALUES, user->lAj);

  /* DMGlobalToLocalBegin(fda, user->ICsi, INSERT_VALUES, user->lICsi); */
  /* DMGlobalToLocalEnd(fda, user->ICsi, INSERT_VALUES, user->lICsi); */

  /* DMGlobalToLocalBegin(fda, user->IEta, INSERT_VALUES, user->lIEta); */
  /* DMGlobalToLocalEnd(fda, user->IEta, INSERT_VALUES, user->lIEta); */

  /* DMGlobalToLocalBegin(fda, user->IZet, INSERT_VALUES, user->lIZet); */
  /* DMGlobalToLocalEnd(fda, user->IZet, INSERT_VALUES, user->lIZet); */

  /* DMGlobalToLocalBegin(fda, user->JCsi, INSERT_VALUES, user->lJCsi); */
  /* DMGlobalToLocalEnd(fda, user->JCsi, INSERT_VALUES, user->lJCsi); */

  /* DMGlobalToLocalBegin(fda, user->JEta, INSERT_VALUES, user->lJEta); */
  /* DMGlobalToLocalEnd(fda, user->JEta, INSERT_VALUES, user->lJEta); */

  /* DMGlobalToLocalBegin(fda, user->JZet, INSERT_VALUES, user->lJZet); */
  /* DMGlobalToLocalEnd(fda, user->JZet, INSERT_VALUES, user->lJZet); */

  /* DMGlobalToLocalBegin(fda, user->KCsi, INSERT_VALUES, user->lKCsi); */
  /* DMGlobalToLocalEnd(fda, user->KCsi, INSERT_VALUES, user->lKCsi); */

  /* DMGlobalToLocalBegin(fda, user->KEta, INSERT_VALUES, user->lKEta); */
  /* DMGlobalToLocalEnd(fda, user->KEta, INSERT_VALUES, user->lKEta); */

  /* DMGlobalToLocalBegin(fda, user->KZet, INSERT_VALUES, user->lKZet); */
  /* DMGlobalToLocalEnd(fda, user->KZet, INSERT_VALUES, user->lKZet); */

  /* DMGlobalToLocalBegin(da, user->IAj, INSERT_VALUES, user->lIAj); */
  /* DMGlobalToLocalEnd(da, user->IAj, INSERT_VALUES, user->lIAj); */

  /* DMGlobalToLocalBegin(da, user->JAj, INSERT_VALUES, user->lJAj); */
  /* DMGlobalToLocalEnd(da, user->JAj, INSERT_VALUES, user->lJAj); */

  /* DMGlobalToLocalBegin(da, user->KAj, INSERT_VALUES, user->lKAj); */
  /* DMGlobalToLocalEnd(da, user->KAj, INSERT_VALUES, user->lKAj); */

  /* DMGlobalToLocalBegin(fda, user->GridSpace, INSERT_VALUES, user->lGridSpace); */
  /* DMGlobalToLocalEnd(fda, user->GridSpace, INSERT_VALUES, user->lGridSpace); */

  /* DMGlobalToLocalBegin(fda, user->Cent, INSERT_VALUES, user->lCent); */
  /* DMGlobalToLocalEnd(fda, user->Cent, INSERT_VALUES, user->lCent); */
  
  PetscReal nCsi, nEta, nZet, nAj, nCoor;
  VecNorm(coords, NORM_INFINITY, &nCoor);
  VecNorm(user->lCsi, NORM_INFINITY, &nCsi);  
  VecNorm(user->lEta, NORM_INFINITY, &nEta);  
  VecNorm(user->lZet, NORM_INFINITY, &nZet);
  VecNorm(user->lAj, NORM_INFINITY, &nAj);
  PetscPrintf(PETSC_COMM_WORLD, "!!norm of %le csi %le %le %le %le\n",nCoor, nCsi, nEta, nZet, nAj);
  PetscBarrier(NULL);

  VecDestroy(&user->Csi);
  VecDestroy(&user->Eta);
  VecDestroy(&user->Zet);
  VecDestroy(&user->Aj);

  /* VecDestroy(&user->ICsi); */
  /* VecDestroy(&user->IEta); */
  /* VecDestroy(&user->IZet); */

  /* VecDestroy(&user->JCsi); */
  /* VecDestroy(&user->JEta); */
  /* VecDestroy(&user->JZet); */

  /* VecDestroy(&user->KCsi); */
  /* VecDestroy(&user->KEta); */
  /* VecDestroy(&user->KZet); */

  /* VecDestroy(&user->IAj); */
  /* VecDestroy(&user->JAj); */
  /* VecDestroy(&user->KAj); */

  return(0);
} 


/* ==================================================================================             */

/* ==================================================================================             */
/* struct Cmpn{ */
/*   PetscScalar x, y, z; */
/* }; */	

Cmpnts DIVC(Cmpnts v1, PetscReal c)
{
  Cmpnts v4;
  v4.x  = v1.x/c;
  v4.y  = v1.y/c;
  v4.z  = v1.z/c;
 
  return(v4);
}

Cmpnts PLUS(Cmpnts v1,Cmpnts v2)
{
  Cmpnts v4;
  v4.x  = v1.x + v2.x ;
  v4.y  = v1.y + v2.y ;
  v4.z  = v1.z + v2.z ;
 
  return(v4);
}

Cmpnts MINUS(Cmpnts v1,Cmpnts v2)
{
  Cmpnts v4;
  v4.x  = v1.x - v2.x ;
  v4.y  = v1.y - v2.y ;
  v4.z  = v1.z - v2.z ;
 
  return(v4);
}

Cmpnts AVERAGE(Cmpnts v1,Cmpnts v2, Cmpnts v3)
{
  Cmpnts v4;
  v4.x  = (v1.x + v2.x +v3.x)/3.;
  v4.y  = (v1.y + v2.y +v3.y)/3.;
  v4.z  = (v1.z + v2.z +v3.z)/3.;
 
  return(v4);
}

Cmpnts AVERAGE4(Cmpnts v1,Cmpnts v2, Cmpnts v3, Cmpnts v5)
{
  Cmpnts v4;
  v4.x  = (v1.x + v2.x +v3.x+v5.x)/4.;
  v4.y  = (v1.y + v2.y +v3.y+v5.y)/4.;
  v4.z  = (v1.z + v2.z +v3.z+v5.z)/4.;
 
  return(v4);
}

Cmpnts Cross(Cmpnts v1,Cmpnts v2)
{
  // output = v1 x v2
  Cmpnts v1cv2;
  v1cv2.x  = v1.y * v2.z - v2.y * v1.z;
  v1cv2.y  =-v1.x * v2.z + v2.x * v1.z;
  v1cv2.z  = v1.x * v2.y - v2.x * v1.y;
 
  return(v1cv2);
}

PetscErrorCode CalculateTrinagleArea(Cmpnts p1, Cmpnts p2, Cmpnts p3, Cmpnts *Area)
{
  Cmpnts edge1, edge2;
  edge1.x=(p2.x-p1.x)/2.;
  edge1.y=(p2.y-p1.y)/2.;
  edge1.z=(p2.z-p1.z)/2.;

  edge2.x=p3.x-p1.x;
  edge2.y=p3.y-p1.y;
  edge2.z=p3.z-p1.z;

  *Area=Cross(edge1,edge2);  

  return 0;
}
/* ==================================================================================             */
/*  Calculating the Moment of Area of each side of the grid

    For reference see Vinokur, JCP, 1989 Eqn (61)
    M1234=M123+M143
    M123= r x n dS

    lMArea the numbering is the same as flux (Ucont) 

    numbering
    p1    p2

    p3    p4

    Synopsis
    PetscErrorCode FormAreaMoment(UserCtx *user, Cmpnts a_c)
    
    user: user context
    a_c : The reference fram location

*/
/* ==================================================================================             */
PetscErrorCode FormAreaMoment(UserCtx *user, Cmpnts a_c)
{
  DM da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  Cmpnts     ***marea, p1,p2,p3,p4, A123, A432, r123,r432;
 
  Vec        Coor;
  Cmpnts     ***coor;
  Cmpnts     ***csi, ***eta, ***zet;

  DMGetCoordinatesLocal(da, &Coor);
  DMDAVecGetArray(fda, Coor, &coor);
  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);

  DMDAVecGetArray(fda, user->lMAreaCsi, &marea);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<lxe; i++) {
	p4 = coor[k  ][j  ][i];
	p2 = coor[k-1][j  ][i];
	p3 = coor[k  ][j-1][i];
	p1 = coor[k-1][j-1][i];

	CalculateTrinagleArea(p1,p2,p3, &A123);
	CalculateTrinagleArea(p4,p3,p2, &A432);
	
	r123=MINUS(AVERAGE(p1,p2,p3),a_c);
	r432=MINUS(AVERAGE(p4,p2,p3),a_c);

	marea[k][j][i]= PLUS(Cross(r123, A123),Cross(r432,A432));

/* 	r123=MINUS(AVERAGE4(p1,p2,p3,p4),a_c); */
/* 	marea[k][j][i]=Cross(r123,csi[k][j][i]); */
	  
      }
    }
  }

  DMDAVecRestoreArray(fda, user->lMAreaCsi, &marea);

  DMDAVecGetArray(fda, user->lMAreaEta, &marea);

  for (k=lzs; k<lze; k++) {
    for (j=ys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	p4 = coor[k  ][j][i  ];
	p3 = coor[k-1][j][i  ];
	p2 = coor[k  ][j][i-1];
	p1 = coor[k-1][j][i-1];

	CalculateTrinagleArea(p1,p2,p3, &A123);
	CalculateTrinagleArea(p4,p3,p2, &A432);
	
	r123=MINUS(AVERAGE(p1,p2,p3),a_c);
	r432=MINUS(AVERAGE(p4,p2,p3),a_c);

	marea[k][j][i]= PLUS(Cross(r123, A123),Cross(r432,A432));

/* 	r123=MINUS(AVERAGE4(p1,p2,p3,p4),a_c); */
/* 	marea[k][j][i]=Cross(r123,eta[k][j][i]); */
      }
    }
  }

  DMDAVecRestoreArray(fda, user->lMAreaEta, &marea);

  DMDAVecGetArray(fda, user->lMAreaZet, &marea);

  for (k=zs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	p4 = coor[k][j  ][i  ];
	p3 = coor[k][j  ][i-1];
	p2 = coor[k][j-1][i  ];
	p1 = coor[k][j-1][i-1];

	CalculateTrinagleArea(p1,p2,p3, &A123);
	CalculateTrinagleArea(p4,p3,p2, &A432);
	
	r123=MINUS(AVERAGE(p1,p2,p3),a_c);
	r432=MINUS(AVERAGE(p4,p2,p3),a_c);

	marea[k][j][i]= PLUS(Cross(r123, A123),Cross(r432,A432));

/* 	r123=MINUS(AVERAGE4(p1,p2,p3,p4),a_c); */
/* 	marea[k][j][i]=Cross(r123,zet[k][j][i]); */
	
      }
    }
  }

  DMDAVecRestoreArray(fda, user->lMAreaZet, &marea);

  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  DMDAVecRestoreArray(fda, Coor, &coor);
  //VecDestroy(&Coor);

  DMLocalToLocalBegin(fda, user->lMAreaCsi, INSERT_VALUES, user->lMAreaCsi);
  DMLocalToLocalEnd(fda, user->lMAreaCsi, INSERT_VALUES, user->lMAreaCsi);

  DMLocalToLocalBegin(fda, user->lMAreaEta, INSERT_VALUES, user->lMAreaEta);
  DMLocalToLocalEnd(fda, user->lMAreaEta, INSERT_VALUES, user->lMAreaEta);

  DMLocalToLocalBegin(fda, user->lMAreaZet, INSERT_VALUES, user->lMAreaZet);
  DMLocalToLocalEnd(fda, user->lMAreaZet, INSERT_VALUES, user->lMAreaZet);

  return 0;
}


PetscErrorCode MomentAreaDivergence(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  Vec		Div;
  PetscReal	***div, ***aj;
  Cmpnts	***csi, ***eta, ***zet;

  PetscReal	maxdiv;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda, user->lMAreaCsi, &csi);
  DMDAVecGetArray(fda, user->lMAreaEta, &eta);
  DMDAVecGetArray(fda, user->lMAreaZet, &zet);

  DMDAVecGetArray(da, user->lAj, &aj);
  VecDuplicate(user->P, &Div);
  VecSet(Div, 0.);
  DMDAVecGetArray(da, Div, &div);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	maxdiv = fabs((csi[k][j][i].x - csi[k][j][i-1].x +
		       eta[k][j][i].x - eta[k][j-1][i].x +
		       zet[k][j][i].x - zet[k-1][j][i].x +
		       csi[k][j][i].y - csi[k][j][i-1].y +
		       eta[k][j][i].y - eta[k][j-1][i].y +
		       zet[k][j][i].y - zet[k-1][j][i].y +
		       csi[k][j][i].z - csi[k][j][i-1].z +
		       eta[k][j][i].z - eta[k][j-1][i].z +
		       zet[k][j][i].z - zet[k-1][j][i].z)*aj[k][j][i]);
	div[k][j][i] = maxdiv;
      }
    }
  }

  DMDAVecRestoreArray(da, Div, &div);  
  
  DMDAVecRestoreArray(fda, user->lMAreaCsi, &csi);
  DMDAVecRestoreArray(fda, user->lMAreaEta, &eta);
  DMDAVecRestoreArray(fda, user->lMAreaZet, &zet);

  DMDAVecRestoreArray(da, user->lAj, &aj);

  VecMax(Div, &i, &maxdiv);
  PetscPrintf(PETSC_COMM_WORLD, "Maxdiv Moment of AREA Metrics %d %e\n", i, maxdiv);

  /* PetscInt mi; */
  /* for (k=zs; k<ze; k++) { */
  /*   for (j=ys; j<ye; j++) { */
  /*     for (mi=xs; mi<xe; mi++) { */
  /* 	if (Gidx(mi,j,k,user) ==i) { */
  /* 	  PetscPrintf(PETSC_COMM_SELF, "MMa %d %d %d\n", mi,j, k); */
  /* 	} */
  /*     } */
  /*   } */
  /* } */

  VecDestroy(&Div);
  return(0);
}


PetscErrorCode MetricsDivergence(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  Vec		Div;
  PetscReal	***div, ***aj;
  Cmpnts	***csi, ***eta, ***zet;

  PetscReal	maxdiv;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);

  DMDAVecGetArray(da, user->lAj, &aj);
  VecDuplicate(user->P, &Div);
  VecSet(Div, 0.);
  DMDAVecGetArray(da, Div, &div);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	maxdiv = fabs((csi[k][j][i].x - csi[k][j][i-1].x +
		       eta[k][j][i].x - eta[k][j-1][i].x +
		       zet[k][j][i].x - zet[k-1][j][i].x +
		       csi[k][j][i].y - csi[k][j][i-1].y +
		       eta[k][j][i].y - eta[k][j-1][i].y +
		       zet[k][j][i].y - zet[k-1][j][i].y +
		       csi[k][j][i].z - csi[k][j][i-1].z +
		       eta[k][j][i].z - eta[k][j-1][i].z +
		       zet[k][j][i].z - zet[k-1][j][i].z)*aj[k][j][i]);
	div[k][j][i] = maxdiv;
      }
    }
  }

  DMDAVecRestoreArray(da, Div, &div);
  
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);

  DMDAVecRestoreArray(da, user->lAj, &aj);

  VecMax(Div, &i, &maxdiv);
  PetscPrintf(PETSC_COMM_WORLD, "Maxdiv Metrics %d %e\n", i, maxdiv);
  /* PetscInt mi; */
  /* for (k=zs; k<ze; k++) { */
  /*   for (j=ys; j<ye; j++) { */
  /*     for (mi=xs; mi<xe; mi++) { */
  /* 	if (Gidx(mi,j,k,user) ==i) { */
  /* 	  PetscPrintf(PETSC_COMM_SELF, "MMa %d %d %d\n", mi,j, k); */
  /* 	} */
  /*     } */
  /*   } */
  /* } */

  VecDestroy(&Div);
  return(0);
}
