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
Cmpnts IndxCmpnt(Cmpnts v1, PetscInt i){
    if (i==0) {return v1.x;}
    else if (i==1) {return v1.y;}
    else {return v1.z;}
}

PetscErrorCode formViscous(UserCtx *user, Vec Visc)
{

  Vec		Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;
  
  
  Cmpnts        ***ucat, ucat_h;
  
  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***alpha[3],t,Ss;
  Vec           Alpha1,Alpha2,Alpha3;
  
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
  PetscReal u_csi[3][3], p_csi[3];// u_csi[0][0]->dudc, u_csi[1][0]->dvdc, u_csi[2][0]->dvdc etc

  if (twoD!=1) {
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs-1; i<lxe; i++) {
	    // ∂ui/∂ξ --> cartesian velocity gradient in curvilinear (IndxCmpnt to Get .x,.y.z components for calc)
        for(PetscInt a = 0; a<3; a++){
        u_csi[a][0] = IndxCmpnt(ucat[k][j][i+1],a)- IndxCmpnt(ucat[k][j][i],a);
        }
    	p_csi[0] = p[k][j][i+1]/q[k][j][i+1].rho - p[k][j][i]/q[k][j][i].rho;
	
	// ∂ui/∂η --> cartesian velocity gradient in curvilinear
	if ((nvert[k][j+1][i  ]> solid && nvert[k][j+1][i  ]<innerblank)  ||
	    (nvert[k][j+1][i+1]> solid && nvert[k][j+1][i+1]<innerblank)) {
        for(a=0; a<3; a++){
          u_csi[a][1] = (IndxCmpnt(ucat[k][j  ][i+1],a) + IndxCmpnt(ucat[k][j  ][i],a) -
		                IndxCmpnt(ucat[k][j-1][i+1],a) - IndxCmpnt(ucat[k][j-1][i],a)) * 0.5;
        }
	    p_csi[1] = (p[k][j  ][i+1]/q[k][j][i+1].rho + p[k][j  ][i]/q[k][j][i].rho -
		            p[k][j-1][i+1]/q[k][j-1][i+1].rho - p[k][j-1][i]/q[k][j-1][i].rho) * 0.5;
	}
	else if  ((nvert[k][j-1][i  ]> solid && nvert[k][j-1][i  ]<innerblank)  ||
		        (nvert[k][j-1][i+1]> solid && nvert[k][j-1][i+1]<innerblank)) {
        for(a=0; a<3; a++){
          u_csi[a][1] = (IndxCmpnt(ucat[k][j+1][i+1],a) + IndxCmpnt(ucat[k][j+1][i],a) -
		                IndxCmpnt(ucat[k][j][i+1],a) - IndxCmpnt(ucat[k][j][i],a)) * 0.5;
        }
	    p_csi[1] = (p[k][j+1][i+1]/q[k][j+1][i+1].rho + p[k][j+1][i]/q[k][j+1][i].rho -
		            p[k][j  ][i+1]/q[k][j  ][i+1].rho - p[k][j  ][i]/q[k][j  ][i].rho) * 0.5;
	}
	else {
        for(a=0; a<3; a++){
          u_csi[a][1] = (IndxCmpnt(ucat[k][j+1][i+1],a) + IndxCmpnt(ucat[k][j+1][i],a) -
		                IndxCmpnt(ucat[k][j-1][i+1],a) - IndxCmpnt(ucat[k][j-1][i],a)) * 0.25;
        }
	        p_csi[1] =  (p[k][j+1][i+1]/q[k][j+1][i+1].rho + p[k][j+1][i]/q[k][j+1][i].rho -
		                p[k][j-1][i+1]/q[k][j-1][i+1].rho - p[k][j-1][i]/q[k][j-1][i].rho) * 0.25; 
	  
	}

	// ∂ui/∂ζ --> cartesian velocity gradient in curvilinear	
	if ((nvert[k+1][j][i  ]> solid && nvert[k+1][j][i  ]<innerblank)||
	    (nvert[k+1][j][i+1]> solid && nvert[k+1][j][i+1]<innerblank)) {
        for(a=0; a<3; a++){
          u_csi[a][2] = (IndxCmpnt(ucat[k][j][i+1],a) + IndxCmpnt(ucat[k][j][i],a) -
		                IndxCmpnt(ucat[k-1][j][i+1],a) - IndxCmpnt(ucat[k-1][j][i],a)) * 0.5;
        }   
	     p_csi[2] = (p[k][j][i+1]/q[k][j][i+1].rho + p[k][j][i]/q[k][j][i].rho -
		            p[k-1][j][i+1]/q[k-1][j][i+1].rho - p[k][j][i]/q[k][j][i].rho) * 0.5;		
	}
	else if ((nvert[k-1][j][i  ]> solid && nvert[k-1][j][i  ]<innerblank) ||
		    (nvert[k-1][j][i+1]> solid && nvert[k-1][j][i+1]<innerblank)) {
	  for(a=0; a<3; a++){
          u_csi[a][2] = (IndxCmpnt(ucat[k+1][j][i+1],a) + IndxCmpnt(ucat[k+1][j][i],a) -
		                IndxCmpnt(ucat[k][j][i+1],a) - IndxCmpnt(ucat[k][j][i],a)) * 0.5;
        }   
	     p_csi[2] = (p[k+1][j][i+1]/q[k+1][j][i+1].rho + p[k+1][j][i]/q[k+1][j][i].rho -
		            p[k][j][i+1]/q[k][j][i+1].rho - p[k][j][i]/q[k][j][i].rho) * 0.5;
	}
	else {
    
      for(a=0; a<3; a++){
          u_csi[a][2] = (IndxCmpnt(ucat[k+1][j][i+1],a) + IndxCmpnt(ucat[k+1][j][i],a) -
		                IndxCmpnt(ucat[k-1][j][i+1],a) - IndxCmpnt(ucat[k-1][j][i],a)) * 0.25;
        }   
	     p_csi[2] = (p[k+1][j][i+1]/q[k+1][j][i+1].rho + p[k+1][j][i]/q[k+1][j][i].rho -
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
					 α = uᵢ ξ𐞥ⱼ ∂ (ρ x τᵢⱼ)
								∂ξ𐞥
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
			}// End if LES CC
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