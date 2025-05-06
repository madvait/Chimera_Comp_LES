#include "variables.h"
#include "petsctime.h" 
//#include "poisson.c"
static char help[] = "CURVIB compressible LES with PETSc/3.12.4-intel-2020a-Python-3.8.2 \n\n";

/* Good programming indicates that we should try to avoid global variables as much as 
   Possible to increase protability and modularity of the code.  */

//PetscInt ti,tistart=0;
//PetscReal	Flux_in = 4.104388e-04, angle = 0;
//PetscInt tiout = 10;
//PetscInt block_number = 1;
//PetscReal FluxInSum, FluxOutSum;
PetscInt tistart = 0;             //Advait: brought back to global scope!
PetscBool rstart_flg = PETSC_FALSE;

/* immeresed value>1 determines the type of correction
   1      constant velocity correction
   2      proportional (needs FluxOutSum>0)
   3      proportional to flux
   4      proportional to normal velocity (flux/area)
*/
PetscInt immersed = 0;
PetscInt wavy = 0,isothermal=0;
PetscInt calc_sigma =0; 
PetscInt invicid = 0;
PetscInt movefsi = 0, rotatefsi=0, sediment=0;
//PetscBool rstart_flg;
//PetscInt implicit = 0;
//PetscInt imp_MAX_IT = 50; 
//PetscInt radi=10;
//PetscInt inletprofile=2; 
PetscReal CMx_c=0., CMy_c=0., CMz_c=0.;
//PetscInt  mg_MAX_IT=30, mg_idx=1, mg_preItr=1, mg_poItr=1;
//PetscReal imp_atol=1e-7, imp_rtol=1e-4, imp_stol=1.e-8;
PetscInt TwoD = 0;
PetscInt STRONG_COUPLING=0;
PetscInt rstart_fsi=0;
PetscInt cop=0, regime=1; // 1 escape regime --- 2 cruise regime
PetscInt fish=0, fish_c=0, eel=0, fishcyl=0, rheology=0,duplicate=0,turbine=0,Pipe=0,pizza=0;
PetscInt wing=0, hydro=0;
PetscReal St_exp=0.5,wavelength=0.95;
PetscInt MHV=0,LV=0;
PetscReal max_angle = -54.*3.1415926/180.;// for MHV; min=0
PetscInt thin=0;
PetscInt dgf_z=0,dgf_y=1,dgf_x=0;
PetscInt dgf_az=0,dgf_ay=0,dgf_ax=1 ;
PetscInt NumberOfBodies=1;
PetscInt moveframe=0,rotateframe=0;//, blank=0;
//PetscReal L_dim=1.;
//PetscReal Mach = 1.5;
PetscInt central=2;
PetscInt les=0, rans=0; //les ==1 Constant, >1 Dynamic Smagorinsky
PetscInt wallfunction=0;
PetscInt averaging=0;
PetscInt mixed=0;
PetscInt clark=0;
PetscInt dynamic_freq=1;
PetscInt periodic=0;
PetscReal max_cs=0.5;
PetscInt grid1d=0;
PetscInt i_periodic=0;
PetscInt j_periodic=0;
PetscInt k_periodic=0;
PetscInt blkpbc=10;
PetscInt pseudo_periodic=0;
PetscInt testfilter_ik=0;
PetscInt testfilter_1d=0;
PetscInt i_homo_filter=0;
PetscInt j_homo_filter=0;
PetscInt k_homo_filter=0;
PetscInt channel=0;
//double poisson_tol=5.e-9;	// relative tolerance
double pi=3.14158692;
//PetscInt  fluxorder=2; // type==1 Godunov, 2 MUSCL, 3 WENO


//IBMNodes	*ibm_ptr;



PetscErrorCode Ucont_P_Binary_Input(UserCtx *user)
{

  PetscViewer	viewer;

  char filen[90];

  
  sprintf(filen, "vfield%5.5d_%1.1d.dat", user->ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  PetscInt N;

  VecGetSize(user->Ucont, &N);
  PetscPrintf(PETSC_COMM_WORLD, "PPP %d\n", N);
  VecLoad((user->Ucont),viewer);
  
  PetscViewerDestroy(&viewer);

  PetscBarrier(NULL);

  PetscOptionsClearValue(NULL, "-vecload_block_size");
  //PetscOptionsClearValue("-vecload_block_size");


  sprintf(filen, "pfield%5.5d_%1.1d.dat", user->ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad((user->P),viewer);
  PetscViewerDestroy(&viewer);

  sprintf(filen, "nvfield%5.5d_%1.1d.dat", user->ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad((user->Nvert_o),viewer);
  PetscViewerDestroy(&viewer);

  DMGlobalToLocalBegin(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);
  DMGlobalToLocalEnd(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);

  if(averaging) {	// Seokkoo Kang
    sprintf(filen, "su0_%06d_%1d.dat",user->ti, user->_this);
    FILE *fp=fopen(filen, "r");
    
    VecSet(user->Ucat_sum, 0);
    VecSet(user->Ucat_cross_sum, 0);
    VecSet(user->Ucat_square_sum, 0);
    VecSet(user->P_sum, 0);
    
    if(fp==NULL) {
      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, setting the statistical quantities to zero and contiues the computation ... ***\n\n", filen);
    }
    else {
      fclose(fp);
      PetscBarrier(NULL);
      sprintf(filen, "su0_%06d_%1d.dat", user->ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->Ucat_sum,viewer);
      PetscViewerDestroy(&viewer);
      
      sprintf(filen, "su1_%06d_%1d.dat", user->ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad(user->Ucat_cross_sum,viewer);
      PetscViewerDestroy(&viewer);
      
      sprintf(filen, "su2_%06d_%1d.dat", user->ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad((user->Ucat_square_sum),viewer);
      PetscViewerDestroy(&viewer);
      
      sprintf(filen, "sp_%06d_%1d.dat", user->ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->P_sum,viewer);
      PetscViewerDestroy(&viewer);
      
      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Read %s, continuing averaging ... ***\n\n", filen);
    }
  }
  
  if(les) {
    Vec Cs;
    
    VecDuplicate(user->P, &Cs);
    
    sprintf(filen, "cs_%06d_%1d.dat", user->ti, user->_this);
    FILE *fp=fopen(filen, "r");
    
    if(fp==NULL) {
      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, setting Cs to 0 and contiues the computation ... ***\n\n", filen);
      VecSet(Cs, 0);
    }
    else {
      fclose(fp);
      
      PetscBarrier(NULL);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( Cs,viewer);
      PetscViewerDestroy(&viewer);
    }
    
    DMGlobalToLocalBegin(user->da, Cs, INSERT_VALUES, user->lCs);
    DMGlobalToLocalEnd(user->da, Cs, INSERT_VALUES, user->lCs);
    
    VecDestroy(&Cs);
  }

  
   if(rans) {
    // K-Omega
    sprintf(filen, "kfield%06d_%1d.dat", user->ti, user->_this);
    FILE *fp=fopen(filen, "r");
    
    if(fp!=NULL) {
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad(user->K_Omega,viewer);
      PetscViewerDestroy(&viewer);
    }
    else {
      //      K_Omega_IC(user);
    }
    
    VecCopy(user->K_Omega, user->K_Omega_o);
    
    DMGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
    
    DMGlobalToLocalBegin(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);    
  }

  return 0;
}

PetscErrorCode Q_P_Binary_Input(UserCtx *user)
{
  PetscViewer	 viewer;  
  char filen2[90];
  //  PetscOptionsClearValue("-vecload_block_size");

  sprintf(filen2, "pfield%5.5d_%1.1d.dat", user->ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &viewer);
  VecLoad((user->P),viewer);
  PetscReal norm;
  VecNorm(user->P, NORM_INFINITY, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "PIn %le\n", norm);
  PetscViewerDestroy(&viewer);

  sprintf(filen2, "nvfield%5.5d_%1.1d.dat", user->ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &viewer);
  VecLoad((user->Nvert),viewer);
  PetscViewerDestroy(&viewer);

  PetscInt N;
  VecGetSize(user->Q, &N);
  PetscPrintf(PETSC_COMM_WORLD, "Input Q size %d\n", N);

  sprintf(filen2, "Qfield%5.5d_%1.1d.dat", user->ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &viewer);
  VecLoad((user->Q),viewer);
  PetscViewerDestroy(&viewer);

  return(0);
}


PetscErrorCode Q_P_Binary_Output(UserCtx *user)
{
  PetscViewer	viewer;
  char filen[80];
  
  sprintf(filen, "Qfield%5.5d_%1.1d.dat", user->ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  VecView(user->Q, viewer);
  PetscViewerDestroy(&viewer);

  if (les){ // if (wavy) {
  sprintf(filen, "sigma%5.5d_%1.1d.dat", user->ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  VecView(user->Sigma, viewer);
  PetscViewerDestroy(&viewer);
  }
 
  sprintf(filen, "pfield%5.5d_%1.1d.dat", user->ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  VecView(user->P, viewer);
  PetscViewerDestroy(&viewer);
 
  sprintf(filen, "nvfield%5.5d_%1.1d.dat", user->ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  VecView(user->Nvert, viewer);
  PetscViewerDestroy(&viewer);

  if (averaging && user->ti != 0){
    sprintf(filen, "Qfield_ave%5.5d_%1.1d.dat", user->ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Q_sum, viewer);
    PetscViewerDestroy(&viewer);
        
    sprintf(filen, "su2_Sq%5.5d_%1.1d.dat", user->ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_square_sum, viewer);
    PetscViewerDestroy(&viewer);

    sprintf(filen, "su_cross%5.5d_%1d.dat", user->ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_cross_sum, viewer);
    PetscViewerDestroy(&viewer);
    
    sprintf(filen, "psum%5.5d_%1.1d.dat", user->ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->P_sum, viewer);
    PetscViewerDestroy(&viewer);   
  }
  
  if(les){
  sprintf(filen, "nu_t%5.5d_%1.1d.dat", user->ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  VecView(user->Nu_t, viewer);
  PetscViewerDestroy(&viewer);
  }

  return 0;
}

PetscErrorCode Q_P_VTK_Input(UserCtx *user)
{
  // seems that PETSC does not have support for vtk read

 PetscViewer viewer;
 char filen[80];
  
 sprintf(filen, "Qfield%5.5d_%1.1d.vts", user->ti, user->_this);

 PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
 PetscViewerSetType(viewer, PETSCVIEWERVTK);
 PetscViewerFileSetMode(viewer, FILE_MODE_READ);
 PetscViewerFileSetName(viewer, filen);

 // PetscViewerVTKOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
 PetscObjectSetName((PetscObject) user->Q, "Q");
 PetscObjectSetOptionsPrefix((PetscObject) user->Q, "Q");
 VecLoad(user->Q, viewer);
 PetscObjectSetName((PetscObject) user->P, "P");
 PetscObjectSetOptionsPrefix((PetscObject) user->P, "P");
 VecLoad(user->P, viewer);
 PetscObjectSetName((PetscObject) user->Nvert, "NV");
 PetscObjectSetOptionsPrefix((PetscObject) user->Nvert, "NV");
 VecLoad(user->Nvert, viewer);

 PetscReal norm;
 VecNorm(user->P, NORM_INFINITY, &norm);
 PetscPrintf(PETSC_COMM_WORLD, "PIn %le\n", norm);

 if (averaging){
   VecLoad(user->Q_sum, viewer);
   VecLoad(user->P_sum, viewer);
   VecLoad(user->Ucat_square_sum, viewer);   
 }

 PetscViewerDestroy(&viewer);

 return 0;
}

PetscErrorCode Q_P_VTK_Output(UserCtx *user)
{
 PetscViewer viewer;
 char filen[80];
  
 sprintf(filen, "Qfield%5.5d_%1.1d.vts", user->ti, user->_this);
 PetscObjectSetName((PetscObject) user->Q, "Q");
 PetscObjectSetOptionsPrefix((PetscObject) user->Q, "Q");
 PetscObjectSetName((PetscObject) user->P, "P");
 PetscObjectSetOptionsPrefix((PetscObject) user->P, "P");
 PetscObjectSetName((PetscObject) user->Nvert, "NV");
 PetscObjectSetOptionsPrefix((PetscObject) user->Nvert, "NV");

 PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
 PetscViewerSetType(viewer, PETSCVIEWERVTK);
 PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
 PetscViewerFileSetName(viewer, filen);   
 VecView(user->Q, viewer);
 VecView(user->P, viewer);
 VecView(user->Nvert, viewer);

 if (averaging){
   // sprintf(filen, "Qfield_ave%5.5d_%1.1d.dat", user->ti, user->_this);
   PetscObjectSetName((PetscObject) user->Q_sum, "Q_sum");
   PetscObjectSetName((PetscObject) user->P_sum, "P_sum");
   PetscObjectSetName((PetscObject) user->Ucat_square_sum, "u_sqr_sum");

   //   PetscViewerFileSetName(viewer, filen);   
   VecView(user->Q_sum, viewer);
   VecView(user->P_sum, viewer);
   VecView(user->Ucat_square_sum, viewer);   
 }

 PetscViewerDestroy(&viewer);

 return 0;
}

PetscErrorCode Ucont_P_Binary_Output(UserCtx *user)
{
  PetscViewer	viewer;
  char filen[80];

  int rank;
  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  sprintf(filen, "vfield%5.5d_%1.1d.dat", user->ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);  
  VecView(user->Ucont, viewer);  
  PetscViewerDestroy(&viewer);
  
  
  sprintf(filen, "ufield%5.5d_%1.1d.dat", user->ti, user->_this);  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);  
  VecView(user->Ucat, viewer);
  PetscViewerDestroy(&viewer);
  
  sprintf(filen, "pfield%5.5d_%1.1d.dat", user->ti, user->_this);  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  VecView(user->P, viewer);
  PetscViewerDestroy(&viewer);

  sprintf(filen, "nvfield%5.5d_%1.1d.dat", user->ti, user->_this);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  VecView(user->Nvert, viewer);
  PetscViewerDestroy(&viewer);
  
  if(averaging && user->ti!=0) {	// Seokkoo Kang
    sprintf(filen, "su0_%06d_%1d.dat", user->ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_sum, viewer);
    PetscViewerDestroy(&viewer);
    sprintf(filen, "su0_%06d_%1d.dat.info", user->ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(NULL);
    
    sprintf(filen, "su1_%06d_%1d.dat",  user->ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_cross_sum, viewer);
    PetscViewerDestroy(&viewer);  
    sprintf(filen, "su1_%06d_%1d.dat.info", user->ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(NULL);
    
    sprintf(filen, "su2_%06d_%1d.dat", user->ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_square_sum, viewer);
    PetscViewerDestroy(&viewer);
    sprintf(filen, "su2_%06d_%1d.dat.info",user->ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(NULL);
    
    sprintf(filen, "sp_%06d_%1d.dat",user->ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->P_sum, viewer);
    PetscViewerDestroy(&viewer);
    sprintf(filen, "sp_%06d_%1d.dat.info",user->ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(NULL);
  }
  
/*  if(les) {
    Vec Cs;
    
    VecDuplicate(user->P, &Cs);
    DMLocalToGlobalBegin(user->da, user->lCs, INSERT_VALUES, Cs);
    DMLocalToGlobalEnd(user->da, user->lCs, INSERT_VALUES, Cs);
    
    sprintf(filen, "cs_%06d_%1d.dat",  user->ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(Cs, viewer);
    PetscViewerDestroy(&viewer);
    
    PetscBarrier(NULL);
    VecDestroy(&Cs);
  }
  */
  if(rans) {
    sprintf(filen, "kfield%06d_%1d.dat", user->ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->K_Omega, viewer);
    PetscViewerDestroy(&viewer);
    
    PetscBarrier(NULL);
  }
  
  return 0;
}

PetscInt Gidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user)

{
  PetscInt nidx;
  DMDALocalInfo	info = user->info;

  PetscInt	mx = info.mx, my = info.my;
  
  AO ao;
  DMDAGetAO(user->da, &ao);
  nidx=i+j*mx+k*mx*my;
  
  AOApplicationToPetsc(ao,1,&nidx);
  
  return (nidx);
}


PetscErrorCode GridDivergence(UserCtx *user)
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
  PetscReal	***div, ***aj, ***nvert;
  Cmpnts	***ucont;
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

  PetscReal    norm;
  VecNorm(user->Vcont, NORM_INFINITY, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "Grid Flux norm  %le\n", norm);

  DMDAVecGetArray(fda,user->lVcont, &ucont);
  DMDAVecGetArray(da, user->lAj, &aj);
  VecDuplicate(user->P, &Div);
  DMDAVecGetArray(da, Div, &div);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	maxdiv = fabs((ucont[k][j][i].x - ucont[k][j][i-1].x +
		       ucont[k][j][i].y - ucont[k][j-1][i].y +
		       ucont[k][j][i].z - ucont[k-1][j][i].z)*aj[k][j][i]);
	if (nvert[k][j][i] + nvert[k+1][j][i] + nvert[k-1][j][i] +
	    nvert[k][j+1][i] + nvert[k][j-1][i] +
	    nvert[k][j][i+1] + nvert[k][j][i-1] > 0.1) maxdiv = 0.;
	div[k][j][i] = maxdiv;
      }
    }
  }

  if (zs==0) {
    k=0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ze == mz) {
    k=mz-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xs==0) {
    i=0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xe==mx) {
    i=mx-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0;
      }
    }
  }

  if (ys==0) {
    j=0;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ye==my) {
    j=my-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }
  DMDAVecRestoreArray(da, Div, &div);
  VecMax(Div, &i, &maxdiv);
  PetscPrintf(PETSC_COMM_WORLD, "Maxdiv Grid %d %d %e\n", user->ti, i, maxdiv);
  PetscInt mi;
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (mi=xs; mi<xe; mi++) {
	if (Gidx(mi,j,k,user) ==i) {
	  PetscPrintf(PETSC_COMM_SELF, "Maxdiv Grid location %d %d %d\n", mi,j, k);
	}
      }
    }
  }
    
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lVcont, &ucont);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  VecDestroy(&Div);
  return(0);
}

PetscErrorCode RestartApplication(UserCtx *user,IBMNodes *ibm, 
				  FSInfo *fsi, Cstart *cstart, 
				  PetscInt tistart)
{
  PetscInt	i, bi, ibi;   

  for (bi=0; bi<user[0].block_number; bi++) {
    user[bi].ti = tistart;

    Q_P_Binary_Input(&(user[bi]));
    //Q_P_VTK_Input(&(user[bi]));

    DMGlobalToLocalBegin(user[bi].da, user[bi].P,
			 INSERT_VALUES, user[bi].lP);
    DMGlobalToLocalEnd(user[bi].da, user[bi].P,
		       INSERT_VALUES, user[bi].lP);

    DMGlobalToLocalBegin(user[bi].da, user[bi].Nvert,
			 INSERT_VALUES, user[bi].lNvert);
    DMGlobalToLocalEnd(user[bi].da, user[bi].Nvert,
		       INSERT_VALUES, user[bi].lNvert);

    DMGlobalToLocalBegin(user[bi].cda, user[bi].Q,
			 INSERT_VALUES, user[bi].lQ);
    DMGlobalToLocalEnd(user[bi].cda, user[bi].Q,
		       INSERT_VALUES, user[bi].lQ);
    VecCopy(user[bi].Q, user[bi].Qold);
    VecCopy(user[bi].lQ, user[bi].lQold);

    if (immersed) {
      VecCopy(user[bi].Nvert, user[bi].Nvert_o);
      DMGlobalToLocalBegin(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
      DMGlobalToLocalEnd(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
    }

    if (rstart_fsi) {
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	
	FSI_DATA_Input(&fsi[ibi],ibi,tistart);
	
	if (movefsi) {
	  if (!moveframe)
	    Elmt_Move_FSI_TRANS(&fsi[ibi], &ibm[ibi],ibi);
	  for (i=0;i<6;i++){
	    fsi[ibi].S_realm1[i]=fsi[ibi].S_real[i];
	    fsi[ibi].S_real[i]=fsi[ibi].S_new[i];
	  }
	  for (i=0; i<ibm[ibi].n_v; i++) {
	    ibm[ibi].u[i].x = fsi[ibi].S_real[1];
	    ibm[ibi].u[i].y = fsi[ibi].S_real[3];
	    ibm[ibi].u[i].z = fsi[ibi].S_real[5];

	    ibm[ibi].uold[i].x = fsi[ibi].S_real[1];
	    ibm[ibi].uold[i].y = fsi[ibi].S_real[3];
	    ibm[ibi].uold[i].z = fsi[ibi].S_real[5];
	  }
	  for (i=0; i<ibm[ibi].n_v; i++) {
	    ibm[ibi].urm1[i].x = fsi[ibi].S_realm1[1];
	    ibm[ibi].urm1[i].y = fsi[ibi].S_realm1[3];
	    ibm[ibi].urm1[i].z = fsi[ibi].S_realm1[5];
	  }
	}

      } // ibi
    } // rstart_fsi


  }// bi
  return(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv)
{
  UserCtx	*user;

  PetscInt	i, bi, ibi;
  
  IBMNodes	*ibm;//, *ibm0, *ibm1;
  //  IBMInfo	*ibm_intp;
  IBMVNodes     *ibmv;
  // Added for fsi
  FSInfo        *fsi;
  PetscBool     DoSCLoop;
  PetscInt      itr_sc;
  Cstart        cstart;

  PetscInt      level;
  UserMG        usermg;
  PetscBool	flg;   // rstart_flg, flg; // Advait
  PetscInt      ti=0, tiout=10, ansys=0;//,tistart = 0;

  PetscInitialize(&argc, &argv, (char *)0, help);
  PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE);
  PetscOptionsGetInt(NULL, NULL, "-tio", &tiout, NULL);
  PetscOptionsGetInt(NULL, NULL, "-imm", &immersed, NULL);
  PetscOptionsGetInt(NULL, NULL, "-ibm_ansys", &ansys, NULL);
  
  PetscOptionsGetInt(NULL, NULL, "-wavy", &wavy, NULL);
  PetscOptionsGetInt(NULL, NULL, "-isothermal", &isothermal, NULL);
  PetscOptionsGetInt(NULL, NULL, "-inv", &invicid, NULL);
  PetscOptionsGetInt(NULL, NULL, "-rstart", &tistart, &rstart_flg);

  //  PetscOptionsGetInt(NULL, NULL, "-imp", &implicit, NULL);
  //  PetscOptionsGetInt(NULL, NULL, "-imp_MAX_IT", &imp_MAX_IT, NULL);
  PetscOptionsGetInt(NULL, NULL, "-fsi", &movefsi, NULL);
  PetscOptionsGetInt(NULL, NULL, "-rfsi", &rotatefsi, NULL);
  //  PetscOptionsGetInt(NULL, NULL, "-radi", &radi, NULL);
  //  PetscOptionsGetInt(NULL, NULL, "-inlet", &inletprofile, NULL);
  PetscOptionsGetInt(NULL, NULL, "-str", &STRONG_COUPLING, NULL);
  PetscOptionsGetInt(NULL, NULL, "-rs_fsi", &rstart_fsi, NULL);
  PetscOptionsGetInt(NULL, NULL, "-cop", &cop, NULL);
  PetscOptionsGetInt(NULL, NULL, "-fish", &fish, NULL);
  PetscOptionsGetInt(NULL, NULL, "-rheology", &rheology, NULL);
  PetscOptionsGetInt(NULL, NULL, "-duplicate", &duplicate, NULL);
  PetscOptionsGetInt(NULL, NULL, "-Pipe", &Pipe, NULL);
  PetscOptionsGetInt(NULL, NULL, "-turbine", &turbine, NULL);
  PetscOptionsGetInt(NULL, NULL, "-fishcyl", &fishcyl, NULL);
  PetscOptionsGetInt(NULL, NULL, "-eel", &eel, NULL);
  PetscOptionsGetInt(NULL, NULL, "-cstart", &fish_c, NULL);
  PetscOptionsGetInt(NULL, NULL, "-wing", &wing, NULL);
  PetscOptionsGetInt(NULL, NULL, "-sediment", &sediment, NULL);
  PetscOptionsGetInt(NULL, NULL, "-mhv", &MHV, NULL);
  PetscOptionsGetInt(NULL, NULL, "-hydro", &hydro, NULL);
  PetscOptionsGetInt(NULL, NULL, "-lv", &LV, NULL);
  PetscOptionsGetInt(NULL, NULL, "-reg", &regime, NULL);
  PetscOptionsGetInt(NULL, NULL, "-TwoD", &TwoD, NULL);
  PetscOptionsGetInt(NULL, NULL, "-thin", &thin, NULL);
  PetscOptionsGetInt(NULL, NULL, "-dgf_z", &dgf_z, NULL);
  PetscOptionsGetInt(NULL, NULL, "-dgf_y", &dgf_y, NULL);
  PetscOptionsGetInt(NULL, NULL, "-dgf_x", &dgf_x, NULL);
  PetscOptionsGetInt(NULL, NULL, "-dgf_az", &dgf_az, NULL);
  PetscOptionsGetInt(NULL, NULL, "-dgf_ay", &dgf_ay, NULL);
  PetscOptionsGetInt(NULL, NULL, "-dgf_ax", &dgf_ax, NULL);
  PetscOptionsGetInt(NULL, NULL, "-body", &NumberOfBodies, NULL);
  PetscOptionsGetInt(NULL, NULL, "-mframe", &moveframe, NULL);
  PetscOptionsGetInt(NULL, NULL, "-rframe", &rotateframe, NULL);
  //  PetscOptionsGetInt(NULL, NULL, "-blk", &blank, NULL);
  //  PetscOptionsGetInt(NULL, NULL, "-init1", &InitialGuessOne, NULL);

  PetscOptionsGetReal(NULL, NULL, "-x_c", &(CMx_c), NULL);
  PetscOptionsGetReal(NULL, NULL, "-y_c", &(CMy_c), NULL);
  PetscOptionsGetReal(NULL, NULL, "-z_c", &(CMz_c), NULL);

  /* PetscOptionsGetReal(NULL, "-imp_atol", &(imp_atol), NULL); */
  /* PetscOptionsGetReal(NULL, "-imp_rtol", &(imp_rtol), NULL); */
  /* PetscOptionsGetReal(NULL, "-imp_stol", &(imp_stol), NULL); */
  /* PetscOptionsGetReal(NULL, "-poisson_tol", &poisson_tol, NULL); */		// Seokkoo Kang: tolerance of implicit matrix free solver. 1.e-4 is enough for most cases.

  PetscOptionsGetInt(NULL, NULL, "-les", &les, NULL); //if 1 Smagorinsky with Cs=0.1, if >=2 Dynamic model
  PetscOptionsGetInt(NULL, NULL, "-rans", &rans, NULL);		       
  //  if (rans) wallfunction=1;
  PetscOptionsGetInt(NULL, NULL, "-wallfunction", &wallfunction, NULL);	//  1 Cabot or 2 Log Law, 3 slip wall
  PetscOptionsGetInt(NULL, NULL, "-mixed", &mixed, NULL);			// Seokkoo Kang: mixed model optiion for LES
  PetscOptionsGetInt(NULL, NULL, "-clark", &clark, NULL);			// Seokkoo Kang: mixed model opton for LES
  PetscOptionsGetInt(NULL, NULL, "-dynamic_freq", &dynamic_freq, NULL);		// Seokkoo Kang: LES dynamic compute frequency 
  PetscOptionsGetInt(NULL, NULL, "-averaging", &averaging, NULL);	// Seokkoo Kang: if 1 do averaging; always begin with -rstart 0

  PetscOptionsGetInt(NULL, NULL, "-grid1d", &grid1d, NULL);
  PetscOptionsGetInt(NULL, NULL, "-i_periodic", &i_periodic, NULL);	
  PetscOptionsGetInt(NULL, NULL, "-j_periodic", &j_periodic, NULL);	
  PetscOptionsGetInt(NULL, NULL, "-k_periodic", &k_periodic, NULL);	
 
  PetscOptionsGetInt(NULL, NULL, "-pbc_domain", &blkpbc, NULL);

  PetscOptionsGetInt(NULL, NULL, "-testfilter_ik", &testfilter_ik, NULL);
  PetscOptionsGetInt(NULL, NULL, "-testfilter_1d", &testfilter_1d, NULL);
  //  PetscOptionsGetInt(NULL, NULL, "-poisson", &poisson, NULL);
  PetscOptionsGetInt(NULL, NULL, "-i_homo_filter", &i_homo_filter, NULL);
  PetscOptionsGetInt(NULL, NULL, "-j_homo_filter", &j_homo_filter, NULL);
  PetscOptionsGetInt(NULL, NULL, "-k_homo_filter", &k_homo_filter, NULL);
  PetscOptionsGetReal(NULL, NULL, "-max_cs", &max_cs, NULL);
  
  PetscOptionsGetReal(NULL, NULL, "-St_exp", &(St_exp), NULL);
  PetscOptionsGetReal(NULL, NULL, "-wlngth", &(wavelength), NULL);
  //  PetscOptionsGetReal(NULL, NULL, "-Mach", &(Mach), NULL); // use -Ma option to set Mach number in FormInitialize in init.c
  PetscOptionsGetInt(NULL, NULL, "-central", &central, NULL); // ==2 second order, ==4 4th order 
  PetscOptionsGetInt(NULL, NULL, "-channel", &channel, NULL);
  PetscOptionsGetInt(NULL, NULL, "-calc_sigma", &calc_sigma, NULL); // 1 if switch function is on
 

    
  //  PetscOptionsGetInt(NULL, NULL, "-fluxorder", &fluxorder, NULL);
  //  PetscPrintf(PETSC_COMM_WORLD, "!!flux order %i\n", fluxorder);

  //  PetscPrintf(PETSC_COMM_WORLD, "tiout %d %le %le thin %d!\n",tiout, imp_atol,imp_rtol,thin);

  if (rstart_flg) 
    ti=tistart;
  
  
  if (immersed) {
    // PetscMalloc(NumberOfBodies*sizeof(IBMNodes), &ibm);
    ibm=(IBMNodes *)calloc(NumberOfBodies,sizeof(IBMNodes));
    // if (rheology) PetscMalloc(NumberOfBodies*sizeof(IBMVNodes), &ibmv);
    if (rheology) ibmv=(IBMVNodes *)calloc(NumberOfBodies,sizeof(IBMVNodes));
    // PetscMalloc(NumberOfBodies*sizeof(FSInfo), &fsi);
    fsi=(FSInfo *)calloc(NumberOfBodies,sizeof(FSInfo));
    for (ibi=0;ibi<NumberOfBodies;ibi++) 
      FsiInitialize(0, &fsi[ibi], ibi);
    //ibm_ptr=ibm;
  }

  MG_Initial(&usermg, ibm);
  
  level = usermg.mglevels-1;
  user = usermg.mgctx[level].user;
  if (immersed) {  
    for (bi=0; bi<user[0].block_number; bi++) {
      user[bi].ibmlist=(IBMList *)malloc(NumberOfBodies*sizeof(IBMList));
      //  PetscMalloc(NumberOfBodies*sizeof(IBMList), &(user[bi].ibmlist));
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	InitIBMList(&(user[bi].ibmlist[ibi]));
      }
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	if (ansys) 
	  ibm_read_Ansys(&ibm[ibi], ibi);
	else
	  ibm_read_ucd(&ibm[ibi], ibi);
	ibm_surface_VTKOut(&ibm[ibi],ibi,0);
      }
      if (fish) { 
	ibi=0;
	fish_init(&(user[bi].dt));
	
	fish_swim(&ibm[ibi],ti, user[bi].dt);
	//ibm_surface_VTKOut(&ibm[ibi],ibi,0);
	
      }

      for (ibi=0;ibi<NumberOfBodies;ibi++){
	ibm_search_advanced(&user[bi],&ibm[ibi], ibi);
    	PetscPrintf(PETSC_COMM_WORLD, "Search ibm %d\n",ibi);
      } 
    }    
  }

  // rstart based on options
  if (rstart_flg) {
    RestartApplication(user, ibm, fsi, &cstart, tistart);
  } // if rstart
  
  PetscInt tisteps = 5;
  PetscOptionsGetInt(NULL, NULL, "-totalsteps", &tisteps, &flg);
  if (tistart==0)
    for (bi=0; bi<user[0].block_number; bi++) {
      user[bi].ti=0;
      SetInitialCondition(&user[bi]);      
      Q_P_Binary_Output(&(user[bi]));
    }
  
  if (immersed) {   
    for (bi=0; bi<user[0].block_number; bi++) {
      for (ibi=0;ibi<NumberOfBodies;ibi++){
	ibm_search_advanced(&user[bi],&ibm[ibi], ibi);
	PetscPrintf(PETSC_COMM_WORLD, "Search ibm ibi %d\n",ibi);
      } 
    }    
  } 
  
  // Pass passive variables through User//
  for (bi=0; bi<user[0].block_number; bi++){
    user[bi].tiout=tiout;
    //user[bi].tistart = tistart;
    //user[bi].rstart_flg = rstart_flg;
  } 

  if (ti==tistart){ 
    tistart++;
    for (bi=0; bi<user[0].block_number; bi++){
      //user[bi].tistart = tistart;
    }
  }
  PetscReal     t_flow_s, t_flow_e, t_struc_e, t_struc_s, t_ti_s, t_ti_e;
  PetscTime(&t_ti_s);
  /* ===================================================================*/
  /*   physical time Step Loop */
  for (ti = tistart; ti<tistart + tisteps; ti++) {
    
    PetscPrintf(PETSC_COMM_WORLD, "Time %d\n", ti);
    for (bi=0; bi<user[0].block_number; bi++) 
      user[bi].ti=ti;
    /* =================================================================*/
    /*     Strong-Coupling (SC) Loop */
    DoSCLoop= PETSC_TRUE ; itr_sc = 0;
    while (DoSCLoop) {

      itr_sc++;
      PetscPrintf(PETSC_COMM_WORLD, "SC LOOP itr # %d\n", itr_sc);
      
      /*     Structral Solver! */
      if (immersed && !turbine) {
	//	PetscGetCPUTime(&t_struc_s);
	Struc_Solver(&usermg,ibm,fsi,&cstart, itr_sc,tistart, &DoSCLoop);
	/*       Struc_predictor(&usermg, ibm, fsi, itr_sc,tistart, &DoSCLoop); */
	// PetscGetCPUTime(&t_struc_e);
      }
      else
	DoSCLoop = PETSC_FALSE;
      
      
      /*     Flow Solver! */
      // PetscGetCPUTime(&t_flow_s);
      Flow_Solver(&usermg, ibm, fsi);
      // PetscGetCPUTime(&t_flow_e);  
      
    }// End of while SC loop
    /* ==================================================================================             */

/*  put the time accuracy coefficient back to 1.5
    after the 1st real-time step */
/*     COEF_TIME_ACCURACY=1.5; */


    PetscInt profiling=0; 
    int rank, size;
    PetscOptionsGetInt(NULL,NULL, "-profiling", &profiling, NULL);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (profiling && !rank){
      FILE *f;
      char filen[80];
      sprintf(filen, "CPU_TIME%4.4d.dat", size);
      f = fopen(filen, "a");
      PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le\n",ti, t_flow_e-t_flow_s, t_struc_e-t_struc_s, t_flow_e-t_flow_s+t_struc_e-t_struc_s);
      fclose(f);    
    }
   
  


/* ==================================================================================             */
/*     Save the old values (at ti) for later */
/*     Vec Error; */
/*     PetscReal error=0.0; */
 

    level = usermg.mglevels-1;
    user = usermg.mgctx[level].user;
    for (bi=0; bi<user[0].block_number; bi++) {

      if (immersed) {
	VecCopy(user[bi].Nvert, user[bi].Nvert_o);
	DMGlobalToLocalBegin(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
	DMGlobalToLocalEnd(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
      }

      // Copy Ucont to Ucont_o for the finest level
      VecCopy(user[bi].Q, user[bi].Qold);
           
      DMGlobalToLocalBegin(user[bi].cda, user[bi].Qold, INSERT_VALUES, user[bi].lQold);
      DMGlobalToLocalEnd(user[bi].cda, user[bi].Qold, INSERT_VALUES, user[bi].lQold);

      if (rans) {
       //		VecCopy(user[bi].K_Omega_o, user[bi].K_Omega_rm1);
       VecCopy(user[bi].K_Omega, user[bi].K_Omega_o);
       
       DMGlobalToLocalBegin(user[bi].fda2, user[bi].K_Omega_o, INSERT_VALUES, user[bi].lK_Omega_o);
       DMGlobalToLocalEnd(user[bi].fda2, user[bi].K_Omega_o, INSERT_VALUES, user[bi].lK_Omega_o);       
      }
    } // for bi

    
    if (immersed) {// && (movefsi || rotatefsi || cop || fish || MHV || LV || fish_c || sediment)){// && ti<tistart+3) {

      for (ibi=0;ibi<NumberOfBodies;ibi++) {

      for (i=0; i<ibm[ibi].n_v; i++) {
	ibm[ibi].x_bp_o[i] = ibm[ibi].x_bp[i];
	ibm[ibi].y_bp_o[i] = ibm[ibi].y_bp[i];
	ibm[ibi].z_bp_o[i] = ibm[ibi].z_bp[i];

	ibm[ibi].urm1[i].x = ibm[ibi].uold[i].x;
	ibm[ibi].urm1[i].y = ibm[ibi].uold[i].y;
	ibm[ibi].urm1[i].z = ibm[ibi].uold[i].z;

	ibm[ibi].uold[i].x = ibm[ibi].u[i].x;
	ibm[ibi].uold[i].y = ibm[ibi].u[i].y;
	ibm[ibi].uold[i].z = ibm[ibi].u[i].z;
      }

      for (i=0;i<6;i++){
	fsi[ibi].S_realm1[i]=fsi[ibi].S_real[i];
	fsi[ibi].S_real[i]=fsi[ibi].S_new[i];

	fsi[ibi].S_ang_rm1[i]=fsi[ibi].S_ang_r[i];
	fsi[ibi].S_ang_r[i]=fsi[ibi].S_ang_n[i];
      }

      for (i=0;i<3;i++){
	fsi[ibi].L_r[i]=fsi[ibi].L_n[i];
      }
      for (i=0;i<4;i++){
	fsi[ibi].q_r[i]=fsi[ibi].q[i];
      }

      fsi[ibi].F_x_real=fsi[ibi].F_x;
      fsi[ibi].F_y_real=fsi[ibi].F_y;
      fsi[ibi].F_z_real=fsi[ibi].F_z;
      
      fsi[ibi].M_x_rm3=fsi[ibi].M_x_rm2;
      fsi[ibi].M_y_rm3=fsi[ibi].M_y_rm2;
      fsi[ibi].M_z_rm3=fsi[ibi].M_z_rm2;

      fsi[ibi].M_x_rm2=fsi[ibi].M_x_real;
      fsi[ibi].M_y_rm2=fsi[ibi].M_y_real;
      fsi[ibi].M_z_rm2=fsi[ibi].M_z_real;

      fsi[ibi].M_x_real=fsi[ibi].M_x;
      fsi[ibi].M_y_real=fsi[ibi].M_y;
      fsi[ibi].M_z_real=fsi[ibi].M_z;



      } //ibi
    } // if immersed

/* ==================================================================================             */
    
  } // ti (physical time) loop
/* ==================================================================================             */
//  nextp:   PetscPrintf(PETSC_COMM_WORLD, "Number of iteration is  %d \n",ti);
  PetscTime(&t_ti_e);
  PetscPrintf(PETSC_COMM_WORLD, "cpu time %le\n",t_ti_e-t_ti_s);

  MG_Finalize(&usermg);
  PetscFinalize();

/* ==================================================================================             */
  return(0);

}

 /*  if (Z[ibi]==2){ */
/* 	      ibm_duplicate(&ibm[NumberOfBodies+particle_added-1],&ibm[ibi], 0.0, 0.0,-l_z); */
	     
/* 	      fsi[NumberOfBodies+particle_added-1].x_c=fsi[ibi].x_c; */
/* 	      fsi[NumberOfBodies+particle_added-1].y_c=fsi[ibi].y_c; */
/* 	      fsi[NumberOfBodies+particle_added-1].z_c=fsi[ibi].z_c-l_z; */

/* 	      fsi[NumberOfBodies+particle_added-1].a_c[0]= fsi[ibi].a_c[0]; */
/* 	      fsi[NumberOfBodies+particle_added-1].a_c[1]= fsi[ibi].a_c[1]; */
/* 	      fsi[NumberOfBodies+particle_added-1].a_c[2]= fsi[ibi].a_c[2]-l_z; */
/* 	    } */
/* 	    else if (Z[ibi]==1){ */
/* 	      ibm_duplicate(&ibm[NumberOfBodies+particle_added-1],&ibm[ibi], 0.0, 0.0,l_z); */
	     
/* 	      fsi[NumberOfBodies+particle_added-1].x_c=fsi[ibi].x_c; */
/* 	      fsi[NumberOfBodies+particle_added-1].y_c=fsi[ibi].y_c; */
/* 	      fsi[NumberOfBodies+particle_added-1].z_c=fsi[ibi].z_c+l_z; */

/* 	      fsi[NumberOfBodies+particle_added-1].a_c[0]= fsi[ibi].a_c[0]; */
/* 	      fsi[NumberOfBodies+particle_added-1].a_c[1]= fsi[ibi].a_c[1]; */
/* 	      fsi[NumberOfBodies+particle_added-1].a_c[2]= fsi[ibi].a_c[2]+l_z; */
/* 	    } */

/* 	    else if (X[ibi]==2){ */
/* 	      ibm_duplicate(&ibm[NumberOfBodies-1],&ibm[ibi], 0.0, 0.0,-l_x); */
	     
/* 	      fsi[NumberOfBodies+particle_added-1].x_c=fsi[ibi].x_c-l_x; */
/* 	      fsi[NumberOfBodies+particle_added-1].y_c=fsi[ibi].y_c; */
/* 	      fsi[NumberOfBodies+particle_added-1].z_c=fsi[ibi].z_c; */

/* 	      fsi[NumberOfBodies+particle_added-1].a_c[0]= fsi[ibi].a_c[0]-l_x; */
/* 	      fsi[NumberOfBodies+particle_added-1].a_c[1]= fsi[ibi].a_c[1]; */
/* 	      fsi[NumberOfBodies+particle_added-1].a_c[2]= fsi[ibi].a_c[2]; */
/* 	    } */
/* 	    else if (X[ibi]==1){ */
/* 	      ibm_duplicate(&ibm[NumberOfBodies-1],&ibm[ibi], 0.0, 0.0,l_x); */
	    
/* 	      fsi[NumberOfBodies+particle_added-1].x_c=fsi[ibi].x_c+l_x; */
/* 	      fsi[NumberOfBodies+particle_added-1].y_c=fsi[ibi].y_c; */
/* 	      fsi[NumberOfBodies+particle_added-1].z_c=fsi[ibi].z_c; */
	      
/* 	      fsi[NumberOfBodies+particle_added-1].a_c[0]= fsi[ibi].a_c[0]+l_x; */
/* 	      fsi[NumberOfBodies+particle_added-1].a_c[1]= fsi[ibi].a_c[1]; */
/* 	      fsi[NumberOfBodies+particle_added-1].a_c[2]= fsi[ibi].a_c[2]; */
	     
/* 	    } */
/* 	    else   PetscPrintf(PETSC_COMM_WORLD, "Particle duplication error!!!!!!!!!!!\n"); */
	   
/* 	    for (i=0;i<3;i++){ */
/* 	      fsi[NumberOfBodies+particle_added-1].pbc[i]=3; */
/* 	      fsi[NumberOfBodies+particle_added-1].S_ang_n[2*i+1]= fsi[ibi].S_ang_n[2*i+1]; */
/* 	      fsi[NumberOfBodies+particle_added-1].alpha[i]=fsi[ibi].alpha[i]; */
/* 	      fsi[NumberOfBodies+particle_added-1].acc[i]=fsi[ibi].acc[i]; */
/* 	      fsi[NumberOfBodies+particle_added-1].S_new[2*i+1]=fsi[ibi].S_new[2*i+1]; */
/* 	      fsi[NumberOfBodies+particle_added-1].S_real[2*i]=fsi[ibi].S_real[2*i]; */
/* 	      fsi[NumberOfBodies+particle_added-1].S_real[2*i+1]=fsi[ibi].S_real[2*i+1]; */
/* 	      fsi[NumberOfBodies+particle_added-1].L_n[i]=fsi[ibi].L_n[i]; */
/* 	    } */
	   

/* 	    fsi[NumberOfBodies+particle_added-1].mu_s=fsi[ibi].mu_s; */
/* 	    fsi[NumberOfBodies+particle_added-1].x_c=fsi[NumberOfBodies+particle_added-1].a_c[0]+fsi[NumberOfBodies+particle_added-1].S_new[0]; */
/* 	    fsi[NumberOfBodies+particle_added-1].y_c=fsi[NumberOfBodies+particle_added-1].a_c[1]+fsi[NumberOfBodies+particle_added-1].S_new[2]; */
/* 	    fsi[NumberOfBodies+particle_added-1].z_c=fsi[NumberOfBodies+particle_added-1].a_c[2]+fsi[NumberOfBodies+particle_added-1].S_new[4]; */
	 
/* 	    for (i=0;i<3;i++){ */
/* 	      for (j=0;j<3;j++){ */
/* 		fsi[NumberOfBodies+particle_added-1].R[i][j]=fsi[ibi].R[i][j]; */
/* 	      } */
/* 	    } */

/* 	    rotate_quarternion(&ibm[NumberOfBodies+particle_added-1],&fsi[NumberOfBodies+particle_added-1],user[bi].dt); */

/* 	    calc_ibm_normal(&ibm[NumberOfBodies+particle_added-1]); */
	   
/* 	    omega_c.x=fsi[NumberOfBodies+particle_added-1].S_ang_n[1]; */
/* 	    omega_c.y=fsi[NumberOfBodies+particle_added-1].S_ang_n[3]; */
/* 	    omega_c.z=fsi[NumberOfBodies+particle_added-1].S_ang_n[5]; */

/* 	    a_c.x=fsi[NumberOfBodies+particle_added-1].x_c; */
/* 	    a_c.y=fsi[NumberOfBodies+particle_added-1].y_c; */
/* 	    a_c.z=fsi[NumberOfBodies+particle_added-1].z_c; */

/* 	    for (i=0; i<ibm[NumberOfBodies+particle_added-1].n_v; i++) { */
	  
/* 	      rx = ibm[NumberOfBodies+particle_added-1].x_bp[i]-a_c.x; */
/* 	      ry = ibm[NumberOfBodies+particle_added-1].y_bp[i]-a_c.y; */
/* 	      rz = ibm[NumberOfBodies+particle_added-1].z_bp[i]-a_c.z; */
/* 	      ibm[NumberOfBodies+particle_added-1].u[i].x =   (rz*omega_c.y-omega_c.z*ry)+fsi[NumberOfBodies+particle_added-1].S_new[1]  ; */
/* 	      ibm[NumberOfBodies+particle_added-1].u[i].y =   (rx*omega_c.z-omega_c.x*rz)+fsi[NumberOfBodies+particle_added-1].S_new[3]  ; */
/* 	      ibm[NumberOfBodies+particle_added-1].u[i].z =   (ry*omega_c.x-omega_c.y*rx)+fsi[NumberOfBodies+particle_added-1].S_new[5]  ; */
	  
/* 	    } */
