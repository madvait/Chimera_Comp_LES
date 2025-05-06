#include "variables.h"
extern PetscInt invicid; 
extern PetscInt radi,fish,TwoD,cop, wing, rheology;
extern PetscReal St_exp,wavelength;
extern PetscInt tiout,STRONG_COUPLING,rheology,NumberOfBodies;
extern PetscInt dgf_z,dgf_y,dgf_x;

PetscErrorCode FsiInitialize(PetscInt n_elmt, FSInfo *fsi,PetscInt ibi)
{
  /*  Note: This routine shouldn't be called before ibm_read */
  PetscInt  i,j;

  for (i=0;i<6;i++) {
    fsi->S_old[i]=0.;
    fsi->S_new[i]=0.;
    fsi->S_realm1[i]=0.;
    fsi->S_real[i]=0.;

    fsi->S_ang_n[i]=0.;
    fsi->S_ang_o[i]=0.;
    fsi->S_ang_r[i]=0.;
    fsi->S_ang_rm1[i]=0.;
  }


  fsi->F_x_old = 0.;
  fsi->F_y_old = 0.;
  fsi->F_z_old = 0.;

  fsi->F_x_real = 0.;
  fsi->F_y_real = 0.;
  fsi->F_z_real = 0.;

  fsi->M_x_old = 0.;
  fsi->M_y_old = 0.;
  fsi->M_z_old = 0.;

  fsi->M_x = 0.;
  fsi->M_y = 0.;
  fsi->M_z = 0.;

  fsi->M_x_real = 0.;
  fsi->M_y_real = 0.;
  fsi->M_z_real = 0.;

  fsi->Mdpdn_x=0.;
  fsi->Mdpdn_y=0.;
  fsi->Mdpdn_z=0.;

  fsi->x_c=0.05;fsi->y_c=6.;fsi->z_c=15.;

  fsi->red_vel=0.52;//1.5;//0.3921;
  fsi->damp=.02;
  fsi->mu_s=500.;//0.025;//0.00568;

/* rheology parameters */
  fsi->clone=0.0;

  fsi->q[0]=1.000;
  fsi->q[1]=0.0;
  fsi->q[2]=0.0;
  fsi->q[3]=0.0;

  fsi->q_r[0]=1.000;
  fsi->q_r[1]=0.0;
  fsi->q_r[2]=0.0;
  fsi->q_r[3]=0.0;


  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      fsi-> R[i][j]=0.0;
      fsi->I_inv[i][j]=0.0;
    }
  }
  
  for (i=0;i<3;i++){
    fsi->alpha[i]=0.0;
    fsi->acc[i]=0.0;
    fsi->pbc[i]=0;
    fsi->L_n[i]=0.0;
    fsi->L_o[i]=0.0;
    fsi->L_r[i]=0.0;
  }

 
  fsi->Max_xbc= 1e23;fsi->Max_ybc= 1e23;fsi->Max_zbc= 1e23;
  fsi->Min_xbc=-1e23;fsi->Min_ybc=-1e23;fsi->Min_zbc=-1e23;
  

  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-red_vel", &(fsi->red_vel), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-damp", &(fsi->damp), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-mu_s", &(fsi->mu_s), PETSC_NULL);

  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-x_c", &(fsi->x_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-y_c", &(fsi->y_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-z_c", &(fsi->z_c), PETSC_NULL);

  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-Max_xbc", &(fsi->Max_xbc), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-Min_xbd", &(fsi->Min_xbc), PETSC_NULL);



  return(0);
}

/* ==================================================================================             */
PetscErrorCode FSI_DATA_Input(FSInfo *FSinf, PetscInt ibi, PetscInt ti)
{
  PetscInt  i;
  PetscReal t;

  FILE *f;
  char filen[80];  
  sprintf(filen, "DATA_FSI%5.5d_%2.2d.dat",ti, ibi);
  f = fopen(filen, "r");
  if (!f) {
    SETERRQ(PETSC_COMM_WORLD,1, "Cannot open FSI DATA file");
    PetscPrintf(PETSC_COMM_WORLD, "FSI_data cannot open file !!!!!!!!!!!!\n");
  }
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input begin %d %s\n",ti,filen);

  fscanf(f, "%le %le %le", &t, &t, &t);	  
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input red vel damp mu %le %le %le \n",FSinf->red_vel,FSinf->damp,FSinf->mu_s);
  fscanf(f, "%le %le %le", &(FSinf->x_c), &(FSinf->y_c), &(FSinf->z_c));	  
  fscanf(f, "%le %le %le \n", &(FSinf->F_x),&(FSinf->F_y), &(FSinf->F_z));	 
  fscanf(f, "%le %le %le \n", &(FSinf->M_x),&(FSinf->M_y), &(FSinf->M_z));	  
 

  for (i=0; i<6; i++) {
    fscanf(f, "%le %le %le %le", &(FSinf->S_new[i]),&(FSinf->S_old[i]), &(FSinf->S_real[i]), &(FSinf->S_realm1[i]));
    fscanf(f, "%le %le %le %le", &(FSinf->S_ang_n[i]),&(FSinf->S_ang_o[i]), &(FSinf->S_ang_r[i]), &(FSinf->S_ang_rm1[i]));
  }
  if (rheology){
    fscanf(f, "%le %le %le ", &FSinf->R[0][0],&FSinf->R[0][1],&FSinf->R[0][2]);
    fscanf(f, "%le %le %le ", &FSinf->R[1][0],&FSinf->R[1][1],&FSinf->R[1][2]);
    fscanf(f, "%le %le %le ", &FSinf->R[2][0],&FSinf->R[2][1],&FSinf->R[2][2]);
    fscanf(f, "%le %le %le ", &FSinf->L_n[0],&FSinf->L_n[1],&FSinf->L_n[2]);
  }
  fclose(f);
  if (rheology)  PetscPrintf(PETSC_COMM_WORLD,"%le %le %le %le %le %le \n",FSinf->R[0][0],FSinf->R[0][1],FSinf->R[0][2],FSinf->R[1][1],FSinf->R[1][2],FSinf->R[2][2]); 
  else{
    PetscPrintf(PETSC_COMM_WORLD, "FSI_data input z, dz/dt  %le %le %le %le\n",FSinf->S_new[4],FSinf->S_new[5],FSinf->red_vel,FSinf->damp);
    PetscPrintf(PETSC_COMM_WORLD, "FSI_data input y, dy/dt  %le %le %le %le\n",FSinf->S_new[2],FSinf->S_new[3],FSinf->red_vel,FSinf->damp);
  }
 
  return(0);
}
/* ==================================================================================             */

/* ==================================================================================             */

PetscErrorCode FSI_position_Output(FSInfo *FSinfo, PetscInt ibi, PetscInt ti)
{
  PetscInt rank, i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscBarrier(PETSC_NULL);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "FSI_position%2.2d",ibi);
    f = fopen(filen, "a");

    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le %le\n",ti, FSinfo->S_new[2],FSinfo->S_new[3],FSinfo->F_y, FSinfo->S_new[4],FSinfo->S_new[5],FSinfo->F_z,FSinfo->S_new[0],FSinfo->S_new[1],FSinfo->F_x, FSinfo->Power);
    fclose(f);

    sprintf(filen, "FSI_Agnle%2.2d",ibi);
    f = fopen(filen, "a");
    if (rheology)  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le \n",ti,FSinfo->S_ang_n[1],FSinfo->M_x,FSinfo->S_ang_n[3],FSinfo->S_ang_n[5],FSinfo->alpha[0],FSinfo->alpha[1],FSinfo->alpha[2]);
    else
      PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, FSinfo->S_ang_n[0],FSinfo->S_ang_n[1],FSinfo->M_x,FSinfo->S_ang_n[2],FSinfo->S_ang_n[3],FSinfo->S_ang_n[4],FSinfo->S_ang_n[5]);
    fclose(f);
  }
  return(0);
}


PetscErrorCode FSI_DATA_Output(FSInfo *FSinfo, PetscInt ibi, PetscInt ti)
{
  PetscInt rank, i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscBarrier(PETSC_NULL);
  if (!rank) {
    FILE *f;
    char filen[80];

    sprintf(filen, "DATA_FSI%5.5d_%2.2d.dat",ti, ibi);
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->red_vel, FSinfo->damp, FSinfo->mu_s);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->x_c, FSinfo->y_c, FSinfo->z_c);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->F_x, FSinfo->F_y, FSinfo->F_z);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->M_x, FSinfo->M_y, FSinfo->M_z);	  
    for (i=0; i<6; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_new[i],FSinfo->S_old[i], FSinfo->S_real[i], FSinfo->S_realm1[i]);
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_ang_n[i],FSinfo->S_ang_o[i], FSinfo->S_ang_r[i], FSinfo->S_ang_rm1[i]);
    }
    if (rheology){
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->R[0][0],FSinfo->R[0][1],FSinfo->R[0][2]);
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->R[1][0],FSinfo->R[1][1],FSinfo->R[1][2]);
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->R[2][0],FSinfo->R[2][1],FSinfo->R[2][2]);
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->L_n[0], FSinfo->L_n[1], FSinfo->L_n[2]);
    }	
    fclose(f);
  }
  return(0);
}

/* ==================================================================================             */
/*    Calculate forces on the immersed body */
/* ==================================================================================             */

PetscErrorCode Calc_forces_SI(UserCtx *user,FSInfo *fsi,
			       IBMNodes *ibm,PetscInt ti, 
			       PetscInt ibi, PetscInt bi)
{
  // This subroutine calculates forces on nvert==1 on the line between nv=1 and nv=3
  DM	da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz; 
  PetscInt      lxs, lxe, lys, lye, lzs, lze;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal     rei= 1./user->ren;
  Vec           Coor;
  PetscInt	i, j, k, elmt;
  PetscReal     sb;

  Cmpnts        ***coor, ***ucat;
  PetscReal     ***p, ***nvert;
  PetscReal     dwdz,dwdy,dwdx;
  PetscReal     dvdz,dvdy,dvdx;
  PetscReal     dudz,dudy,dudx;
  PetscReal     Txx,Tyy,Tzz;
  PetscReal     Tzy,Tzx,Tyx;
  PetscReal     nfx,nfy,nfz;
  PetscReal     nsx,nsy,nsz;
  PetscReal     ntx,nty,ntz;
  PetscReal     F_px,F_py,F_pz,Ap_x,Ap_y,Ap_z; //Forces and Area
  PetscReal     F_nx,F_ny,F_nz,An_x,An_y,An_z; 
  //PetscReal     Cp_x,Cp_y,Cp_z; //Pressure Forces
  //PetscReal     Cs_x,Cs_y,Cs_z; //Surface Forces
  PetscReal     Cp_nx,Cp_ny,Cp_nz; //Pressure Forces - side
  PetscReal     Cs_nx,Cs_ny,Cs_nz; //Surface Forces - side
  PetscReal     Cp_px,Cp_py,Cp_pz; //Pressure Forces + side
  PetscReal     Cs_px,Cs_py,Cs_pz; //Surface Forces + side
  PetscReal     F_xSum,F_ySum,F_zSum,A_totSum; //Surface Force
  PetscReal     F_pxSum,F_pySum,F_pzSum; //Surface Force
  PetscReal     F_nxSum,F_nySum,F_nzSum; //Surface Force
  PetscReal     Ap_xSum,Ap_ySum,Ap_zSum; // + side
  PetscReal     An_xSum,An_ySum,An_zSum; // - side
  PetscReal     A_xSum,A_ySum,A_zSum;
  PetscReal     Cp_xSum,Cp_ySum,Cp_zSum; //Pressure Force
  PetscReal     Cp_pxSum,Cp_pySum,Cp_pzSum; //Pressure Force
  PetscReal     Cp_nxSum,Cp_nySum,Cp_nzSum; //Pressure Force

  // Moments
  PetscReal      MF_px,MF_py,MF_pz; //Forces for Moment Calc
  PetscReal      MF_nx,MF_ny,MF_nz; //Forces for Moment Calc
  PetscReal      M_nx,M_ny,M_nz;   //Moments
  PetscReal      M_px,M_py,M_pz;   //Moments
  PetscReal      r_x,r_y,r_z;   //Anchor dist
  PetscReal      x, y, z;       //cell coord
  PetscReal      X_c,Y_c,Z_c;   //center of rotation coord
  PetscReal      M_xSum,M_ySum,M_zSum; //Surface Mom on all processors
  PetscReal      M_pxSum,M_pySum,M_pzSum; //Surface Mom on all processors
  PetscReal      M_nxSum,M_nySum,M_nzSum; //Surface Mom on all processors
  PetscReal      Iap_x,Iap_y,Iap_z;
  PetscReal      Ian_x,Ian_y,Ian_z;
  PetscReal      Iap_xSum,Iap_ySum,Iap_zSum; // + side
  PetscReal      Ian_xSum,Ian_ySum,Ian_zSum; // - side
  PetscReal      Ia_xSum,Ia_ySum,Ia_zSum;

  PetscReal      A_x,A_y,A_z;
  Cmpnts         ***csi,***eta,***zet;
  PetscReal      csi1,csi2,csi3;
  PetscReal      eta1,eta2,eta3;
  PetscReal      zet1,zet2,zet3;
  PetscReal      ***aj;//***iaj,***jaj,***kaj;

  PetscReal      MFdpdn_px,MFdpdn_py,MFdpdn_pz;
  PetscReal      Mdpdn_px,Mdpdn_py,Mdpdn_pz;
  PetscReal      Mdpdn_pxSum,Mdpdn_pySum,Mdpdn_pzSum;
  PetscReal      MFdpdn_nx,MFdpdn_ny,MFdpdn_nz;
  PetscReal      Mdpdn_nx,Mdpdn_ny,Mdpdn_nz;
  PetscReal      Mdpdn_nxSum,Mdpdn_nySum,Mdpdn_nzSum;
  PetscReal      Mdpdn_xSum,Mdpdn_ySum,Mdpdn_zSum; 
  PetscReal      val=1.9;

  IBMInfo        *ibminfo;
  IBMListNode    *current;

/* ==================================================================================             */
/*   Init var */
  F_px=0.;F_py=0.;F_pz=0.;
  F_nx=0.;F_ny=0.;F_nz=0.;
  Ap_x=0.;Ap_y=0.;Ap_z=0.;
  An_x=0.;An_y=0.;An_z=0.;
  Cp_px=0.;Cp_py=0.;Cp_pz=0.;
  Cs_px=0.;Cs_py=0.;Cs_pz=0.;
  Cp_nx=0.;Cp_ny=0.;Cp_nz=0.;
  Cs_nx=0.;Cs_ny=0.;Cs_nz=0.;

  M_px=0.;M_py=0.;M_pz=0.;
  M_nx=0.;M_ny=0.;M_nz=0.;
  Iap_x=0.;Iap_y=0.;Iap_z=0.;
  Ian_x=0.;Ian_y=0.;Ian_z=0.;

  Mdpdn_px=0.;Mdpdn_py=0.;Mdpdn_pz=0.;
  Mdpdn_nx=0.;Mdpdn_ny=0.;Mdpdn_nz=0.;
  MFdpdn_px=0.;MFdpdn_py=0.;MFdpdn_pz=0.;
  MFdpdn_nx=0.;MFdpdn_ny=0.;MFdpdn_nz=0.;

  /*   Check Later Y_c !!!!!!!!!!!!!!!!!!!!!!!! */
  X_c=fsi->x_c; Y_c=fsi->y_c; Z_c=fsi->z_c;

  PetscPrintf(PETSC_COMM_WORLD, "RE in calc_force  %le X_c %le %le %le \n",rei, X_c,Y_c,Z_c);

/* ==================================================================================             */
/*   Get Working arrays */
  DMGetCoordinatesLocal(da, &Coor);
  DMDAVecGetArray(fda, Coor, &coor);
  //  DMDAVecGetArray(fda, user->lCent, &cent);
  DMDAVecGetArray(fda, user->lUcat, &ucat);
  DMDAVecGetArray(da, user->lP, &p);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lAj, &aj);
  /* DMDAVecGetArray(da, user->lIAj, &iaj); */
  /* DMDAVecGetArray(da, user->lJAj, &jaj); */
  /* DMDAVecGetArray(da, user->lKAj, &kaj); */

/* ==================================================================================             */
/* Loop around all ibm nodes */
  current = user->ibmlist[ibi].head;
  while (current) {
    ibminfo = &current->ibm_intp;
    current = current->next;	
    i = ibminfo->ni; j= ibminfo->nj; k = ibminfo->nk;
    elmt = ibminfo->cell;
    sb = ibminfo->d_s;
    //sb = ibminfo->d_i;// |dn| from node to bndry


    // normal 
    nfx=ibm->nf_x[elmt];
    nfy=ibm->nf_y[elmt];
    nfz=ibm->nf_z[elmt];

    // 1st bi-normal of nf
    nsx=ibm->ns_x[elmt];
    nsy=ibm->ns_y[elmt];
    nsz=ibm->ns_z[elmt];

    // 2nd bi-normal of nf
    ntx=ibm->nt_x[elmt];
    nty=ibm->nt_y[elmt];
    ntz=ibm->nt_z[elmt];

    MF_px=0.;MF_py=0.;MF_pz=0.;
    MF_nx=0.;MF_ny=0.;MF_nz=0.;
    MFdpdn_px=0.;MFdpdn_py=0.;MFdpdn_pz=0.;
    MFdpdn_nx=0.;MFdpdn_ny=0.;MFdpdn_nz=0.;
    
   
    if (i>=xs && i<xe && j>=ys && j<ye && k>=zs && k<ze) {	    

/* ==================================================================================             */
/*       Shear 
Stress Tensor at ib node (nvert==1) (2nd & 1st order)  */
      if (!invicid) {
      if (nvert[k+1][j][i]<2.5 && nvert[k-1][j][i]<2.5 && k+1<mz-1 && k-1>0) {
	zet1 = 0.25*((zet[k][j][i].x+zet[k-1][j][i].x)*
		     ( aj[k][j][i]  + aj[k-1][j][i]));
	zet2 = 0.25*((zet[k][j][i].y+zet[k-1][j][i].y)*
		     ( aj[k][j][i]  + aj[k-1][j][i]));
	zet3 = 0.25*((zet[k][j][i].z+zet[k-1][j][i].z)*
		     ( aj[k][j][i]  + aj[k-1][j][i]));
	dwdz = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet3;
	dvdz = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet3;
	dudz = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet3;

	dwdy = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet2;
	dvdy = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet2;
	dudy = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet2;

	dwdx = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet1;
	dvdx = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet1;
	dudx = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet1;

      } else if (nvert[k+1][j][i]<2.5 && k+1<mz-1) {
	zet1 = (zet[k][j][i].x)*aj[k][j][i];
	zet2 = (zet[k][j][i].y)*aj[k][j][i];
	zet3 = (zet[k][j][i].z)*aj[k][j][i];

	dwdz = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet3;
	dvdz = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet3;
	dudz = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet3;

	dwdy = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet2;
	dvdy = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet2;
	dudy = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet2;

	dwdx = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet1;
	dvdx = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet1;
	dudx = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet1;

      } else if (nvert[k-1][j][i]<2.5 && k-1>0){
	zet1 = (zet[k-1][j][i].x)*aj[k-1][j][i];
	zet2 = (zet[k-1][j][i].y)*aj[k-1][j][i];
	zet3 = (zet[k-1][j][i].z)*aj[k-1][j][i];

	dwdz = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet3;
	dvdz = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet3;
	dudz = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet3;

	dwdy = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet2;
	dvdy = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet2;
	dudy = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet2;

	dwdx = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet1;
	dvdx = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet1;
	dudx = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet1;

      } else {
	dwdz = 0.;
	dvdz = 0.;
	dudz = 0.;

	dwdy = 0.;
	dvdy = 0.;
	dudy = 0.;

	dwdx = 0.;
	dvdx = 0.;
	dudx = 0.;
      }


      if (nvert[k][j+1][i]<2.5 && j+1<my-1 && nvert[k][j-1][i]<2.5 && j-1>0) {
	eta1 = 0.25*(eta[k][j][i].x+eta[k][j-1][i].x)*
	              (aj[k][j][i]+aj[k][j-1][i]);
	eta2 = 0.25*(eta[k][j][i].y+eta[k][j-1][i].y)*
	              (aj[k][j][i]+aj[k][j-1][i]);
	eta3 = 0.25*(eta[k][j][i].z+eta[k][j-1][i].z)*
	              (aj[k][j][i]+aj[k][j-1][i]);
	dwdz += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta3;
	dvdz += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta3;
	dudz += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta3;

	dwdy += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta2;
	dvdy += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta2;
	dudy += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta2;

	dwdx += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta1;
	dvdx += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta1;
	dudx += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta1;

      } else if (nvert[k][j+1][i]<2.5 && j+1<my-1) {
	eta1 = eta[k][j][i].x*aj[k][j][i];
	eta2 = eta[k][j][i].y*aj[k][j][i];
	eta3 = eta[k][j][i].z*aj[k][j][i];

	dwdz += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta3;
	dvdz += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta3;
	dudz += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta3;

	dwdy += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta2;
	dvdy += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta2;
	dudy += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta2;

	dwdx += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta1;
	dvdx += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta1;
	dudx += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta1;

      } else if (nvert[k][j-1][i]<2.5 && j-1>0){
	eta1 = eta[k][j-1][i].x*aj[k][j-1][i];
	eta2 = eta[k][j-1][i].y*aj[k][j-1][i];
	eta3 = eta[k][j-1][i].z*aj[k][j-1][i];

	dwdz += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta3;
	dvdz += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta3;
	dudz += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta3;

	dwdy += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta2;
	dvdy += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta2;
	dudy += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta2;

	dwdx += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta1;
	dvdx += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta1;
	dudx += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta1;
      } 

      if (nvert[k][j][i+1]<2.5 && i+1<mx-1 && nvert[k][j][i-1]<2.5 && i-1>0) {
	csi1 = 0.25*(csi[k][j][i].x+csi[k][j][i-1].x)*
	              (aj[k][j][i]+aj[k][j][i-1]);
	csi2 = 0.25*(csi[k][j][i].y+csi[k][j][i-1].y)*
	              (aj[k][j][i]+aj[k][j][i-1]);
	csi3 = 0.25*(csi[k][j][i].z+csi[k][j][i-1].z)*
	              (aj[k][j][i]+aj[k][j][i-1]);

	dwdz += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi3;
	dvdz += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi3;
	dudz += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi3;

	dwdy += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi2;
	dvdy += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi2;
	dudy += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi2;

	dwdx += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi1;
	dvdx += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi1;
	dudx += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi1;

      } else if (nvert[k][j][i+1]<2.5 && i+1<mx-1) {
	csi1 = csi[k][j][i].x*aj[k][j][i];
	csi2 = csi[k][j][i].y*aj[k][j][i];
	csi3 = csi[k][j][i].z*aj[k][j][i];

	dwdz += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi3;
	dvdz += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi3;
	dudz += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi3;

	dwdy += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi2;
	dvdy += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi2;
	dudy += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi2;

	dwdx += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi1;
	dvdx += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi1;
	dudx += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi1;

      } else if (nvert[k][j][i-1]<2.5 && i-1>0){
	csi1 = csi[k][j][i-1].x*aj[k][j][i-1];
	csi2 = csi[k][j][i-1].y*aj[k][j][i-1];
	csi3 = csi[k][j][i-1].z*aj[k][j][i-1];

	dwdz += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi3;
	dvdz += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi3;
	dudz += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi3;

	dwdy += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi2;
	dvdy += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi2;
	dudy += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi2;

	dwdx += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi1;
	dvdx += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi1;
	dudx += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi1;
      }
 
      Tzz = rei * (dwdz + dwdz);
      Tyy = rei * (dvdy + dvdy);
      Txx = rei * (dudx + dudx);
      Tzy = rei * (dwdy + dvdz);
      Tzx = rei * (dwdx + dudz);
      Tyx = rei * (dvdx + dudy);
      } else { // if inviscid
	Tzz = 0.;
	Tyy = 0.;
	Txx = 0.;
	Tzy = 0.;
	Tzx = 0.;
	Tyx = 0.;
      }
/* ==================================================================================             */
/*       Pressure Force */

      /* if (i==1 && j==100 && k==119) { */
      /* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d P %le %le %le %le\n", i, j, k, p[k][j][i],p[k+1][j][i],p[k-1][j][i],p[k-2][j][i]); */
      /* 	//	flagprint =1; */
      /* } */

      //dr = sqrt(r_x*r_x+r_y*r_y+r_z*r_z);

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

      // k+ calculate Fz at k+1/2   
      if (k<mz-1)
      if (nvert[k+1][j][i]>val ){
	A_x=-0.5*(zet[k][j][i].x+zet[k+1][j][i].x);//zet[k  ][j][i].x;
	A_y=-0.5*(zet[k][j][i].y+zet[k+1][j][i].y);//zet[k  ][j][i].y;
	A_z=-0.5*(zet[k][j][i].z+zet[k+1][j][i].z);//zet[k  ][j][i].z;

	x = 0.5*(coor[k][j][i].x+coor[k+1][j][i].x);
	y = 0.5*(coor[k][j][i].y+coor[k+1][j][i].y);
	z = 0.5*(coor[k][j][i].z+coor[k+1][j][i].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

/* 	Cp_pz += -p[k+1][j][i]*A_z; */
/* 	MF_pz += -p[k+1][j][i]*A_z; */
	Cp_pz += -p[k][j][i]*A_z;
	MF_pz += -p[k][j][i]*A_z;
	Ap_z  += A_z;
	Iap_x -=   r_y*A_z;
	Iap_y -= - r_x*A_z;
	MFdpdn_pz += sb*(nfy*r_z-nfz*r_y)*A_z;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

/* 	Cp_py += -p[k+1][j][i]*A_y; */
/* 	MF_py += -p[k+1][j][i]*A_y; */
	Cp_py += -p[k][j][i]*A_y;
	MF_py += -p[k][j][i]*A_y;
	Ap_y  += A_y;
	Iap_x -=  - r_z*A_y;
	Iap_z -=    r_x*A_y;
	MFdpdn_py += sb*(nfy*r_z-nfz*r_y)*A_y;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

/* 	Cp_px += -p[k+1][j][i]*A_x; */
/* 	MF_px += -p[k+1][j][i]*A_x; */
	Cp_px += -p[k][j][i]*A_x;
	MF_px += -p[k][j][i]*A_x;
	Ap_x  += A_x;
	Iap_z -=   - r_y*A_x;
	Iap_y -= -(- r_z*A_x);
	MFdpdn_px += sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_pz += Tzz*A_z ;
	Cs_py += Tzy*A_z ;
	Cs_px += Tzx*A_z ;
	  
	MF_pz += Tzz*A_z ;
	MF_py += Tzy*A_z ;
	MF_px += Tzx*A_z ;

	Cs_pz += Tzy*A_y ;
	Cs_py += Tyy*A_y ;
	Cs_px += Tyx*A_y ;
	  
	MF_pz += Tzy*A_y ;
	MF_py += Tyy*A_y ;
	MF_px += Tyx*A_y ;

	Cs_pz += Tzx*A_x;
	Cs_py += Tyx*A_x;
	Cs_px += Txx*A_x;
	  
	MF_pz += Tzx*A_x;
	MF_py += Tyx*A_x;
	MF_px += Txx*A_x;

	M_px -=   r_y*MF_pz - r_z*MF_py;//
	M_py -= -(r_x*MF_pz - r_z*MF_px);
	M_pz -=   r_x*MF_py - r_y*MF_px;
      }   

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;
   
      //k- calculate Fz at k-1/2
      if (k>0)
      if (nvert[k-1][j][i]>val) {
	A_x=-0.5*(zet[k][j][i].x+zet[k-1][j][i].x);//zet[k-1][j][i].x;
	A_y=-0.5*(zet[k][j][i].y+zet[k-1][j][i].y);//zet[k-1][j][i].y;
	A_z=-0.5*(zet[k][j][i].z+zet[k-1][j][i].z);//zet[k-1][j][i].z;

	x = 0.5*(coor[k][j][i].x+coor[k-1][j][i].x);
	y = 0.5*(coor[k][j][i].y+coor[k-1][j][i].y);
	z = 0.5*(coor[k][j][i].z+coor[k-1][j][i].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

/* 	Cp_nz +=  p[k-1][j][i]*A_z; */
/* 	MF_nz +=  p[k-1][j][i]*A_z; */
	Cp_nz +=  p[k][j][i]*A_z;
	MF_nz +=  p[k][j][i]*A_z;
	An_z  += A_z;
	Ian_x -=   r_y*A_z;
	Ian_y -= - r_x*A_z;
	MFdpdn_nz += -sb*(nfy*r_z-nfz*r_y)*A_z;

/* 	Cp_ny +=  p[k-1][j][i]*A_y; */
/* 	MF_ny +=  p[k-1][j][i]*A_y; */
	Cp_ny +=  p[k][j][i]*A_y;
	MF_ny +=  p[k][j][i]*A_y;
	An_y  += A_y;
	Ian_x -=  - r_z*A_y;
	Ian_z -=    r_x*A_y;
	MFdpdn_ny += -sb*(nfy*r_z-nfz*r_y)*A_y;

/* 	Cp_nx +=  p[k-1][j][i]*A_x; */
/* 	MF_nx +=  p[k-1][j][i]*A_x; */
	Cp_nx +=  p[k][j][i]*A_x;
	MF_nx +=  p[k][j][i]*A_x;
	An_x  += A_x;
	Ian_z -=   - r_y*A_x;
	Ian_y -= -(- r_z*A_x);
	MFdpdn_nx += -sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_nz -= Tzz*A_z ;
	Cs_ny -= Tzy*A_z ;
	Cs_nx -= Tzx*A_z ;
	
	MF_nz -= Tzz*A_z ;
	MF_ny -= Tzy*A_z ;
	MF_nx -= Tzx*A_z ;
	
	Cs_nz -= Tzy*A_y ;
	Cs_ny -= Tyy*A_y ;
	Cs_nx -= Tyx*A_y ;
	
	MF_nz -= Tzy*A_y ;
	MF_ny -= Tyy*A_y ;
	MF_nx -= Tyx*A_y ;
	
	Cs_nz -= Tzx*A_x;
	Cs_ny -= Tyx*A_x;
	Cs_nx -= Txx*A_x;
	
	MF_nz -= Tzx*A_x;
	MF_ny -= Tyx*A_x;
	MF_nx -= Txx*A_x;

	M_nx -=   r_y*MF_nz - r_z*MF_ny;//
	M_ny -= -(r_x*MF_nz - r_z*MF_nx);
	M_nz -=   r_x*MF_ny - r_y*MF_nx;
      }
     
      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

      //j+
      if (j< my-1)
      if (nvert[k][j+1][i]>val ){
	A_x=-0.5*(eta[k][j][i].x+eta[k][j+1][i].x);//eta[k][j  ][i].x;
	A_y=-0.5*(eta[k][j][i].y+eta[k][j+1][i].y);//eta[k][j  ][i].y;
	A_z=-0.5*(eta[k][j][i].z+eta[k][j+1][i].z);//eta[k][j  ][i].z;

	x = 0.5*(coor[k][j][i].x+coor[k][j+1][i].x);
	y = 0.5*(coor[k][j][i].y+coor[k][j+1][i].y);
	z = 0.5*(coor[k][j][i].z+coor[k][j+1][i].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

/* 	Cp_pz += -p[k][j+1][i]*A_z; */
/* 	MF_pz += -p[k][j+1][i]*A_z; */
	Cp_pz += -p[k][j][i]*A_z;
	MF_pz += -p[k][j][i]*A_z;
	Ap_z  += A_z;
	Iap_x -=   r_y*A_z;
	Iap_y -= - r_x*A_z;
	MFdpdn_pz += sb*(nfy*r_z-nfz*r_y)*A_z;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

/* 	Cp_py += -p[k][j+1][i]*A_y; */
/* 	MF_py += -p[k][j+1][i]*A_y; */
	Cp_py += -p[k][j][i]*A_y;
	MF_py += -p[k][j][i]*A_y;
	Ap_y  += A_y;
	Iap_x -=  - r_z*A_y;
	Iap_z -=    r_x*A_y;
	MFdpdn_py += sb*(nfy*r_z-nfz*r_y)*A_y;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

/* 	Cp_px += -p[k][j+1][i]*A_x; */
/* 	MF_px += -p[k][j+1][i]*A_x; */
	Cp_px += -p[k][j][i]*A_x;
	MF_px += -p[k][j][i]*A_x;
	Ap_x  += A_x;
	Iap_z -=   - r_y*A_x;
	Iap_y -= -(- r_z*A_x);
	MFdpdn_px += sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_pz += Tzz*A_z ;
	Cs_py += Tzy*A_z ;
	Cs_px += Tzx*A_z ;
	  
	MF_pz += Tzz*A_z ;
	MF_py += Tzy*A_z ;
	MF_px += Tzx*A_z ;

	Cs_pz += Tzy*A_y ;
	Cs_py += Tyy*A_y ;
	Cs_px += Tyx*A_y ;
	  
	MF_pz += Tzy*A_y ;
	MF_py += Tyy*A_y ;
	MF_px += Tyx*A_y ;

	Cs_pz += Tzx*A_x;
	Cs_py += Tyx*A_x;
	Cs_px += Txx*A_x;
	  
	MF_pz += Tzx*A_x;
	MF_py += Tyx*A_x;
	MF_px += Txx*A_x;

	M_px -=   r_y*MF_pz - r_z*MF_py;//
	M_py -= -(r_x*MF_pz - r_z*MF_px);
	M_pz -=   r_x*MF_py - r_y*MF_px;

      }

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

      //j-
      if (j>0)
      if (nvert[k][j-1][i]>val) {
	A_x=-0.5*(eta[k][j][i].x+eta[k][j-1][i].x);//eta[k][j-1][i].x;
	A_y=-0.5*(eta[k][j][i].y+eta[k][j-1][i].y);//eta[k][j-1][i].y;
	A_z=-0.5*(eta[k][j][i].z+eta[k][j-1][i].z);//eta[k][j-1][i].z;

	x = 0.5*(coor[k][j][i].x+coor[k][j-1][i].x);
	y = 0.5*(coor[k][j][i].y+coor[k][j-1][i].y);
	z = 0.5*(coor[k][j][i].z+coor[k][j-1][i].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

/* 	Cp_nz +=  p[k][j-1][i]*A_z; */
/* 	MF_nz +=  p[k][j-1][i]*A_z; */
	Cp_nz +=  p[k][j][i]*A_z;
	MF_nz +=  p[k][j][i]*A_z;
	An_z  += A_z;
	Ian_x -=   r_y*A_z;
	Ian_y -= - r_x*A_z;
	MFdpdn_nz += -sb*(nfy*r_z-nfz*r_y)*A_z;

/* 	Cp_ny +=  p[k][j-1][i]*A_y; */
/* 	MF_ny +=  p[k][j-1][i]*A_y; */
	Cp_ny +=  p[k][j][i]*A_y;
	MF_ny +=  p[k][j][i]*A_y;	     
	An_y  += A_y;
	Ian_x -=  - r_z*A_y;
	Ian_z -=    r_x*A_y;
	MFdpdn_ny += -sb*(nfy*r_z-nfz*r_y)*A_y;

/* 	Cp_nx +=  p[k][j-1][i]*A_x; */
/* 	MF_nx +=  p[k][j-1][i]*A_x; */
	Cp_nx +=  p[k][j][i]*A_x;
	MF_nx +=  p[k][j][i]*A_x;
	An_x  += A_x;
	Ian_z -=   - r_y*A_x;
	Ian_y -= -(- r_z*A_x);
	MFdpdn_nx += -sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_nz -= Tzz*A_z ;
	Cs_ny -= Tzy*A_z ;
	Cs_nx -= Tzx*A_z ;
	
	MF_nz -= Tzz*A_z ;
	MF_ny -= Tzy*A_z ;
	MF_nx -= Tzx*A_z ;
	
	Cs_nz -= Tzy*A_y ;
	Cs_ny -= Tyy*A_y ;
	Cs_nx -= Tyx*A_y ;
	
	MF_nz -= Tzy*A_y ;
	MF_ny -= Tyy*A_y ;
	MF_nx -= Tyx*A_y ;
	
	Cs_nz -= Tzx*A_x;
	Cs_ny -= Tyx*A_x;
	Cs_nx -= Txx*A_x;
	
	MF_nz -= Tzx*A_x;
	MF_ny -= Tyx*A_x;
	MF_nx -= Txx*A_x;

	M_nx -=   r_y*MF_nz - r_z*MF_ny;//
	M_ny -= -(r_x*MF_nz - r_z*MF_nx);
	M_nz -=   r_x*MF_ny - r_y*MF_nx;
      }

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

      //i+
      if (i<mx-1)
      if (nvert[k][j][i+1]>val){
	A_x=-0.5*(csi[k][j][i].x+csi[k][j][i+1].x);//csi[k][j][i].x;
	A_y=-0.5*(csi[k][j][i].y+csi[k][j][i+1].y);//csi[k][j][i].y;
	A_z=-0.5*(csi[k][j][i].z+csi[k][j][i+1].z);//csi[k][j][i].z;

	x = 0.5*(coor[k][j][i].x+coor[k][j][i+1].x);
	y = 0.5*(coor[k][j][i].y+coor[k][j][i+1].y);
	z = 0.5*(coor[k][j][i].z+coor[k][j][i+1].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

/* 	Cp_pz += -p[k][j][i+1]*A_z; */
/* 	MF_pz += -p[k][j][i+1]*A_z; */
	Cp_pz += -p[k][j][i]*A_z;
	MF_pz += -p[k][j][i]*A_z;
	Ap_z  += A_z;
	Iap_x -=   r_y*A_z;
	Iap_y -= - r_x*A_z;
	MFdpdn_pz += sb*(nfy*r_z-nfz*r_y)*A_z;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

/* 	Cp_py += -p[k][j][i+1]*A_y; */
/* 	MF_py += -p[k][j][i+1]*A_y; */
	Cp_py += -p[k][j][i]*A_y;
	MF_py += -p[k][j][i]*A_y;
	Ap_y  += A_y;
	Iap_x -=  - r_z*A_y;
	Iap_z -=    r_x*A_y;
	MFdpdn_py += sb*(nfy*r_z-nfz*r_y)*A_y;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

/* 	Cp_px += -p[k][j][i+1]*A_x; */
/* 	MF_px += -p[k][j][i+1]*A_x; */
	Cp_px += -p[k][j][i]*A_x;
	MF_px += -p[k][j][i]*A_x;
	Ap_x  += A_x;
	Iap_z -=   - r_y*A_x;
	Iap_y -= -(- r_z*A_x);
	MFdpdn_px += sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_pz += Tzz*A_z ;
	Cs_py += Tzy*A_z ;
	Cs_px += Tzx*A_z ;
	  
	MF_pz += Tzz*A_z ;
	MF_py += Tzy*A_z ;
	MF_px += Tzx*A_z ;

	Cs_pz += Tzy*A_y ;
	Cs_py += Tyy*A_y ;
	Cs_px += Tyx*A_y ;
	  
	MF_pz += Tzy*A_y ;
	MF_py += Tyy*A_y ;
	MF_px += Tyx*A_y ;

	Cs_pz += Tzx*A_x;
	Cs_py += Tyx*A_x;
	Cs_px += Txx*A_x;
	  
	MF_pz += Tzx*A_x;
	MF_py += Tyx*A_x;
	MF_px += Txx*A_x;

	M_px -=   r_y*MF_pz - r_z*MF_py;//
	M_py -= -(r_x*MF_pz - r_z*MF_px);
	M_pz -=   r_x*MF_py - r_y*MF_px;
      }

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

      //i-
      if (i>0)
      if  (nvert[k][j][i-1]>val) {
	A_x=-0.5*(csi[k][j][i].x+csi[k][j][i-1].x);//csi[k][j][i-1].x;
	A_y=-0.5*(csi[k][j][i].y+csi[k][j][i-1].y);//csi[k][j][i-1].y;
	A_z=-0.5*(csi[k][j][i].z+csi[k][j][i-1].z);//csi[k][j][i-1].z;

	x = 0.5*(coor[k][j][i].x+coor[k][j][i-1].x);
	y = 0.5*(coor[k][j][i].y+coor[k][j][i-1].y);
	z = 0.5*(coor[k][j][i].z+coor[k][j][i-1].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

/* 	Cp_nz +=  p[k][j][i-1]*A_z; */
/* 	MF_nz +=  p[k][j][i-1]*A_z; */
	Cp_nz +=  p[k][j][i]*A_z;
	MF_nz +=  p[k][j][i]*A_z;
	An_z  += A_z;
	Ian_x -=   r_y*A_z;
	Ian_y -= - r_x*A_z;
	MFdpdn_nz += -sb*(nfy*r_z-nfz*r_y)*A_z;

/* 	Cp_ny +=  p[k][j][i-1]*A_y; */
/* 	MF_ny +=  p[k][j][i-1]*A_y; */
	Cp_ny +=  p[k][j][i]*A_y;
	MF_ny +=  p[k][j][i]*A_y;	     
	An_y  += A_y;
	Ian_x -=  - r_z*A_y;
	Ian_z -=    r_x*A_y;
	MFdpdn_ny += -sb*(nfy*r_z-nfz*r_y)*A_y;

/* 	Cp_nx +=  p[k][j][i-1]*A_x; */
/* 	MF_nx +=  p[k][j][i-1]*A_x; */
	Cp_nx +=  p[k][j][i]*A_x;
	MF_nx +=  p[k][j][i]*A_x;
	An_x  += A_x;
	Ian_z -=   - r_y*A_x;
	Ian_y -= -(- r_z*A_x);
	MFdpdn_nx += -sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_nz -= Tzz*A_z ;
	Cs_ny -= Tzy*A_z ;
	Cs_nx -= Tzx*A_z ;
	
	MF_nz -= Tzz*A_z ;
	MF_ny -= Tzy*A_z ;
	MF_nx -= Tzx*A_z ;
	
	Cs_nz -= Tzy*A_y ;
	Cs_ny -= Tyy*A_y ;
	Cs_nx -= Tyx*A_y ;
	
	MF_nz -= Tzy*A_y ;
	MF_ny -= Tyy*A_y ;
	MF_nx -= Tyx*A_y ;
	
	Cs_nz -= Tzx*A_x;
	Cs_ny -= Tyx*A_x;
	Cs_nx -= Txx*A_x;
	
	MF_nz -= Tzx*A_x;
	MF_ny -= Tyx*A_x;
	MF_nx -= Txx*A_x;
	
	M_nx -=   r_y*MF_nz - r_z*MF_ny;//
	M_ny -= -(r_x*MF_nz - r_z*MF_nx);
	M_nz -=   r_x*MF_ny - r_y*MF_nx;

      }
    
/*       M_px -=   r_y*MF_pz - r_z*MF_py;// */
/*       M_py -= -(r_x*MF_pz - r_z*MF_px); */
/*       M_pz -=   r_x*MF_py - r_y*MF_px; */

/*       M_nx -=   r_y*MF_nz - r_z*MF_ny;// */
/*       M_ny -= -(r_x*MF_nz - r_z*MF_nx); */
/*       M_nz -=   r_x*MF_ny - r_y*MF_nx; */

      Mdpdn_px -=   r_y*MFdpdn_pz - r_z*MFdpdn_py;//
      Mdpdn_py -= -(r_x*MFdpdn_pz - r_z*MFdpdn_px);
      Mdpdn_pz -=   r_x*MFdpdn_py - r_y*MFdpdn_px;

      Mdpdn_nx -=   r_y*MFdpdn_nz - r_z*MFdpdn_ny;//
      Mdpdn_ny -= -(r_x*MFdpdn_nz - r_z*MFdpdn_nx);
      Mdpdn_nz -=   r_x*MFdpdn_ny - r_y*MFdpdn_nx;
      
/*       if (ibi==1 & i==2) */
/* 	PetscPrintf(PETSC_COMM_SELF, "MF %d %d %le %le %le %le %le %le %le %le\n",j,k,MF_py,MF_ny,MF_pz,MF_nz,M_px,M_nx,r_y,r_z); */

    }
/* ==================================================================================             */
/*     End of Loop ibm nodes */
  }

/* ==================================================================================             */
/*   Total Force on each processor */
  F_px = Cp_px + Cs_px; 
  F_py = Cp_py + Cs_py;
  F_pz = Cp_pz + Cs_pz;

  F_nx = Cp_nx + Cs_nx; 
  F_ny = Cp_ny + Cs_ny;
  F_nz = Cp_nz + Cs_nz;

/* ==================================================================================             */  
/*   Global Sum */
  MPI_Allreduce(&F_px, &F_pxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&F_py, &F_pySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&F_pz, &F_pzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&F_nx, &F_nxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&F_ny, &F_nySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&F_nz, &F_nzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Ap_x, &Ap_xSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Ap_y, &Ap_ySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Ap_z, &Ap_zSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&An_x, &An_xSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&An_y, &An_ySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&An_z, &An_zSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Cp_nx, &Cp_nxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Cp_ny, &Cp_nySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Cp_nz, &Cp_nzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Cp_px, &Cp_pxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Cp_py, &Cp_pySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Cp_pz, &Cp_pzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&M_px, &M_pxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&M_py, &M_pySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&M_pz, &M_pzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&M_nx, &M_nxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&M_ny, &M_nySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&M_nz, &M_nzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Iap_x, &Iap_xSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Iap_y, &Iap_ySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Iap_z, &Iap_zSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Ian_x, &Ian_xSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Ian_y, &Ian_ySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Ian_z, &Ian_zSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Mdpdn_px, &Mdpdn_pxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Mdpdn_py, &Mdpdn_pySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Mdpdn_pz, &Mdpdn_pzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Mdpdn_nx, &Mdpdn_nxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Mdpdn_ny, &Mdpdn_nySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Mdpdn_nz, &Mdpdn_nzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  
/* ==================================================================================             */
/*   Scale Check later !!!!! */

  A_xSum = 0.5 * (Ap_xSum + An_xSum);
  A_ySum = 0.5 * (Ap_ySum + An_ySum);
  A_zSum = 0.5 * (Ap_zSum + An_zSum);

/*   if (fabs(Ap_xSum)>1e-6) */
/*     F_pxSum=F_pxSum*A_xSum/Ap_xSum; */
/*   if (fabs(Ap_ySum)>1e-6) */
/*     F_pySum=F_pySum*A_ySum/Ap_ySum; */
/*   if (fabs(Ap_zSum)>1e-6) */
/*     F_pzSum=F_pzSum*A_zSum/Ap_zSum; */

/*   if (fabs(An_xSum)>1e-6) */
/*     F_nxSum=F_nxSum*A_xSum/An_xSum; */
/*   if (fabs(An_ySum)>1e-6) */
/*     F_nySum=F_nySum*A_ySum/An_ySum; */
/*   if (fabs(An_zSum)>1e-6) */
/*     F_nzSum=F_nzSum*A_zSum/An_zSum; */

  F_xSum = F_pxSum + F_nxSum;
  F_ySum = F_pySum + F_nySum;
  F_zSum = F_pzSum + F_nzSum;

  //if (!fish && !cop) {
  /* if (fabs(A_xSum)>1e-6) */
  /*   F_xSum=F_xSum/A_xSum*2.; */
  /* if (fabs(A_ySum)>1e-6) */
  /*   F_ySum=F_ySum/A_ySum*2.; */
  /* if (fabs(A_zSum)>1e-6) */
  /*   F_zSum=F_zSum/A_zSum*2.; */
  //  }

/*   if (fabs(Ap_xSum)>1e-6) */
/*     Cp_pxSum=Cp_pxSum*A_xSum/Ap_xSum; */
/*   if (fabs(Ap_ySum)>1e-6) */
/*     Cp_pySum=Cp_pySum*A_ySum/Ap_ySum; */
/*   if (fabs(Ap_zSum)>1e-6) */
/*     Cp_pzSum=Cp_pzSum*A_zSum/Ap_zSum; */

/*   if (fabs(An_xSum)>1e-6) */
/*     Cp_nxSum=Cp_nxSum*A_xSum/An_xSum; */
/*   if (fabs(An_ySum)>1e-6) */
/*     Cp_nySum=Cp_nySum*A_ySum/An_ySum; */
/*   if (fabs(An_zSum)>1e-6) */
/*     Cp_nzSum=Cp_nzSum*A_zSum/An_zSum; */

  Cp_xSum = Cp_pxSum + Cp_nxSum;
  Cp_ySum = Cp_pySum + Cp_nySum;
  Cp_zSum = Cp_pzSum + Cp_nzSum;

  /* if (fabs(A_xSum)>1e-6) */
  /*   Cp_xSum=Cp_xSum/A_xSum*2.; */
  /* if (fabs(A_ySum)>1e-6) */
  /*   Cp_ySum=Cp_ySum/A_ySum*2.; */
  /* if (fabs(A_zSum)>1e-6) */
  /*   Cp_zSum=Cp_zSum/A_zSum*2.; */    

  Ia_xSum = 0.5 * (Iap_xSum - Ian_xSum);
  Ia_ySum = 0.5 * (Iap_ySum - Ian_ySum);
  Ia_zSum = 0.5 * (Iap_zSum - Ian_zSum);

/*   if (fabs(Iap_xSum)>1e-6) */
/*     M_pxSum=M_pxSum*Ia_xSum/Iap_xSum; */
/*   if (fabs(Iap_ySum)>1e-6) */
/*     M_pySum=M_pySum*Ia_ySum/Iap_ySum; */
/*   if (fabs(Iap_zSum)>1e-6) */
/*     M_pzSum=M_pzSum*Ia_zSum/Iap_zSum ; */

/*   if (fabs(Ian_xSum)>1e-6) */
/*     M_nxSum=M_nxSum*Ia_xSum/Ian_xSum; */
/*   if (fabs(Ian_ySum)>1e-6) */
/*     M_nySum=M_nySum*Ia_ySum/Ian_ySum; */
/*   if (fabs(Ian_zSum)>1e-6) */
/*     M_nzSum=M_nzSum*Ia_zSum/Ian_zSum; */

  M_xSum = M_pxSum + M_nxSum;
  M_ySum = M_pySum + M_nySum;
  M_zSum = M_pzSum + M_nzSum;

/*   if (fabs(Iap_xSum)>1e-6) */
/*     Mdpdn_pxSum=Mdpdn_pxSum*Ia_xSum/Iap_xSum; */
/*   if (fabs(Iap_ySum)>1e-6) */
/*     Mdpdn_pySum=Mdpdn_pySum*Ia_ySum/Iap_ySum; */
/*   if (fabs(Iap_zSum)>1e-6) */
/*     Mdpdn_pzSum=Mdpdn_pzSum*Ia_zSum/Iap_zSum ; */

/*   if (fabs(Ian_xSum)>1e-6) */
/*     Mdpdn_nxSum=Mdpdn_nxSum*Ia_xSum/Ian_xSum; */
/*   if (fabs(Ian_ySum)>1e-6) */
/*     Mdpdn_nySum=Mdpdn_nySum*Ia_ySum/Ian_ySum; */
/*   if (fabs(Ian_zSum)>1e-6) */
/*     Mdpdn_nzSum=Mdpdn_nzSum*Ia_zSum/Ian_zSum; */

  Mdpdn_xSum = Mdpdn_pxSum + Mdpdn_nxSum;
  Mdpdn_ySum = Mdpdn_pySum + Mdpdn_nySum;
  Mdpdn_zSum = Mdpdn_pzSum + Mdpdn_nzSum;

/*   if (fabs(Ap_zSum)>1e-6) */
/*     M_xSum=M_xSum/A_zSum*2.; */
/*   if (fabs(Ap_ySum)>1e-6) */
/*     M_ySum=M_ySum/A_ySum*2.; */
/*   if (fabs(Ap_zSum)>1e-6) */
/*     M_zSum=M_zSum/A_zSum*2.; */

  A_totSum = Ap_xSum + Ap_ySum + Ap_zSum;

/* ==================================================================================             */
/*   store results in fsi */
  fsi->F_x = F_xSum; fsi->F_y = F_ySum; fsi->F_z = F_zSum;
  fsi->A_tot = A_totSum;
  fsi->M_x = M_xSum; fsi->M_y = M_ySum; fsi->M_z = M_zSum;
  fsi->Mdpdn_x = Mdpdn_xSum; fsi->Mdpdn_y = Mdpdn_ySum; fsi->Mdpdn_z = Mdpdn_zSum;

/* ==================================================================================             */
/*   output values */
  PetscPrintf(PETSC_COMM_WORLD, "F_x,F_y,F_z SI, %le %le %le Az %le %le Ay %le %le\n",F_xSum,F_ySum,F_zSum,Ap_zSum,An_zSum,Ap_ySum,An_ySum);
  PetscPrintf(PETSC_COMM_WORLD, "M_x,M_y,M_z SI, %le %le %le Ia_x %le %le Ip_y %le %le\n",M_xSum,M_ySum,M_zSum,Iap_xSum,Ian_xSum,Iap_ySum,Ian_ySum);
  PetscPrintf(PETSC_COMM_WORLD, "Mdpdn_x,Mdpdn_y,Mdpdn_z SI, %le %le %le\n",Mdpdn_xSum,Mdpdn_ySum,Mdpdn_zSum);

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Force_Coeff_SI%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le\n",ti, F_xSum, F_ySum, F_zSum,Cp_xSum,Cp_ySum,Cp_zSum, A_xSum,A_ySum,A_zSum);
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le %le\n",ti,F_ySum, F_zSum,Cp_ySum,Cp_zSum,F_ySum-Cp_ySum,F_zSum-Cp_zSum,A_ySum,F_xSum, M_xSum); */
    fclose(f);

    sprintf(filen, "Momt_Coeff_SI%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le\n",ti, M_xSum, M_ySum, M_zSum,fsi->M_x_old,fsi->M_x_real,fsi->Mdpdn_x);
    fclose(f);

  }

/* ==================================================================================             */
/*   Restore Working arrays */
  //DMDAVecRestoreArray(fda, user->lCent, &cent);
  DMDAVecRestoreArray(fda, Coor, &coor);
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
  DMDAVecRestoreArray(da, user->lP, &p);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  /* DMDAVecRestoreArray(da, user->lIAj, &iaj); */
  /* DMDAVecRestoreArray(da, user->lJAj, &jaj); */
  /* DMDAVecRestoreArray(da, user->lKAj, &kaj); */
  // VecDestroy(&Coor);
  return(0);
}

PetscErrorCode Calc_forces_SI2(UserCtx *user,FSInfo *fsi,
			       IBMNodes *ibm,PetscInt ti, 
			       PetscInt ibi, PetscInt bi)
{
  // This subroutine calculates forces on nvert==1
  DM	da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz; 
  PetscInt      lxs, lxe, lys, lye, lzs, lze;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal     rei= 1./user->ren;
  Vec           Coor;
  PetscInt	i, j, k, elmt;
  PetscReal     sb;

  Cmpnts        ***coor, ***ucat;
  PetscReal     ***p, ***nvert;
  PetscReal     dwdz,dwdy,dwdx;
  PetscReal     dvdz,dvdy,dvdx;
  PetscReal     dudz,dudy,dudx;
  PetscReal     Txx,Tyy,Tzz;
  PetscReal     Tzy,Tzx,Tyx;
  PetscReal     nfx,nfy,nfz;
  PetscReal     nsx,nsy,nsz;
  PetscReal     ntx,nty,ntz;
  PetscReal     F_px,F_py,F_pz,Ap_x,Ap_y,Ap_z; //Forces and Area
  PetscReal     F_nx,F_ny,F_nz,An_x,An_y,An_z; 
  //PetscReal     Cp_x,Cp_y,Cp_z; //Pressure Forces
  //PetscReal     Cs_x,Cs_y,Cs_z; //Surface Forces
  PetscReal     Cp_nx,Cp_ny,Cp_nz; //Pressure Forces - side
  PetscReal     Cs_nx,Cs_ny,Cs_nz; //Surface Forces - side
  PetscReal     Cp_px,Cp_py,Cp_pz; //Pressure Forces + side
  PetscReal     Cs_px,Cs_py,Cs_pz; //Surface Forces + side
  PetscReal     F_xSum,F_ySum,F_zSum,A_totSum; //Surface Force
  PetscReal     F_pxSum,F_pySum,F_pzSum; //Surface Force
  PetscReal     F_nxSum,F_nySum,F_nzSum; //Surface Force
  PetscReal     Ap_xSum,Ap_ySum,Ap_zSum; // + side
  PetscReal     An_xSum,An_ySum,An_zSum; // - side
  PetscReal     A_xSum,A_ySum,A_zSum;
  PetscReal     Cp_xSum,Cp_ySum,Cp_zSum; //Pressure Force
  PetscReal     Cp_pxSum,Cp_pySum,Cp_pzSum; //Pressure Force
  PetscReal     Cp_nxSum,Cp_nySum,Cp_nzSum; //Pressure Force

  // Moments
  PetscReal      MF_px,MF_py,MF_pz; //Forces for Moment Calc
  PetscReal      MF_nx,MF_ny,MF_nz; //Forces for Moment Calc
  PetscReal      M_nx,M_ny,M_nz;   //Moments
  PetscReal      M_px,M_py,M_pz;   //Moments
  PetscReal      r_x,r_y,r_z;   //Anchor dist
  PetscReal      x, y, z;       //cell coord
  PetscReal      X_c,Y_c,Z_c;   //center of rotation coord
  PetscReal      M_xSum,M_ySum,M_zSum; //Surface Mom on all processors
  PetscReal      M_pxSum,M_pySum,M_pzSum; //Surface Mom on all processors
  PetscReal      M_nxSum,M_nySum,M_nzSum; //Surface Mom on all processors
  PetscReal      Iap_x,Iap_y,Iap_z;
  PetscReal      Ian_x,Ian_y,Ian_z;
  PetscReal      Iap_xSum,Iap_ySum,Iap_zSum; // + side
  PetscReal      Ian_xSum,Ian_ySum,Ian_zSum; // - side
  PetscReal      Ia_xSum,Ia_ySum,Ia_zSum;

  PetscReal      A_x,A_y,A_z;
  Cmpnts         ***csi,***eta,***zet;
  PetscReal      csi1,csi2,csi3;
  PetscReal      eta1,eta2,eta3;
  PetscReal      zet1,zet2,zet3;
  PetscReal      ***aj;//***iaj,***jaj,***kaj;

  PetscReal      MFdpdn_px,MFdpdn_py,MFdpdn_pz;
  PetscReal      Mdpdn_px,Mdpdn_py,Mdpdn_pz;
  PetscReal      Mdpdn_pxSum,Mdpdn_pySum,Mdpdn_pzSum;
  PetscReal      MFdpdn_nx,MFdpdn_ny,MFdpdn_nz;
  PetscReal      Mdpdn_nx,Mdpdn_ny,Mdpdn_nz;
  PetscReal      Mdpdn_nxSum,Mdpdn_nySum,Mdpdn_nzSum;
  PetscReal      Mdpdn_xSum,Mdpdn_ySum,Mdpdn_zSum; 

  IBMInfo        *ibminfo;
  IBMListNode    *current;

/* ==================================================================================             */
/*   Init var */
  F_px=0.;F_py=0.;F_pz=0.;
  F_nx=0.;F_ny=0.;F_nz=0.;
  Ap_x=0.;Ap_y=0.;Ap_z=0.;
  An_x=0.;An_y=0.;An_z=0.;
  Cp_px=0.;Cp_py=0.;Cp_pz=0.;
  Cs_px=0.;Cs_py=0.;Cs_pz=0.;
  Cp_nx=0.;Cp_ny=0.;Cp_nz=0.;
  Cs_nx=0.;Cs_ny=0.;Cs_nz=0.;

  M_px=0.;M_py=0.;M_pz=0.;
  M_nx=0.;M_ny=0.;M_nz=0.;
  Iap_x=0.;Iap_y=0.;Iap_z=0.;
  Ian_x=0.;Ian_y=0.;Ian_z=0.;

  Mdpdn_px=0.;Mdpdn_py=0.;Mdpdn_pz=0.;
  Mdpdn_nx=0.;Mdpdn_ny=0.;Mdpdn_nz=0.;
  MFdpdn_px=0.;MFdpdn_py=0.;MFdpdn_pz=0.;
  MFdpdn_nx=0.;MFdpdn_ny=0.;MFdpdn_nz=0.;

  /*   Check Later Y_c !!!!!!!!!!!!!!!!!!!!!!!! */
  X_c=fsi->x_c; Y_c=fsi->y_c; Z_c=fsi->z_c;

  PetscPrintf(PETSC_COMM_WORLD, "RE in calc_force  %le X_c %le %le %le \n",rei, X_c,Y_c,Z_c);

/* ==================================================================================             */
/*   Get Working arrays */
  DMGetCoordinatesLocal(da, &Coor);
  DMDAVecGetArray(fda, Coor, &coor);
  //  DMDAVecGetArray(fda, user->lCent, &cent);
  DMDAVecGetArray(fda, user->lUcat, &ucat);
  DMDAVecGetArray(da, user->lP, &p);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lAj, &aj);
  /* DMDAVecGetArray(da, user->lIAj, &iaj); */
  /* DMDAVecGetArray(da, user->lJAj, &jaj); */
  /* DMDAVecGetArray(da, user->lKAj, &kaj); */

/* ==================================================================================             */
/* Loop around all ibm nodes */
  current = user->ibmlist[ibi].head;
  while (current) {
    ibminfo = &current->ibm_intp;
    current = current->next;	
    i = ibminfo->ni; j= ibminfo->nj; k = ibminfo->nk;
    elmt = ibminfo->cell;
    sb = ibminfo->d_s;
    //sb = ibminfo->d_i;// |dn| from node to bndry


    // normal 
    nfx=ibm->nf_x[elmt];
    nfy=ibm->nf_y[elmt];
    nfz=ibm->nf_z[elmt];

    // 1st bi-normal of nf
    nsx=ibm->ns_x[elmt];
    nsy=ibm->ns_y[elmt];
    nsz=ibm->ns_z[elmt];

    // 2nd bi-normal of nf
    ntx=ibm->nt_x[elmt];
    nty=ibm->nt_y[elmt];
    ntz=ibm->nt_z[elmt];
    if (i>=lxs && i<lxe && j>=ys && j<lye && k>=lzs && k<lze) {	        

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;
      MFdpdn_px=0.;MFdpdn_py=0.;MFdpdn_pz=0.;
      MFdpdn_nx=0.;MFdpdn_ny=0.;MFdpdn_nz=0.;


/* ==================================================================================             */
/*       Shear 
Stress Tensor at ib node (nvert==1) (2nd & 1st order)  */
      if (!invicid) {
      if (nvert[k+1][j][i]<2.5 && nvert[k-1][j][i]<2.5 && k+1<mz-1 && k-1>0) {
	zet1 = 0.25*((zet[k][j][i].x+zet[k-1][j][i].x)*
		     ( aj[k][j][i]  + aj[k-1][j][i]));
	zet2 = 0.25*((zet[k][j][i].y+zet[k-1][j][i].y)*
		     ( aj[k][j][i]  + aj[k-1][j][i]));
	zet3 = 0.25*((zet[k][j][i].z+zet[k-1][j][i].z)*
		     ( aj[k][j][i]  + aj[k-1][j][i]));
	dwdz = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet3;
	dvdz = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet3;
	dudz = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet3;

	dwdy = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet2;
	dvdy = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet2;
	dudy = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet2;

	dwdx = (ucat[k+1][j][i].z - ucat[k-1][j][i].z)/2.*zet1;
	dvdx = (ucat[k+1][j][i].y - ucat[k-1][j][i].y)/2.*zet1;
	dudx = (ucat[k+1][j][i].x - ucat[k-1][j][i].x)/2.*zet1;

      } else if (nvert[k+1][j][i]<2.5 && k+1<mz-1) {
	zet1 = (zet[k][j][i].x)*aj[k][j][i];
	zet2 = (zet[k][j][i].y)*aj[k][j][i];
	zet3 = (zet[k][j][i].z)*aj[k][j][i];

	dwdz = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet3;
	dvdz = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet3;
	dudz = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet3;

	dwdy = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet2;
	dvdy = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet2;
	dudy = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet2;

	dwdx = (ucat[k+1][j][i].z - ucat[k][j][i].z)*zet1;
	dvdx = (ucat[k+1][j][i].y - ucat[k][j][i].y)*zet1;
	dudx = (ucat[k+1][j][i].x - ucat[k][j][i].x)*zet1;

      } else if (nvert[k-1][j][i]<2.5 && k-1>0){
	zet1 = (zet[k-1][j][i].x)*aj[k-1][j][i];
	zet2 = (zet[k-1][j][i].y)*aj[k-1][j][i];
	zet3 = (zet[k-1][j][i].z)*aj[k-1][j][i];

	dwdz = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet3;
	dvdz = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet3;
	dudz = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet3;

	dwdy = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet2;
	dvdy = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet2;
	dudy = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet2;

	dwdx = (ucat[k][j][i].z - ucat[k-1][j][i].z)*zet1;
	dvdx = (ucat[k][j][i].y - ucat[k-1][j][i].y)*zet1;
	dudx = (ucat[k][j][i].x - ucat[k-1][j][i].x)*zet1;

      } else {
	dwdz = 0.;
	dvdz = 0.;
	dudz = 0.;

	dwdy = 0.;
	dvdy = 0.;
	dudy = 0.;

	dwdx = 0.;
	dvdx = 0.;
	dudx = 0.;
      }


      if (nvert[k][j+1][i]<2.5 && j+1<my-1 && nvert[k][j-1][i]<2.5 && j-1>0) {
	eta1 = 0.25*(eta[k][j][i].x+eta[k][j-1][i].x)*
	              (aj[k][j][i]+aj[k][j-1][i]);
	eta2 = 0.25*(eta[k][j][i].y+eta[k][j-1][i].y)*
	              (aj[k][j][i]+aj[k][j-1][i]);
	eta3 = 0.25*(eta[k][j][i].z+eta[k][j-1][i].z)*
	              (aj[k][j][i]+aj[k][j-1][i]);
	dwdz += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta3;
	dvdz += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta3;
	dudz += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta3;

	dwdy += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta2;
	dvdy += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta2;
	dudy += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta2;

	dwdx += (ucat[k][j+1][i].z - ucat[k][j-1][i].z)/2.*eta1;
	dvdx += (ucat[k][j+1][i].y - ucat[k][j-1][i].y)/2.*eta1;
	dudx += (ucat[k][j+1][i].x - ucat[k][j-1][i].x)/2.*eta1;

      } else if (nvert[k][j+1][i]<2.5 && j+1<my-1) {
	eta1 = eta[k][j][i].x*aj[k][j][i];
	eta2 = eta[k][j][i].y*aj[k][j][i];
	eta3 = eta[k][j][i].z*aj[k][j][i];

	dwdz += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta3;
	dvdz += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta3;
	dudz += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta3;

	dwdy += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta2;
	dvdy += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta2;
	dudy += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta2;

	dwdx += (ucat[k][j+1][i].z - ucat[k][j][i].z)*eta1;
	dvdx += (ucat[k][j+1][i].y - ucat[k][j][i].y)*eta1;
	dudx += (ucat[k][j+1][i].x - ucat[k][j][i].x)*eta1;

      } else if (nvert[k][j-1][i]<2.5 && j-1>0){
	eta1 = eta[k][j-1][i].x*aj[k][j-1][i];
	eta2 = eta[k][j-1][i].y*aj[k][j-1][i];
	eta3 = eta[k][j-1][i].z*aj[k][j-1][i];

	dwdz += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta3;
	dvdz += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta3;
	dudz += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta3;

	dwdy += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta2;
	dvdy += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta2;
	dudy += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta2;

	dwdx += (ucat[k][j][i].z - ucat[k][j-1][i].z)*eta1;
	dvdx += (ucat[k][j][i].y - ucat[k][j-1][i].y)*eta1;
	dudx += (ucat[k][j][i].x - ucat[k][j-1][i].x)*eta1;
      } 

      if (nvert[k][j][i+1]<2.5 && i+1<mx-1 && nvert[k][j][i-1]<2.5 && i-1>0) {
	csi1 = 0.25*(csi[k][j][i].x+csi[k][j][i-1].x)*
	              (aj[k][j][i]+aj[k][j][i-1]);
	csi2 = 0.25*(csi[k][j][i].y+csi[k][j][i-1].y)*
	              (aj[k][j][i]+aj[k][j][i-1]);
	csi3 = 0.25*(csi[k][j][i].z+csi[k][j][i-1].z)*
	              (aj[k][j][i]+aj[k][j][i-1]);

	dwdz += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi3;
	dvdz += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi3;
	dudz += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi3;

	dwdy += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi2;
	dvdy += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi2;
	dudy += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi2;

	dwdx += (ucat[k][j][i+1].z - ucat[k][j][i-1].z)/2.*csi1;
	dvdx += (ucat[k][j][i+1].y - ucat[k][j][i-1].y)/2.*csi1;
	dudx += (ucat[k][j][i+1].x - ucat[k][j][i-1].x)/2.*csi1;

      } else if (nvert[k][j][i+1]<2.5 && i+1<mx-1) {
	csi1 = csi[k][j][i].x*aj[k][j][i];
	csi2 = csi[k][j][i].y*aj[k][j][i];
	csi3 = csi[k][j][i].z*aj[k][j][i];

	dwdz += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi3;
	dvdz += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi3;
	dudz += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi3;

	dwdy += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi2;
	dvdy += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi2;
	dudy += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi2;

	dwdx += (ucat[k][j][i+1].z - ucat[k][j][i].z)*csi1;
	dvdx += (ucat[k][j][i+1].y - ucat[k][j][i].y)*csi1;
	dudx += (ucat[k][j][i+1].x - ucat[k][j][i].x)*csi1;

      } else if (nvert[k][j][i-1]<2.5 && i-1>0){
	csi1 = csi[k][j][i-1].x*aj[k][j][i-1];
	csi2 = csi[k][j][i-1].y*aj[k][j][i-1];
	csi3 = csi[k][j][i-1].z*aj[k][j][i-1];

	dwdz += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi3;
	dvdz += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi3;
	dudz += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi3;

	dwdy += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi2;
	dvdy += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi2;
	dudy += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi2;

	dwdx += (ucat[k][j][i].z - ucat[k][j][i-1].z)*csi1;
	dvdx += (ucat[k][j][i].y - ucat[k][j][i-1].y)*csi1;
	dudx += (ucat[k][j][i].x - ucat[k][j][i-1].x)*csi1;
      }
 
      Tzz = rei * (dwdz + dwdz);
      Tyy = rei * (dvdy + dvdy);
      Txx = rei * (dudx + dudx);
      Tzy = rei * (dwdy + dvdz);
      Tzx = rei * (dwdx + dudz);
      Tyx = rei * (dvdx + dudy);
      } else { // if inviscid
	Tzz = 0.;
	Tyy = 0.;
	Txx = 0.;
	Tzy = 0.;
	Tzx = 0.;
	Tyx = 0.;
      }
/* ==================================================================================             */
/*       Pressure Force */

      /* if (i==1 && j==100 && k==119) { */
      /* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d P %le %le %le %le\n", i, j, k, p[k][j][i],p[k+1][j][i],p[k-1][j][i],p[k-2][j][i]); */
      /* 	//	flagprint =1; */
      /* } */

      //dr = sqrt(r_x*r_x+r_y*r_y+r_z*r_z);

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

      // k+ calculate Fz at k+1/2      
      if (nvert[k+1][j][i]<0.9 ){
	A_x=0.5*(zet[k][j][i].x+zet[k+1][j][i].x);//zet[k  ][j][i].x;
	A_y=0.5*(zet[k][j][i].y+zet[k+1][j][i].y);//zet[k  ][j][i].y;
	A_z=0.5*(zet[k][j][i].z+zet[k+1][j][i].z);//zet[k  ][j][i].z;

	x = 0.5*(coor[k][j][i].x+coor[k+1][j][i].x);
	y = 0.5*(coor[k][j][i].y+coor[k+1][j][i].y);
	z = 0.5*(coor[k][j][i].z+coor[k+1][j][i].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

/* 	Cp_pz += -p[k+1][j][i]*A_z; */
/* 	MF_pz += -p[k+1][j][i]*A_z; */
	Cp_pz += -p[k][j][i]*A_z;
	MF_pz += -p[k][j][i]*A_z;
	Ap_z  += A_z;
	Iap_x -=   r_y*A_z;
	Iap_y -= - r_x*A_z;
	MFdpdn_pz += sb*(nfy*r_z-nfz*r_y)*A_z;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

/* 	Cp_py += -p[k+1][j][i]*A_y; */
/* 	MF_py += -p[k+1][j][i]*A_y; */
	Cp_py += -p[k][j][i]*A_y;
	MF_py += -p[k][j][i]*A_y;
	Ap_y  += A_y;
	Iap_x -=  - r_z*A_y;
	Iap_z -=    r_x*A_y;
	MFdpdn_py += sb*(nfy*r_z-nfz*r_y)*A_y;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

/* 	Cp_px += -p[k+1][j][i]*A_x; */
/* 	MF_px += -p[k+1][j][i]*A_x; */
	Cp_px += -p[k][j][i]*A_x;
	MF_px += -p[k][j][i]*A_x;
	Ap_x  += A_x;
	Iap_z -=   - r_y*A_x;
	Iap_y -= -(- r_z*A_x);
	MFdpdn_px += sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_pz += Tzz*A_z ;
	Cs_py += Tzy*A_z ;
	Cs_px += Tzx*A_z ;
	  
	MF_pz += Tzz*A_z ;
	MF_py += Tzy*A_z ;
	MF_px += Tzx*A_z ;

	Cs_pz += Tzy*A_y ;
	Cs_py += Tyy*A_y ;
	Cs_px += Tyx*A_y ;
	  
	MF_pz += Tzy*A_y ;
	MF_py += Tyy*A_y ;
	MF_px += Tyx*A_y ;

	Cs_pz += Tzx*A_x;
	Cs_py += Tyx*A_x;
	Cs_px += Txx*A_x;
	  
	MF_pz += Tzx*A_x;
	MF_py += Tyx*A_x;
	MF_px += Txx*A_x;

	M_px -=   r_y*MF_pz - r_z*MF_py;//
	M_py -= -(r_x*MF_pz - r_z*MF_px);
	M_pz -=   r_x*MF_py - r_y*MF_px;
      }   

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;
   
      //k- calculate Fz at k-1/2
      if (nvert[k-1][j][i]<0.9) {
	A_x=0.5*(zet[k][j][i].x+zet[k-1][j][i].x);//zet[k-1][j][i].x;
	A_y=0.5*(zet[k][j][i].y+zet[k-1][j][i].y);//zet[k-1][j][i].y;
	A_z=0.5*(zet[k][j][i].z+zet[k-1][j][i].z);//zet[k-1][j][i].z;

	x = 0.5*(coor[k][j][i].x+coor[k-1][j][i].x);
	y = 0.5*(coor[k][j][i].y+coor[k-1][j][i].y);
	z = 0.5*(coor[k][j][i].z+coor[k-1][j][i].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

/* 	Cp_nz +=  p[k-1][j][i]*A_z; */
/* 	MF_nz +=  p[k-1][j][i]*A_z; */
	Cp_nz +=  p[k][j][i]*A_z;
	MF_nz +=  p[k][j][i]*A_z;
	An_z  += A_z;
	Ian_x -=   r_y*A_z;
	Ian_y -= - r_x*A_z;
	MFdpdn_nz += -sb*(nfy*r_z-nfz*r_y)*A_z;

/* 	Cp_ny +=  p[k-1][j][i]*A_y; */
/* 	MF_ny +=  p[k-1][j][i]*A_y; */
	Cp_ny +=  p[k][j][i]*A_y;
	MF_ny +=  p[k][j][i]*A_y;
	An_y  += A_y;
	Ian_x -=  - r_z*A_y;
	Ian_z -=    r_x*A_y;
	MFdpdn_ny += -sb*(nfy*r_z-nfz*r_y)*A_y;

/* 	Cp_nx +=  p[k-1][j][i]*A_x; */
/* 	MF_nx +=  p[k-1][j][i]*A_x; */
	Cp_nx +=  p[k][j][i]*A_x;
	MF_nx +=  p[k][j][i]*A_x;
	An_x  += A_x;
	Ian_z -=   - r_y*A_x;
	Ian_y -= -(- r_z*A_x);
	MFdpdn_nx += -sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_nz -= Tzz*A_z ;
	Cs_ny -= Tzy*A_z ;
	Cs_nx -= Tzx*A_z ;
	
	MF_nz -= Tzz*A_z ;
	MF_ny -= Tzy*A_z ;
	MF_nx -= Tzx*A_z ;
	
	Cs_nz -= Tzy*A_y ;
	Cs_ny -= Tyy*A_y ;
	Cs_nx -= Tyx*A_y ;
	
	MF_nz -= Tzy*A_y ;
	MF_ny -= Tyy*A_y ;
	MF_nx -= Tyx*A_y ;
	
	Cs_nz -= Tzx*A_x;
	Cs_ny -= Tyx*A_x;
	Cs_nx -= Txx*A_x;
	
	MF_nz -= Tzx*A_x;
	MF_ny -= Tyx*A_x;
	MF_nx -= Txx*A_x;

	M_nx -=   r_y*MF_nz - r_z*MF_ny;//
	M_ny -= -(r_x*MF_nz - r_z*MF_nx);
	M_nz -=   r_x*MF_ny - r_y*MF_nx;
      }
     
      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

      //j+
      if (nvert[k][j+1][i]<0.9  ){
	A_x=0.5*(eta[k][j][i].x+eta[k][j+1][i].x);//eta[k][j  ][i].x;
	A_y=0.5*(eta[k][j][i].y+eta[k][j+1][i].y);//eta[k][j  ][i].y;
	A_z=0.5*(eta[k][j][i].z+eta[k][j+1][i].z);//eta[k][j  ][i].z;

	x = 0.5*(coor[k][j][i].x+coor[k][j+1][i].x);
	y = 0.5*(coor[k][j][i].y+coor[k][j+1][i].y);
	z = 0.5*(coor[k][j][i].z+coor[k][j+1][i].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

/* 	Cp_pz += -p[k][j+1][i]*A_z; */
/* 	MF_pz += -p[k][j+1][i]*A_z; */
	Cp_pz += -p[k][j][i]*A_z;
	MF_pz += -p[k][j][i]*A_z;
	Ap_z  += A_z;
	Iap_x -=   r_y*A_z;
	Iap_y -= - r_x*A_z;
	MFdpdn_pz += sb*(nfy*r_z-nfz*r_y)*A_z;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

/* 	Cp_py += -p[k][j+1][i]*A_y; */
/* 	MF_py += -p[k][j+1][i]*A_y; */
	Cp_py += -p[k][j][i]*A_y;
	MF_py += -p[k][j][i]*A_y;
	Ap_y  += A_y;
	Iap_x -=  - r_z*A_y;
	Iap_z -=    r_x*A_y;
	MFdpdn_py += sb*(nfy*r_z-nfz*r_y)*A_y;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

/* 	Cp_px += -p[k][j+1][i]*A_x; */
/* 	MF_px += -p[k][j+1][i]*A_x; */
	Cp_px += -p[k][j][i]*A_x;
	MF_px += -p[k][j][i]*A_x;
	Ap_x  += A_x;
	Iap_z -=   - r_y*A_x;
	Iap_y -= -(- r_z*A_x);
	MFdpdn_px += sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_pz += Tzz*A_z ;
	Cs_py += Tzy*A_z ;
	Cs_px += Tzx*A_z ;
	  
	MF_pz += Tzz*A_z ;
	MF_py += Tzy*A_z ;
	MF_px += Tzx*A_z ;

	Cs_pz += Tzy*A_y ;
	Cs_py += Tyy*A_y ;
	Cs_px += Tyx*A_y ;
	  
	MF_pz += Tzy*A_y ;
	MF_py += Tyy*A_y ;
	MF_px += Tyx*A_y ;

	Cs_pz += Tzx*A_x;
	Cs_py += Tyx*A_x;
	Cs_px += Txx*A_x;
	  
	MF_pz += Tzx*A_x;
	MF_py += Tyx*A_x;
	MF_px += Txx*A_x;

	M_px -=   r_y*MF_pz - r_z*MF_py;//
	M_py -= -(r_x*MF_pz - r_z*MF_px);
	M_pz -=   r_x*MF_py - r_y*MF_px;

      }

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

      //j-
      if (j==0) {
	A_x=eta[k][j][i].x;//eta[k][j-1][i].x;
	A_y=eta[k][j][i].y;//eta[k][j-1][i].y;
	A_z=eta[k][j][i].z;//eta[k][j-1][i].z;

	x = coor[k][j][i].x;//+coor[k][j-1][i].x);
	y = coor[k][j][i].y;//+coor[k][j-1][i].y);
	z = coor[k][j][i].z;//+coor[k][j-1][i].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

/* 	Cp_nz +=  p[k][j-1][i]*A_z; */
/* 	MF_nz +=  p[k][j-1][i]*A_z; */
	Cp_nz +=  p[k][j][i]*A_z;
	MF_nz +=  p[k][j][i]*A_z;
	An_z  += A_z;
	Ian_x -=   r_y*A_z;
	Ian_y -= - r_x*A_z;
	MFdpdn_nz += -sb*(nfy*r_z-nfz*r_y)*A_z;

/* 	Cp_ny +=  p[k][j-1][i]*A_y; */
/* 	MF_ny +=  p[k][j-1][i]*A_y; */
	Cp_ny +=  p[k][j][i]*A_y;
	MF_ny +=  p[k][j][i]*A_y;	     
	An_y  += A_y;
	Ian_x -=  - r_z*A_y;
	Ian_z -=    r_x*A_y;
	MFdpdn_ny += -sb*(nfy*r_z-nfz*r_y)*A_y;

/* 	Cp_nx +=  p[k][j-1][i]*A_x; */
/* 	MF_nx +=  p[k][j-1][i]*A_x; */
	Cp_nx +=  p[k][j][i]*A_x;
	MF_nx +=  p[k][j][i]*A_x;
	An_x  += A_x;
	Ian_z -=   - r_y*A_x;
	Ian_y -= -(- r_z*A_x);
	MFdpdn_nx += -sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_nz -= Tzz*A_z ;
	Cs_ny -= Tzy*A_z ;
	Cs_nx -= Tzx*A_z ;
	
	MF_nz -= Tzz*A_z ;
	MF_ny -= Tzy*A_z ;
	MF_nx -= Tzx*A_z ;
	
	Cs_nz -= Tzy*A_y ;
	Cs_ny -= Tyy*A_y ;
	Cs_nx -= Tyx*A_y ;
	
	MF_nz -= Tzy*A_y ;
	MF_ny -= Tyy*A_y ;
	MF_nx -= Tyx*A_y ;
	
	Cs_nz -= Tzx*A_x;
	Cs_ny -= Tyx*A_x;
	Cs_nx -= Txx*A_x;
	
	MF_nz -= Tzx*A_x;
	MF_ny -= Tyx*A_x;
	MF_nx -= Txx*A_x;

	M_nx -=   r_y*MF_nz - r_z*MF_ny;//
	M_ny -= -(r_x*MF_nz - r_z*MF_nx);
	M_nz -=   r_x*MF_ny - r_y*MF_nx;

      } else if (nvert[k][j-1][i]<0.9) {
	A_x=0.5*(eta[k][j][i].x+eta[k][j-1][i].x);//eta[k][j-1][i].x;
	A_y=0.5*(eta[k][j][i].y+eta[k][j-1][i].y);//eta[k][j-1][i].y;
	A_z=0.5*(eta[k][j][i].z+eta[k][j-1][i].z);//eta[k][j-1][i].z;

	x = 0.5*(coor[k][j][i].x+coor[k][j-1][i].x);
	y = 0.5*(coor[k][j][i].y+coor[k][j-1][i].y);
	z = 0.5*(coor[k][j][i].z+coor[k][j-1][i].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

/* 	Cp_nz +=  p[k][j-1][i]*A_z; */
/* 	MF_nz +=  p[k][j-1][i]*A_z; */
	Cp_nz +=  p[k][j][i]*A_z;
	MF_nz +=  p[k][j][i]*A_z;
	An_z  += A_z;
	Ian_x -=   r_y*A_z;
	Ian_y -= - r_x*A_z;
	MFdpdn_nz += -sb*(nfy*r_z-nfz*r_y)*A_z;

/* 	Cp_ny +=  p[k][j-1][i]*A_y; */
/* 	MF_ny +=  p[k][j-1][i]*A_y; */
	Cp_ny +=  p[k][j][i]*A_y;
	MF_ny +=  p[k][j][i]*A_y;	     
	An_y  += A_y;
	Ian_x -=  - r_z*A_y;
	Ian_z -=    r_x*A_y;
	MFdpdn_ny += -sb*(nfy*r_z-nfz*r_y)*A_y;

/* 	Cp_nx +=  p[k][j-1][i]*A_x; */
/* 	MF_nx +=  p[k][j-1][i]*A_x; */
	Cp_nx +=  p[k][j][i]*A_x;
	MF_nx +=  p[k][j][i]*A_x;
	An_x  += A_x;
	Ian_z -=   - r_y*A_x;
	Ian_y -= -(- r_z*A_x);
	MFdpdn_nx += -sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_nz -= Tzz*A_z ;
	Cs_ny -= Tzy*A_z ;
	Cs_nx -= Tzx*A_z ;
	
	MF_nz -= Tzz*A_z ;
	MF_ny -= Tzy*A_z ;
	MF_nx -= Tzx*A_z ;
	
	Cs_nz -= Tzy*A_y ;
	Cs_ny -= Tyy*A_y ;
	Cs_nx -= Tyx*A_y ;
	
	MF_nz -= Tzy*A_y ;
	MF_ny -= Tyy*A_y ;
	MF_nx -= Tyx*A_y ;
	
	Cs_nz -= Tzx*A_x;
	Cs_ny -= Tyx*A_x;
	Cs_nx -= Txx*A_x;
	
	MF_nz -= Tzx*A_x;
	MF_ny -= Tyx*A_x;
	MF_nx -= Txx*A_x;

	M_nx -=   r_y*MF_nz - r_z*MF_ny;//
	M_ny -= -(r_x*MF_nz - r_z*MF_nx);
	M_nz -=   r_x*MF_ny - r_y*MF_nx;
      }

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

      //i+
      if (nvert[k][j][i+1]<0.9){
	A_x=0.5*(csi[k][j][i].x+csi[k][j][i+1].x);//csi[k][j][i].x;
	A_y=0.5*(csi[k][j][i].y+csi[k][j][i+1].y);//csi[k][j][i].y;
	A_z=0.5*(csi[k][j][i].z+csi[k][j][i+1].z);//csi[k][j][i].z;

	x = 0.5*(coor[k][j][i].x+coor[k][j][i+1].x);
	y = 0.5*(coor[k][j][i].y+coor[k][j][i+1].y);
	z = 0.5*(coor[k][j][i].z+coor[k][j][i+1].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

/* 	Cp_pz += -p[k][j][i+1]*A_z; */
/* 	MF_pz += -p[k][j][i+1]*A_z; */
	Cp_pz += -p[k][j][i]*A_z;
	MF_pz += -p[k][j][i]*A_z;
	Ap_z  += A_z;
	Iap_x -=   r_y*A_z;
	Iap_y -= - r_x*A_z;
	MFdpdn_pz += sb*(nfy*r_z-nfz*r_y)*A_z;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

/* 	Cp_py += -p[k][j][i+1]*A_y; */
/* 	MF_py += -p[k][j][i+1]*A_y; */
	Cp_py += -p[k][j][i]*A_y;
	MF_py += -p[k][j][i]*A_y;
	Ap_y  += A_y;
	Iap_x -=  - r_z*A_y;
	Iap_z -=    r_x*A_y;
	MFdpdn_py += sb*(nfy*r_z-nfz*r_y)*A_y;

/*       if (i==2 && ibi==2) { */
/* 	PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d Ia %le %le %le %le %le\n", i, j, k,Iap_x, r_y, A_z, r_z,A_y) ; */
/* 	//	flagprint =1; */
/*       } */

/* 	Cp_px += -p[k][j][i+1]*A_x; */
/* 	MF_px += -p[k][j][i+1]*A_x; */
	Cp_px += -p[k][j][i]*A_x;
	MF_px += -p[k][j][i]*A_x;
	Ap_x  += A_x;
	Iap_z -=   - r_y*A_x;
	Iap_y -= -(- r_z*A_x);
	MFdpdn_px += sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_pz += Tzz*A_z ;
	Cs_py += Tzy*A_z ;
	Cs_px += Tzx*A_z ;
	  
	MF_pz += Tzz*A_z ;
	MF_py += Tzy*A_z ;
	MF_px += Tzx*A_z ;

	Cs_pz += Tzy*A_y ;
	Cs_py += Tyy*A_y ;
	Cs_px += Tyx*A_y ;
	  
	MF_pz += Tzy*A_y ;
	MF_py += Tyy*A_y ;
	MF_px += Tyx*A_y ;

	Cs_pz += Tzx*A_x;
	Cs_py += Tyx*A_x;
	Cs_px += Txx*A_x;
	  
	MF_pz += Tzx*A_x;
	MF_py += Tyx*A_x;
	MF_px += Txx*A_x;

	M_px -=   r_y*MF_pz - r_z*MF_py;//
	M_py -= -(r_x*MF_pz - r_z*MF_px);
	M_pz -=   r_x*MF_py - r_y*MF_px;
      }

      MF_px=0.;MF_py=0.;MF_pz=0.;
      MF_nx=0.;MF_ny=0.;MF_nz=0.;

      //i-
      if  (nvert[k][j][i-1]<0.9) {
	A_x=0.5*(csi[k][j][i].x+csi[k][j][i-1].x);//csi[k][j][i-1].x;
	A_y=0.5*(csi[k][j][i].y+csi[k][j][i-1].y);//csi[k][j][i-1].y;
	A_z=0.5*(csi[k][j][i].z+csi[k][j][i-1].z);//csi[k][j][i-1].z;

	x = 0.5*(coor[k][j][i].x+coor[k][j][i-1].x);
	y = 0.5*(coor[k][j][i].y+coor[k][j][i-1].y);
	z = 0.5*(coor[k][j][i].z+coor[k][j][i-1].z);

	r_x = x-X_c;
	r_y = y-Y_c;
	r_z = z-Z_c;

/* 	Cp_nz +=  p[k][j][i-1]*A_z; */
/* 	MF_nz +=  p[k][j][i-1]*A_z; */
	Cp_nz +=  p[k][j][i]*A_z;
	MF_nz +=  p[k][j][i]*A_z;
	An_z  += A_z;
	Ian_x -=   r_y*A_z;
	Ian_y -= - r_x*A_z;
	MFdpdn_nz += -sb*(nfy*r_z-nfz*r_y)*A_z;

/* 	Cp_ny +=  p[k][j][i-1]*A_y; */
/* 	MF_ny +=  p[k][j][i-1]*A_y; */
	Cp_ny +=  p[k][j][i]*A_y;
	MF_ny +=  p[k][j][i]*A_y;	     
	An_y  += A_y;
	Ian_x -=  - r_z*A_y;
	Ian_z -=    r_x*A_y;
	MFdpdn_ny += -sb*(nfy*r_z-nfz*r_y)*A_y;

/* 	Cp_nx +=  p[k][j][i-1]*A_x; */
/* 	MF_nx +=  p[k][j][i-1]*A_x; */
	Cp_nx +=  p[k][j][i]*A_x;
	MF_nx +=  p[k][j][i]*A_x;
	An_x  += A_x;
	Ian_z -=   - r_y*A_x;
	Ian_y -= -(- r_z*A_x);
	MFdpdn_nx += -sb*(nfy*r_z-nfz*r_y)*A_x;

	Cs_nz -= Tzz*A_z ;
	Cs_ny -= Tzy*A_z ;
	Cs_nx -= Tzx*A_z ;
	
	MF_nz -= Tzz*A_z ;
	MF_ny -= Tzy*A_z ;
	MF_nx -= Tzx*A_z ;
	
	Cs_nz -= Tzy*A_y ;
	Cs_ny -= Tyy*A_y ;
	Cs_nx -= Tyx*A_y ;
	
	MF_nz -= Tzy*A_y ;
	MF_ny -= Tyy*A_y ;
	MF_nx -= Tyx*A_y ;
	
	Cs_nz -= Tzx*A_x;
	Cs_ny -= Tyx*A_x;
	Cs_nx -= Txx*A_x;
	
	MF_nz -= Tzx*A_x;
	MF_ny -= Tyx*A_x;
	MF_nx -= Txx*A_x;
	
	M_nx -=   r_y*MF_nz - r_z*MF_ny;//
	M_ny -= -(r_x*MF_nz - r_z*MF_nx);
	M_nz -=   r_x*MF_ny - r_y*MF_nx;

      }
    
/*       M_px -=   r_y*MF_pz - r_z*MF_py;// */
/*       M_py -= -(r_x*MF_pz - r_z*MF_px); */
/*       M_pz -=   r_x*MF_py - r_y*MF_px; */

/*       M_nx -=   r_y*MF_nz - r_z*MF_ny;// */
/*       M_ny -= -(r_x*MF_nz - r_z*MF_nx); */
/*       M_nz -=   r_x*MF_ny - r_y*MF_nx; */

      Mdpdn_px -=   r_y*MFdpdn_pz - r_z*MFdpdn_py;//
      Mdpdn_py -= -(r_x*MFdpdn_pz - r_z*MFdpdn_px);
      Mdpdn_pz -=   r_x*MFdpdn_py - r_y*MFdpdn_px;

      Mdpdn_nx -=   r_y*MFdpdn_nz - r_z*MFdpdn_ny;//
      Mdpdn_ny -= -(r_x*MFdpdn_nz - r_z*MFdpdn_nx);
      Mdpdn_nz -=   r_x*MFdpdn_ny - r_y*MFdpdn_nx;
      
/*       if (ibi==1 & i==2) */
/* 	PetscPrintf(PETSC_COMM_SELF, "MF %d %d %le %le %le %le %le %le %le %le\n",j,k,MF_py,MF_ny,MF_pz,MF_nz,M_px,M_nx,r_y,r_z); */

    }
/* ==================================================================================             */
/*     End of Loop ibm nodes */
  }

/* ==================================================================================             */
/*   Total Force on each processor */
  F_px = Cp_px + Cs_px; 
  F_py = Cp_py + Cs_py;
  F_pz = Cp_pz + Cs_pz;

  F_nx = Cp_nx + Cs_nx; 
  F_ny = Cp_ny + Cs_ny;
  F_nz = Cp_nz + Cs_nz;

/* ==================================================================================             */  
/*   Global Sum */
  MPI_Allreduce(&F_px, &F_pxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&F_py, &F_pySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&F_pz, &F_pzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&F_nx, &F_nxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&F_ny, &F_nySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&F_nz, &F_nzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Ap_x, &Ap_xSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Ap_y, &Ap_ySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Ap_z, &Ap_zSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&An_x, &An_xSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&An_y, &An_ySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&An_z, &An_zSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Cp_nx, &Cp_nxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Cp_ny, &Cp_nySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Cp_nz, &Cp_nzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Cp_px, &Cp_pxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Cp_py, &Cp_pySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Cp_pz, &Cp_pzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&M_px, &M_pxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&M_py, &M_pySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&M_pz, &M_pzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&M_nx, &M_nxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&M_ny, &M_nySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&M_nz, &M_nzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Iap_x, &Iap_xSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Iap_y, &Iap_ySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Iap_z, &Iap_zSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Ian_x, &Ian_xSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Ian_y, &Ian_ySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Ian_z, &Ian_zSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Mdpdn_px, &Mdpdn_pxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Mdpdn_py, &Mdpdn_pySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Mdpdn_pz, &Mdpdn_pzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  MPI_Allreduce(&Mdpdn_nx, &Mdpdn_nxSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Mdpdn_ny, &Mdpdn_nySum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&Mdpdn_nz, &Mdpdn_nzSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

  
/* ==================================================================================             */
/*   Scale Check later !!!!! */

  A_xSum = 0.5 * (Ap_xSum + An_xSum);
  A_ySum = 0.5 * (Ap_ySum + An_ySum);
  A_zSum = 0.5 * (Ap_zSum + An_zSum);

/*   if (fabs(Ap_xSum)>1e-6) */
/*     F_pxSum=F_pxSum*A_xSum/Ap_xSum; */
/*   if (fabs(Ap_ySum)>1e-6) */
/*     F_pySum=F_pySum*A_ySum/Ap_ySum; */
/*   if (fabs(Ap_zSum)>1e-6) */
/*     F_pzSum=F_pzSum*A_zSum/Ap_zSum; */

/*   if (fabs(An_xSum)>1e-6) */
/*     F_nxSum=F_nxSum*A_xSum/An_xSum; */
/*   if (fabs(An_ySum)>1e-6) */
/*     F_nySum=F_nySum*A_ySum/An_ySum; */
/*   if (fabs(An_zSum)>1e-6) */
/*     F_nzSum=F_nzSum*A_zSum/An_zSum; */

  F_xSum = F_pxSum + F_nxSum;
  F_ySum = F_pySum + F_nySum;
  F_zSum = F_pzSum + F_nzSum;

  //if (!fish && !cop) {
  /* if (fabs(A_xSum)>1e-6) */
  /*   F_xSum=F_xSum/A_xSum*2.; */
  /* if (fabs(A_ySum)>1e-6) */
  /*   F_ySum=F_ySum/A_ySum*2.; */
  /* if (fabs(A_zSum)>1e-6) */
  /*   F_zSum=F_zSum/A_zSum*2.; */
  //  }

/*   if (fabs(Ap_xSum)>1e-6) */
/*     Cp_pxSum=Cp_pxSum*A_xSum/Ap_xSum; */
/*   if (fabs(Ap_ySum)>1e-6) */
/*     Cp_pySum=Cp_pySum*A_ySum/Ap_ySum; */
/*   if (fabs(Ap_zSum)>1e-6) */
/*     Cp_pzSum=Cp_pzSum*A_zSum/Ap_zSum; */

/*   if (fabs(An_xSum)>1e-6) */
/*     Cp_nxSum=Cp_nxSum*A_xSum/An_xSum; */
/*   if (fabs(An_ySum)>1e-6) */
/*     Cp_nySum=Cp_nySum*A_ySum/An_ySum; */
/*   if (fabs(An_zSum)>1e-6) */
/*     Cp_nzSum=Cp_nzSum*A_zSum/An_zSum; */

  Cp_xSum = Cp_pxSum + Cp_nxSum;
  Cp_ySum = Cp_pySum + Cp_nySum;
  Cp_zSum = Cp_pzSum + Cp_nzSum;

  /* if (fabs(A_xSum)>1e-6) */
  /*   Cp_xSum=Cp_xSum/A_xSum*2.; */
  /* if (fabs(A_ySum)>1e-6) */
  /*   Cp_ySum=Cp_ySum/A_ySum*2.; */
  /* if (fabs(A_zSum)>1e-6) */
  /*   Cp_zSum=Cp_zSum/A_zSum*2.; */    

  Ia_xSum = 0.5 * (Iap_xSum - Ian_xSum);
  Ia_ySum = 0.5 * (Iap_ySum - Ian_ySum);
  Ia_zSum = 0.5 * (Iap_zSum - Ian_zSum);

/*   if (fabs(Iap_xSum)>1e-6) */
/*     M_pxSum=M_pxSum*Ia_xSum/Iap_xSum; */
/*   if (fabs(Iap_ySum)>1e-6) */
/*     M_pySum=M_pySum*Ia_ySum/Iap_ySum; */
/*   if (fabs(Iap_zSum)>1e-6) */
/*     M_pzSum=M_pzSum*Ia_zSum/Iap_zSum ; */

/*   if (fabs(Ian_xSum)>1e-6) */
/*     M_nxSum=M_nxSum*Ia_xSum/Ian_xSum; */
/*   if (fabs(Ian_ySum)>1e-6) */
/*     M_nySum=M_nySum*Ia_ySum/Ian_ySum; */
/*   if (fabs(Ian_zSum)>1e-6) */
/*     M_nzSum=M_nzSum*Ia_zSum/Ian_zSum; */

  M_xSum = M_pxSum + M_nxSum;
  M_ySum = M_pySum + M_nySum;
  M_zSum = M_pzSum + M_nzSum;

/*   if (fabs(Iap_xSum)>1e-6) */
/*     Mdpdn_pxSum=Mdpdn_pxSum*Ia_xSum/Iap_xSum; */
/*   if (fabs(Iap_ySum)>1e-6) */
/*     Mdpdn_pySum=Mdpdn_pySum*Ia_ySum/Iap_ySum; */
/*   if (fabs(Iap_zSum)>1e-6) */
/*     Mdpdn_pzSum=Mdpdn_pzSum*Ia_zSum/Iap_zSum ; */

/*   if (fabs(Ian_xSum)>1e-6) */
/*     Mdpdn_nxSum=Mdpdn_nxSum*Ia_xSum/Ian_xSum; */
/*   if (fabs(Ian_ySum)>1e-6) */
/*     Mdpdn_nySum=Mdpdn_nySum*Ia_ySum/Ian_ySum; */
/*   if (fabs(Ian_zSum)>1e-6) */
/*     Mdpdn_nzSum=Mdpdn_nzSum*Ia_zSum/Ian_zSum; */

  Mdpdn_xSum = Mdpdn_pxSum + Mdpdn_nxSum;
  Mdpdn_ySum = Mdpdn_pySum + Mdpdn_nySum;
  Mdpdn_zSum = Mdpdn_pzSum + Mdpdn_nzSum;

/*   if (fabs(Ap_zSum)>1e-6) */
/*     M_xSum=M_xSum/A_zSum*2.; */
/*   if (fabs(Ap_ySum)>1e-6) */
/*     M_ySum=M_ySum/A_ySum*2.; */
/*   if (fabs(Ap_zSum)>1e-6) */
/*     M_zSum=M_zSum/A_zSum*2.; */

  A_totSum = Ap_xSum + Ap_ySum + Ap_zSum;

/* ==================================================================================             */
/*   store results in fsi */
  fsi->F_x = F_xSum; fsi->F_y = F_ySum; fsi->F_z = F_zSum;
  fsi->A_tot = A_totSum;
  fsi->M_x = M_xSum; fsi->M_y = M_ySum; fsi->M_z = M_zSum;
  fsi->Mdpdn_x = Mdpdn_xSum; fsi->Mdpdn_y = Mdpdn_ySum; fsi->Mdpdn_z = Mdpdn_zSum;

/* ==================================================================================             */
/*   output values */
  PetscPrintf(PETSC_COMM_WORLD, "F_x,F_y,F_z SI, %le %le %le Az %le %le Ay %le %le\n",F_xSum,F_ySum,F_zSum,Ap_zSum,An_zSum,Ap_ySum,An_ySum);
  PetscPrintf(PETSC_COMM_WORLD, "M_x,M_y,M_z SI, %le %le %le Ia_x %le %le Ip_y %le %le\n",M_xSum,M_ySum,M_zSum,Iap_xSum,Ian_xSum,Iap_ySum,Ian_ySum);
  PetscPrintf(PETSC_COMM_WORLD, "Mdpdn_x,Mdpdn_y,Mdpdn_z SI, %le %le %le\n",Mdpdn_xSum,Mdpdn_ySum,Mdpdn_zSum);

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Force_Coeff_SI%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le\n",ti, F_xSum, F_ySum, F_zSum,Cp_xSum,Cp_ySum,Cp_zSum, A_xSum,A_ySum,A_zSum);
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le %le\n",ti,F_ySum, F_zSum,Cp_ySum,Cp_zSum,F_ySum-Cp_ySum,F_zSum-Cp_zSum,A_ySum,F_xSum, M_xSum); */
    fclose(f);

    sprintf(filen, "Momt_Coeff_SI%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le\n",ti, M_xSum, M_ySum, M_zSum,fsi->M_x_old,fsi->M_x_real,fsi->Mdpdn_x);
    fclose(f);

  }

/* ==================================================================================             */
/*   Restore Working arrays */
  //DMDAVecRestoreArray(fda, user->lCent, &cent);
  DMDAVecRestoreArray(fda, Coor, &coor);
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
  DMDAVecRestoreArray(da, user->lP, &p);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  /* DMDAVecRestoreArray(da, user->lIAj, &iaj); */
  /* DMDAVecRestoreArray(da, user->lJAj, &jaj); */
  /* DMDAVecRestoreArray(da, user->lKAj, &kaj); */
  // VecDestroy(&Coor);
  return(0);
}

/* ==================================================================================             */
/*   subroutines related to projection of forces onto ibm! */
/* ==================================================================================             */
PetscInt IB_node_total_each_cpu(UserCtx *user, IBMNodes *ibm)
{
 
  IBMListNode  *current;
  PetscInt      number_of_IBnodes_at_this_cpu;
  // point to the first node at each cpu
  //IBMListNode * current = user->ibmlist[0].head;
  current = user->ibmlist[0].head;
  number_of_IBnodes_at_this_cpu = 0;
  while (current)
    {
      current = current->next;
      number_of_IBnodes_at_this_cpu++;
    }
  return number_of_IBnodes_at_this_cpu;
} 

PetscErrorCode Initalize_Projecting(IBMNodes * ibm ) {
  PetscInt	n_elmt = ibm->n_elmt;
  PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->tau0));
  PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->tauN));
  PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->pres));  

  PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Bvel_u));  
  PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Bvel_v));  
  PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->Bvel_w));  
  return(0);
}

PetscErrorCode Projecting( UserCtx * user, IBMNodes * ibm )
{
  IBMInfo          *ibminfo;
  //IBMListNode      *current;
  PetscInt         *current_element;
  DM               da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

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

  PetscReal         *NIB;
  PetscInt         processors,rank;
  PetscInt         number_of_IBnodes=0, total_IBnodes;
  PetscInt         *rec, *dis, *stride, start;
   
  PetscReal        uu, vv, ww, tau0, tauN, pres;
  PetscReal        *current_tau0,  *current_tauN,  *current_pres, *current_uvel, *current_vvel, *current_wvel;       
  PetscReal        *tau0_buffer,*tauN_buffer, *pres_buffer, *uvel_buffer, *vvel_buffer, *wvel_buffer;
  
  PetscInt         i,j,k,elmt;
  PetscInt         *element_buffer;
  MPI_Datatype     stype, ltype;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &processors);

  // Creating the buffer to send and recieve information
  PetscMalloc(processors*sizeof(PetscInt), &rec);
  PetscMalloc(processors*sizeof(PetscInt), &dis);
  PetscMalloc(processors*sizeof(PetscInt), &stride);

  // Finding the total number of IB nodes on each cpu

  IB_node_total_each_cpu(user, ibm);
  number_of_IBnodes = IB_node_total_each_cpu(user, ibm);

  // Finding total number of IB nodes in all cpu's to all cpu's knowing the total
  MPI_Allreduce(&number_of_IBnodes,&total_IBnodes,1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // All cpu's sending thier number of IB nodes to the root cpu naming "stride" in root
  MPI_Gather(&number_of_IBnodes, 1, MPI_INT, stride, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Now root Bcasts nomber of each cpu's IB node as stride[i] to all other processors
  MPI_Bcast(stride, processors, MPI_INT, 0, MPI_COMM_WORLD);

   PetscPrintf(PETSC_COMM_SELF, "Ibm nodes total %d on each cpu %d!\n",total_IBnodes, number_of_IBnodes); 
  // Allocate the rec and dis (placment) to each and all processors
  start = 0;
  for (i=0;i<processors;i++) { 
    dis[i] = start;
    rec[i] = stride[i];
    start = start + stride[i];
  }
  
 
  PetscMalloc(number_of_IBnodes*sizeof(PetscInt), &current_element);
  PetscMalloc(number_of_IBnodes*sizeof(PetscReal), &current_tau0);
  PetscMalloc(number_of_IBnodes*sizeof(PetscReal), &current_tauN);
  PetscMalloc(number_of_IBnodes*sizeof(PetscReal), &current_pres);
  PetscMalloc(number_of_IBnodes*sizeof(PetscReal), &current_uvel);
  PetscMalloc(number_of_IBnodes*sizeof(PetscReal), &current_vvel);
  PetscMalloc(number_of_IBnodes*sizeof(PetscReal), &current_wvel);
  

  PetscMalloc(total_IBnodes*sizeof(PetscInt), &element_buffer);
  PetscMalloc(total_IBnodes*sizeof(PetscReal), &tau0_buffer);
  PetscMalloc(total_IBnodes*sizeof(PetscReal), &tauN_buffer);
  PetscMalloc(total_IBnodes*sizeof(PetscReal), &pres_buffer);
  PetscMalloc(total_IBnodes*sizeof(PetscReal), &uvel_buffer);
  PetscMalloc(total_IBnodes*sizeof(PetscReal), &vvel_buffer);
  PetscMalloc(total_IBnodes*sizeof(PetscReal), &wvel_buffer);

  // PetscPrintf(PETSC_COMM_SELF, "Ibm nodes total %d on each cpu %d rank %d!\n",total_IBnodes, number_of_IBnodes, rank); 
  // Now putting the IB infor into a local current array to be sent to the root process later
  // point to the first IB node..
  PetscReal     ***ustar, ***p;
  Cmpnts	***ucat,***coor;
  extern PetscInt wallfunction, invicid;
  DMDAVecGetArray(da, user->P, &p);
  DMDAVecGetArray(fda, user->Ucat, &ucat);
  Vec       coords;
  DMGetCoordinatesLocal(da, &coords);
  DMDAVecGetArray(fda, coords, &coor); // change for compressible
  //  DMDAVecGetArray(fda, user->Cent, &coor);
  
  if (wallfunction) DMDAVecGetArray(da, user->lUstar, &ustar);
  IBMListNode * current = user->ibmlist[0].head;
  //current = user->ibmlist.head;
  PetscInt ii= 0;
  while(current) {
    ibminfo = &current->ibm_intp;
    current=current->next;
    i = ibminfo->ni, j= ibminfo->nj, k = ibminfo->nk;
    // point to the related neighbour bed cell
    elmt = ibminfo->cell;
    double cs1 = ibminfo->cs1, cs2 = ibminfo->cs2, cs3 = ibminfo->cs3;
    double nfx , nfy, nfz;	
    double sb = ibminfo->d_s, Ua_N;
    Cmpnts Ua, Ub;
    nfx=coor[k][j][i].x-ibminfo->pmin.x;
    nfy=coor[k][j][i].y-ibminfo->pmin.y;
    nfz=coor[k][j][i].z-ibminfo->pmin.z;
    double dr=sqrt(nfx*nfx+nfy*nfy+nfz*nfz);
    nfx=nfx/dr;
    nfy=nfy/dr;
    nfz=nfz/dr;

    if (elmt>=0) {
      Ua.x = ibm->u[ibm->nv1[elmt]].x * cs1 + ibm->u[ibm->nv2[elmt]].x * cs2 + ibm->u[ibm->nv3[elmt]].x * cs3;
      Ua.y = ibm->u[ibm->nv1[elmt]].y * cs1 + ibm->u[ibm->nv2[elmt]].y * cs2 + ibm->u[ibm->nv3[elmt]].y * cs3;
      Ua.z = ibm->u[ibm->nv1[elmt]].z * cs1 + ibm->u[ibm->nv2[elmt]].z * cs2 + ibm->u[ibm->nv3[elmt]].z * cs3;

      Ua_N =  Ua.x*nfx + Ua.y*nfy + Ua.z*nfz;
    }  else {
      Ua.x = Ua.y = Ua.z = 0;
      Ua_N=0.;
    }
    if (ii<number_of_IBnodes) {
      current_element[ii] = elmt;
      if( i >= lxs && i < lxe && j >= lys && j < lye && k >= lzs && k < lze ) {
	current_pres[ii] = p[k][j][i];
	current_uvel[ii] = ucat[k][j][i].x;
	current_vvel[ii] = ucat[k][j][i].y;
	current_wvel[ii] = ucat[k][j][i].z;

	Ub.x = ucat[k][j][i].x;
	Ub.y = ucat[k][j][i].y;
	Ub.z = ucat[k][j][i].z;
	// calcuate shear
	if (wallfunction) 
	  current_tau0[ii] = ustar[k][j][i]*ustar[k][j][i];
	else if (invicid) current_tau0[ii]=0.;
	else {
	  Cmpnts Udef;
	  PetscReal Udef_n;
	  Udef.x= (Ub.x-Ua.x);
	  Udef.y= (Ub.y-Ua.y);
	  Udef.z= (Ub.z-Ua.z);
	  Udef_n= Udef.x*nfx + Udef.y*nfy + Udef.z*nfz;
	  //tangential defect
	  Udef.x -= Udef_n * nfx;
	  Udef.y -= Udef_n * nfy;
	  Udef.z -= Udef_n * nfz;

	  current_tau0[ii]=sqrt(Udef.x*Udef.x +Udef.y*Udef.y +Udef.z*Udef.z)/sb/user->ren;
	  current_tauN[ii]=2*( (Ub.x-Ua.x)*nfx + (Ub.y-Ua.y)*nfy + (Ub.z-Ua.z)*nfz )/sb/user->ren;

	}
      }
    }
    ii++;
  } // End of local IB node on each process ( end of IBM list nodes)
  // DMDAVecRestoreArray(fda, user->Cent, &coor);
  DMDAVecRestoreArray(fda, coords,&coor);
  DMDAVecRestoreArray(da, user->P, &p);
  DMDAVecRestoreArray(fda, user->Ucat, &ucat);
  if (wallfunction) DMDAVecRestoreArray(da, user->lUstar, &ustar);

  
  // first sending the number of elements to the ROOT
  // then sending the local info from processes to ROOT
  MPI_Type_contiguous(number_of_IBnodes, MPI_INT, &stype);
  MPI_Type_commit(&stype);
  MPI_Gatherv(current_element, 1, stype, element_buffer, rec, dis, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Type_contiguous(number_of_IBnodes, MPI_DOUBLE, &ltype);
  MPI_Type_commit(&ltype);
  MPI_Gatherv(current_tau0, 1, ltype, tau0_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(current_tauN, 1, ltype, tauN_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(current_pres, 1, ltype, pres_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(current_uvel, 1, ltype, uvel_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(current_vvel, 1, ltype, vvel_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(current_wvel, 1, ltype, wvel_buffer, rec, dis, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  //  PetscBarrier(PETSC_NULL);
  // Now all info are at ROOT, then ROOt will postprocess the info to assemble the info on the bed surface

  if(rank == 0)
    {
      PetscMalloc(ibm->n_elmt*sizeof(PetscReal), &NIB);
      // clear former info
      for(i=0; i<ibm->n_elmt;i++){
	ibm->Bvel_u[i] = 0.0;
	ibm->Bvel_v[i] = 0.0;
	ibm->Bvel_w[i] = 0.0;
	ibm->pres[i] = 0.0;
	ibm->tau0[i] = 0.0;
	ibm->tauN[i] = 0.0;
	NIB[i]=0.;
      }
      
      
      // Projecting onto the body surface
      for(i=0; i<total_IBnodes; i++) {
   
	elmt = element_buffer[i];
	if(elmt < ibm->n_elmt && elmt>0) { 
	  tau0 = tau0_buffer[i];
	  tauN = tauN_buffer[i];
	  pres = pres_buffer[i];
	  uu = uvel_buffer[i];
	  vv = vvel_buffer[i];
	  ww = wvel_buffer[i];
	  
	  // transforming IB node's info into the body surface element 
	  if(fabs(uu) + fabs(vv) + fabs(ww) > 0.0)   {
	    ibm->Bvel_u[elmt] += uu;
	    ibm->Bvel_v[elmt] += vv;
	    ibm->Bvel_w[elmt] += ww;
	    ibm->pres[elmt] += pres;
	    ibm->tau0[elmt] += tau0;
	    ibm->tauN[elmt] += tauN;
	    NIB[elmt] += 1.;
	  }
	}
      } // End of all IB nodes  
 
      // Now averaging the variables on each triangle cell of body surface
      for(i=0; i<ibm->n_elmt; i++){
	if(NIB[i]>0.1) {
	  ibm->Bvel_u[i] /= NIB[i];
	  ibm->Bvel_v[i] /= NIB[i];
	  ibm->Bvel_w[i] /= NIB[i];
	  ibm->pres[i] /= NIB[i];
	  ibm->tau0[i] /= NIB[i];
	  ibm->tauN[i] /= NIB[i];
	}
      }
      
      PetscInt	nelmt, nii;
      PetscReal riter;
      // assign & smooth values for elmts that dont have a value
      for(i=0;i<3;i++){
      for( nelmt = 0; nelmt < ibm->n_elmt; nelmt++ )    {
	if(NIB[nelmt]<0.5)     {
	  riter = 0.;
	  for(nii = 0; nii < ibm->n_elmt; nii++)    {
	    // if different elements & on the same side n1.n2>0
	    if(nelmt!=nii  && ((ibm->nf_x[nelmt]*ibm->nf_x[nii] +
				ibm->nf_y[nelmt]*ibm->nf_y[nii] +
				ibm->nf_z[nelmt]*ibm->nf_z[nii])>0.))    {
	      //---------- interpolation between all around vertics (more than 8 neighbor elemts) ************
	      if(ibm->nv1[ nelmt ] == ibm->nv1[ nii ] ||
		 ibm->nv1[ nelmt ] == ibm->nv2[ nii ] ||
		 ibm->nv1[ nelmt ] == ibm->nv3[ nii ] ||
		 ibm->nv2[ nelmt ] == ibm->nv1[ nii ] ||
		 ibm->nv2[ nelmt ] == ibm->nv2[ nii ] ||
		 ibm->nv2[ nelmt ] == ibm->nv3[ nii ] ||
		 ibm->nv3[ nelmt ] == ibm->nv1[ nii ] ||
		 ibm->nv3[ nelmt ] == ibm->nv2[ nii ] ||
		 ibm->nv3[ nelmt ] == ibm->nv3[ nii ])       {
		riter = riter + 1.;
		ibm->pres[ nelmt ] = ( ibm->pres[ nii ] + ibm->pres[ nelmt ] * ( riter - 1. ) ) / riter;
		ibm->tau0[ nelmt ] = ( ibm->tau0[ nii ] + ibm->tau0[ nelmt ] * ( riter - 1. ) ) / riter;
		ibm->tauN[ nelmt ] = ( ibm->tauN[ nii ] + ibm->tauN[ nelmt ] * ( riter - 1. ) ) / riter;
		ibm->Bvel_u[ nelmt] = ( ibm->Bvel_u[ nii ] + ibm->Bvel_u[ nelmt ] * ( riter - 1. ) ) / riter;
		ibm->Bvel_v[ nelmt] = ( ibm->Bvel_v[ nii ] + ibm->Bvel_v[ nelmt ] * ( riter - 1. ) ) / riter;
		ibm->Bvel_w[ nelmt] = ( ibm->Bvel_w[ nii ] + ibm->Bvel_w[ nelmt ] * ( riter - 1. ) ) / riter;
	      }
	    }
	  }
	}
      }
      }
  
      
      PetscFree(NIB);

      // Now ROOT Bcast the info to all other processes
      MPI_Bcast(ibm->Bvel_u, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast(ibm->Bvel_v, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast(ibm->Bvel_w, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast(ibm->pres, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast(ibm->tau0, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast(ibm->tauN, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);

    }  // End of ROOT process work (end of rank == 0)


  // other cpu's getting the info from ROOT here
  else  {
      MPI_Bcast(ibm->Bvel_u, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast(ibm->Bvel_v, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast(ibm->Bvel_w, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast(ibm->pres, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast(ibm->tau0, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast(ibm->tauN, ibm->n_elmt, MPIU_REAL, 0, MPI_COMM_WORLD);
     
    } // End of recieving info from ROOT 
 
  PetscBarrier(PETSC_NULL);

  PetscFree(rec);
  PetscFree(dis);
  PetscFree(stride);
  
  PetscFree(current_element);
  PetscFree(current_tau0);
  PetscFree(current_tauN);
  PetscFree(current_pres);
  PetscFree(current_uvel);
  PetscFree(current_vvel);
  PetscFree(current_wvel);
  PetscFree(element_buffer);
  PetscFree(tau0_buffer);
  PetscFree(tauN_buffer);
  PetscFree(pres_buffer);
  PetscFree(uvel_buffer);
  PetscFree(vvel_buffer);
  PetscFree(wvel_buffer);

  return 0;
}


PetscErrorCode Finalizing_Projecting( UserCtx * user, IBMNodes * ibm )
{


  PetscReal       riter;
  PetscInt	  i,nii,nelmt;
  

  
  /*for(elmt=0; elmt<n_elmt; elmt++)
    {
    if(ibm->nf_z[elmt] < 1.e-7 || ibm->elmt_depth[nelmt] > 0.2)
    {
    ibm->Shvel[elmt] = 0.;
    ibm->Bvel[elmt].x = 0.;
    ibm->Bvel[elmt].y = 0.;
    ibm->Bvel[elmt].z = 0.;
    }
    }*/	   
  

  for( nelmt = 0; nelmt < ibm->n_elmt; nelmt++ )    {      
    double test = fabs(ibm->pres[nelmt])+fabs(ibm->tauN[nelmt])+fabs(ibm->tau0[nelmt])+
                  fabs(ibm->Bvel_u[nelmt])+fabs(ibm->Bvel_v[nelmt])+fabs(ibm->Bvel_w[nelmt]);
    if(test < 1.e-7)     {  
      riter = 0.;
      for(nii = 0; nii < ibm->n_elmt; nii++)    {
	double testi = fabs(ibm->pres[nii])+fabs(ibm->tauN[nii])+fabs(ibm->tau0[nii])
	  +fabs(ibm->Bvel_u[nii])+fabs(ibm->Bvel_v[nii])+fabs(ibm->Bvel_w[nii]);
	if(nelmt!=nii && testi > 1.e-7 )    {			  
	  //--------- interpolation between around elment (max 3 around elemts)***************
	  /*
	    if(((ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || ibm->nv1[ nelmt ] == ibm->nv3[ nii ])
	    &&(ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || ibm->nv2[ nelmt ] == ibm->nv3[ nii ]))
	    
	    || ((ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || ibm->nv1[ nelmt ] == ibm->nv3[ nii ])
	    &&(ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || ibm->nv3[ nelmt ] == ibm->nv2[ nii ] || ibm->nv3[ nelmt ] == ibm->nv3[ nii ]))
	    
	    
	    || ((ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || ibm->nv2[ nelmt ] == ibm->nv3[ nii ])
	    &&(ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || ibm->nv3[ nelmt ] == ibm->nv2[ nii ] || ibm->nv3[ nelmt ] == ibm->nv3[ nii ])))
	    
	  */
	  //---------- interpolation between all around vertics (more than 8 neighbor elemts) ************
	  if(ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || 
	     ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || 
	     ibm->nv1[ nelmt ] == ibm->nv3[ nii ] || 
	     ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || 
	     ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || 
	     ibm->nv2[ nelmt ] == ibm->nv3[ nii ] || 
	     ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || 
	     ibm->nv3[ nelmt ] == ibm->nv2[ nii ] || 
	     ibm->nv3[ nelmt ] == ibm->nv3[ nii ])       {
	    riter = riter + 1.;
	    ibm->pres[ nelmt ] = ( ibm->pres[ nii ] + ibm->pres[ nelmt ] * ( riter - 1. ) ) / riter;
	    ibm->tau0[ nelmt ] = ( ibm->tau0[ nii ] + ibm->tau0[ nelmt ] * ( riter - 1. ) ) / riter;
	    ibm->tauN[ nelmt ] = ( ibm->tauN[ nii ] + ibm->tauN[ nelmt ] * ( riter - 1. ) ) / riter;
	    ibm->Bvel_u[ nelmt] = ( ibm->Bvel_u[ nii ] + ibm->Bvel_u[ nelmt ] * ( riter - 1. ) ) / riter;
	    ibm->Bvel_v[ nelmt] = ( ibm->Bvel_v[ nii ] + ibm->Bvel_v[ nelmt ] * ( riter - 1. ) ) / riter;
	    ibm->Bvel_w[ nelmt] = ( ibm->Bvel_w[ nii ] + ibm->Bvel_w[ nelmt ] * ( riter - 1. ) ) / riter;
	  }
	}
      }
    }
  }

  //Smoothing the variables at all cells to correct the high gradient regions
  for(i=0;i<2;i++){
    for( nelmt = 0; nelmt < ibm->n_elmt; nelmt++ )    {
      riter = 0.;
      for(nii = 0; nii < ibm->n_elmt; nii++)	{
	if(nii!=nelmt)  {	  	  
	  //---------- interpolation between all around vertics (more than 8 neighbor elemts) ************
	  if(ibm->nv1[ nelmt ] == ibm->nv1[ nii ] || 
	     ibm->nv1[ nelmt ] == ibm->nv2[ nii ] || 
	     ibm->nv1[ nelmt ] == ibm->nv3[ nii ] || 
	     ibm->nv2[ nelmt ] == ibm->nv1[ nii ] || 
	     ibm->nv2[ nelmt ] == ibm->nv2[ nii ] || 
	     ibm->nv2[ nelmt ] == ibm->nv3[ nii ] || 
	     ibm->nv3[ nelmt ] == ibm->nv1[ nii ] || 
	     ibm->nv3[ nelmt ] == ibm->nv2[ nii ] || 
	     ibm->nv3[ nelmt ] == ibm->nv3[ nii ])   	  {

	    riter = riter + 1.;
	    ibm->pres[ nelmt ] = ( ibm->pres[ nii ] + ibm->pres[ nelmt ] * ( riter - 1. ) ) / riter;
	    ibm->tau0[ nelmt ] = ( ibm->tau0[ nii ] + ibm->tau0[ nelmt ] * ( riter - 1. ) ) / riter;
	    ibm->tauN[ nelmt ] = ( ibm->tauN[ nii ] + ibm->tauN[ nelmt ] * ( riter - 1. ) ) / riter;
	    ibm->Bvel_u[ nelmt] = ( ibm->Bvel_u[ nii ] + ibm->Bvel_u[ nelmt ] * ( riter - 1. ) ) / riter;
	    ibm->Bvel_v[ nelmt] = ( ibm->Bvel_v[ nii ] + ibm->Bvel_v[ nelmt ] * ( riter - 1. ) ) / riter;
	    ibm->Bvel_w[ nelmt] = ( ibm->Bvel_w[ nii ] + ibm->Bvel_w[ nelmt ] * ( riter - 1. ) ) / riter;

	  }	  
	}
      }      
    }
  }        
  
/*   // finalizing the projection of velocity vector and bed shear stress on body  surface  */
/*   PetscInt   nv1,nv2,nv3; */
/*   PetscReal  M_z, pr, dA; */
/*   PetscReal  F_z,F_x, F_y, r_x, r_y; */
/*   // calcualte moment as well */
/*   M_z= 0.; */
/*   F_z= 0.; */
/*   for( elmt = 0; elmt < n_elmt; elmt++ )   { */
/*     nfx = ibm->nf_x[ elmt ]; */
/*     nfy = ibm->nf_y[ elmt ]; */
/*     nfz = ibm->nf_z[ elmt ]; */
    
/*     ucx = ibm->Bvel_u[ elmt ]; */
/*     ucy = ibm->Bvel_v[ elmt ]; */
/*     ucz = ibm->Bvel_w[ elmt ]; */
    
    
/*     ibm->Bvel_u[ elmt ] = ucx - (ucx * nfx + ucy * nfy + ucz * nfz) * nfx; */
/*     ibm->Bvel_v[ elmt ] = ucy - (ucx * nfx + ucy * nfy + ucz * nfz) * nfy; */
/*     ibm->Bvel_w[ elmt ] = ucz - (ucx * nfx + ucy * nfy + ucz * nfz) * nfz; */
/*     //ibm->Shvel[ elmt ] = ibm->Shvel[elmt];  */
    
/*     nv1=ibm->nv1[elmt]; */
/*     nv2=ibm->nv2[elmt]; */
/*     nv3=ibm->nv3[elmt]; */
    
/*     r_x= ibm->cent_x[elmt]; */
/*     r_y= ibm->cent_y[elmt]; */
    
/*     dA = ibm->dA[elmt]; */

/*     nfx= ibm->nf_x[elmt]; */
/*     nfy= ibm->nf_y[elmt]; */
/*     nfz= ibm->nf_z[elmt]; */

/*     pr = (ibm->pres[elmt]); */

/*     F_x = -pr*dA*nfx; */
/*     F_y = -pr*dA*nfy; */

/*     ibm->nt_z[elmt] =  (r_x*F_y - r_y*F_x)/dA; */
/*     ibm->ns_x[elmt] =  (r_x*F_y - r_y*F_x); */

/*     F_z += -pr*dA*nfz; */
/*     M_z +=   r_x*F_y - r_y*F_x; */
/*   } */
  
/*   PetscPrintf(PETSC_COMM_WORLD, "The Moment %le Force %le\n", M_z, F_z);  */

  return ( 0 );
}

PetscErrorCode Output_Projecting( IBMNodes * ibm )
{
  FILE            *f;
  char            filen[80];
  PetscInt        i, rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) {
    PetscReal       length_scale,specific_weight=1., Re,Vel=1.,temp;
    
    PetscBool dyn=PETSC_FALSE,Hg=PETSC_FALSE,OSI=PETSC_FALSE;
    
    PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-real_chact_leng", &length_scale, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-specific_weight", &specific_weight, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-ren", &Re, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-real_chact_vel", &Vel, PETSC_NULL);
    PetscOptionsGetBool(PETSC_NULL, PETSC_NULL,"-dynes",&dyn,PETSC_NULL);
    PetscOptionsGetBool(PETSC_NULL, PETSC_NULL,"-Hg",&Hg,PETSC_NULL);
    PetscOptionsGetBool(PETSC_NULL, PETSC_NULL,"-osi",&OSI,PETSC_NULL);
    
    sprintf(filen, "stress.dat");
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z, Pressure, tau0, tauN,u,v,w, M, rxp\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE VARLOCATION=([4-11]=CELLCENTERED)\n", ibm->n_v, ibm->n_elmt);
    PetscPrintf(PETSC_COMM_WORLD, "stress write \n");
    // x - component
    for (i=0; i<ibm->n_v; i++) {
      
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e \n", ibm->x_bp[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    // y - component
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e \n", ibm->y_bp[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    // z - component
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%e \n", ibm->z_bp[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
    
    
    //Write out the Pressure
    for (i=0; i<ibm->n_elmt; i++) 
      {
	//1000 is specific weight of blood = appox water
	temp = ibm->pres[i] * specific_weight * Vel * Vel;
	if (Hg)
	  temp =temp*760/101325;
	
	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", temp);
	
      }
    for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->tau0[i]);
    }

    for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->tauN[i]);
    }

    for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->Bvel_u[i]);
    }

    for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->Bvel_v[i]);
    }

    for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->Bvel_w[i]);
    }

    for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->nt_z[i]);
    }

    for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%le\n", ibm->ns_x[i]);
    }
	
    PetscFPrintf(PETSC_COMM_WORLD, f, "\n");

    
    //Write out the link nodes
    for (i=0; i<ibm->n_elmt; i++) 
      {
	
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1);
      }   
    
   
    fclose(f);

  }  
  
  return(0);  
}

PetscErrorCode Projecting_VTKOut(IBMNodes *ibm, PetscInt ibi, PetscInt ti)
{
    // vtk file name
  PetscInt n_cells=3;
  PetscInt rank,i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "surface%2.2d_%5.5d.vtk", ibi,ti);
    f = fopen(filen, "w"); // open file

    PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 2.0\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "Surface Grid\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "ASCII\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET UNSTRUCTURED_GRID\n");
   
    //    PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  5523993 double\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  %d float\n",ibm->n_v);
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->x_bp0[i],ibm->y_bp0[i],ibm->z_bp0[i]);
    }

    PetscFPrintf(PETSC_COMM_WORLD, f, "CELLS %d %d\n",ibm->n_elmt, (n_cells+1)*ibm->n_elmt);
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD,f, "%d  %d %d %d\n",n_cells, ibm->nv1[i], ibm->nv2[i], ibm->nv3[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_TYPES %d\n",ibm->n_elmt);
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD,f, "%d\n",5);
    }

    PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n", ibm->n_v);
    PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS u float\n");
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->u[i].x,ibm->u[i].y,ibm->u[i].z);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_DATA %d\n", ibm->n_elmt);
    PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS nf float\n");
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->nf_x[i],ibm->nf_y[i],ibm->nf_z[i]);
    }
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS nt float\n"); */
/*     for (i=0; i<ibm->n_elmt; i++) { */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->nt_x[i],ibm->nt_y[i],ibm->nt_z[i]); */
/*     } */
    PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS p double\n");    
    PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n");
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", ibm->pres[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS tau_w double\n");    
    PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n");
    for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", ibm->tau0[i]);
    }

    PetscFPrintf(PETSC_COMM_WORLD, f, "SCALARS tau_n double\n");    
    PetscFPrintf(PETSC_COMM_WORLD, f, "LOOKUP_TABLE default\n");
    for (i=0; i<ibm->n_elmt; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n", ibm->tauN[i]);
    }

    PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS u_proj float\n");
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->Bvel_u[i],ibm->Bvel_v[i],ibm->Bvel_w[i]);
    }
  }
  return(0);
}

PetscErrorCode Projecting_Power(IBMNodes *ibm, FSInfo *fsi,
				PetscInt ti)
{

  PetscInt   elmt, n_elmt=ibm->n_elmt;
  PetscInt   nv1,nv2,nv3;
  PetscReal  M_z, p, dA;
  PetscReal  n_x,n_y,n_z, F_x, F_y,F_z,Pw,Pw_y,r_x, r_y;
  PetscReal  nt_x,nt_y,nt_z;
  PetscReal  u_x,u_y,u_z;
  PetscReal  F_xSum,F_ySum,F_zSum,tauN, tau0; //Surface Force
  Cmpnts     Udef;
  PetscReal  Udef_n;

  M_z= 0.; Pw=0;  Pw_y=0.;
  F_xSum=0.; F_ySum=0.; F_zSum=0.;
  PetscInt rank=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  

  for (elmt=0; elmt<n_elmt; elmt++) {
    nv1=ibm->nv1[elmt];
    nv2=ibm->nv2[elmt];
    nv3=ibm->nv3[elmt];
    
    r_x= ibm->cent_x[elmt];
    r_y= ibm->cent_y[elmt];
    
    dA = ibm->dA[elmt];

    n_x= ibm->nf_x[elmt];
    n_y= ibm->nf_y[elmt];
    n_z= ibm->nf_z[elmt];

    u_x = (ibm->u[nv1].x + ibm->u[nv2].x + ibm->u[nv3].x)/3.;
    u_y = (ibm->u[nv1].y + ibm->u[nv2].y + ibm->u[nv3].y)/3.;
    u_z = (ibm->u[nv1].z + ibm->u[nv2].z + ibm->u[nv3].z)/3.;

    Udef.x= (ibm->Bvel_u[elmt]-u_x);
    Udef.y= (ibm->Bvel_v[elmt]-u_y);
    Udef.z= (ibm->Bvel_w[elmt]-u_z);
    Udef_n= Udef.x*n_x + Udef.y*n_y + Udef.z*n_z;
    //tangential defect
    Udef.x -= Udef_n * n_x;
    Udef.y -= Udef_n * n_y;
    Udef.z -= Udef_n * n_z;
    
    nt_x= Udef.x/(sqrt(Udef.x*Udef.x+Udef.y*Udef.y+Udef.z*Udef.z));
    nt_y= Udef.y/(sqrt(Udef.x*Udef.x+Udef.y*Udef.y+Udef.z*Udef.z));
    nt_z= Udef.z/(sqrt(Udef.x*Udef.x+Udef.y*Udef.y+Udef.z*Udef.z));
    
    p = (ibm->pres[elmt]);
    tauN=ibm->tauN[elmt] ;
    tau0=ibm->tau0[elmt] ;

    F_x = (-p+tauN)*dA*n_x + tau0*dA*nt_x;
    F_y = (-p+tauN)*dA*n_y + tau0*dA*nt_y ;
    F_z = (-p+tauN)*dA*n_z + tau0*dA*nt_z ;

    F_xSum += F_x;
    F_ySum += F_y;
    F_zSum += F_z;

    Pw += F_x*u_x + F_y*u_y + F_z*u_z;    
    Pw_y += F_y*u_y;    
  }
  
  fsi->Power= Pw;

  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "power_projected.dat");
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le\n", ti,Pw, Pw_y,F_xSum, F_ySum, F_zSum);
    fclose(f);
  }
  PetscPrintf(PETSC_COMM_WORLD, "The Power %le forces %le %le %le\n", Pw, F_xSum, F_ySum, F_zSum); 

  return(0); 
}

/* ==================================================================================             */
/* For strong-coupling */
/* ==================================================================================             */
PetscErrorCode Calc_FSI_pos_SC(FSInfo *FSinfo,IBMNodes *ibm, 
			       PetscReal dt, PetscReal dtime, 
			       PetscReal Re) 
{ 
  PetscInt     i,j;
  PetscInt     itr=23;
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z; //Forces and Area
  
  // Calc Forces
  //PetscPrintf(PETSC_COMM_WORLD, "n_elmt in calc_pos  %d\n", ibm->n_elmt);
  //Calc_forces(&FSinfo,&ibm,Re,ti);
  
  // init values
  for (i=0;i<6;i++) {
    S_new[i]=FSinfo->S_real[i];//FSinfo->S_new[i];
    S_real[i]=FSinfo->S_real[i];
    S_realm1[i]=FSinfo->S_realm1[i];
    //FSinfo->S_old[i]=FSinfo->S_new[i];
    //S_old[i]=S_new[i];
  }
  
  red_vel=FSinfo->red_vel;
  damp   =FSinfo->damp   ;
  mu_s   =FSinfo->mu_s   ;

  F_x = FSinfo->F_x;
  F_y = FSinfo->F_y;
  F_z = FSinfo->F_z;
/*   F_x = 0.5*(FSinfo->F_x + FSinfo->F_x_real); */
/*   F_y = 0.5*(FSinfo->F_y + FSinfo->F_y_real); */
/*   F_z = 0.5*(FSinfo->F_z + FSinfo->F_z_real); */

  PetscPrintf(PETSC_COMM_WORLD, "FSI  %le %le %le %le %le %le %le %le %le\n",red_vel,damp,mu_s,dt,dtime,F_y,F_z, S_real[2],S_realm1[2] );

  // solve lin mom equ
  for (i=0; i< itr;i++) { 

    for (j=0;j<6;j++) {
      S_old[j]=S_new[j];
    }

    if (dgf_x) {
    S_new[0]=S_new[0]-dt/2./dtime*
      (3.*S_new[0]-4.*S_real[0]+S_realm1[0])+S_new[1]*dt; // x
    S_new[1]=S_new[1]-dt/2./dtime*
         (3.*S_new[1]-4.*S_real[1]+S_realm1[1])+
      dt*(-2.*damp*(red_vel)*S_new[1]
	  -(red_vel*red_vel)*S_old[0]
	  + mu_s*F_x); // dx/dt
    }
    if (dgf_y) {
    S_new[2]=S_new[2]-dt/2./dtime*
      (3.*S_new[2]-4.*S_real[2]+S_realm1[2])+S_new[3]*dt; // y    
    // dy/dt
    S_new[3]=S_new[3]-dt/2./dtime*
         (3.*S_new[3]-4.*S_real[3]+S_realm1[3])+
      dt*(-2.*damp*(red_vel)*S_new[3]
	  -(red_vel*red_vel)*S_old[2]
	  + mu_s*F_y);
    }
    if (dgf_z) {
    S_new[4]=S_new[4]-dt/2./dtime*
      (3*S_new[4]-4*S_real[4]+S_realm1[4])+S_new[5]*dt; //z
    S_new[5]=S_new[5]-dt/2./dtime*(3*S_new[5]-4*S_real[5]+S_realm1[5])
              +dt*(-2.*damp*(red_vel)*S_new[5]
		   -(red_vel*red_vel)*(S_old[4])
		   + mu_s*F_z); //dz/dt
    }

    // FSI convergence
    //PetscPrintf(PETSC_COMM_SELF, "FSI convergence z: %le  u_z:%le\n", S_new[4]-S_old[4],S_new[5]-S_old[5]);
    //PetscPrintf(PETSC_COMM_WORLD, "FSI convergence y: %le  u_y:%le\n", S_new[2]-S_old[2],S_new[3]-S_old[3]);

  }
    
  // store results
  for (i=0;i<6;i++){
    FSinfo->S_new[i]=S_new[i];
  }

  // output values
  PetscPrintf(PETSC_COMM_WORLD, "z, dz/dt %le %le %le\n",S_new[4],S_new[5], F_z);
  PetscPrintf(PETSC_COMM_WORLD, "y, dy/dt %le %le %le\n",S_new[2],S_new[3], F_y);
  PetscPrintf(PETSC_COMM_WORLD, "x, dx/dt %le %le %le\n",S_new[0],S_new[1], F_x);

  return(0);
}

/* ==================================================================================             */
/*  integral equation For strong/weak-coupling */
/* ==================================================================================             */
PetscErrorCode Calc_FSI_pos_intg(FSInfo *fsi,IBMNodes *ibm,
				 PetscReal dt)				  
{ 
  PetscInt     i;
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];  
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z, F_mg=0.;//0.123; //Forces and Area
  PetscReal    w=0.75;


  // init values
  for (i=0;i<6;i++) {
    S_new[i]=0.;//fsi->S_real[i];//fsi->S_new[i];
    S_real[i]=fsi->S_real[i];
    S_realm1[i]=fsi->S_realm1[i];
    S_old[i]=fsi->S_old[i];    
  }
  
  red_vel=fsi->red_vel;
  damp   =fsi->damp   ;
  mu_s   =fsi->mu_s   ;

  F_mg = 0.;//8./6. * mu_s;

  F_x = 0.5*(fsi->F_x + fsi->F_x_real);
  F_y = 0.5*(fsi->F_y + fsi->F_y_real);
  F_z = 0.5*(fsi->F_z + fsi->F_z_real);
  
  PetscPrintf(PETSC_COMM_WORLD, "mass %le F_x %le F_y %le F_z %le \n",mu_s,F_x,F_y,F_z);
  
  PetscReal     L_x,L_y,L_z; 
  PetscReal	*x_bp =ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;
  PetscReal	xbp_min, ybp_min, zbp_min, xbp_max, ybp_max, zbp_max;
  PetscInt      n_v,mode=0;

  // solve lin mom equ.
  if (dgf_x) {
    S_new[1]=S_real[1]+ dt*(mu_s*F_x); // u=u_r + int(F/mdt)

    S_new[0]=S_real[0]+0.5*(S_new[1]+S_real[1])*dt; //x=x_r+u_avedt
  }
  if (dgf_y && mode<1) {
    S_new[3]=S_real[3]+ dt*(mu_s*F_y); // u=u_r + int(F/mdt)

    S_new[2]=S_real[2]+0.5*(S_new[3]+S_real[3])*dt; //y=y_r+u_avedt
  }
  if (dgf_z) {
    S_new[5]=S_real[5]+ dt*(mu_s*F_z-F_mg); // u=u_r + int(F/mdt)

    S_new[4]=S_real[4]+0.5*(S_new[5]+S_real[5])*dt; //x=x_r+u_avedt
  }

  // Relaxation
  if (STRONG_COUPLING) {
    fsi->atk=0.;
    for (i=1;i<6;i+=2){
      fsi->dS[i]=S_old[i]-S_new[i];
    
      if (fabs(fsi->dS[i]-fsi->dS_o[i])>1e-8  &&
     	  fsi->atk_o!=0.3) {
	fsi->atk+=(fsi->dS[i])/
	  (fsi->dS_o[i]-fsi->dS[i]);
      }
    }
    fsi->atk=fsi->atk_o+(fsi->atk_o-1)*fsi->atk;
    if (fsi->atk>.9) fsi->atk=.9;
    if (fsi->atk<-.2) fsi->atk=-0.2;
    
    w=(1.-fsi->atk);
    //if (rheology)  PetscOptionsGetReal(PETSC_NULL, "-w_str", &w, PETSC_NULL);
    PetscPrintf(PETSC_COMM_WORLD, "under relaxation coefficient %le \n",w);
    
    for (i=1;i<6;i+=2){
      S_new[i]=w*S_new[i]+(1.-w)*S_old[i];
      S_new[i-1]=S_real[i-1]+0.5*(S_new[i]+S_real[i])*dt;
    }
  }

  // if (rheology){
  fsi->x_c=fsi->a_c[0]+S_new[0];
  fsi->y_c=fsi->a_c[1]+S_new[2];
  fsi->z_c=fsi->a_c[2]+S_new[4];
  //  }
  
  // store results
  for (i=0;i<6;i++){
    fsi->S_new[i]=S_new[i];
  }
  
 for (i=0;i<3;i++){
   fsi->acc[i]=(S_new[2*i+1]-S_real[2*i+1])/dt;
  }

  // output values
 PetscPrintf(PETSC_COMM_WORLD, "z, dz/dt a_z %le %le %le \n",S_new[4],S_new[5],fsi->acc[2]);
 PetscPrintf(PETSC_COMM_WORLD, "y, dy/dt a_y %le %le %le \n",S_new[2],S_new[3],fsi->acc[1]);
 PetscPrintf(PETSC_COMM_WORLD, "x, dx/dt a_x %le %le %le \n",S_new[0],S_new[1],fsi->acc[0]);

  return(0);
}


PetscErrorCode calc_quarternion(FSInfo *FSinfo,PetscReal dt,PetscInt ibi)
{
  PetscInt     i,j,k;
 
  PetscReal    w_r[3],w_o[3],w_n[3];
  PetscReal    q_r[4],q_o[4],q_n[4];
  PetscReal    M_x,M_y,M_z; //Torque
  PetscReal    L_n[3],L_r[3],L_o[3]; //Angular momentum
  PetscReal    residue=0.0,w;
  PetscReal    dtaw=0.5*dt;
  PetscReal    rhs[4],R[3][3],RT[3][3],X[3][3],I_n[3][3];

  /* initialize values  */
  for (i=0;i<3;i++){
    w_r[i]=FSinfo->S_ang_r[i];
    w_o[i]=FSinfo->S_ang_r[i];
    w_n[i]=0.0;
    L_n[i]=0.0;
    L_o[i]=0.0;
    L_r[i]=FSinfo->L_r[i];
  }
  for (i=0;i<4;i++){
    q_r[i]=FSinfo->q_r[i];
    q_o[i]=FSinfo->q_r[i];
    q_n[i]=0.0;
  }

  
 // Trapezoidal rule
  M_x = 0.5*(FSinfo->M_x + FSinfo->M_x_real);
  M_y = 0.5*(FSinfo->M_y + FSinfo->M_y_real);
  M_z = 0.5*(FSinfo->M_z + FSinfo->M_z_real);

/*   PetscPrintf(PETSC_COMM_WORLD, "q[0] %le q[1] %le q[2] %le q[3] %le \n",FSinfo->q[0],FSinfo->q[1],FSinfo->q[2],FSinfo->q[3]); */
/*   PetscPrintf(PETSC_COMM_WORLD, "M_x %le M_y %le M_z %le \n",M_x,M_y,M_z); */
/*   PetscPrintf(PETSC_COMM_WORLD, "w_x_r %le w_y_r %le w_z_r %le \n",FSinfo->S_ang_r[1],FSinfo->S_ang_r[3],FSinfo->S_ang_r[5]); */


  
  //  PetscPrintf(PETSC_COMM_WORLD, "L_r[0] %le L_r[1] %le L_r[2] %le \n",L_r[0],L_r[1],L_r[2]);
  
  L_n[0]=dt*M_x+L_r[0];
  L_n[1]=dt*M_y+L_r[1];
  L_n[2]=dt*M_z+L_r[2];

  // Relaxation
  if (STRONG_COUPLING) {
   /*  if (fabs(FSinfo->M_x_old)<1.0e-10) FSinfo->M_x_old= FSinfo->M_x; */
/*     if (fabs(FSinfo->M_y_old)<1.0e-10) FSinfo->M_y_old= FSinfo->M_y; */
/*     if (fabs(FSinfo->M_z_old)<1.0e-10) FSinfo->M_z_old= FSinfo->M_z; */

    M_x = 0.5*((0.5*FSinfo->M_x+0.5*FSinfo->M_x_old) + FSinfo->M_x_real);
    M_y = 0.5*((0.5*FSinfo->M_y+0.5*FSinfo->M_y_old) + FSinfo->M_y_real);
    M_z = 0.5*((0.5*FSinfo->M_z+0.5*FSinfo->M_z_old) + FSinfo->M_z_real);

    //  PetscPrintf(PETSC_COMM_WORLD, "Stroung_Coupling works for quaternion \n");

    for (i=0;i<3;i++){
      L_r[i] = FSinfo->L_r[i];
      L_o[i] = FSinfo->L_o[i];
    }

    FSinfo->atk=0.;
    for (i=1;i<3;i++){
      FSinfo->dS[i]=L_o[i]-L_n[i];
    
      if (fabs(FSinfo->dS[i]-FSinfo->dS_o[i])>1e-8  &&  FSinfo->atk_o!=0.9) {
	FSinfo->atk+=(FSinfo->dS[i])/(FSinfo->dS_o[i]-FSinfo->dS[i]);
      }
    }
    FSinfo->atk=FSinfo->atk_o+(FSinfo->atk_o-1)*FSinfo->atk;
    if (FSinfo->atk>.9) FSinfo->atk=.9;
    if (FSinfo->atk<0.1) FSinfo->atk=0.1;
    w=1.-FSinfo->atk;
    PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-w_str", &w, PETSC_NULL);
    //   w=0.1;
    PetscPrintf(PETSC_COMM_WORLD, "under relaxation coefficient %le \n",w);
 
    for (i=0;i<3;i++){
      L_n[i]=w*L_n[i]+(1.-w)*L_o[i];
    }
  
  } // strong-coupling

  //  PetscPrintf(PETSC_COMM_WORLD, "L_x %le L_y  %le L_z %le  \n",L_n[0],L_n[1],L_n[2]);

  for (k=0;k<50;k++){ /*  while (residue>0.001){ */
    rhs[0]=0.25*(-w_o[0]*q_o[1]-w_o[1]*q_o[2]-w_o[2]*q_o[3])+0.25*(-w_r[0]*q_r[1]-w_r[1]*q_r[2]-w_r[2]*q_r[3]);
    rhs[1]=0.25*( w_o[0]*q_o[0]+w_o[1]*q_o[3]-w_o[2]*q_o[2])+0.25*( w_r[0]*q_r[0]+w_r[1]*q_r[3]-w_r[2]*q_r[2]);
    rhs[2]=0.25*(-w_o[0]*q_o[3]+w_o[1]*q_o[0]+w_o[2]*q_o[1])+0.25*(-w_r[0]*q_r[3]+w_r[1]*q_r[0]+w_r[2]*q_r[1]);
    rhs[3]=0.25*( w_o[0]*q_o[2]-w_o[1]*q_o[1]+w_o[2]*q_o[0])+0.25*( w_r[0]*q_r[2]-w_r[1]*q_r[1]+w_r[2]*q_r[0]);
    
    for (i=0;i<4;i++){
      q_n[i]=q_o[i]+dtaw*(rhs[i]+(q_r[i]-q_o[i])/dt);
    }
    //  PetscPrintf(PETSC_COMM_WORLD, "q[0] %le q[1] %le q[2] %le q[3] %le \n",q_n[0],q_n[1],q_n[2],q_n[3]);
    R[0][0]=1.0-2.0*q_n[2]*q_n[2]-2.0*q_n[3]*q_n[3];
    R[0][1]=2.0*(q_n[1]*q_n[2]-q_n[0]*q_n[3]);
    R[0][2]=2.0*(q_n[1]*q_n[3]+q_n[0]*q_n[2]);
    R[1][0]=2.0*(q_n[1]*q_n[2]+q_n[0]*q_n[3]);
    R[1][1]=1.0-2.0*q_n[1]*q_n[1]-2.0*q_n[3]*q_n[3];
    R[1][2]=2.0*(q_n[2]*q_n[3]-q_n[0]*q_n[1]);
    R[2][0]=2.0*(q_n[1]*q_n[3]-q_n[0]*q_n[2]);
    R[2][1]=2.0*(q_n[2]*q_n[3]+q_n[0]*q_n[1]);
    R[2][2]=1.0-2.0*q_n[1]*q_n[1]-2.0*q_n[2]*q_n[2];
    //  PetscPrintf(PETSC_COMM_WORLD, "R: R[0][0] %le R[0][1] %le R[0][2] %le R[1][1] %le R[1][2] %le R[2][2] %le \n",R[0][0],R[0][1],R[0][2],R[1][1],R[1][2],R[2][2]);

    for (i=0;i<3;i++){
      for (j=0;j<3;j++){
	if (i==j) RT[i][j]=R[i][j];
	else      RT[i][j]=R[j][i];
      }
    }
    // PetscPrintf(PETSC_COMM_WORLD, "I_inv: I_n[0][0] %le I_inv[0][1] %le I_inv[0][2] %le I_inv[1][1] %le I_inv[1][2] %le I_inv[2][2] %le \n",FSinfo->I_inv[0][0],FSinfo->I_inv[0][1],FSinfo->I_inv[0][2],FSinfo->I_inv[1][1],FSinfo->I_inv[1][2],FSinfo->I_inv[2][2]);
    //  Mat_Maltiply(X,R,RT);
    for (i=0;i<3;i++){
      for (j=0;j<3;j++){
	X[i][j]=FSinfo->I_inv[i][0]*RT[0][j]+FSinfo->I_inv[i][1]*RT[1][j]+FSinfo->I_inv[i][2]*RT[2][j];
      }
    }
    //  Mat_Maltiply(X,FSinfo->I_inv,RT);

    for (i=0;i<3;i++){
      for (j=0;j<3;j++){
	I_n[i][j]=R[i][0]*X[0][j]+R[i][1]*X[1][j]+R[i][2]*X[2][j];
      }
    }
    //    Mat_Maltiply(I_n,R,X);
    //  PetscPrintf(PETSC_COMM_WORLD, "I_n: I_n[0][0] %le I_n[0][1] %le I_n[0][2] %le I_n[1][1] %le I_n[1][2] %le I_n[2][2] %le \n",I_n[0][0],I_n[0][1],I_n[0][2],I_n[1][1],I_n[1][2],I_n[2][2]);
    w_n[0]=I_n[0][0]*L_n[0]+I_n[0][1]*L_n[1]+I_n[0][2]*L_n[2];
    w_n[1]=I_n[1][0]*L_n[0]+I_n[1][1]*L_n[1]+I_n[1][2]*L_n[2];
    w_n[2]=I_n[2][0]*L_n[0]+I_n[2][1]*L_n[1]+I_n[2][2]*L_n[2];

    residue=0.0;

    for (i=0;i<4;i++){
      residue +=(q_n[i]-q_o[i])*(q_n[i]-q_o[i]);
    }

    // PetscPrintf(PETSC_COMM_WORLD, "w_n[0] %le w_n[1]  %le w_n[2] %le  \n",w_n[0],w_n[1],w_n[2]);
    residue=sqrt(residue/4);
    //  PetscPrintf(PETSC_COMM_WORLD, "residue is  %le  \n",residue);
    for (i=0;i<4;i++){
      q_o[i]=q_n[i];
    }
    for (i=0;i<3;i++){
      w_o[i]=w_n[i];
    }
  }

  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      FSinfo->R[i][j]=R[i][j];
    }
  }
 
  //  PetscPrintf(PETSC_COMM_WORLD, "Rotation Matrix: R[0][0] %le R[0][1] %le R[0][2] %le R[1][1] %le R[1][2] %le R[2][2] %le \n",FSinfo->R[0][0],FSinfo->R[0][1],FSinfo->R[0][2],FSinfo->R[1][1],FSinfo->R[1][2],FSinfo->R[2][2]);

  // PetscPrintf(PETSC_COMM_WORLD, "Inverse of Moment of Inertia: I_n[0][0]  %le I_n[0][1] %le I_n[0][2] %le I_n[1][1] %le I_n[1][2] %le I_n[2][2] %le \n",I_n[0][0],I_n[0][1],I_n[0][2],I_n[1][1],I_n[1][2],I_n[2][2]);


  for (i=0;i<4;i++){
    FSinfo->q[i]=q_n[i];
  }
  
 /*  FSinfo->S_ang_n[1]=w_r[0]; */
/*   FSinfo->S_ang_n[3]=w_r[1]; */
/*   FSinfo->S_ang_n[5]=w_r[2]; */

  FSinfo->S_ang_n[1]=w_n[0];
  FSinfo->S_ang_n[3]=w_n[1];
  FSinfo->S_ang_n[5]=w_n[2];

  // FSinfo->S_ang_n[0]=FSinfo->S_ang_r[0]+0.5*(FSinfo->S_ang_n[1]+FSinfo->S_ang_r[1])*dt;
  
  FSinfo->alpha[0]=(FSinfo->S_ang_n[1]-FSinfo->S_ang_r[1])/dt;
  FSinfo->alpha[1]=(FSinfo->S_ang_n[3]-FSinfo->S_ang_r[3])/dt;
  FSinfo->alpha[2]=(FSinfo->S_ang_n[5]-FSinfo->S_ang_r[5])/dt;

  for (i=0;i<3;i++){
    FSinfo->L_n[i]=L_n[i];
  }
  
  // output values
  PetscPrintf(PETSC_COMM_WORLD, "w_x %le w_y  %le w_z %le  \n",FSinfo->S_ang_n[1],FSinfo->S_ang_n[3],FSinfo->S_ang_n[5]);
  PetscPrintf(PETSC_COMM_WORLD, "a_x %le a_y  %le a_z %le  \n",FSinfo->alpha[0],FSinfo->alpha[1],FSinfo->alpha[2]);
  PetscPrintf(PETSC_COMM_WORLD, "q[0] %le q[1] %le q[2] %le q[3] %le \n",FSinfo->q[0],FSinfo->q[1],FSinfo->q[2],FSinfo->q[3]);
  return 0;
}


/* ==================================================================================             */
/*   Moving the IBM body   */
/* ==================================================================================             */
PetscErrorCode Elmt_Move_FSI_TRANS(FSInfo *FSinfo, IBMNodes *ibm,
				   PetscInt ibi)
{
  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;
  PetscInt i;

  PetscPrintf(PETSC_COMM_WORLD, "MOVE BODY x: %le  y:%le z:%le\n", FSinfo->S_new[0],FSinfo->S_new[2],FSinfo->S_new[4]);

  for (i=0; i<n_v; i++) {
    // change for stat case 4/9/06 iman
    ibm->x_bp[i] = ibm->x_bp0[i]+(FSinfo->S_new[0]);//-FSinfo->S_real[0]);
    ibm->y_bp[i] = ibm->y_bp0[i]+(FSinfo->S_new[2]);//-FSinfo->S_real[2]);
    ibm->z_bp[i] = ibm->z_bp0[i]+(FSinfo->S_new[4]);//-FSinfo->S_real[4]);
  }
  
  return(0);
}

/* ==================================================================================             */
PetscErrorCode Elmt_Move_FSI_ROT(FSInfo *FSinfo, IBMNodes *ibm, 
				 PetscReal dt, PetscInt ibi)
{
  PetscInt n_v = ibm->n_v, n_elmt = ibm->n_elmt;
  PetscReal  x_c=FSinfo->x_c, y_c=FSinfo->y_c, z_c=FSinfo->z_c;
  PetscReal  rot_x=FSinfo->S_ang_n[0];//-FSinfo->S_ang_o[0];
  PetscInt i;
  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  PetscReal rx,ry,rz;
  PetscReal wx=FSinfo->S_ang_n[1], wy=FSinfo->S_ang_n[3], wz=FSinfo->S_ang_n[5];
  
  // only rot around x
  wy=0.;wz=0.;

  PetscPrintf(PETSC_COMM_WORLD, "Body Rotate: %le Center %le %le %le\n",rot_x, x_c,y_c,z_c);
  
  for (i=0; i<n_v; i++) {      
    ibm->x_bp[i] = ibm->x_bp0[i];
    ibm->z_bp[i] = z_c + (ibm->z_bp0[i]-z_c)*cos(rot_x) - (ibm->y_bp0[i]-y_c)*sin(rot_x);
    ibm->y_bp[i] = y_c + (ibm->z_bp0[i]-z_c)*sin(rot_x) + (ibm->y_bp0[i]-y_c)*cos(rot_x);
  }
  
 
  for (i=0; i<n_elmt; i++) {
    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
    dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e]; 
    dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e]; 
    dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e]; 
    
    dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e]; 
    dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e]; 
    dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e]; 
    
    ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
    ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
    ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;

    dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] + 
	      ibm->nf_z[i]*ibm->nf_z[i]);
    
    ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;
    
    /*       PetscPrintf(PETSC_COMM_WORLD, "NFZ %d %d %d %d %e\n", i, ibm->nv1[i], ibm->nv2[i], ibm->nv3[i], ibm->nf_z[i]); */
    //      PetscPrintf(PETSC_COMM_WORLD, "%le %le %le %le %le %le\n", x_bp[n1e], y_bp[n1e], ibm->x_bp0[n1e], ibm->y_bp0[n1e], x_bp[n3e], y_bp[n3e]);
    // Addedd 5/7/06 iman
    // ns = nf x k
    if ((((1.-ibm->nf_z[i])<=1e-6 )&&((-1.+ibm->nf_z[i])<1e-6))||
	(((ibm->nf_z[i]+1.)<=1e-6 )&&((-1.-ibm->nf_z[i])<1e-6))) {
      ibm->ns_x[i] = 1.;
      ibm->ns_y[i] = 0.;
      ibm->ns_z[i] = 0 ;

      // nt = ns x nf
      ibm->nt_x[i] = 0.;
      ibm->nt_y[i] = 1.;
      ibm->nt_z[i] = 0.;
      } else {
      ibm->ns_x[i] =  ibm->nf_y[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->ns_y[i] = -ibm->nf_x[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->ns_z[i] = 0 ;

      // nt = ns x nf
      ibm->nt_x[i] = -ibm->nf_x[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->nt_y[i] = -ibm->nf_y[i]*ibm->nf_z[i]/ sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      ibm->nt_z[i] = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i]);
      }

  }

  // 1st order approx. not good!
  for (i=0; i<n_v; i++) {
/*       ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / dt; */
/*       ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / dt; */
/*       ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / dt; */
    rx = ibm->x_bp[i]-x_c;
    ry = ibm->y_bp[i]-y_c;
    rz = ibm->z_bp[i]-z_c;      
    ibm->u[i].x =   ry*wz-wy*rz  ;
    ibm->u[i].y =-( rx*wz-wx*rz );
    ibm->u[i].z =   rx*wy-wx*ry  ;      
  }
  
  return(0);
}


PetscErrorCode Elmt_Move_FSI_TRANSpROT_quarternion(IBMNodes *ibm, FSInfo *fsi)
{
  PetscInt n_v = ibm->n_v;
  PetscInt i,j;
  Cmpnts   p,a_c;
  PetscReal Rot[3][3];

   
  a_c.x=fsi->a_c[0];
  a_c.y=fsi->a_c[1];
  a_c.z=fsi->a_c[2];

  PetscPrintf(PETSC_COMM_WORLD, " initial center of mass: a_c.x %le a_c.y %le a_c.z %le \n",a_c.x,a_c.y,a_c.z);
 
  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      Rot[i][j]=fsi->R[i][j];
      //  PetscPrintf(PETSC_COMM_WORLD, "R[%d][%d]= %le \n",i,j,Rot[i][j]);
    }
  }
 
  for (i=0; i<n_v; i++) {
   
    p.x = ibm->x_bp0[i]-fsi->a_c[0];
    p.y = ibm->y_bp0[i]-fsi->a_c[1];
    p.z = ibm->z_bp0[i]-fsi->a_c[2];
  
    ibm->x_bp[i] =Rot[0][0]*p.x+Rot[0][1]*p.y+Rot[0][2]*p.z+fsi->x_c;
    ibm->y_bp[i] =Rot[1][0]*p.x+Rot[1][1]*p.y+Rot[1][2]*p.z+fsi->y_c;
    ibm->z_bp[i] =Rot[2][0]*p.x+Rot[2][1]*p.y+Rot[2][2]*p.z+fsi->z_c;
  }


  //  PetscPrintf(PETSC_COMM_WORLD, "fsi[%d].S_new[0] %le \n",NumberOfBodies+particle_added-1,fsi[NumberOfBodies+particle_added-1].S_new[0]);


  return(0);
}

/* ==================================================================================             */
/*  Set velocity of the ibm  */
/* ==================================================================================             */
PetscErrorCode Set_ibm_velocity_to_fsi(FSInfo *fsi, IBMNodes *ibm,
				       PetscInt ibi,
				       PetscInt addtranslation,
				       PetscInt addrotation)
{
  PetscInt n_v = ibm->n_v;
  PetscInt i;
  Cmpnts	omega_c, a_c;
  PetscReal     rx,ry,rz;

  for (i=0; i<n_v; i++) {
    if (addtranslation) {
      ibm->u[i].x = fsi->S_new[1];
      ibm->u[i].y = fsi->S_new[3];
      ibm->u[i].z = fsi->S_new[5];
    } else {
      ibm->u[i].x = 0.;
      ibm->u[i].y = 0.;
      ibm->u[i].z = 0.;
    }
  }

  if (addrotation){
    omega_c.x=fsi->S_ang_n[1];
    omega_c.y=fsi->S_ang_n[3];
    omega_c.z=fsi->S_ang_n[5];
    
    // axis of frame rotation
    a_c.x=fsi->x_c;
    a_c.y=fsi->y_c;
    a_c.z=fsi->z_c;
    
    for (i=0; i<n_v; i++) {
      rx = ibm->x_bp[i]-a_c.x;
      ry = ibm->y_bp[i]-a_c.y;
      rz = ibm->z_bp[i]-a_c.z;
      ibm->u[i].x +=   (rz*omega_c.y-omega_c.z*ry)  ;
      ibm->u[i].y +=   (rx*omega_c.z-omega_c.x*rz)  ;
      ibm->u[i].z +=   (ry*omega_c.x-omega_c.y*rx)  ;   
    }
  }
  
  return(0);
}

