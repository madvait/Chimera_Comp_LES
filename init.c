#include "variables.h"
//extern PetscInt block_number, inletprofile, blank;
//extern PetscReal L_dim;
extern PetscInt  moveframe,rotateframe;
extern PetscInt  les;
extern int averaging;

//PetscErrorCode InflowWaveFormRead(UserCtx *user);
/* PetscErrorCode GridRestriction(PetscInt i, PetscInt j, PetscInt k, */
/* 			       PetscInt *ih, PetscInt *jh, PetscInt *kh, */
/* 			       UserCtx *user); */
PetscErrorCode FormMetrics(UserCtx *user);
//PetscErrorCode Update_Metrics_PBC(UserCtx *user);
PetscErrorCode MetricsDivergence(UserCtx *user);

PetscErrorCode GridRestriction(PetscInt i, PetscInt j, PetscInt k,
			       PetscInt *ih, PetscInt *jh, PetscInt *kh,
			       UserCtx *user)
{
  if ((user->isc)) {
    *ih = i;
  }
  else {
    *ih = 2 * i;
  }

  if ((user->jsc)) {
    *jh = j;
  }
  else {
    *jh = 2 * j;
  }

  if ((user->ksc)) {
    *kh = k;
  }
  else {
    *kh = 2 * k;
  }

  return 0;
}


PetscErrorCode FormInitialize(UserCtx *user)
{
  user->assignedA = PETSC_FALSE;
  user->multinullspace=PETSC_FALSE;

  VecSet(user->Nvert,0.);
  VecSet(user->lNvert, 0.);


  if (user->thislevel==user->mglevels-1) {
  VecSet(user->P,0.);
  VecSet(user->lP,0.);
  VecSet(user->Nvert_o,0.);
  VecSet(user->lNvert_o, 0.);

/*   VecSet(user->Ucont,0.); */
/*   VecSet(user->lUcont,0.); */
/*   VecSet(user->Ucat,0.); */
/*   VecSet(user->Bcs.Ubcs,0.); */
/*   VecSet(user->DUold, 0.); */
/*   if (moveframe || rotateframe) { */
/*     VecSet(user->Vcont,0.); */
/*     VecSet(user->lVcont,0.); */
/*   }	   */
/*   if (rotateframe) { */
/*     VecSet(user->Wcat,0.); */
/*     VecSet(user->lWcat,0.); */
/*   } */
  } 

  user->FluxIntpSum=0.;
  user->FluxInSum=0.;
  user->FluxIntfcSum=0.;
  user->FluxOutSum=0.;
  user->ratio=0.;

  user->ren = 3000.;
  user->Ma = 2.81;
  user->gamma=1.4;
  user->Pr=0.72;
  user->dt = 0.01;
  user->cfl=0.2;
  user->vnn=0.2;
 

  PetscOptionsGetReal(NULL, PETSC_NULL, "-ren", &user->ren, PETSC_NULL);
  PetscOptionsGetReal(NULL, PETSC_NULL, "-Ma", &user->Ma, PETSC_NULL);
  PetscOptionsGetReal(NULL, PETSC_NULL, "-Pr", &user->Pr, PETSC_NULL);
  PetscOptionsGetReal(NULL, PETSC_NULL, "-gamma", &user->gamma, PETSC_NULL);
  PetscOptionsGetReal(NULL, PETSC_NULL, "-dt", &user->dt, PETSC_NULL);
  PetscOptionsGetReal(NULL, PETSC_NULL, "-cfl", &user->cfl, PETSC_NULL);
  PetscOptionsGetReal(NULL, PETSC_NULL, "-vnn", &user->vnn, PETSC_NULL);

  user->st = 1.;//0.038406145;
  
  PetscPrintf(PETSC_COMM_WORLD, "Re %le St %le Ma %le Pr %le gamma %\le dt %le \n",user->ren,user->st,user->Ma, user->Pr, user->gamma, user->dt);

  return(0);
}

PetscErrorCode MGDACreate(UserMG *usermg, PetscInt bi)
{
  DMBoundaryType  xperiod,yperiod,zperiod;

  MGCtx *mgctx = usermg->mgctx;

  UserCtx *user, *user_high;

  PetscErrorCode ierr;

  PetscInt l;
  
  // claculalte grid size for coarser levels
  for (l=usermg->mglevels-2; l>=0; l--) {
    user = mgctx[l].user;
    user_high = mgctx[l+1].user;
         
    if (user[bi].isc) {	 
      user[bi].IM = user_high[bi].IM;
    }
    else {
      user[bi].IM = (user_high[bi].IM + 1) / 2;
    }
    if (user[bi].jsc) {
      user[bi].JM = user_high[bi].JM;
    }
    else {
      user[bi].JM = (user_high[bi].JM + 1) / 2;
    }

    if (user[bi].ksc) {
      user[bi].KM = user_high[bi].KM;
    }
    else {
      user[bi].KM = (user_high[bi].KM + 1) / 2;
    }
    PetscPrintf(PETSC_COMM_WORLD, "Grid dimension in i,j and k direction for level %d are %d %d %d :\n",l,user[bi].IM, user[bi].JM, user[bi].KM);
    if (user[bi].IM*(2- (user[bi].isc))-(user_high[bi].IM+1-(user[bi].isc)) +
	user[bi].JM*(2- (user[bi].jsc))-(user_high[bi].JM+1-(user[bi].jsc)) +
	user[bi].KM*(2- (user[bi].ksc))-(user_high[bi].KM+1-(user[bi].ksc))) {
      PetscPrintf(PETSC_COMM_WORLD, "Grid at level %d can't be further restricted!\n", l);
    } //end if
  } // end for l

  l = 0;
  user = mgctx[l].user;

  PetscInt M,N,P,m,n,p,s,MM,NN,PP, *lx, *ly, *lz;
  PetscInt *lxSum, *lySum, *lzSum;

  //stencil width
  s=2;
 
  if(user[bi].bctype[0]==7  && user[bi].bctype[2]==7  && user[bi].bctype[4]==7) {
    xperiod = DM_BOUNDARY_PERIODIC;
    yperiod = DM_BOUNDARY_PERIODIC;
    zperiod = DM_BOUNDARY_PERIODIC;
  }
  
  else if(user[bi].bctype[0]==7  && user[bi].bctype[2]==7) {
    xperiod = DM_BOUNDARY_PERIODIC;
    yperiod = DM_BOUNDARY_PERIODIC;
    zperiod = DM_BOUNDARY_NONE;
  }
  else if(user[bi].bctype[2]==7  && user[bi].bctype[4]==7) {
    xperiod = DM_BOUNDARY_NONE;
    yperiod = DM_BOUNDARY_PERIODIC;
    zperiod = DM_BOUNDARY_PERIODIC;
  }
  else if(user[bi].bctype[0]==7  && user[bi].bctype[4]==7) {
    xperiod = DM_BOUNDARY_PERIODIC;
    yperiod = DM_BOUNDARY_NONE;
    zperiod = DM_BOUNDARY_PERIODIC;
  }
  else if(user[bi].bctype[0]==7){
    xperiod = DM_BOUNDARY_PERIODIC;
    yperiod = DM_BOUNDARY_NONE;
    zperiod = DM_BOUNDARY_NONE;
  }
  else if(user[bi].bctype[2]==7) {
    xperiod = DM_BOUNDARY_NONE;
    yperiod = DM_BOUNDARY_PERIODIC;
    zperiod = DM_BOUNDARY_NONE;
  }
  else if(user[bi].bctype[4]==7 ) {
    xperiod = DM_BOUNDARY_NONE;
    yperiod = DM_BOUNDARY_NONE; 
    zperiod = DM_BOUNDARY_PERIODIC;
  }
  else {
    s=2;
    xperiod = DM_BOUNDARY_NONE;
    yperiod = DM_BOUNDARY_NONE;
    zperiod = DM_BOUNDARY_NONE;
  }

  // global dimension in each direction of the array
  M=user[bi].IM;
  N=user[bi].JM;
  P=user[bi].KM;

  ierr=DMDACreate3d(PETSC_COMM_WORLD,xperiod,yperiod,zperiod,DMDA_STENCIL_BOX,M,N,P,
		    PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,s,PETSC_NULL, PETSC_NULL,
		    PETSC_NULL,&(user[bi].da)); //Mohsen Jan 12  
  PetscPrintf(PETSC_COMM_WORLD, "DA is Created for coarsest level and global dimension in i,j and k direction are %d %d %d :\n",M,N,P); 
  DMSetFromOptions(user[bi].da);
  DMSetUp(user[bi].da);

  DMDAGetInfo(user[bi].da, PETSC_NULL, &MM,&NN,&PP,&m, &n, &p, PETSC_NULL, PETSC_NULL,PETSC_NULL, PETSC_NULL,PETSC_NULL,PETSC_NULL);
  if (ierr) 
    SETERRQ1(PETSC_COMM_WORLD,1, "problem creating DA %d",ierr);
  PetscInt compressiblevarNum = 5;
  ierr=DMDACreate3d(PETSC_COMM_WORLD,xperiod,yperiod,zperiod,DMDA_STENCIL_BOX,M, N, P,
		    m, n, p,compressiblevarNum,s, PETSC_NULL, PETSC_NULL, PETSC_NULL,&(user[bi].cda));
  DMSetFromOptions(user[bi].cda);
  DMSetUp(user[bi].cda);
  PetscPrintf(PETSC_COMM_WORLD, "DA is Created for coarsest level and global dimension in i,j and k direction are %d %d %d and dof %d proc  %d %d %d:\n",MM,NN,PP,compressiblevarNum,m,n,p); 
 
  if (ierr) 
    SETERRQ1(PETSC_COMM_WORLD,1, "problem creating DA %d",ierr);
  user[bi].aotopetsc = PETSC_FALSE;
  DMDASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  DMGetCoordinateDM(user[bi].da, &(user[bi].fda));
  DMDAGetLocalInfo(user[bi].da, &(user[bi].info));

  /* DMDALocalInfo	info = user[bi].info; */
  /* PetscInt	xs = info.xs, xe = info.xs + info.xm; */
  /* PetscInt      ys = info.ys, ye = info.ys + info.ym; */
  /* PetscInt	zs = info.zs, ze = info.zs + info.zm; */
  /* PetscInt	mx = info.mx, my = info.my, mz = info.mz; */
  /* PetscPrintf(PETSC_COMM_WORLD, "READ info 1 error xe,ye,ze mx,my,mz: %d %d %d %d  %d %d %d %d %d %d\n",ierr,xs,ys,zs, xe,ye,ze, mx,my,mz); */

  
  if (rans && usermg->mglevels-1==0) { // fda2 for k-omega at the finest level
    DMDAGetInfo(user[bi].da, PETSC_NULL, &MM,&NN,&PP,&m, &n, &p, PETSC_NULL, PETSC_NULL,PETSC_NULL, PETSC_NULL,PETSC_NULL,PETSC_NULL);
 
    ierr=DMDACreate3d(PETSC_COMM_WORLD,xperiod,yperiod,zperiod,DMDA_STENCIL_BOX,user[bi].IM, user[bi].JM, user[bi].KM,
		    m,n,p,2,s, PETSC_NULL, PETSC_NULL, PETSC_NULL,&(user[bi].fda2));
    if (ierr) SETERRQ1(PETSC_COMM_SELF,1, "problem creating FDA2 for RANS%d",ierr);
  }

  // for other levels
  for (l=1; l<usermg->mglevels; l++) {
    user = mgctx[l].user;
    
    // Get info about the coarse DA
    DMDAGetInfo(mgctx[l-1].user[bi].da, PETSC_NULL, &MM,&NN,&PP,&m,&n,&p,PETSC_NULL,PETSC_NULL,PETSC_NULL, PETSC_NULL,PETSC_NULL,PETSC_NULL);
    if (l==1) PetscPrintf(PETSC_COMM_WORLD, "Corresponing numberof procs in i,j and k direction are %d %d %d :\n",m,n,p);
    // Find the distribution of the coarse DA
    PetscMalloc(m*sizeof(PetscInt), &lx);
    PetscMalloc(n*sizeof(PetscInt), &ly);
    PetscMalloc(p*sizeof(PetscInt), &lz);
    PetscMalloc(m*sizeof(PetscInt), &lxSum);
    PetscMalloc(n*sizeof(PetscInt), &lySum);
    PetscMalloc(p*sizeof(PetscInt), &lzSum);

    PetscInt rank;
    for (rank=0; rank<m; rank++) lx[rank]=0;
    for (rank=0; rank<n; rank++) ly[rank]=0;
    for (rank=0; rank<p; rank++) lz[rank]=0;

    DMDALocalInfo	info = mgctx[l-1].user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt    ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
    // Calculationg the distribution of the fine DA 
    // such that the refined portion of the coarse DA is on the
    // same processor as the coarse DA
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if ( user[bi].isc) 
      lx[rank % m]= (xe-xs);
    else {
      if (m==1)
	lx[rank % m]= user[bi].IM+1;
      else if (xs==0) 
	lx[rank % m]= 2*xe-1;
      else if (xe==mx)
	lx[rank % m]= user[bi].IM+1-2*xs+1;
      else
	lx[rank % m]= (xe-xs)*2;
    }
    if ( user[bi].jsc) 
      ly[(rank % (m*n))/m]= (ye-ys);
    else {
      if (n==1)
	ly[(rank % (m*n))/m]= user[bi].JM+1;
      else if (ys==0)
	ly[(rank % (m*n))/m]= 2*ye-1;
      else if (ye==my) 
	ly[(rank % (m*n))/m]= user[bi].JM+1-2*ys+1;
      else
	ly[(rank % (m*n))/m]= (ye-ys)*2;
    }
    if ( user[bi].ksc) 
      lz[rank/(m*n)]=ze-zs;
    else{
      if (p==1)
	lz[rank/(m*n)]= user[bi].KM+1;
      else if (zs==0)
	lz[rank/(m*n)]= 2*ze-1;
      else if (ze==mz)
	lz[rank/(m*n)]= user[bi].KM+1-2*zs+1;
      else
	lz[rank/(m*n)]= (ze-zs)*2;
    }
    
    // All_reduce to all processes
    MM=0;
    for (rank=0; rank<m; rank++) {
      MPI_Allreduce(&lx[rank],&lxSum[rank],1,MPIU_INT,MPI_MAX,PETSC_COMM_WORLD);
      MM+=lxSum[rank];
    }
    NN=0;
    for (rank=0; rank<n; rank++) {
      MPI_Allreduce(&ly[rank],&lySum[rank],1,MPIU_INT,MPI_MAX,PETSC_COMM_WORLD);      
      NN+=lySum[rank];
    }
    PP=0;
    for (rank=0; rank<p; rank++) {
      MPI_Allreduce(&lz[rank],&lzSum[rank],1,MPIU_INT,MPI_MAX,PETSC_COMM_WORLD);      
      PP+=lzSum[rank];
    }
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    PetscBarrier(PETSC_NULL);
    
    // Create the refined DA based on the new distribution
   
    ierr=DMDACreate3d(PETSC_COMM_WORLD,xperiod,yperiod,zperiod,DMDA_STENCIL_BOX,user[bi].IM,user[bi].JM,user[bi].KM,m,n,p,1,s,lxSum,lySum,lzSum,&(user[bi].da));
    PetscPrintf(PETSC_COMM_WORLD, "Global dimension of DA in i,j and k direction for level %d are %d %d %d :\n",l,user[bi].IM,user[bi].JM,user[bi].KM);
    ierr=DMDACreate3d(PETSC_COMM_WORLD,xperiod,yperiod,zperiod,DMDA_STENCIL_BOX,user[bi].IM,user[bi].JM,user[bi].KM,
		      m,n,p,5,s,lxSum,lySum,lzSum,&(user[bi].cda));
 
    if (ierr) SETERRQ1(PETSC_COMM_SELF,1, "problem creating DA %d",ierr);
    
    if (rans &&  l == usermg->mglevels-1) { // fda2 for k-omega at the finest level
      ierr=DMDACreate3d(PETSC_COMM_WORLD,xperiod,yperiod,zperiod,DMDA_STENCIL_BOX,
		      user[bi].IM+1, user[bi].JM+1, user[bi].KM+1,
		      m, n, p,
		      2, s, lxSum, lySum, lzSum,
		      &(user[bi].fda2));
      if (ierr) SETERRQ1(PETSC_COMM_SELF,1, "problem creating FDA2 for RANS%d",ierr);
    }
    user[bi].aotopetsc = PETSC_FALSE;
    DMDASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    DMGetCoordinateDM(user[bi].da, &(user[bi].fda));
    DMDAGetLocalInfo(user[bi].da, &(user[bi].info));
    
    PetscFree(lx);PetscFree(ly);PetscFree(lz);  
    PetscFree(lxSum);PetscFree(lySum);PetscFree(lzSum);  
  }
  return 0;
}

PetscErrorCode MG_Initial(UserMG *usermg, IBMNodes *ibm)
{
  MGCtx *mgctx;

  PetscInt level, block_number, blank=0;

  PetscInt rank;
  PetscInt IM, JM, KM;

  PetscInt i, j, k;
  PetscInt ih, jh, kh;
  PetscInt bi,sb;

  PetscErrorCode ierr; 

  Vec Coor, gCoor, Coor_high;
  Cmpnts ***coor, ***coor_high;
  PetscReal *gc;
  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscReal cl = 1., L_dim=1;
  PetscOptionsGetReal(NULL, PETSC_NULL, "-chact_leng", &cl, PETSC_NULL);
  PetscOptionsGetInt(NULL, PETSC_NULL, "-blk", &blank, PETSC_NULL);

  /* How many MG levels, the default is 1 */
  /* if (poisson)  */
  /*   usermg->mglevels = 3; */
  /* else */
  usermg->mglevels = 1;
  PetscOptionsGetInt(NULL, PETSC_NULL, "-mg_level", &usermg->mglevels, PETSC_NULL);
  
  PetscPrintf(PETSC_COMM_WORLD, "generate grid! mglevels %d \n", usermg->mglevels);
  usermg->ksc = PETSC_FALSE;
  usermg->jsc = PETSC_FALSE;
  usermg->isc = PETSC_FALSE;
  
/*   PetscOptionsGetTruth(NULL, PETSC_NULL, "-mg_k_semi", &usermg->ksc, PETSC_NULL); */
/*   PetscOptionsGetTruth(NULL, PETSC_NULL, "-mg_j_semi", &usermg->jsc, PETSC_NULL); */
/*   PetscOptionsGetTruth(NULL, PETSC_NULL, "-mg_i_semi", &usermg->isc, PETSC_NULL); */
  

  PetscMalloc(usermg->mglevels*sizeof(MGCtx), &(usermg->mgctx));
  mgctx = usermg->mgctx;
  FILE *fd, *fd1;
  
  /* Read in number of blocks and allocate memory for UserCtx */
  PetscInt generate_grid=0, grid1d=0;
  PetscOptionsGetInt(NULL, PETSC_NULL, "-grid", &generate_grid, PETSC_NULL);
  PetscOptionsGetInt(NULL, PETSC_NULL, "-grid1d", &grid1d, PETSC_NULL);
  if (generate_grid) {
    block_number=1;
    PetscOptionsGetInt(NULL, PETSC_NULL, "-block_number", &block_number, PETSC_NULL);
    PetscPrintf(PETSC_COMM_WORLD, "generate grid!\n");
  } else {
    if (!rank) {
    fd = fopen("grid.dat", "r");
    
    fscanf(fd, "%i\n", &block_number);
    MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
/*     fclose(fd); */
    }
    else {
      MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
  }
  PetscBool set; 
  PetscInt ksc[block_number], jsc[block_number], isc[block_number], nblk=block_number;
  PetscOptionsGetIntArray(NULL, PETSC_NULL, "-mg_k_semi", ksc, &nblk, &set);
  PetscOptionsGetIntArray(NULL, PETSC_NULL, "-mg_j_semi", jsc, &nblk, &set);
  PetscOptionsGetIntArray(NULL, PETSC_NULL, "-mg_i_semi", isc, &nblk, &set);
 
  PetscInt cgrid[block_number];
  
  PetscOptionsGetIntArray(NULL, PETSC_NULL, "-cgrid",cgrid, &nblk, &set);
 
  for (level=0; level<usermg->mglevels; level++) {
    PetscMalloc(block_number*sizeof(UserCtx), &mgctx[level].user);
    for (bi=0; bi<block_number; bi++) {
      mgctx[level].user[bi].ibm = ibm;
      mgctx[level].user[bi].isc = isc[bi];//&usermg->isc;
      mgctx[level].user[bi].jsc = jsc[bi];//&usermg->jsc;
      mgctx[level].user[bi].ksc = ksc[bi];//&usermg->ksc;
      mgctx[level].user[bi].cgrid = cgrid[bi];
      if (level==0) PetscPrintf(PETSC_COMM_WORLD, "SEMI_COARSENING for block %d in i=%d j=%d k=%d\n", bi, isc[bi], jsc[bi], ksc[bi]);
      if (level==0) PetscPrintf(PETSC_COMM_WORLD, "C-grid for block %d is %d\n", bi,cgrid[bi]);
    }
  }
  
  /* Read BC type and interface control file at the 
     finest grid level*/
  level = usermg->mglevels-1;
  UserCtx *user;

  user = mgctx[level].user;
  if (!rank) {
    
    if (block_number>1) {
      fd1 = fopen("interface.dat", "r"); if (!fd1) SETERRQ(PETSC_COMM_SELF,1, "Cannot open interface.dat to write!!");
      for (bi=0; bi<block_number; bi++) {
	if (blank && bi==0) {    
	  for (sb=1; sb<block_number; sb++) { 
	    fscanf(fd1, "%i %i %i %i %i %i %i\n", &user[bi].ip[sb],&user[bi].jp[sb],&user[bi].kp[sb],
		   &user[bi].dispx[sb],&user[bi].dispy[sb],&user[bi].dispz[sb],&user[bi].dispnn[sb]);
	    PetscPrintf(PETSC_COMM_WORLD,"BLANK ip,jp,kp,dispy,dispz,dispnn - %i,%i,%i,%i,%i,%i\n", user[bi].ip[sb], user[bi].jp[sb], 
			user[bi].kp[sb], user[bi].dispy[sb], user[bi].dispz[sb], user[bi].dispnn[sb]);
	  }
	  MPI_Bcast(&user[bi].ip, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].jp, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].kp, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispx, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispy, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispz, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispnn, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  
	}   
	fscanf(fd1, "%i\n", &(user[bi].itfcptsnumber));
	MPI_Bcast(&(user[bi].itfcptsnumber), 1, MPI_INT, 0, PETSC_COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD, "itfcp number # %d\n",user[bi].itfcptsnumber);

	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfcI));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfcJ));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfcK));
	
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostI));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostJ));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostK));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostB));
	
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscReal),
		    &(user[bi].itfchostx));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscReal),
		    &(user[bi].itfchosty));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscReal),
		    &(user[bi].itfchostz));
      
	for (i=0; i<user[bi].itfcptsnumber; i++) {
	  fscanf(fd1, "%i %i %i\n", &(user[bi].itfcI[i]), &(user[bi].itfcJ[i]),
		 &(user[bi].itfcK[i]));
	  fscanf(fd1, "%i %i %i %i\n", &(user[bi].itfchostI[i]),
		 &(user[bi].itfchostJ[i]), &(user[bi].itfchostK[i]),
		 &(user[bi].itfchostB[i]));
	  fscanf(fd1, "%le %le %le\n", &(user[bi].itfchostx[i]),
		 &(user[bi].itfchosty[i]), &(user[bi].itfchostz[i]));
	}

	MPI_Bcast(user[bi].itfcI, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfcJ, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfcK, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
      
	MPI_Bcast(user[bi].itfchostI, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostJ, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostK, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostB, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
      
	MPI_Bcast(user[bi].itfchostx, user[bi].itfcptsnumber, MPIU_REAL, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchosty, user[bi].itfcptsnumber, MPIU_REAL, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostz, user[bi].itfcptsnumber, MPIU_REAL, 0,
		  PETSC_COMM_WORLD);
      }
      fclose(fd1);
      
      /* Read the interface control file */
      fd1 = fopen("intcontrol.dat", "r"); if (!fd1) SETERRQ(PETSC_COMM_SELF,1, "Cannot open intcontrol.dat to write!!");
      for (bi=0; bi<block_number; bi++) {
	fscanf(fd1, "%i %i %i %i %i %i %i\n", &(user[bi].inttype[0]),
	       &(user[bi].inttype[1]), &(user[bi].inttype[2]),
	       &(user[bi].inttype[3]), &(user[bi].inttype[4]),
	       &(user[bi].inttype[5]), &(user[bi].inttype[6]));
	MPI_Bcast(&(user[bi].inttype[0]), 7, MPI_INT, 0, PETSC_COMM_WORLD);
      }
      fclose(fd1);
    } // if block_number>1
    
    /* Read in bcs.dat for boundary conditions at 6 boundary surfaces 
     First put the data onto the finest level and restrict to the coarser
     levels */
    
    fd1 = fopen("bcs.dat", "r"); if (!fd1) SETERRQ(PETSC_COMM_SELF,1, "Cannot open bcs.dat to write!!");

    for (bi=0; bi<block_number; bi++) {
      fscanf(fd1, "%i %i %i %i %i %i\n", &(user[bi].bctype[0]),
	     &(user[bi].bctype[1]), &(user[bi].bctype[2]),
	     &(user[bi].bctype[3]), &(user[bi].bctype[4]),
	     &(user[bi].bctype[5]));
      if (user[bi].bctype[0]==7 || user[bi].bctype[1]==7) i_periodic=1;
      if (user[bi].bctype[2]==7 || user[bi].bctype[3]==7) j_periodic=1;
      if (user[bi].bctype[4]==7 || user[bi].bctype[5]==7) k_periodic=1;

      PetscPrintf(PETSC_COMM_WORLD, "BC type x %d %d y %d %d z %d %d\n", user[bi].bctype[0], user[bi].bctype[1], user[bi].bctype[2], user[bi].bctype[3], user[bi].bctype[4], user[bi].bctype[5]);

      MPI_Bcast(&(user[bi].bctype[0]), 6, MPI_INT, 0, PETSC_COMM_WORLD);
     
    }
    fclose(fd1);

  }  else { // if !rank
    if (block_number>1) {  
      for (bi=0; bi<block_number; bi++) {
  	if (blank && bi==0) {     
	  MPI_Bcast(&user[bi].ip, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].jp, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].kp, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispx, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispy, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispz, 10, MPI_INT, 0, PETSC_COMM_WORLD);
	  MPI_Bcast(&user[bi].dispnn, 10, MPI_INT, 0, PETSC_COMM_WORLD);	 
	}    

	MPI_Bcast(&(user[bi].itfcptsnumber), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfcI));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfcJ));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfcK));
      
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostI));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostJ));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostK));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscInt),
		    &(user[bi].itfchostB));
      
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscReal),
		    &(user[bi].itfchostx));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscReal),
		    &(user[bi].itfchosty));
	PetscMalloc((user[bi].itfcptsnumber)*sizeof(PetscReal),
		    &(user[bi].itfchostz));

      
	MPI_Bcast(user[bi].itfcI, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfcJ, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfcK, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
      
	MPI_Bcast(user[bi].itfchostI, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostJ, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostK, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostB, user[bi].itfcptsnumber, MPI_INT, 0,
		  PETSC_COMM_WORLD);
      
	MPI_Bcast(user[bi].itfchostx, user[bi].itfcptsnumber, MPIU_REAL, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchosty, user[bi].itfcptsnumber, MPIU_REAL, 0,
		  PETSC_COMM_WORLD);
	MPI_Bcast(user[bi].itfchostz, user[bi].itfcptsnumber, MPIU_REAL, 0,
		  PETSC_COMM_WORLD);
      }

      for (bi=0; bi<block_number; bi++) {
	MPI_Bcast(&(user[bi].inttype[0]), 7, MPI_INT, 0, PETSC_COMM_WORLD);
      }
    }

    for (bi=0; bi<block_number; bi++) {
      MPI_Bcast(&(user[bi].bctype[0]), 6, MPI_INT, 0, PETSC_COMM_WORLD);
      if (user[bi].bctype[0]==7 || user[bi].bctype[1]==7) i_periodic=1;
      if (user[bi].bctype[2]==7 || user[bi].bctype[3]==7) j_periodic=1;
      if (user[bi].bctype[4]==7 || user[bi].bctype[5]==7) k_periodic=1;
    }
  }
  
  PetscPrintf(PETSC_COMM_WORLD, "i_periodic is %d  , j_periodic is %d and k_periodic is %d \n", i_periodic, j_periodic, k_periodic);
 
  /* The bcs.dat for boundary conditions at 6 boundary surfaces 
     on the finest level is restricted to the coarser
     levels */

  for (level = usermg->mglevels-2; level>=0; level--) {
      for (bi=0; bi<block_number; bi++) {
      mgctx[level].user[bi].bctype[0] = mgctx[level+1].user[bi].bctype[0];
      mgctx[level].user[bi].bctype[1] = mgctx[level+1].user[bi].bctype[1];
      mgctx[level].user[bi].bctype[2] = mgctx[level+1].user[bi].bctype[2];
      mgctx[level].user[bi].bctype[3] = mgctx[level+1].user[bi].bctype[3];
      mgctx[level].user[bi].bctype[4] = mgctx[level+1].user[bi].bctype[4];
      mgctx[level].user[bi].bctype[5] = mgctx[level+1].user[bi].bctype[5];
    }
  }

  /* Read from grid.dat the number of grid points along I, J, K directions
     and the detail coordinates of grid nodes */
  level = usermg->mglevels-1;
  user = mgctx[level].user;

/*   if (inletprofile==3) { */
/*   PetscPrintf(PETSC_COMM_WORLD, "READ INFLOW WAVE FORM!!!! mg_init\n"); */
/*   for (bi=0; bi<block_number; bi++) { */
/*     InflowWaveFormRead(&user[bi]); */
/*   } */
/*   PetscPrintf(PETSC_COMM_WORLD, "READ INFLOW WAVE FORM!!!! mg_init\n"); */
/*   } */

  PetscPrintf(PETSC_COMM_WORLD, "Rank %d \n", rank);

  PetscInt imm[block_number], kmm[block_number], jmm[block_number];
  if (generate_grid) {
    PetscOptionsGetIntArray(NULL, PETSC_NULL, "-im", imm, &nblk, &set);
    PetscOptionsGetIntArray(NULL, PETSC_NULL, "-jm", jmm, &nblk, &set);
    PetscOptionsGetIntArray(NULL, PETSC_NULL, "-km", kmm, &nblk, &set);
  }
  
  for (bi=0; bi<block_number; bi++) {

    if (!rank) {
      if (!generate_grid)
	fscanf(fd, "%i %i %i\n", &(user[bi].IM),
	     &(user[bi].JM), &(user[bi].KM));
      else {
	user[bi].IM=imm[bi];
	user[bi].JM=jmm[bi];
	user[bi].KM=kmm[bi];
      }	
      IM = user[bi].IM;
      JM = user[bi].JM;
      KM = user[bi].KM;
      
      MPI_Bcast(&(IM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(JM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(KM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(&(IM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(JM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(KM), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      user[bi].IM = IM;
      user[bi].JM = JM;
      user[bi].KM = KM;
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "READ GRID.dat %le I/j/KM %d %d %d\n", L_dim/cl, user[bi].IM,user[bi].JM,user[bi].KM);

    // Create DAs
    MGDACreate(usermg, bi);

    level = usermg->mglevels-1;
    user = mgctx[level].user;

    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    //    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    PetscPrintf(PETSC_COMM_WORLD, "READ info xe,ye,ze %d %d %d\n", xe,ye,ze);
    
    if (grid1d) PetscMalloc((IM+JM+KM)*sizeof(PetscReal), &gc);
    else        PetscMalloc(3*(IM*JM*KM)*sizeof(PetscReal), &gc);

    PetscReal L_x,L_y,L_z;

    if (generate_grid) {
       PetscOptionsGetReal(NULL, PETSC_NULL, "-L_x", &L_x, PETSC_NULL);
       PetscOptionsGetReal(NULL, PETSC_NULL, "-L_y", &L_y, PETSC_NULL);
       PetscOptionsGetReal(NULL, PETSC_NULL, "-L_z", &L_z, PETSC_NULL);
       PetscPrintf(PETSC_COMM_WORLD,"ksi eta zeta length  are %le %le %le \n",L_x,L_y,L_z);
    }
  
    // DMDAGetGhostedCoordinates(user[bi].da, &Coor);
    DMGetCoordinatesLocal(user[bi].da, &Coor);
    DMGetCoordinates(user[bi].da, &gCoor);
    VecSet(gCoor,0.0);
    
    DMDAVecGetArray(user[bi].fda, gCoor, &coor);
  
    if (!rank) {
      if (grid1d) {
	PetscReal xx;
	// read i
	for (i=0; i<IM; i++) 
	  fscanf(fd, "%le %le %le\n",&gc[i],&xx,&xx);
	// read j
	for (j=0; j<JM; j++) 
	  fscanf(fd, "%le %le %le\n",&xx,&gc[IM+j],&xx);
	// read k
	for (i=0; i<KM; i++) 
	  fscanf(fd, "%le %le %le\n",&xx,&xx,&gc[IM+JM+i]);
	
	MPI_Bcast(gc, (IM+JM+KM), MPIU_REAL, 0, PETSC_COMM_WORLD);
      
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    for (i=xs; i<xe; i++) {
	      if (k<KM && j<JM && i<IM) {
		coor[k][j][i].x = *(gc + i)/cl*L_dim;
		coor[k][j][i].y = *(gc + IM + j)/cl*L_dim;
		coor[k][j][i].z = *(gc + IM + JM + k)/cl*L_dim;
	      }
	    }
	  }
	}
      
      } else { // if 3d gridgen file
	
	for (k=0; k<KM; k++) {
	  for (j=0; j<JM; j++) {
	    for (i=0; i<IM; i++) {
	      if (!generate_grid)
		fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3);
	      else
		*(gc+(k*JM*IM + j*IM + i)*3) = L_x/(IM-1.) * i;	
	    }
	  }
	}
	
	for (k=0; k<KM; k++) {
	  for (j=0; j<JM; j++) {
	    for (i=0; i<IM; i++) {
	      if (!generate_grid)
		fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 1);
	      else
		*(gc+(k*JM*IM + j*IM + i)*3+1) = L_y/(JM-1.) * j;
	    }
	  }
	}
	
	for (k=0; k<KM; k++) {
	  for (j=0; j<JM; j++) {
	    for (i=0; i<IM; i++) {
	      if (!generate_grid)
		fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 2);
	      else
		*(gc+(k*JM*IM + j*IM + i)*3+2) = L_z/(KM-1.) * k ;
	    }
	  }
	}
	
	MPI_Bcast(gc, 3*(IM*JM*KM), MPIU_REAL, 0, PETSC_COMM_WORLD);
	
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    for (i=xs; i<xe; i++) {
	      if (k<KM && j<JM && i<IM) {
		coor[k][j][i].x = *(gc + (k * (IM*JM) + j * IM + i) * 3  )/cl*L_dim;
		coor[k][j][i].y = *(gc + (k * (IM*JM) + j * IM + i) * 3+1)/cl*L_dim;
		coor[k][j][i].z = *(gc + (k * (IM*JM) + j * IM + i) * 3+2)/cl*L_dim;
		//		PetscPrintf(PETSC_COMM_WORLD, "!!coor %le %le %le \n",coor[k][j][i].x, coor[k][j][i].y, coor[k][j][i].z);
	      }
	    }
	  }
	}
	
      } // grid-1d
    } // rank
    else {
      if (grid1d) {
	MPI_Bcast(gc, (IM+JM+KM), MPIU_REAL, 0, PETSC_COMM_WORLD);

	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    for (i=xs; i<xe; i++) {
	      if (k<KM && j<JM && i<IM) {
		coor[k][j][i].x = *(gc + i)/cl*L_dim;
		coor[k][j][i].y = *(gc + IM + j)/cl*L_dim;
		coor[k][j][i].z = *(gc + IM + JM + k)/cl*L_dim;
	      }
	    }
	  }
	}
	
      } else { // if 3d gridgen file
	
	MPI_Bcast(gc, 3*(IM*JM*KM), MPIU_REAL, 0, PETSC_COMM_WORLD);
	
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    for (i=xs; i<xe; i++) {
	      if (k<KM && j<JM && i<IM) {
		coor[k][j][i].x = *(gc + (k * (IM*JM) + j * IM + i) * 3  )/cl*L_dim;
		coor[k][j][i].y = *(gc + (k * (IM*JM) + j * IM + i) * 3+1)/cl*L_dim;
		coor[k][j][i].z = *(gc + (k * (IM*JM) + j * IM + i) * 3+2)/cl*L_dim;
	      }
	    }
	  }
	}
      } // grid-1d
    }
    PetscPrintf(PETSC_COMM_WORLD, "!! gc %le \n",gc[100]);
    PetscFree(gc);
    DMDAVecRestoreArray(user[bi].fda, gCoor, &coor);
    
    PetscReal nCoor;
    VecNorm(gCoor, NORM_INFINITY, &nCoor);
    PetscPrintf(PETSC_COMM_WORLD, "!!norm of Coor %le \n",nCoor);


    /* DMLocalToGlobalBegin(user[bi].fda, Coor, INSERT_VALUES, gCoor); */
    /* DMLocalToGlobalEnd(user[bi].fda, Coor, INSERT_VALUES, gCoor); */
    
    DMGlobalToLocalBegin(user[bi].fda, gCoor, INSERT_VALUES, Coor);
    DMGlobalToLocalEnd(user[bi].fda, gCoor, INSERT_VALUES, Coor);
    
    // VecDestroy(&Coor);
    // VecDestroy(&gCoor);
  } //bi

  UserCtx *user_high;
  for (level=usermg->mglevels-2; level>-1; level--) {
    user = mgctx[level].user;
    user_high = mgctx[level+1].user;
    for (bi = 0; bi<block_number; bi++) {
      
      DMDALocalInfo	info = user[bi].info;
      PetscInt	xs = info.xs, xe = info.xs + info.xm;
      PetscInt  ys = info.ys, ye = info.ys + info.ym;
      PetscInt	zs = info.zs, ze = info.zs + info.zm;
      PetscInt	mx = info.mx, my = info.my, mz = info.mz;
      
      // DMDAGetGhostedCoordinates(user_high[bi].da, &Coor_high);
      DMGetCoordinatesLocal(user_high[bi].da, &Coor_high);
      DMGetCoordinates(user[bi].da, &Coor);

      DMDAVecGetArray(user_high[bi].fda,Coor_high, &coor_high);
      DMDAVecGetArray(user[bi].fda, Coor, &coor);

      
      if (xe==mx) xe--;
      if (ye==my) ye--;
      if (ze==mz) ze--;
     
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {

	    GridRestriction(i, j, k, &ih, &jh, &kh, &user[bi]);
	    coor[k][j][i].x = coor_high[kh][jh][ih].x;
	    coor[k][j][i].y = coor_high[kh][jh][ih].y;
	    coor[k][j][i].z = coor_high[kh][jh][ih].z;
	  }
	}
      }
      
      DMDAVecRestoreArray(user[bi].fda, Coor, &coor);
      DMDAVecRestoreArray(user_high[bi].fda, Coor_high, &coor_high);

      // DMDAGetGhostedCoordinates(user[bi].da, &gCoor)
      DMGetCoordinatesLocal(user[bi].da, &gCoor);

      DMGlobalToLocalBegin(user[bi].fda, Coor, INSERT_VALUES, gCoor);
      DMGlobalToLocalEnd(user[bi].fda, Coor, INSERT_VALUES, gCoor);

      // VecDestroy(&Coor);
      // VecDestroy(&Coor_high);
      //  VecDestroy(&gCoor);
      
    }
  }

  if (!rank && !generate_grid) fclose(fd);


  /* Create the global and local vectors at all grid levels */

  for (level=usermg->mglevels-1; level>=0; level--) {
    user = mgctx[level].user;
    
    for (bi=0; bi<block_number; bi++) {

      user[bi].thislevel = level;
      user[bi]._this = bi;
      user[bi].block_number = block_number;
      user[bi].mglevels = usermg->mglevels;
      /* if (level > 0) { */
      /* 	user[bi].da_c = &mgctx[level-1].user[bi].da; */
      /* 	user[bi].lNvert_c = &mgctx[level-1].user[bi].lNvert; */
      /* 	user[bi].user_c = &mgctx[level-1].user[bi]; */
      /* } */
      /* if (level < usermg->mglevels-1) { */
      /* 	user[bi].da_f = &mgctx[level+1].user[bi].da; */
      /* 	user[bi].user_f = &mgctx[level+1].user[bi]; */
      /* } */

      ierr = DMCreateGlobalVector(user[bi].cda, &(user[bi].Q));
      VecDuplicate(user[bi].Q, &(user[bi].Qold));
      VecDuplicate(user[bi].Q, &(user[bi].Rhs));
      ierr = DMCreateLocalVector(user[bi].cda, &(user[bi].lQ));
      VecDuplicate(user[bi].lQ, &(user[bi].lQold));

      //Destroyed in metrics.c from here !!!1
      ierr = DMCreateGlobalVector(user[bi].fda, &(user[bi].Csi));
      ierr = VecDuplicate(user[bi].Csi, &(user[bi].Eta));
      ierr = VecDuplicate(user[bi].Csi, &(user[bi].Zet));

      /* VecDuplicate(user[bi].Csi, &(user[bi].ICsi)); */
      /* VecDuplicate(user[bi].Csi, &(user[bi].IEta)); */
      /* VecDuplicate(user[bi].Csi, &(user[bi].IZet)); */
      /* VecDuplicate(user[bi].Csi, &(user[bi].JCsi)); */
      /* VecDuplicate(user[bi].Csi, &(user[bi].JEta)); */
      /* VecDuplicate(user[bi].Csi, &(user[bi].JZet)); */
      /* VecDuplicate(user[bi].Csi, &(user[bi].KCsi)); */
      /* VecDuplicate(user[bi].Csi, &(user[bi].KEta)); */
      /* VecDuplicate(user[bi].Csi, &(user[bi].KZet)); */
      //Destroyed in metrics.c to here!!!

      /* VecDuplicate(user[bi].Csi, &(user[bi].Cent)); */
      
      /* VecDuplicate(user[bi].Csi, &(user[bi].GridSpace)); */

      if(level==usermg->mglevels-1) {		     				     
	/* VecDuplicate(user[bi].Csi, &(user[bi].Ucont)); */
	/* VecDuplicate(user[bi].Csi, &(user[bi].Ucont_o)); */
	/* VecDuplicate(user[bi].Csi, &(user[bi].Ucont_rm1)); */
	VecDuplicate(user[bi].Csi, &(user[bi].Ucat));
	/* VecDuplicate(user[bi].Csi, &(user[bi].Ucat_o)); */
	/* VecDuplicate(user[bi].Csi, &(user[bi].DUold)); */
	/* VecDuplicate(user[bi].Csi, &(user[bi].Bcs.Ubcs)); */
	/* VecDuplicate(user[bi].Csi, &(user[bi].psuedot)); */
	
	/* VecDuplicate(user[bi].Csi, &(user[bi].Vcont)); */
	/* if (rotateframe) */
	/*   VecDuplicate(user[bi].Csi, &(user[bi].Wcat)); */

	/* VecDuplicate(user[bi].Csi, &(user[bi].Itfc)); */
      }
      //      VecDuplicate(user[bi].Csi, &(user[bi].Rhs));
/*       if (level < usermg->mglevels-1) { */
/* 	VecDuplicate(user[bi].Csi, &(user[bi].Forcing)); */
/* 	VecDuplicate(user[bi].Csi, &(user[bi].Ucont_MG)); */
/*       } */
    
      ierr = DMCreateGlobalVector(user[bi].da, &(user[bi].Aj)); CHKERRQ(ierr); //Destroyed in metrics.c
      /* VecDuplicate(user[bi].Aj, &(user[bi].IAj));//Destroyed in metrics.c */
      /* VecDuplicate(user[bi].Aj, &(user[bi].JAj));//Destroyed in metrics.c */
      /* VecDuplicate(user[bi].Aj, &(user[bi].KAj));//Destroyed in metrics.c */

      VecDuplicate(user[bi].Aj, &(user[bi].Nvert));
/*       VecDuplicate(user[bi].Aj, &(user[bi].Phi)); */

      if(level==usermg->mglevels-1) {
      VecDuplicate(user[bi].Aj, &(user[bi].P));
/*       VecDuplicate(user[bi].Aj, &(user[bi].P_o)); */
      VecDuplicate(user[bi].Aj, &(user[bi].Nvert_o));
      }
    
      if(averaging && level==usermg->mglevels-1) {// New for Averaging
	VecDuplicate(user[bi].Csi, &user[bi].Ucat_sum);
	VecDuplicate(user[bi].Csi, &user[bi].Ucat_cross_sum);
	VecDuplicate(user[bi].Csi, &user[bi].Ucat_square_sum);
	VecDuplicate(user[bi].Q, &user[bi].Q_sum);

	VecSet(user[bi].Ucat_sum,0);
	VecSet(user[bi].Ucat_cross_sum,0);
	VecSet(user[bi].Ucat_square_sum,0);
	  
	VecDuplicate(user[bi].Aj, &user[bi].P_sum);
	VecSet(user[bi].P_sum,0);
      }
    
      //      VecDuplicate(user[bi].Aj, &(user[bi].Volume));

      DMCreateLocalVector(user[bi].fda, &(user[bi].lCsi));
      VecDuplicate(user[bi].lCsi, &(user[bi].lEta));
      VecDuplicate(user[bi].lCsi, &(user[bi].lZet));
    
      /* VecDuplicate(user[bi].lCsi, &(user[bi].lICsi)); */
      /* VecDuplicate(user[bi].lCsi, &(user[bi].lIEta)); */
      /* VecDuplicate(user[bi].lCsi, &(user[bi].lIZet)); */
      /* VecDuplicate(user[bi].lCsi, &(user[bi].lJCsi)); */
      /* VecDuplicate(user[bi].lCsi, &(user[bi].lJEta)); */
      /* VecDuplicate(user[bi].lCsi, &(user[bi].lJZet)); */
      /* VecDuplicate(user[bi].lCsi, &(user[bi].lKCsi)); */
      /* VecDuplicate(user[bi].lCsi, &(user[bi].lKEta)); */
      /* VecDuplicate(user[bi].lCsi, &(user[bi].lKZet)); */
      /* VecDuplicate(user[bi].lCsi, &(user[bi].lGridSpace)); */
      /* VecDuplicate(user[bi].lCsi, &(user[bi].lCent)); */
/*       VecDuplicate(user[bi].lCsi, &(user[bi].Centx)); */
/*       VecDuplicate(user[bi].lCsi, &(user[bi].Centy)); */
/*       VecDuplicate(user[bi].lCsi, &(user[bi].Centz));	 */
/*       if(level==usermg->mglevels-1) {			  */
/*       VecDuplicate(user[bi].lCsi, &(user[bi].lUcont)); */
      VecDuplicate(user[bi].lCsi, &(user[bi].lUcat));

/*       VecDuplicate(user[bi].lCsi, &(user[bi].lItfc)); */

/*       VecDuplicate(user[bi].lCsi, &(user[bi].lVcont)); */
/*       if (rotateframe) { */
/*       VecDuplicate(user[bi].lCsi, &(user[bi].lWcat)); */
/*       VecDuplicate(user[bi].lCsi, &(user[bi].lMAreaCsi)); */
/*       VecDuplicate(user[bi].lCsi, &(user[bi].lMAreaEta)); */
/*       VecDuplicate(user[bi].lCsi, &(user[bi].lMAreaZet)); */
/*       //   FormAreaMoment(&user[bi]); */
/*       } */
/*       VecDuplicate(user[bi].lCsi, &(user[bi].lUcont_o)); */
//* *       VecDuplicate(user[bi].lCsi, &(user[bi].lUcont_rm1)); *\/ */
/*       //      VecDuplicate(user[bi].lCsi, &(user[bi].lArea)); */
/*       } */
      DMCreateLocalVector(user[bi].da, &(user[bi].lAj));
      /* VecDuplicate(user[bi].lAj, &(user[bi].lIAj)); */
      /* VecDuplicate(user[bi].lAj, &(user[bi].lJAj)); */
      /* VecDuplicate(user[bi].lAj, &(user[bi].lKAj)); */
      //    VecDuplicate(user[bi].lAj, &(user[bi].lItfcP));
    
      VecDuplicate(user[bi].lAj, &(user[bi].lNvert));

      if (les) {
      VecDuplicate(user[bi].lCsi, &(user[bi].lBeta));
	
      VecDuplicate(user[bi].Csi, &(user[bi].Sigma));
      VecDuplicate(user[bi].lCsi, &(user[bi].lSigma));        

      VecDuplicate(user[bi].lAj, &(user[bi].lNu_t));		
      VecDuplicate(user[bi].Aj, &(user[bi].Nu_t));		
      VecDuplicate(user[bi].lAj, &(user[bi].lCS));		
      VecDuplicate(user[bi].Aj, &(user[bi].CS));		
      VecDuplicate(user[bi].lAj, &(user[bi].lBBeta));
      VecSet(user[bi].lNu_t, 0);			
      }
    /*   VecDuplicate(user[bi].lAj, &(user[bi].lPhi)); */
 
      if(level==usermg->mglevels-1) {
      VecDuplicate(user[bi].lAj, &(user[bi].lP));
      VecDuplicate(user[bi].lAj, &(user[bi].lNvert_o));
      /* VecDuplicate(user[bi].lAj, &(user[bi].lCS));		 */
      /* VecDuplicate(user[bi].Aj, &(user[bi].CS)); */			
      /* VecDuplicate(user[bi].lAj, &(user[bi].lBBeta));	 */		
	
      }
      //      VecDuplicate(user[bi].lAj, &(user[bi].lVolume));
    
      if(rans  && level==usermg->mglevels-1) {
	DMCreateLocalVector(user[bi].fda2, &user[bi].lK_Omega);	
	VecSet(user[bi].lK_Omega, 0);	// seokkoo
	DMCreateLocalVector(user[bi].fda2, &user[bi].lK_Omega_o);	
	VecSet(user[bi].lK_Omega_o, 0);	// seokkoo	
	VecDuplicate(user[bi].lP, &user[bi].Distance);
	VecSet(user[bi].Distance, 0);

	DMCreateGlobalVector(user[bi].fda2, &user[bi].K_Omega);	
	VecSet(user[bi].K_Omega, 0);	// seokkoo
	VecDuplicate(user[bi].K_Omega, &(user[bi].K_Omega_o));	
	VecSet(user[bi].K_Omega_o, 0);// seokkoo
	//VecDuplicate(user[bi].K_Omega, &(user[bi].K_Omega_rm1));
	//VecDuplicate(user[bi].lP, &(user[bi].lSrans));		VecSet(user[bi].lSrans, 0);// seokkoo
	VecDuplicate(user[bi].lP, &(user[bi].lNu_t));		
	VecSet(user[bi].lNu_t, 0);// seokkoo
        VecDuplicate(user[bi].P, &(user[bi].Nu_t));		
	VecSet(user[bi].Nu_t, 0);// seokkoo


	if(rans==3) {
	  VecDuplicate(user[bi].lP, &(user[bi].lF1));
	  VecSet(user[bi].lF1, 0);
	}
      }
      
      // New for LES
      /* if(les && level==usermg->mglevels-1) { */
      /* 	VecDuplicate(user[bi].lAj, &user[bi].lCs);	 */
      /* 	VecDuplicate(user[bi].lAj, &(user[bi].lNu_t));		 */
      /* 	VecDuplicate(user[bi].Aj, &user[bi].Cs);	 */
      /* 	VecDuplicate(user[bi].Aj, &(user[bi].Nu_t));	 */
      /* 	VecSet(user[bi].lNu_t, 0); */
      /* } */
      if ((wallfunction || rans) && level==usermg->mglevels-1) 
	VecDuplicate(user[bi].lP, &user[bi].lUstar);

      PetscPrintf(PETSC_COMM_WORLD,"global and local vectors at all level %d  are created\n",level);

      FormInitialize(&(user[bi]));

      PetscPrintf(PETSC_COMM_WORLD,"Initialization for level %d is done \n",level);

      // read Nvert for blanking in overset grids
      if (block_number>1 && bi==0 && blank && level==usermg->mglevels-1) {
	PetscViewer	viewer;
	char filen[80];
	sprintf(filen, "nvfield_blank.dat");
	
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
	VecLoad((user[bi].Nvert_o),viewer);
	PetscViewerDestroy(&viewer);
	
	DMGlobalToLocalBegin(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
	DMGlobalToLocalEnd(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
	
	VecCopy(user[bi].Nvert_o, user[bi].Nvert);
	VecCopy(user[bi].lNvert_o, user[bi].lNvert);    
      }
    
      FormMetrics(&(user[bi])); 
      
      //      if (user[bi].bctype[0]==7 || user[bi].bctype[2]==7 || user[bi].bctype[4]==7) Update_Metrics_PBC(&(user[bi])); 
    
      PetscPrintf(PETSC_COMM_WORLD,"Metrics formation for level %d is done \n",level);

      if (level==usermg->mglevels-1) MetricsDivergence(&(user[bi]));
    }
    
  }

  return(0);
}


PetscErrorCode MG_Finalize(UserMG *usermg)
{

  MGCtx *mgctx;

  PetscInt level, bi;

  UserCtx *user;  

  mgctx = usermg->mgctx;
  for (level=usermg->mglevels-1; level>=0; level--) {
    user=mgctx[level].user;
    for (bi=0; bi<user[0].block_number; bi++) {

      if(level==usermg->mglevels-1 && averaging) {
	VecDestroy(&user[bi].Ucat_square_sum);
	VecDestroy(&user[bi].Ucat_cross_sum);
	VecDestroy(&user[bi].Ucat_sum);
      }

      if ((wallfunction || rans) && level==usermg->mglevels-1) 
	VecDestroy(&user[bi].lUstar);

      /* VecDestroy(&user[bi].Cent); */
      /* VecDestroy(&user[bi].GridSpace); */
/*       VecDestroy(&user[bi].Centx); */
/*       VecDestroy(&user[bi].Centy); */
/*       VecDestroy(&user[bi].Centz); */
      
      VecDestroy(&user[bi].Q);
      VecDestroy(&user[bi].lQ);
      VecDestroy(&user[bi].Qold);
      VecDestroy(&user[bi].lQold);
      VecDestroy(&(user[bi].Rhs));
      
      // Csi, Eta, Zet, Aj destroyed in metrics.c
      /* VecDestroy(&user[bi].Csi); */
      /* VecDestroy(&user[bi].Eta); */
      /* VecDestroy(&user[bi].Zet); */
      /* VecDestroy(&user[bi].Aj); */

      if(level==usermg->mglevels-1) {
      /* VecDestroy(&user[bi].Ucont); */
      /* VecDestroy(&user[bi].Ucont_o); */
      /* VecDestroy(&user[bi].Ucont_rm1); */
      VecDestroy(&user[bi].Ucat);
      /* VecDestroy(&user[bi].Ucat_o); */
      /* VecDestroy(&user[bi].DUold); */
      /* VecDestroy(&user[bi].Bcs.Ubcs); */
      /* VecDestroy(&user[bi].psuedot); */

      /* VecDestroy(&user[bi].Vcont); */
      /* if (rotateframe) */
      /* VecDestroy(&user[bi].Wcat);  */

      /* VecDestroy(&user[bi].Itfc); */
      }
     
   /*    if (level < usermg->mglevels-1) { */
/* 	VecDestroy(&user[bi].Forcing); */
/* 	VecDestroy(&user[bi].Ucont_MG); */
/*       } */

      VecDestroy(&user[bi].Nvert);
      //      VecDestroy(&user[bi].Phi);

      if(level==usermg->mglevels-1) {
      VecDestroy(&user[bi].P);
      /* VecDestroy(&user[bi].P_o); */
      VecDestroy(&user[bi].Nvert_o);
      }      

      if(averaging && level==usermg->mglevels-1) {// New for Averaging
	VecDestroy(&user[bi].Ucat_sum);
	VecDestroy(&user[bi].Ucat_cross_sum);
	VecDestroy(&user[bi].Ucat_square_sum);
	VecDestroy(&user[bi].Q_sum);

	VecDestroy(&user[bi].P_sum);
      }

      VecDestroy(&user[bi].lCsi);
      VecDestroy(&user[bi].lEta);
      VecDestroy(&user[bi].lZet);

      /* VecDestroy(&user[bi].lICsi); */
      /* VecDestroy(&user[bi].lIEta); */
      /* VecDestroy(&user[bi].lIZet); */
      /* VecDestroy(&user[bi].lJCsi); */
      /* VecDestroy(&user[bi].lJEta); */
      /* VecDestroy(&user[bi].lJZet); */
      /* VecDestroy(&user[bi].lKCsi); */
      /* VecDestroy(&user[bi].lKEta); */
      /* VecDestroy(&user[bi].lKZet); */
      /* VecDestroy(&user[bi].lGridSpace); */
      /* VecDestroy(&user[bi].lCent); */

      if(level==usermg->mglevels-1) {
/*       VecDestroy(&user[bi].lUcont); */
      VecDestroy(&user[bi].lUcat);
/*       VecDestroy(&user[bi].lUcont_o); */
/*       VecDestroy(&user[bi].lUcont_rm1); */
      /* VecDestroy(&user[bi].Bcs.Ubcs); */
      /* VecDestroy(&user[bi].psuedot); */

      /* VecDestroy(&user[bi].lVcont); */
      /* if (rotateframe) */
      /* VecDestroy(&user[bi].lWcat); */

      /* VecDestroy(&user[bi].lItfc); */
      }
     
      VecDestroy(&user[bi].lAj);
      /* VecDestroy(&user[bi].lIAj); */
      /* VecDestroy(&user[bi].lJAj); */
      /* VecDestroy(&user[bi].lKAj); */
  
      VecDestroy(&user[bi].lNvert);
      
      if (les) {
      VecDestroy(&user[bi].Sigma);
      VecDestroy(&user[bi].lSigma);

      VecDestroy(&user[bi].lNu_t);
      VecDestroy(&user[bi].Nu_t);
      VecDestroy(&user[bi].lCS);
      VecDestroy(&user[bi].CS);
      VecDestroy(&user[bi].lBBeta);
      VecDestroy(&user[bi].lBeta);
      }
      /*  VecDestroy(&user[bi].lPhi); */

      if(level==usermg->mglevels-1) {
      VecDestroy(&user[bi].lP);
      VecDestroy(&user[bi].lNvert_o);
      }         

      if(rans  && level==usermg->mglevels-1) {
	VecDestroy(&user[bi].lK_Omega);	
	VecDestroy(&user[bi].lK_Omega_o);	
	VecDestroy(&user[bi].Distance);

	VecDestroy(&user[bi].K_Omega);	
	VecDestroy(&user[bi].K_Omega_o);	
	VecDestroy(&user[bi].lNu_t);		
	
	if(rans==3) {
	  VecDestroy(&user[bi].lF1);
	}
	//	DMDestroy(&user[bi].fda2);
      }

      if ((wallfunction || rans) && level==usermg->mglevels-1) 
	VecDestroy(&user[bi].lUstar);	

      //    DMDestroy(&user[bi].fda);
      DMDestroy(&user[bi].da);
      DMDestroy(&user[bi].cda);
    }//for bi
    
    /* PetscInt implicit; */
    /* PetscOptionsGetInt(NULL, PETSC_NULL, "-imp", &implicit, PETSC_NULL); */
    /* if(level==usermg->mglevels-1 && implicit==5) { */
/*       DMCompositeDestroy(usermg->mgctx[level].packer); */
/*       DMMGDestroy(usermg->mgctx[level].dmmg); */
/*     } */
    PetscFree(user);
  }//for level
  PetscFree(usermg->mgctx);
  //PetscFree(usermg);
  return 0;
}
