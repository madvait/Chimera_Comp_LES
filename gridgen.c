//#include "variables.h"
//#include "poisson.c"
static char help[] = "CURVIB compressible LES with PETSc/3.12.4-intel-2020a \n\n";

#include "petsctime.h" 

/* #undef __FUNCT__ */
/* #define __FUNCT__ "main" */

int main(int argc, char **argv)
{

  PetscInt i, j, k;  
  PetscReal L1=0.6, L=4.4, dx=0.01, r=1.1;
  PetscInt IM=2001, JM=3001, KM=5001, jh;

 

  PetscReal x[1000],x_all[IM],xx=-1;
  PetscReal y[2000],y_all[JM];
  PetscReal z[1000],z_all[KM];

  PetscInt grid_case=2;
  /*
    This code generates a uniform mesh (size dx) from -L1 to L1 and 
    then stretches it from L1 to L. 
    Most of the parameters need to be approximated by hand. 
    - use n=log(1-(L-L1)/dx(1-r))/log(r) to approximate the number of points
    n required for the stretched region L-L1. here it is about 40 based on 
    the above prameters
   */
  
  PetscInitialize(&argc, &argv, (char *)0, help);
  
 
  /* case == 1 lanl shock particle 
             2 LES cyl 2D
	     3 LES cyl 3D
   */

  //  PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE);          
  PetscOptionsGetInt(NULL, PETSC_NULL, "-g_case", &grid_case, PETSC_NULL);
  PetscOptionsGetInt(NULL, PETSC_NULL, "-im", &IM,PETSC_NULL);
  PetscOptionsGetInt(NULL, PETSC_NULL, "-jm", &JM,PETSC_NULL);
  PetscOptionsGetInt(NULL, PETSC_NULL, "-km", &KM,PETSC_NULL);

  if (grid_case==1) {
  //for i=1:62
  // uniform from 0 to L1
    IM=201; JM=201;
    
  for (i=0; i<63; i++) 
    x[i]=i*dx;
  
  //for i=63:100
  // stretch from L1 to L
  for (i=63; i<100; i++) 
    x[i]=x[i-1]+(x[i-1]-x[i-2])*r;
  
  x[100]=5.;

  j=100;
  //  for i=1:101
  // mirror the grid 
  for (i=0; i<101; i++) {
    x_all[i]=-x[j];
    j=j-1;
  }


  //for i=102:201
  for (i=101; i<201; i++)  x_all[i]=x[i-100];

  }  else if (grid_case == 2) {
    // dy=0.0025 from 0 to 1.5 ==> n=1.5/dx=600
    IM=3;
   
    //    KM=3201;
    dx=0.0025;
    r=1.1;

    for (i=0; i<IM; i++) x_all[i]=dx*4*i;

    /// z direction 
    // 0-1 uniform 0.01
    for (i=0; i<101; i++) z_all[i]=4*dx*i;
    // stretch to 0.0024
    // need to go from 4dx to dx with r=1.1 ==> 14 pts
    for (i=101; i<115; i++)
      z_all[i]=z_all[i-1]+(z_all[i-1]-z_all[i-2])/r;

    PetscPrintf(PETSC_COMM_WORLD, "!location z %f z_i-z_i-1 %f\n",z_all[114], z_all[114]-z_all[113] );
    jh=115+8/dx-1;
    // uniform with dx for L=8, n=8/dx=3200
    for (i=115; i<jh; i++) z_all[i]=dx*(i-114)+z_all[114];

    PetscPrintf(PETSC_COMM_WORLD, "!location 9 z %f\n",z_all[jh-1]);

    // stretch back to 4dx
    for (i=jh; i<(jh+14); i++)
      z_all[i]=z_all[i-1]+(z_all[i-1]-z_all[i-2])*r;

    PetscPrintf(PETSC_COMM_WORLD, "!location z %f z_i-z_i-1 %f\n",z_all[jh+13], z_all[jh+13]-z_all[jh+12] );

    for (i=(jh+14); i<KM; i++) z_all[i]=4*dx*(i-jh-13)+z_all[jh+13];

    PetscPrintf(PETSC_COMM_WORLD, "!location end z %f\n",z_all[KM-1]);

    ///// y direction
    for (i=0; i<801; i++) 
      y[i]=i*dx;

    PetscPrintf(PETSC_COMM_WORLD, "!location 1.5 y %f\n",y[800]);
    // need to go from dx to 4dx with r=1.1 ==> 14 pts
    for (i=801; i<815; i++)
      y[i]=y[i-1]+(y[i-1]-y[i-2])*r;

    PetscPrintf(PETSC_COMM_WORLD, "!location y %f y_i-y_i-1 %f\n",y[814], y[814]-y[813] );
    // the streched length is 0.07. The total lenght is L1=4.25-1.5-0.07=2.68
    // with 4dx=0.01 need 269 pts
    jh=1081;
    for (i=815; i<(jh+1); i++) 
       y[i]=(i-814)*4.*dx+y[814];

     PetscPrintf(PETSC_COMM_WORLD, "!location 4.25 y %f\n",y[jh]);

     JM=jh*2+1;

     j=jh;
     //  for i=1:101
     // mirror the grid 
     for (i=0; i<jh; i++) {
       y_all[i]=-y[j];
       j=j-1;
     }

     PetscPrintf(PETSC_COMM_WORLD, "!location JM %i\n",JM);
     for (i=jh; i<JM; i++)  y_all[i]=y[i-jh];

     PetscPrintf(PETSC_COMM_WORLD, "!location JM %i\n",JM);
  } else if (grid_case == 3) {
    // dy=0.005 from 0 to 1.5 ==> n=1.5/dx=300
    IM=101;
    PetscInt j1;
    //    KM=3201;
    dx=0.005;
    r=1.1;

    for (i=0; i<IM; i++) x_all[i]=dx*2*i;

    /// z direction 
    // 0-1 uniform 0.01
    for (i=0; i<101; i++) z_all[i]=2*dx*i;
    // stretch to 0.0024
    // need to go from 2dx to dx with r=1.1 ==> 7 pts
    j1=101+7;
    for (i=101; i<j1+1; i++)
      z_all[i]=z_all[i-1]+(z_all[i-1]-z_all[i-2])/r;

    PetscPrintf(PETSC_COMM_WORLD, "!location z %f z_i-z_i-1 %f\n",z_all[j1], z_all[j1]-z_all[j1-1] );
    jh=j1+8/dx-1;
    // uniform with dx for L=8, n=8/dx=1600
    for (i=j1; i<jh; i++) z_all[i]=dx*(i-j1)+z_all[j1];

    PetscPrintf(PETSC_COMM_WORLD, "!location 9 z %f\n",z_all[jh-1]);

    // stretch back to 2dx
    for (i=jh; i<(jh+7); i++)
      z_all[i]=z_all[i-1]+(z_all[i-1]-z_all[i-2])*r;

    PetscPrintf(PETSC_COMM_WORLD, "!location z %f z_i-z_i-1 %f\n",z_all[jh+6], z_all[jh+6]-z_all[jh+5] );

    for (i=(jh+7); i<KM; i++) 
      if (z_all[i-1]< 15.) 
	z_all[i]=2.*dx*(i-jh-6)+z_all[jh+6];
      else
	z_all[i]=z_all[i-1]+(z_all[i-1]-z_all[i-2])*r;

    
    PetscPrintf(PETSC_COMM_WORLD, "!location end z %f\n",z_all[KM-1]);

    ///// y direction
    for (i=0; i<401; i++) 
      y[i]=i*dx;

    PetscPrintf(PETSC_COMM_WORLD, "!location 1.5 y %f\n",y[400]);
    // need to go from dx to 4dx with r=1.1 ==> 7 pts
    for (i=401; i<408; i++)
      y[i]=y[i-1]+(y[i-1]-y[i-2])*r;

    PetscPrintf(PETSC_COMM_WORLD, "!location y %f y_i-y_i-1 %f\n",y[407], y[407]-y[406] );
    // the streched length is 0.07. The total lenght is L1=4.25-1.5-0.07=2.68
    // with 4dx=0.01 need 269 pts
    jh=674;
    for (i=408; i<(jh+1); i++) 
       y[i]=(i-407)*2.*dx+y[407];

     PetscPrintf(PETSC_COMM_WORLD, "!location 4.75 y %f\n",y[jh]);

     JM=jh*2+1;

     j=jh;
     //  for i=1:101
     // mirror the grid 
     for (i=0; i<jh; i++) {
       y_all[i]=-y[j];
       j=j-1;
     }

     PetscPrintf(PETSC_COMM_WORLD, "!location JM %i\n",JM);
     for (i=jh; i<JM; i++)  y_all[i]=y[i-jh];

     PetscPrintf(PETSC_COMM_WORLD, "!location KM %i IM %i\n",KM, IM);
  }



  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "grid_output.dat");
    f = fopen(filen, "w");
    // write i
    for (i=0; i<IM; i++) 
      PetscFPrintf(PETSC_COMM_WORLD,f, "%le %le %le\n",x_all[i],xx,xx);
    // read j
    for (j=0; j<JM; j++) 
      PetscFPrintf(PETSC_COMM_WORLD,f, "%le %le %le\n",xx,y_all[j],xx);
    // read k
    for (k=0; k<KM; k++)
      PetscFPrintf(PETSC_COMM_WORLD,f, "%le %le %le\n",xx,xx,z_all[k]);
    //      fscanf(fd, "%le %le %le\n",&xx,&xx,&gc[IM+JM+i]);
    fclose(f);
  }
  
  PetscFinalize();
  return(0);

}
