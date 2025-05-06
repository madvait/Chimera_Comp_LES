#include "variables.h"
extern PetscReal CMx_c,CMy_c,CMz_c;
//extern PetscInt ti,tistart;
extern PetscInt eel,pizza,moveframe;
extern PetscReal St_exp,wavelength;
PetscReal ampl,V_slip,omega, kwave,alpha;
PetscReal T_period_fish;
PetscReal c_0,c_1,c_2;
PetscInt  N_period_fish,forward=0;
//amir
PetscReal Y,Z,yy,zz,wn,z2,theta,absomg;
PetscInt  travel,airfoil;

PetscErrorCode fish_init(PetscReal *delti)
{
  PetscReal  pi = 3.141592653589793, St,Wl;
  //  PetscOptionsInsertFile(NULL, PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  PetscOptionsGetInt(NULL, PETSC_NULL, "-N_period_fish", &N_period_fish, PETSC_NULL);
  PetscPrintf(PETSC_COMM_WORLD, "N_period_fish is : %d \n",N_period_fish);
  //amir 
  theta=0.0;
  PetscOptionsGetReal(NULL, PETSC_NULL, "-theta", &theta, PETSC_NULL);
  PetscOptionsGetReal(NULL, PETSC_NULL, "-Croty", &yy, PETSC_NULL);
  PetscOptionsGetReal(NULL, PETSC_NULL, "-Crotz", &zz, PETSC_NULL);
  PetscOptionsGetReal(NULL, PETSC_NULL, "-Wavelength", &Wl, PETSC_NULL); 
  PetscOptionsGetReal(NULL, PETSC_NULL, "-Amplitude", &ampl, PETSC_NULL);
  travel=0;
  PetscOptionsGetInt(NULL, PETSC_NULL, "-travel", &travel, PETSC_NULL);
  PetscOptionsGetInt(NULL, PETSC_NULL, "-airfoil", &airfoil, PETSC_NULL);
  theta=-pi*theta/180.0;
  wn=2*pi/Wl;
  //ampl=.003;
  z2=zz-0.4;

  PetscInt  tistart;
  PetscBool rstart_flg;
  PetscOptionsGetInt(NULL, PETSC_NULL, "-rstart", &tistart, &rstart_flg);

  PetscPrintf(PETSC_COMM_WORLD, "travel %d tistart %d omega %le z2 %le  wavenumber %le Amplitude %le\n",travel,tistart,omega,z2,wn,ampl);
  //*amir
  V_slip=.7;

  //alpha = 2.18;
  alpha=pi/90.0;
  if (eel) 
    ampl=0.01;//0.089;
  else
  c_0=0.02;
  c_1=-0.08;
  c_2=0.16;
/*   c_0=0.02; */
/*   c_2=(0.3*(ampl-c_0)+c_0)/0.291;  */
/*   c_1=ampl - c_0 - c_2 ;//  ! Ampl. is defined at the end of tail */

  //St=St_exp/(2.*ampl);  // non-dim  frequency St=fL/U
  //omega=2*pi*St;//2.*pi; //non-dimensionalized w
  omega=2.*pi*St_exp;
  //kwave=2.*pi*St*V_slip; //non-dimensionalized k
  //wavelength = 0.95;
  if (wavelength>1.e6) kwave=0.0; 
  else  kwave= 2*pi/wavelength;//8.307767;
  V_slip= kwave*0.5/pi/St;

  T_period_fish=2.*pi/omega;
  *delti = T_period_fish/N_period_fish;  

  PetscPrintf(PETSC_COMM_WORLD, "fish init: St_exp %le Amplitude %le coeff a0 a1 a2 %le %le %le\n",St_exp,ampl,c_0,c_1,c_2);
  PetscPrintf(PETSC_COMM_WORLD, "fish init: dt %le w %le f %le k %le lambda %le T %le N %d  V_slip %le\n",*delti,omega,St,kwave,wavelength,T_period_fish,N_period_fish,V_slip);
 
  return(0);
}

PetscErrorCode fish_swim(IBMNodes *ibm, PetscInt ti
			,PetscReal delti)
{
  PetscReal time;
  time=ti*delti;
  PetscInt tistart_o = 0.;
  PetscPrintf(PETSC_COMM_WORLD, "fish swim");
  //  PetscOptionsInsertFile(NULL, PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  PetscOptionsGetInt(NULL, PETSC_NULL, "-pizza", &pizza, PETSC_NULL);
  PetscPrintf(PETSC_COMM_WORLD, "N_period_fish is : %d \n",N_period_fish);
  PetscOptionsGetInt(NULL, PETSC_NULL, "-forward", &forward, PETSC_NULL);
  PetscOptionsGetInt(NULL, PETSC_NULL, "-tistart_o", &tistart_o, PETSC_NULL);
 
  PetscReal  pi = 3.141592653589793;
  PetscReal  h0,h1,h2,h3,h4;

  PetscInt   i, n_v=ibm->n_v, n_elmt=ibm->n_elmt;
  PetscReal  x0,x1,y0,y1,z0,z1;
  PetscReal  az0,dz;

  y0=1.e+5;
  y1=-1.e+5;
  x0=y0;
  x1=y1;
  z0=y0;
  z1=y1;
  for (i=0; i<n_v; i++) {
    x0=PetscMin(x0,(ibm->x_bp0[i]));
    x1=PetscMax(x1,(ibm->x_bp0[i]));
    y0=PetscMin(y0,(ibm->y_bp0[i]));
    y1=PetscMax(y1,(ibm->y_bp0[i]));
    z0=PetscMin(z0,(ibm->z_bp0[i]));
    z1=PetscMax(z1,(ibm->z_bp0[i]));
  }

  PetscPrintf(PETSC_COMM_WORLD, "MAX fish: %d %le %le %le %le %le %le\n",ti,z1,z0,y1,y0,x1,x0);
  if(pizza){
    kwave=0.0;
    h0=exp(-1.0)*ampl*sin(-omega*time);
    h1=exp(-0.75)*ampl*sin(kwave*0.25-omega*time);
    h2=exp(-0.5)*ampl*sin(kwave*0.5-omega*time);
    h3=exp(-0.25)*ampl*sin(kwave*0.75-omega*time);
    h4=ampl*sin(kwave*1.0-omega*time);
    PetscPrintf(PETSC_COMM_WORLD, " t %le  h0 %le  h1 %le h2 %le h3 %le h4 %le\n",time,h0,h1,h2,h3,h4);
  }
  /*   oscillation for Mackerel in x-dir */
/*   for (i=0; i<n_v; i++) { */
/*     dz=(ibm->z_bp0[i]-z0);//+(13.3086/15.0);//for the tail simulation only */
/*     if(pizza){ */
/*       if(dz <=0.25)        ibm->x_bp[i]=ibm->x_bp0[i]+h1+(dz-0.25)*(h1-h0)/(0.25); */
/*       else  if(dz <=0.50)  ibm->x_bp[i]=ibm->x_bp0[i]+h2+(dz-0.5)*(h2-h1)/(0.25); */
/*       else  if(dz <=0.75)  ibm->x_bp[i]=ibm->x_bp0[i]+h3+(dz-0.75)*(h3-h2)/(0.25); */
/*       else  if(dz <=1.0)   ibm->x_bp[i]=ibm->x_bp0[i]+h4+(dz-1.0)*(h4-h3)/(0.25); */
/*       else                 ibm->x_bp[i]=ibm->x_bp0[i]+0.0; */
/*     }else{ */
/*       if (eel) // Anguiliform */
/* 	az0=ampl*exp(alpha*(dz-1.)); */
/*       else    // Carangiform: */
/* 	az0=c_0+c_1*dz+c_2*dz*dz; */

/*       ibm->x_bp[i]=ibm->x_bp0[i] + */

/* 	az0*sin(kwave*dz-omega*time); */
    
/*     } */
/*   } */
  PetscReal Maxaz=0.0;
  if (travel==1){
    PetscPrintf(PETSC_COMM_WORLD, "z2, travel %d: %le theta %le\n",z2,travel,omega); 
    for (i=0; i<n_v; i++) {	 
    
    
    //wave on y direction
    
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = ibm->z_bp0[i];
    ibm->x_bp[i] = ibm->x_bp0[i];
//     if ( ibm->z_bp0[i]<6.8){
      dz=ibm->z_bp0[i];
      
      az0=ampl;
    //   if (ti-tistart <100)
   //   az0=az0*(ti-tistart)/100;

   time=(ti-tistart_o)*delti; 
      absomg=sin(wn*dz-omega*time);
      ibm->y_bp[i]=ibm->y_bp0[i]+az0*absomg;
      ibm->u[i].y=-az0*omega*cos(wn*dz-omega*time);
     if (forward){ 
      absomg=sin(wn*dz+omega*time);
      ibm->y_bp[i]=ibm->y_bp0[i]+az0*absomg;
      ibm->u[i].y=az0*omega*cos(wn*dz+omega*time);
    }

      if (az0>Maxaz) Maxaz=az0;
      
//     }
    }
    //PetscPrintf(PETSC_COMM_WORLD, "omega %le, Maxaz %le time %le az0 %le\n",omega,Maxaz,time,az0);
    //  rotation
/*     for (i=0; i<n_v; i++) { */
/*       Y=ibm->y_bp[i]; */
/*       Z=ibm->z_bp[i]; */
/*       ibm->y_bp[i]=yy+(Z-zz)*sin(theta)+(Y-yy)*cos(theta); */
/*       ibm->z_bp[i]=zz+(Z-zz)*cos(theta)-(Y-yy)*sin(theta); */
/*       ibm->u[i].z=-ibm->u[i].y*sin(theta); */
/*       ibm->u[i].y=ibm->u[i].y*cos(theta); */
/*     } */
  }

  if (airfoil==1){
    for (i=0; i<n_v; i++) {
 
      Y=ibm->y_bp0[i];
      Z=ibm->z_bp0[i];
      theta=alpha*sin(omega*time);
      ibm->y_bp[i]=yy+(Z-z2)*sin(theta)+(Y-yy)*cos(theta);
      ibm->z_bp[i]=z2+(Z-z2)*cos(theta)-(Y-yy)*sin(theta);
    }
  
  }

    /*    Calculate the new normal & velcity */
    PetscInt n1e, n2e, n3e;
    PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  
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
    }
  
  PetscReal  v_max=0.;
  PetscInt   i_vmax=0;
  if (ti==0) {
    for (i=0; i<n_v; i++) {
      ibm->u[i].x = 0.0;
      ibm->u[i].y = 0.0;
      ibm->u[i].z = 0.0;
    } 

  } if (ti!=0) {
   for (i=0; i<n_v; i++) {
 /*         ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / delti;
          ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / delti;
         ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / delti;
     
   */ 
     if (v_max<fabs(ibm->u[i].z)) {
       i_vmax=i;
       v_max= fabs(ibm->u[i].z);
     }
     if (v_max<fabs(ibm->u[i].y)) {
       i_vmax=i;
       v_max= fabs(ibm->u[i].y);
     }
     if (v_max<fabs(ibm->u[i].x)) {
       i_vmax=i;
       v_max= fabs(ibm->u[i].x);
     }
   }
  }


  PetscPrintf(PETSC_COMM_WORLD, "MAX fish Velocity: %le %le %le %le %le %le %le %le\n",time, v_max, ibm->x_bp[i_vmax],ibm->y_bp[i_vmax],ibm->z_bp[i_vmax],ibm->u[i_vmax].x,ibm->u[i_vmax].y,ibm->u[i_vmax].z);
  
/*   PetscInt rank; */
/*   MPI_Comm_rank(PETSC_COMM_WORLD, &rank); */
/*   if (!rank) { */
/*     if (ti == (ti/tiout)*tiout) { */
/*       FILE *f; */
/*       char filen[80]; */
/*       sprintf(filen, "surface%3.3d.dat",ti); */
/*       f = fopen(filen, "w"); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z\n"); */
/*       PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", n_v, n_elmt); */
/*       for (i=0; i<n_v; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->x_bp[i]); */
/*       } */
/*       for (i=0; i<n_v; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->y_bp[i]); */
/*       } */
/*       for (i=0; i<n_v; i++) {	 */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->z_bp[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_x[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_y[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nf_z[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_x[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_y[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->nt_z[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_x[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_y[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e\n", ibm->ns_z[i]); */
/*       } */
/*       for (i=0; i<n_elmt; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1); */
/*       } */
/*       fclose(f); */

/* /\*       PetscFPrintf(PETSC_COMM_WORLD, f, "Variables=x,y,z\n"); *\/ */
/* /\*       PetscFPrintf(PETSC_COMM_WORLD, f, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n", n_v, n_elmt); *\/ */
/* /\*       for (i=0; i<n_v; i++) { *\/ */
/* /\* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%e %e %e\n", ibm->x_bp[i], ibm->y_bp[i], ibm->z_bp[i]); *\/ */
/* /\*       } *\/ */
/* /\*       for (i=0; i<n_elmt; i++) { *\/ */
/* /\* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%d %d %d\n", ibm->nv1[i]+1, ibm->nv2[i]+1, ibm->nv3[i]+1); *\/ */
/* /\*       } *\/ */
/* /\*       fclose(f); *\/ */

/*     } */
/*   } */

  return(0);
}


