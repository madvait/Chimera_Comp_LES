#include "variables.h"

extern PetscInt  NumberOfBodies,blank, immersed, invicid;
//extern PetscInt  ti, tiout;
extern PetscInt  fish,averaging;
extern PetscInt  movefsi, rotatefsi, moveframe, STRONG_COUPLING;


PetscErrorCode Flow_Solver(UserMG *usermg,IBMNodes *ibm, 
			   FSInfo *fsi)
{

  UserCtx	*user;
  PetscInt	bi, ibi, level;
  //  PetscReal     tm_s, tm_e, tp_e, tpr_e;

  level = usermg->mglevels-1;
  user = usermg->mgctx[level].user;

/* ==================================================================================             */
/*   Turbulence Models! */
/* ==================================================================================             */
  
 if(les && !invicid && les!=10){
    for (bi=0; bi<user[0].block_number; bi++) {
      PetscPrintf(PETSC_COMM_WORLD, "LES\n");
      //Compute_Smagorinsky_Constant_1(&user[bi], user[bi].lQ);
      Compute_Smagorinsky_Constant_1(&user[bi]);
      Compute_eddy_viscosity_LES(&user[bi]);
    }
  }



/* ==================================================================================             */
/*   Compressible Solver! */
/* ==================================================================================             */
  RungeKutta(user, ibm, fsi);

    
/* ==================================================================================             */
/*        BC!        */
/* ==================================================================================             */


 /* ==================================================================================             */
    /*     OUTPUT Values! */
/* ==================================================================================             */
    
  for (bi=0; bi<user[0].block_number; bi++){ 
    if (averaging)
      Do_averaging(&(user[bi]));
    if (user[bi].ti == (user[bi].ti/user[bi].tiout) * user[bi].tiout)// || ti<10) 
      //  Q_P_VTK_Output(&(user[bi]));
      Q_P_Binary_Output(&(user[bi]));
  }
  
  return(0);

}

PetscErrorCode Struc_Solver(UserMG *usermg,IBMNodes *ibm, 
			    FSInfo *fsi, Cstart *cstart,
			    PetscInt itr_sc,
			    PetscInt tistart, 
			    PetscBool *DoSCLoop)
{
  PetscReal     dS_sc, dS_MIN=1e-5, dSmax,dS_rel=0.0;
  UserCtx	*user;
  PetscInt	i,j,bi,ibi, level, Add_dUndt=1,MHV_stuck=0,P_No ;
  Cmpnts	u_c,omega_c, a_c, rr;
  PetscReal     rx,ry,rz;

  level = usermg->mglevels-1;
  user = usermg->mgctx[level].user;

/* ==================================================================================             */
/*     Store old values to determine SC convergence */
/* ==================================================================================             */

  if (STRONG_COUPLING) {
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      for (i=0;i<6;i++){
	fsi[ibi].S_old[i] = fsi[ibi].S_new[i];
	fsi[ibi].S_ang_o[i]=fsi[ibi].S_ang_n[i];      
	if (itr_sc==1) {
	  fsi[ibi].dS[i]=0.;
	  fsi[ibi].atk=0.3;
	}
	fsi[ibi].dS_o[i]=fsi[ibi].dS[i];
	fsi[ibi].atk_o=fsi[ibi].atk;
      }
      if (itr_sc==2){
	fsi[ibi].atk_o=0.298;
      }
      
      fsi[ibi].F_x_old=fsi[ibi].F_x;
      fsi[ibi].F_y_old=fsi[ibi].F_y;
      fsi[ibi].F_z_old=fsi[ibi].F_z;
      
      fsi[ibi].M_x_old=fsi[ibi].M_x;
      fsi[ibi].M_y_old=fsi[ibi].M_y;
      fsi[ibi].M_z_old=fsi[ibi].M_z;
      
      fsi[ibi].L_o[0]=fsi[ibi].L_n[0];
      fsi[ibi].L_o[1]=fsi[ibi].L_n[1];
      fsi[ibi].L_o[2]=fsi[ibi].L_n[2];
      
    } // bi
  } // immersed

/* ==================================================================================             */
/*     Calculating Forces! */
/* ==================================================================================             */
  if (immersed) {   
    for (bi=0; bi<user[0].block_number; bi++) {      
      for (ibi=0;ibi<NumberOfBodies;ibi++) {

	//Calc_forces_SI(&user[bi],&fsi[ibi],&ibm[ibi], ti, ibi, bi);
	Calc_forces_SI2(&user[bi],&fsi[ibi],&ibm[ibi], user[bi].ti, ibi, bi);

      }
    } //bi
  } //immersed

/* ==================================================================================             */
/*     Find The new Position & Move the BODY */
/* ==================================================================================             */
  if (movefsi || rotatefsi || fish){// && ti>tistart+3) {
    for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--) {
      user = usermg->mgctx[level].user;
      for (bi=0; bi<user[0].block_number; bi++) {
	if (immersed) {
	  	  
	  for (ibi=0;ibi<NumberOfBodies;ibi++) {
	    /* calling rotation and tranlatin should be in this order!!!
	       1) call rotation to compute rotation with the same axis 
	       as the moments 
	       2) call translation and translate the axis to the new position */
	     //	    if (rotatefsi) calc_quarternion(&fsi[ibi],user[bi].dt,ibi);
	    if (movefsi) Calc_FSI_pos_intg(&fsi[ibi], &ibm[ibi], user[bi].dt) ;
	    //	  Calc_FSI_pos_SC(&fsi[ibi], &ibm[ibi], 0.5*(user->dt), user->dt) ;
	    //Forced_Motion(&fsi[ibi], 0.5,user->dt);
	    if (fish && ibi==0){
	      PetscPrintf(PETSC_COMM_WORLD, "fish swim! \n");
	      fish_swim(&ibm[ibi],(user[bi].ti)*user[bi].dt, user[bi].dt);
	    }
	  }	
	  
	  // CollisionDetectionOfCylinders(fsi,ibm,NumberOfBodies);

	  if (!moveframe) 
	    for (ibi=0;ibi<NumberOfBodies;ibi++) {
	      if (rotatefsi) {
	      	Elmt_Move_FSI_TRANSpROT_quarternion(&ibm[ibi],&fsi[ibi]);
	      	calc_ibm_normal(&ibm[ibi]);
	      } else if (movefsi)
		Elmt_Move_FSI_TRANS(&fsi[ibi], &ibm[ibi],ibi);
	    }
	      
	  VecSet(user[bi].Nvert,0.);
	  VecSet(user[bi].lNvert,0.);
	  for (ibi=0;ibi<NumberOfBodies;ibi++) {
	    if (!fish)	
	      Set_ibm_velocity_to_fsi(&fsi[ibi], &ibm[ibi],ibi, movefsi, rotatefsi);
	    
	    PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
	    ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
	    if (user[bi].ti == (user[bi].ti/user[bi].tiout)*user[bi].tiout)  ibm_surface_VTKOut(&ibm[ibi],ibi,user[bi].ti);
	  }
	  
	}
      }
    }  
  }
  
/* ==================================================================================             */
/*   Convergence of the SC Loop */
/* ==================================================================================             */
  *DoSCLoop = PETSC_FALSE;
  dSmax=1e-10;
  bi=0;
  for (ibi=0;ibi<NumberOfBodies;ibi++) {
    if (movefsi && STRONG_COUPLING) {
      for (i=0;i<6;i++) {
	dS_sc = fabs(fsi[ibi].S_new[i]-fsi[ibi].S_old[i]);
	if (dS_sc > dS_MIN) *DoSCLoop = PETSC_TRUE;
	if (dS_sc > dSmax) dSmax=dS_sc;
      }
    }
    
    PetscPrintf(PETSC_COMM_WORLD, "S-C Convergence %d %le %le %le SC %d\n", itr_sc, dSmax,fsi[0].S_new[5],fsi[0].S_old[5], *DoSCLoop);
    if (itr_sc>10) *DoSCLoop = PETSC_FALSE;
    
    if ((movefsi || rotatefsi) && !(*DoSCLoop)) {
      FSI_position_Output(&fsi[ibi], ibi, user[bi].ti);
      if (user[bi].ti == (user[bi].ti/user[bi].tiout)*user[bi].tiout) FSI_DATA_Output(&fsi[ibi], ibi, user[bi].ti);
    }
  }
/* ==================================================================================             */   
  return(0);
}
/* ==================================================================================             */



/* ==================================================================================             */
/*   Time advancement! */
/* ==================================================================================             */

PetscErrorCode RungeKutta(UserCtx *user, IBMNodes *ibm, 
			  FSInfo *fsi)
{
  
  PetscInt istage;//, istage0=0;
  PetscReal alfa[4], beta[4];
  PetscInt	bi, ibi; 
  PetscReal normF_bk[user[0].block_number], normQ[user[0].block_number];
  //alfa[0] = 0.25; alfa[1] = 1/3.; alfa[2] = 0.5; alfa[3] = 1.;
  // u(l)= beta u(n) + alpha (u(l)+dt* rhs(l))
  alfa[0]=1.; alfa[1] = 1/4.; alfa[2] = 2/3.;
  beta[0]=1.; beta[1] = 3/4.; beta[2] = 1/3.;
  
   // use first order if not WENO
  PetscInt  fluxorder=3; // type==1 Godunov, 2 MUSCL, 3 WENOe
  PetscOptionsGetInt(NULL, PETSC_NULL, "-fluxorder", &fluxorder, PETSC_NULL);
   if (fluxorder<3) {
     alfa[1]=1/2.; beta[1]=1/2.;
   }

    PetscInt pseudot;
    // pseudo time iteration
    for(pseudot=0; pseudot<1;pseudot++) {
      for (bi=0; bi<user[0].block_number; bi++) {
	for (istage=0; istage<fluxorder; istage++) {

	  formRHS(&(user[bi]), user[bi].Rhs, istage);

	  VecNorm(user[bi].Rhs, NORM_INFINITY, &normF_bk[bi]);
	  VecNorm(user[bi].Q, NORM_INFINITY, &normQ[bi]);
	  PetscPrintf(PETSC_COMM_WORLD, "!!norm of RHS %le Q %le %le\n", normF_bk[bi], normQ[bi], user[bi].dt);
       	  // Advanced in time using RK scheme
	  if (istage==0) 
	    VecWAXPY(user[bi].Q, alfa[istage] * user[bi].dt * user[bi].st, user[bi].Rhs, user[bi].Qold);
	  else {
	    /* Computes y = alpha x + beta y
	       VecAXPBY(Vec y,PetscScalar alpha,PetscScalar beta,Vec x) 
	       compute Q=alfa[istate]*Q+alfa[istage]*dt*RHS */
	    VecAXPBY(user[bi].Q,  alfa[istage] * user[bi].dt * user[bi].st, alfa[istage], user[bi].Rhs);
	    /* Computes y = alpha x + y.
	       VecAXPY(Vec y,PetscScalar alpha,Vec x) 
	       Compute Q=Q+ beta * Qn */
	    VecAXPY(user[bi].Q, beta[istage], user[bi].Qold);
	  }
	  //PetscPrintf(PETSC_COMM_WORLD, "!!FormBCS \n");
       	  FormBCS(&user[bi], alfa[istage]);
	  PetscPrintf(PETSC_COMM_WORLD, "!!FormBCS \n");
       	  
	  DMGlobalToLocalBegin(user[bi].cda, user[bi].Q, INSERT_VALUES, user[bi].lQ);
	  DMGlobalToLocalEnd(user[bi].cda, user[bi].Q, INSERT_VALUES, user[bi].lQ);
	  
	  EqOfState(&user[bi]);

	}//istage
	if (immersed) {
	  for (ibi=0;ibi<NumberOfBodies;ibi++) {
	    ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
	  }
	  // Call form BCs to fix the BCs after ibm intp - do not advance iso wall alfa=0
	  FormBCS(&user[bi], 0.);
	}

      }
    }

   return(0);
}



//!--------------------------------------------------------------------
PetscErrorCode SolveQuadratic(PetscReal II[4],PetscReal *x) 
//!--------------------------------------------------------------------  
{
  PetscReal delta;

  // delta = b^2 - 4ac
  delta=II[1]*II[1] - 4.*II[2]*II[0];
  if (delta<0.){
    PetscPrintf(PETSC_COMM_WORLD, "Imaginary Roots!\n");    
    return(1);
  } else {
    x[0]=(-II[1]+sqrt(delta))/(2.*II[2]);
    x[1]=(-II[1]-sqrt(delta))/(2.*II[2]);
    x[2]=0.;
  }
  return (0); 
}


//!--------------------------------------------------------------------
PetscErrorCode SolveCubic(PetscReal II[4],PetscReal *x) 
//!--------------------------------------------------------------------  
/*
Solution: To find roots to the following cubic equation where a, b, c, and d are real.

   a*x3 + b*x2 + c*x + d = 0

Formula:

  Step 1: Calculate p and q
          p = ( 3*c/a - (b/a)^2 ) / 3
          q = ( 2*(b/a)^3 - 9*b*c/a/a + 27*d/a ) / 27
  Step 2: Calculate discriminant D
          D = (p/3)^3 + (q/2)^2
  Step 3: Depending on the sign of D, you follow different strategy.
          If D<0, three distinct real roots.
          If D=0, three real roots of which at least two are equal.
          If D>0, one real and two complex roots.
  Step 3a: For D>0 and D=0
          Calculate u and v
          u = cubic_root(-q/2 + sqrt(D))
          v = cubic_root(-q/2 - sqrt(D))
          Find the three transformed roots
          y1 = u + v
          y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2
          y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
  Step 3b: Alternately, for D<0, a trigonometric formulation is more convenient
          y1 =  2 * sqrt(|p|/3) * cos(phi/3)
          y2 = -2 * sqrt(|p|/3) * cos((phi+pi)/3)
          y3 = -2 * sqrt(|p|/3) * cos((phi-pi)/3)
          where phi = acos(-q/2/sqrt(|p|^3/27))
                pi  = 3.141592654...
  Step 4  Finally, find the three roots
          x = y - b/a/3

  Things to watch out for:
    1. Make sure you know what is integer, real, and complex.
    2. FORTRAN's SQRT (square root) function takes only non-negative arguments.
    3. FORTRAN's exponentiation is "**", not "^".
    4. There is no "cubic_root" function in FORTRAN; you do it by raising a
       number to a factional power, e.g., 0.333... for cubic root.
    5. Your can only raise a non-negative number to a fractional power.
       e.g., (-2.)**(1./3.) will not work; you need to issue -(2.)**(1./3.)
       And do not forget the decimal point.  For example, in computing 2.**(1/3),
       FORTRAN evaluates 1/3 first, which is 0, not 0.33... as you might expect.
    6  There is no "|" in FORTRAN.  The absolute value function is "ABS".
    7. FORTRAN is case-insensitive; do not simply use "d" (for coefficient)
       and "D" (for discriminant), as FORTRAN thinks both variables are equivalent.
*/

{
  PetscReal p,q,D,u,y1,y2,y3;
  PetscReal phi, pi  = 3.141592654;
  
  p = ( 3*II[2]/II[0] - pow((II[1]/II[0]), 2) ) / 3.;
  q = ( 2*pow((II[1]/II[0]),3) - 9*II[1]*II[2]/II[0]/II[0] + 27*II[3]/II[0] ) / 27.;

  D = pow(p/3,3) + pow(q/2,2);

  if (D<0.) {
    phi = acos(-q/2/sqrt(pow(fabs(p),3)/27));
    y1 =  2 * sqrt(fabs(p)/3) * cos(phi/3);
    y2 = -2 * sqrt(fabs(p)/3) * cos((phi+pi)/3);
    y3 = -2 * sqrt(fabs(p)/3) * cos((phi-pi)/3);
  } else if (fabs(D)<1e-6) {

    u = -pow(fabs(q)/2,1./3.);
    if (q<0) u=-u;

   
    y1 = 2*u ;
    y2 = -u;
    y3 = -u;
  } else {
    
    y1=0;
    y2=0.;
    y3=0.;
    PetscPrintf(PETSC_COMM_WORLD, "Imaginary Roots!\n");    
    return(1);
  }

 

  x[0] = y1 - II[1]/II[0]/3;

  x[1] = y2 - II[1]/II[0]/3;

  x[2] = y3 - II[1]/II[0]/3;

 
  return(0);
}
