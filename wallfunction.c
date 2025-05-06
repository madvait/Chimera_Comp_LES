#include "variables.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* double integrate_F(double nu, double utau, double ya, double yb); */
/* double integrate_Fy(double nu, double utau, double ya, double yb); */
/* double sign(double a); */

double kappa=0.41, B=5.5;

void noslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, 
	     Cmpnts *Ub, double nx, double ny, double nz)
{
	double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
	
	(*Ub).x = sb/sc * u_c;
	(*Ub).y = sb/sc * v_c;
	(*Ub).z = sb/sc * w_c;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
}

void freeslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, 
	       Cmpnts *Ub, double nx, double ny, double nz)
{
  //  in the Normal direction interpolation 
  double Uan = Ua.x * nx + Ua.y * ny + Ua.z * nz;
  double Ucn = Uc.x * nx + Uc.y * ny + Uc.z * nz;
  double Ubn = Uan + (Ucn - Uan) * sb/sc;
  //   in the Tangential direction Ub_t = Uc_t
  double ut = Uc.x - Ucn * nx, vt = Uc.y - Ucn * ny, wt = Uc.z - Ucn * nz;
  
	(*Ub).x = ut;
	(*Ub).y = vt;
	(*Ub).z = wt;
	
	(*Ub).x += Ubn*nx;
	(*Ub).y += Ubn*ny;
	(*Ub).z += Ubn*nz;
}

/* ali added the Wall function for rough and smooth walls using Log-low   */ 
double E_coeff (double utau, double ks, double nu)
{
  double kplus=utau*ks/nu, dB;
  if(kplus<=2.25) dB = 0.0;
  else if(kplus>2.25 && kplus<90.) dB = (B-8.5+1./kappa*log(kplus))*sin(0.4258*(log(kplus)-0.811));
  else if(kplus>=90.) dB = B-8.5+1./kappa*log(kplus);
  return exp(kappa*(B-dB));
}

double u_hydset_roughness(double nu, double y, double utau, double ks)
{
  double y0plus=11.53, y1plus=300,f;
  double yplus = utau*y/nu;
 
  if(yplus<=y0plus){f= utau * yplus;}
  else if (yplus<=y1plus) {f= utau/kappa*log(E_coeff(utau,ks,nu)*yplus);}
  else {f=-1.;}
  return f;
}

double f_hydset(double nu, double u, double y, double utau0, double ks) 
{
double y0plus=11.53, f;
double yplus=utau0*y/nu;
 if (yplus<=y0plus) {f= utau0*yplus-u;}
 else {f= utau0*(1./kappa*log(E_coeff(utau0,ks,nu)*yplus))-u;}
return f;
}

double df_hydset (double nu, double u, double y, double utau0, double ks)
{
double eps=1.e-7;
return (f_hydset(nu, u, y, utau0 + eps, ks) - f_hydset(nu, u, y, utau0 - eps, ks))/(2.*eps);
}


double find_utau_hydset(double nu,double u, double y, double utau_guess, double ks)
{
 double utau,utau0 = utau_guess;  
 int ii;
 for(ii=0;ii<30;ii++){
 utau=utau0-f_hydset(nu,u,y,utau0,ks)/df_hydset(nu,u,y,utau0,ks);
 if (fabs(utau-utau0)<1.e-7)break;
 utau0=utau;
  }
 if (ii==30) PetscPrintf(PETSC_COMM_SELF, "iter utau diff %d %le %le  %le\n",ii,utau,utau0,fabs(utau-utau0));
 //printf("iter utau diff %d %le %le  %le\n",ii,utau,utau0,fabs(utau-utau0));
return utau;
};

//printf("iter  yplus  utau diff %d %le %le  %le\n",i,yplus,utau, fabs(u-uc));


double nu_t(double yplus)// in fact, nu_t/nu
{
	return kappa * yplus * pow ( 1. - exp( - yplus / 19. ) , 2.0 );
};

int pre_integrate_flag=0;
int n_yp=0;
int interval_yp=2;
int max_yp=500000;//5.e5
double *integration_buffer;

void pre_integrate ()
{
	if(pre_integrate_flag) return;
	else pre_integrate_flag=1;
	
	n_yp = ( max_yp / interval_yp ) ;
	
	//	double	integration_buffer  [ n_yp + 1 ];
	double	test_buffer  [ n_yp + 1 ];
	integration_buffer = &test_buffer;

	integration_buffer[0] = 0.;
	PetscInt i,k;
	for(i=1; i<=n_yp; i++) {
		int N=24;
		double ya_p = (double)(i-1) * interval_yp;
		double yb_p = (double)(i) * interval_yp;
		double ydiff = yb_p - ya_p, dy_p = ydiff / (double)N;
		double E[N+1];
		double val=0, ybegin_p=ya_p;
		
		for(k=0; k<=N; k++) E[k] =  1. / ( 1. + nu_t(ybegin_p + dy_p*k ) );
		
		for(k=0; k<N; k++) {
			double F[4];
			F[0] = E[k];
			F[1] = 1./ ( 1. + nu_t (ybegin_p+dy_p*1./3.) );
			F[2] = 1./ ( 1. + nu_t (ybegin_p+dy_p*2./3.) );
			F[3] = E [k+1];
			val += dy_p  / 3.* ( 3 * F[0] + 9 * F[1] + 9 * F[2] + 3 *F[3] ) / 8.;
			ybegin_p += dy_p;
		}
		integration_buffer[i] = integration_buffer[i-1] + val;
	}
};

double integrate_F1 (double nu, double utau, double yb)
{
	double val=0;
	
	pre_integrate ();
	
	double ya_plus = 0 * utau / nu;
	double yb_plus = yb * utau / nu;
	
	
	int ib = (int) ( yb_plus / (double) interval_yp );
	int N=4;
	
	if ( yb_plus <= (double) max_yp ) {
	  // PetscPrintf(PETSC_COMM_WORLD," ingrate buffer yb+ %f ib %d \n", yb_plus, ib);
		double int_a = 0;
		double int_b = ( integration_buffer[ib+1] - integration_buffer[ib] ) / (double) (interval_yp) * ( yb_plus - (double)(ib) * interval_yp ) + integration_buffer[ib];
		
	
		return ( int_b - int_a ) / utau;
	}
	else {
		val = integration_buffer [ n_yp ];
		ya_plus = max_yp;
	}
	
	
	double ydiff = yb_plus - ya_plus, dy = ydiff / (double)N;
	//	std::vector<double> E(N+1);
	double E[N+1];
	double ybegin=ya_plus;
	PetscInt i;
	for(i=0; i<=N; i++) {
		E[i] =  1./ ( 1. + nu_t(ybegin + dy*i ) );
	}
	
	for(i=0; i<N; i++) {
	        double F[4];
	        F[0] = E[i];
		F[1] = 1./ ( 1. + nu_t(ybegin+dy*1./3.) );
		F[2] = 1./ ( 1. + nu_t(ybegin+dy*2./3.) );
		F[3] = E[i+1];
		val += dy / 3.* ( 3 * F[0] + 9 * F[1] + 9 * F[2] + 3 *F[3] ) / 8.;
		ybegin += dy;
	}
	
	val /= utau;
	return val;
};

// ks: roughness size
double integrate_F(double nu, double utau, double yb, double ks)
{
  double ks_plus = utau * ks / nu;

	double yp_shift = 0.9*(sqrt(ks_plus)-(ks_plus)*exp(-ks_plus/6.));
	double y_shift = yp_shift * nu / utau;
	
	return integrate_F1 (nu, utau, yb + y_shift) - integrate_F1 (nu, utau, y_shift);
};


double u_Werner(double nu, double y, double utau)
{
	double yplus = utau * y / nu; 
	double A=8.3, B=1./7.;
	
	if(yplus<=11.81) return yplus*utau;
	else return A * pow(yplus, B) * utau;
};

double f_Werner(double nu, double u, double y, double utau)
{
	double ym=11.81*nu/utau;
	double A=8.3, B=1./7.;
	
	if( fabs(u) <= nu/(2*ym) * pow(A, 2./(1.-B) ) ) return utau*utau - u/y*nu;
	else return utau*utau -  u/fabs(u) * pow( 0.5*(1-B) * pow(A, (1+B)/(1-B)) * pow(nu/y, 1+B) + (1+B)/A * pow(nu/y, B) * fabs(u), 2/(1+B) );
}

double df_Werner(double nu, double u, double y, double utau)
{
	double eps=1.e-7;
	return ( f_Werner(nu, u, y, utau+eps) - f_Werner(nu, u, y, utau-eps) ) / ( 2*eps ) ;
}

double u_Cabot(double nu, double y, double utau, double dpdn)
{
	return utau * utau * integrate_F1( nu, utau, y );// + dpdn * integrate_Fy( nu, utau, 0, y );
};

double u_Cabot_roughness(double nu, double y, double utau, double dpdn, double ks_plus)
{
  return utau * utau * integrate_F( nu, utau, y, ks_plus );// + dpdn * integrate_Fy( nu, utau, 0, y );
};

double f_Cabot(double nu, double u, double y, double utau, double dpdn)
{
	return utau * utau * integrate_F1( nu, utau, y ) - u;
}

double f_Cabot_roughness(double nu, double u, double y, double utau, double dpdn, double ks)
{
	return utau * utau * integrate_F( nu, utau, y, ks ) - u;
}

double df_Cabot(double nu, double u, double y, double utau, double dpdn)
{
	double eps=1.e-7;
	return ( f_Cabot(nu, u, y, utau+eps, dpdn) - f_Cabot(nu, u, y, utau-eps, dpdn) ) / ( 2*eps ) ;
}

double df_Cabot_roughness(double nu, double u, double y, double utau, double dpdn, double ks)
{
	double eps=1.e-7;
	return ( f_Cabot_roughness(nu, u, y, utau+eps, dpdn, ks) - f_Cabot_roughness (nu, u, y, utau-eps, dpdn, ks) ) / ( 2*eps ) ;
}

/*
	integrate F dy 	= integrate [ 1/ (nut + nu) ] dy
				= integrate [ 1/ (nu) / {1 + k*yp( 1- exp(-yp/A) )^2 }  ] dy
				= integrate [ 1/ utau / {1 + k*yp( 1- exp(-yp/A) )^2 }  ] d(yp)
				= 1/ utau * integrate [  1 / {1 + k*yp( 1- exp(-yp/A)^2 ) }  ] d(yp)

	integrate Fy dy	= y/ utau * integrate [  1 / {1 + k*yp( 1- exp(-yp/A)^2 ) }  ] d(yp)
				= nu/utau^2 * integrate [  yp / {1 + k*yp( 1- exp(-yp/A)^2 ) }  ] d(yp)
*/


double find_utau_Cabot(double nu, double u, double y, double guess, double dpdn)
{
	double x, x0=guess;
	
	int i;
	
	for(i=0; i<30; i++) {
		x = x0 - f_Cabot(nu, u, y, x0, dpdn)/df_Cabot(nu, u, y, x0, dpdn);
		if( fabs(x0 - x) < 1.e-10) break;
		x0 = x;
	}
	
	if( fabs(x0 - x) > 1.e-5 && i>=29 ) printf("\n!!!!!!!! Iteration Failed !!!!!!!!\n");
		
	return x;
};

double find_utau_Cabot_roughness(double nu, double u, double y, double guess, double dpdn, double ks)
{
	double x, x0=guess;
	
	int i;
	
	for(i=0; i<30; i++) {
		x = x0 - f_Cabot_roughness(nu, u, y, x0, dpdn, ks)/df_Cabot_roughness(nu, u, y, x0, dpdn, ks);
		if( fabs(x0 - x) < 1.e-10) break;
		x0 = x;
	}
	
	if( fabs(x0 - x) > 1.e-5 && i>=29 ) printf("\n!!!!!!!! Iteration Failed !!!!!!!!\n");
		
	return x;
};

double find_utau_Werner(double nu, double u, double y, double guess)
{
	double x, x0=guess;
	
	int i;
	
	for(i=0; i<20; i++) {
		x = x0 - f_Werner(nu, u, y, x0)/df_Werner(nu, u, y, x0);
		if( fabs(x0 - x) < 1.e-7 ) break;
		x0 = x;
	}
	
	if( fabs(x0 - x) > 1.e-5 && i>=19 ) printf("\n!!!!!!!! Iteration Failed !!!!!!!!\n");
		
	return x;
};

double sign(double a)
{
	if(a>0) return 1;
	else if(a<0) return -1;
	else return 0;
}


double u_loglaw(double y, double utau, double roughness)
{
	return utau *  1./kappa * log( (roughness+y) / roughness ) ;
}

double find_utau_loglaw(double u, double y, double roughness)
{
	return kappa * u / log( (y+roughness) / roughness);
};

void wall_function (UserCtx *user, double sc, double sb, 
		    Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, 
		    double nx, double ny, double nz)
{
	double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
	*ustar = find_utau_Cabot(1./user->ren, ut_mag, sc, 0.01, 0);
	//double ut_mag_modeled = u_Werner(1./user->ren, sb, utau);
	double ut_mag_modeled = u_Cabot(1./user->ren, sb, *ustar, 0);
	
	if(ut_mag>1.e-10) {
		ut *= ut_mag_modeled/ut_mag;
		vt *= ut_mag_modeled/ut_mag;
		wt *= ut_mag_modeled/ut_mag;
	}
	else ut=vt=wt=0;
					
	// u = ut + (u.n)n
	(*Ub).x = ut + sb/sc * un * nx;
	(*Ub).y = vt + sb/sc * un * ny;
	(*Ub).z = wt + sb/sc * un * nz;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
}

void wall_function_loglaw (UserCtx *user, double ks, double sc, double sb, 
			   Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, 
			   double nx, double ny, double nz)
{
	double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	//	double ua_n= Ua.x*nx + Ua.y* ny + Ua.z* nz;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
	//*ustar = find_utau_Cabot_roughness(1./user->ren, ut_mag, sc, 0.01, 0, ks);
        *ustar = find_utau_hydset(1./user->ren, ut_mag, sc, 0.001, ks);
	//double ut_mag_modeled = u_Cabot_roughness(1./user->ren, sb, *ustar, 0, ks);
	double ut_mag_modeled = u_hydset_roughness(1./user->ren, sb, *ustar, ks);
	if (ut_mag_modeled<0.) {
	  //	  freeslip(user, sc,sb,Ua, Uc,Ub,nx,ny,nz);
	  ut_mag_modeled=ut_mag;
	  //*ustar =0.;
	  //  PetscPrintf(PETSC_COMM_SELF, "Y+>300!!!!!!!!!!!!!!! %le ",*ustar);
	}
	
	//PetscReal yp=sc* (*ustar)*user->ren;
	//if (yp>300) PetscPrintf(PETSC_COMM_SELF, "  Y+>300!!!!!!!!!!!!!!!  %le ",yp);
	
	if(ut_mag>1.e-10) {
		ut *= ut_mag_modeled/ut_mag;
		vt *= ut_mag_modeled/ut_mag;
		wt *= ut_mag_modeled/ut_mag;
	}
	else ut=vt=wt=0;
					
	// u = ut + (u.n)n
	(*Ub).x = ut + (sb/sc * un ) * nx;
	(*Ub).y = vt + (sb/sc * un ) * ny;
	(*Ub).z = wt + (sb/sc * un ) * nz;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
	
}



/* // old  */

/* double u_Werner(double nu, double y, double utau) */
/* { */
/* 	double yplus = utau * y / nu;  */
/* 	double A=8.3, B=1./7.; */
	
/* 	if(yplus<=11.81) return yplus*utau; */
/* 	else return A * pow(yplus, B) * utau; */
/* }; */

/* double f_Werner(double nu, double u, double y, double utau) */
/* { */
/* 	double ym=11.81*nu/utau; */
/* 	double A=8.3, B=1./7.; */
	
/* 	if( fabs(u) <= nu/(2*ym) * pow(A, 2./(1.-B) ) ) return utau*utau - u/y*nu; */
/* 	else return utau*utau -  u/fabs(u) * pow( 0.5*(1-B) * pow(A, (1+B)/(1-B)) * pow(nu/y, 1+B) + (1+B)/A * pow(nu/y, B) * fabs(u), 2/(1+B) ); */
/* } */

/* double df_Werner(double nu, double u, double y, double utau) */
/* { */
/* 	double eps=1.e-7; */
/* 	return ( f_Werner(nu, u, y, utau+eps) - f_Werner(nu, u, y, utau-eps) ) / ( 2*eps ) ; */
/* } */

/* double u_Cabot(double nu, double y, double utau, double dpdn) */
/* { */
/* 	//return utau * ( 1./kappa * log (utau*y/nu) + B ); */
/* 	return utau * utau * integrate_F( nu, utau, 0, y );// + dpdn * integrate_Fy( nu, utau, 0, y ); */
/* }; */

/* double u_Cabot2(double nu, double y, double utau_total, double utau_i, double dpdn) */
/* { */
/* 	return utau_i * utau_i * integrate_F( nu, utau_total, 0, y );// + dpdn * integrate_Fy( nu, utau_total, 0, y ); */
/* }; */

/* double f_Cabot(double nu, double u, double y, double utau, double dpdn) */
/* { */
/* 	//return utau * ( 1./kappa * log (utau*y/nu) + B ) - u; */
/* 	return utau * utau * integrate_F( nu, utau, 0, y ) - u;// + dpdn * integrate_Fy( nu, utau, 0, y ); */
/* } */

/* double df_Cabot(double nu, double u, double y, double utau, double dpdn) */
/* { */
/* 	double eps=1.e-7; */
/* 	return ( f_Cabot(nu, u, y, utau+eps, dpdn) - f_Cabot(nu, u, y, utau-eps, dpdn) ) / ( 2*eps ) ; */
/* } */

/* double nu_t(double nu, double utau, double y) */
/* { */
/* 	double yplus = utau * y / nu;  */
/* 	return nu * kappa * yplus * pow ( 1. - exp( - yplus / 19. ) , 2.0 ); */
/* }; */

/* /\* */
/* 	integrate F dy 	= integrate [ 1/ (nut + nu) ] dy */
/* 				= integrate [ 1/ (nu) / {1 + k*yp( 1- exp(-yp/A) )^2 }  ] dy */
/* 				= integrate [ 1/ utau / {1 + k*yp( 1- exp(-yp/A) )^2 }  ] d(yp) */
/* 				= 1/ utau * integrate [  1 / {1 + k*yp( 1- exp(-yp/A)^2 ) }  ] d(yp) */

/* 	integrate Fy dy	= y/ utau * integrate [  1 / {1 + k*yp( 1- exp(-yp/A)^2 ) }  ] d(yp) */
/* 				= nu/utau^2 * integrate [  yp / {1 + k*yp( 1- exp(-yp/A)^2 ) }  ] d(yp) */
/* *\/ */

/* double integrate_F(double nu, double utau, double ya, double yb) */
/* { */
/*   int N=20;//12; */
/* 	double ydiff = yb - ya, dy = ydiff / (double)N; */
/* 	double F[4]; */
/* 	double val=0, ybegin=ya; */
/* 	int i; */
	
/* 	for(i=0; i<N; i++) { */
/* 		F[0] = 1./ ( nu + nu_t(nu, utau, ybegin) ); */
/* 		F[1] = 1./ ( nu + nu_t(nu, utau, ybegin+dy*1/3.) ); */
/* 		F[2] = 1./ ( nu + nu_t(nu, utau, ybegin+dy*2./3.) ); */
/* 		F[3] = 1./ ( nu + nu_t(nu, utau, ybegin+dy) ); */
/* 		val += dy  / 3.* ( 3 * F[0] + 9 * F[1] + 9 * F[2] + 3 *F[3] ) / 8.; */
/* 		ybegin += dy; */
/* 	} */
/* 	return val; */
	
/* }; */
/* /\* */
/* double integrate_Fy(double nu, double utau, double ya, double yb) */
/* { */
	
/* 	if(!preintegrated) { */
/* 		preintegrated=1; */
/* 		pre_integrate(); */
/* 	} */
	
/* 	double yp = utau * yb / nu; */

/* 	int i; */
/* 	double K; */
/* 	i=1; */
/* 	do { */
/* 		if(X[i]>yp && X[i-1]<yp) break; */
/* 		i++; */
/* 	} while(1); */
/* 	K = (Fyp[i] - Fyp[i-1])/(X[i] - X[i-1])*(yp - X[i-1]) + Fyp[i-1]; */
	
/* 	//printf("%d, %e %e %e \n", i, yp, Fyp[i], Fyp[i-1]); */
/* 	//printf("%f\n", K * nu / (utau*utau)); */
/* 	return K * nu / (utau*utau); */
/* } */
/* *\/ */
/* double find_utau_Cabot(double nu,  double u, double y, double guess, double dpdn) */
/* { */
/* 	double x, x0=guess; */
	
/* 	int i; */
	
/* 	for(i=0; i<20; i++) { */
/* 		x = x0 - f_Cabot(nu, u, y, x0, dpdn)/df_Cabot(nu, u, y, x0, dpdn); */
/* 		if( fabs(x0 - x) < 1.e-7) break; */
/* 		x0 = x; */
/* 	} */
	
/* 	if( fabs(x0 - x) > 1.e-5 && i>=19 ) printf("\n!!!!!!!! Iteration Failed !!!!!!!!\n"); */
		
/* 	return x; */
/* }; */

/* double find_utau_Werner(double nu, double u, double y, double guess) */
/* { */
/* 	double x, x0=guess; */
	
/* 	int i; */
	
/* 	for(i=0; i<20; i++) { */
/* 		x = x0 - f_Werner(nu, u, y, x0)/df_Werner(nu, u, y, x0); */
/* 		if( fabs(x0 - x) < 1.e-7 ) break; */
/* 		x0 = x; */
/* 	} */
	
/* 	if( fabs(x0 - x) > 1.e-5 && i>=19 ) printf("\n!!!!!!!! Iteration Failed !!!!!!!!\n"); */
		
/* 	return x; */
/* }; */


/* double sign(double a) */
/* { */
/* 	if(a>0) return 1; */
/* 	else if(a<0) return -1; */
/* 	else return 0; */
/* } */

/* void wall_function (UserCtx *user, double sc, double sb,  */
/* 		    Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub,  */
/* 		    PetscReal *ustar, double nx, double ny, double nz) */
/* { */
/* 	double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z; */
/* 	double un = u_c * nx + v_c * ny + w_c * nz; */
/* 	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz; */
/* 	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt ); */
/* 	*ustar = (PetscReal) find_utau_Cabot(1./user->ren, ut_mag, sc, 0.01, 0); */
/* 	//double ut_mag_modeled = u_Werner(1./user->ren, sb, utau); */
/* 	double ut_mag_modeled = u_Cabot(1./user->ren, sb, *ustar, 0); */
	
/* 	if(ut_mag>1.e-10) { */
/* 		ut *= ut_mag_modeled/ut_mag; */
/* 		vt *= ut_mag_modeled/ut_mag; */
/* 		wt *= ut_mag_modeled/ut_mag; */
/* 	} */
/* 	else ut=vt=wt=0; */
					
/* 	// u = ut + (u.n)n */
/* 	(*Ub).x = ut + sb/sc * un * nx; */
/* 	(*Ub).y = vt + sb/sc * un * ny; */
/* 	(*Ub).z = wt + sb/sc * un * nz; */
	
/* 	(*Ub).x += Ua.x; */
/* 	(*Ub).y += Ua.y; */
/* 	(*Ub).z += Ua.z; */
/* } */
