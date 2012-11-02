/* 
 * File:   Helium_Solver.cpp
 * Author: henriasv
 * 
 * Created on 30. oktober 2012, 21:09
 */

#include "Helium_Solver.h"


Helium_Solver::Helium_Solver(int _N) {
    	alpha = 2;
	N = _N;
}

Helium_Solver::~Helium_Solver() {
}

double Helium_Solver::helium_function(double x1, double x2, double y1, double y2, double z1, double z2) 
{
	double alpha = 2;
	double r1r2 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
	double exp1 = sqrt(x1*x1 + y1*y1 + z1*z1);
	double exp2 = sqrt(x2*x2 + y2*y2 + z2*z2);
	if (r1r2>1e-6) {
	return exp(-2*alpha*(exp1+exp2))/r1r2;
	} else return 0;
}


#define   ZERO       1.0E-10
void Helium_Solver::gauleg(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359; 
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
	   ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
	 p2 =0.0;

   	   /*
	   ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

	 for(j = 1; j <= n; j++) {
	    p3 = p2;
	    p2 = p1;
	    p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
	 }

	   /*
	   ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */
 
	 pp = n * (z * p1 - p2)/(z * z - 1.0);
	 z1 = z;
	 z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /* 
	  ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
}
#undef ZERO
