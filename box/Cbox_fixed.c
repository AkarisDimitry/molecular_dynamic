/* Cbox_fixed.c */

void box_fixed_bouce(float *x, float *v, float *x0, float *xf, int n) {
	float lx, ly, lz;
	int i;
	lx = xf[0]-x0[0];
	ly = xf[0]-x0[0];
     lz = xf[0]-x0[0];

	for (i = 0; i < n; i++)	{
	     if (x[3*i] < x0[0] ) {
		     x[3*i] = 2*x0[0] - x[3*i]; 
 	          v[3*i] = -v[3*i] ;
          }
          if (x[3*i+1] < x0[1] ) {
               x[3*i+1] = 2*x0[1] - x[3*i+1];
               v[3*i+1] = -v[3*i+1] ;
          }
          if (x[3*i+2] < x0[2] ) {
               x[3*i+2] = 2*x0[2] - x[3*i+2];
               v[3*i+2] = -v[3*i+2] ;
          }

          if (x[3*i] > xf[0] ) {
               x[3*i] = 2*xf[0] - x[3*i];
               v[3*i] = -v[3*i] ;
          }
          if (x[3*i+1] > xf[1] ) {
               x[3*i+1] = 2*xf[1] - x[3*i+1];
               v[3*i+1] = -v[3*i+1] ;
          }
          if (x[3*i+2] > xf[2] ) {
               x[3*i+2] = 2*xf[2] - x[3*i+2];
               v[3*i+2] = -v[3*i+2] ;
          }
       }	
}



void box_fixed(float *x, float *v, float *x0, float *xf, int n) {
     int i;

     for (i = 0; i < n; i++)  {
          if (x[3*i] < x0[0]) {
               if (v[3*i] < 0) { v[3*i] = -v[3*i];}
          } else if(x[3*i] > xf[0]) {
               if (v[3*i] > 0) { v[3*i] = -v[3*i];}
          }

          if (x[3*i+1] < x0[1]) {
               if (v[3*i+1] < 0) { v[3*i+1] = -v[3*i+1];}
          } else if(x[3*i+1] > xf[1]) {
               if (v[3*i+1] > 0) { v[3*i+1] = -v[3*i+1];}
          }

          if (x[3*i+2] < x0[2]) {
               if (v[3*i+2] < 0) { v[3*i+2] = -v[3*i+2];}
          } else if(x[3*i+2] > xf[2]) {
               if (v[3*i+2] > 0) { v[3*i+2] = -v[3*i+2];}
          }

     }
}




