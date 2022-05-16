#include <math.h>





void forces(float *x, float *v, float *F, float *E, int n, float *rcut, float *box) {
    int i, j;
    float Eshift, rcut2, rcut6;
    float bc[3], ds[3], r2[1];

    rcut2 = rcut[0]*rcut[0];
    rcut6 = rcut2*rcut2*rcut2;
    Eshift = 4/rcut6*rcut6 - 4/rcut6;

    for (i = 0; i < n; i++) {
	for (j = i + 1; j < n; j++) {
	    if ( (x[3*i]-x[3*j]) > rcut[0]) { bc[0]=box[0]; }	
	    else if ( (x[3*j]-x[3*i]) > rcut[0]) { bc[0]=-box[0]; }
	    else { bc[0]=0; }

	    if ( (x[3*i+1]-x[3*j+1]) > rcut[0]) { bc[1]=box[1]; }
	    else if ( (x[3*j+1]-x[3*i+1]) > rcut[0]) { bc[1]=-box[1]; }
	    else { bc[1]=0; }

	    if ( (x[3*i+2]-x[3*j+2]) > rcut[0]) { bc[2]=box[2]; }
	    else if ( (x[3*j+2]-x[3*i+2]) > rcut[0]) { bc[2]=-box[2]; }
	    else { bc[2]=0; }

		ds[0] = x[3*i]   - x[3*j]   - bc[0];
		ds[1] = x[3*i+1] - x[3*j+1] - bc[1];
		ds[2] = x[3*i+2] - x[3*j+2] - bc[2];

	    r2[0] = (ds[0])*(ds[0]) + (ds[1])*(ds[1]) + (ds[2])*(ds[2]);
	
	    pair_force( ds, r2, F+3*i, F+3*j, rcut[0]);
   	    pair_energ( r2, E, rcut[0], Eshift); 
   	    
	}
    }
}



