/* VelVerlet.c */
#include <math.h>

void pair_force(float *ds, float *r, float *f1, float *f2, int rcut, float *c1, float *c2) {
	float r3[1], c1c2[1];
	float fg1, fg2, fg3;

	r3[0] = r[0]*r[0]*r[0];
	c1c2[0] = c1[0] * c2[0];

	if (r[0] < rcut) {       
	    fg1 = c1c2[0] * ds[0] / (r3[0]);
        fg2 = c1c2[0] * ds[1] / (r3[0]);
	    fg3 = c1c2[0] * ds[2] / (r3[0]);

	    f1[0] += fg1;
	    f1[1] += fg2;
       	f1[2] += fg3;
            
	    f2[0] -= fg1;
	    f2[1] -= fg2;
	    f2[2] -= fg3;
	}	
}

void pair_energ( float *r, float *E, float rcut, float *Eshift, float *c1, float *c2)  {
    if ( r[0] < rcut) {
		E[0] += c1[0] * c2[0] / r[0];
    }
    else { E[0] += Eshift[0]; }
}

void forces_periodic(float *x, float *v, float *charge, float *F, float *E, int n, float *rcut, float *box) {
    int i, j;
    float Eshift[1];
    float bc[3], ds[3], r[1];

    Eshift[0] = 1/rcut[0];

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

	    r[0] = sqrt((ds[0])*(ds[0]) + (ds[1])*(ds[1]) + (ds[2])*(ds[2]));
	
	    pair_force( ds, r, F+3*i, F+3*j, rcut[0], charge+i, charge+j);
   	    pair_energ( r, E, rcut[0], Eshift, charge+i, charge+j); 
   	    
	}
    }
}

void forces(float *x, float *v, float *charge, float *F, float *E, int n, float *rcut) {
    int i, j;
    float Eshift[1];
    float ds[3], r[1];

    Eshift[0] = 1/rcut[0];

    for (i = 0; i < n; i++) {
	for (j = i + 1; j < n; j++) {

		ds[0] = x[3*i]   - x[3*j]   ;
		ds[1] = x[3*i+1] - x[3*j+1] ;
		ds[2] = x[3*i+2] - x[3*j+2] ;

	    r[0] = sqrt((ds[0])*(ds[0]) + (ds[1])*(ds[1]) + (ds[2])*(ds[2]));
	
	    pair_force( ds, r, F+3*i, F+3*j, rcut[0], charge+i, charge+j);
   	    pair_energ( r, E, rcut[0], Eshift, charge+i, charge+j); 
   	    
	}
    }
}





