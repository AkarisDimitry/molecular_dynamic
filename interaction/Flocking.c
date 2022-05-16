/* VelVerlet.c */
#include <math.h>
#include <stdio.h>

void force(float *ds, float *r, float *f1, float *f2, float *coef) {
	float f[3];

	f[0] = coef[0] + coef[1]*ds[0] + coef[2]*r[0];
	f[1] = coef[0] + coef[1]*ds[1] + coef[2]*r[1];
	f[2] = coef[0] + coef[1]*ds[2] + coef[2]*r[2];

	f1[0] += f[0];
	f1[1] += f[1];
	f1[2] += f[2];

	f2[0] -= f[0];
	f2[1] -= f[1];
	f2[2] -= f[2];
}

void pair_force(float *ds, float *r, float *dv, float *f1, float *f2, int rcut, float *coef1, float *coef2, float *coef3) {
	if (r[0] < rcut) {       
	    force(ds, r, f1, f2, coef1);
	    force(dv, r, f1, f2, coef2);
	    force(ds, r, f1, f2, coef3);
	}	
}

void pair_energ( float *r, float *E, float rcut, float *Eshift, float *c1, float *c2)  {
    if ( r[0] < rcut) {
		E[0] += c1[0] * c2[0] / r[0];
    }
    else { E[0] += Eshift[0]; }
}

void norm(float *mat, int n) {
	int i;
	float norm[0];

	for (i = 0; i < n; i++) {
		norm[0] = sqrt(mat[3*i]*mat[3*i] + mat[3*i+1]*mat[3*i+1] + mat[3*i+2]*mat[3*i+2]);
		if ( norm[0] > 0) { 
			mat[3*i]   = mat[3*i]/norm[0];
			mat[3*i+1] = mat[3*i+1]/norm[0];
			mat[3*i+2] = mat[3*i+2]/norm[0];
		}
	}
}

void forces(float *x, float *v, float *charge, float *F, float *E, int n, float *rcut, float *box,  float *coef1,  float *coef2,  float *coef3) {
    int i, j;
    float Eshift[1];
    float bc[3], ds[3], dv[3], r[1];

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

		dv[0] = v[3*i]   - v[3*j]   ;
		dv[1] = v[3*i+1] - v[3*j+1] ;
		dv[2] = v[3*i+2] - v[3*j+2] ;

	    r[0] = sqrt((ds[0])*(ds[0]) + (ds[1])*(ds[1]) + (ds[2])*(ds[2]));
		
	    pair_force( ds, r, dv, F+3*i, F+3*j, rcut[0], coef1, coef2, coef3);
   	    pair_energ( r, E, rcut[0], Eshift, charge+i, charge+j); 
   	    
	}
    }
    norm(F, n);
}





