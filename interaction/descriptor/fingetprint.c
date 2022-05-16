#include <math.h>

#define PI 3.141592653
void Cutoff(float *r2, float *rcut, float *rcut2, float *fc) {

	if ( r2 > rcut2) { fc[0] = 0; }
	else { fc[0] = 0.5*(1+cos(PI*sqrt(r2[0])*rcut[0] )); }
}

/*
void G4() {
	float fc[3];
	Cutoff(r2, rcut, rcut2, fc);
	Cutoff(r2+1, rcut, rcut2, fc+1);
	Cutoff(r2+2, rcut, rcut2, fc+2);
	F[f_index[0]] += pow(1+cst[0], zetta[0]) * exp(-eta[0]*(r2[0]+r2[0]+r2[0])/rcut2[0])*fc[0]*fc[1]*fc[2]; 

}*/

void G2(float *r2, float *rcut, float *rcut2, float *eta, float *F, int *f_index) {
	float fc[1];
	Cutoff(r2, rcut, rcut2, fc);
	F[f_index[0]] += exp(-eta[0]*r2[0]/rcut2[0])*fc[0]; 
}	

void fingerprint(	float *x, int *t, int *t_n, float *F, int *atom_n, float *rcut, float *box, 
					float *G2_eta, float *nG2_eta, float *G2_n, 
					float *G4_n) {
	int i, j, k;
    float rcut2[1], bc[3], ds[3], r2[atom_n[0]*atom_n[0]];
    int f_index[1];

    rcut2[0] = rcut[0]*rcut[0];

    for (i = 0; i < atom_n[0]; i++) {
		for (j = i + 1; j < atom_n[0]; j++) {
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

		    r2[i+j*atom_n[0]] = (ds[0])*(ds[0]) + (ds[1])*(ds[1]) + (ds[2])*(ds[2]);

		    for (k = 0; k < nG2_eta[0]; k++) {
		    	f_index[0] = t[i]*(G2_n[0]+G4_n[0])+k*t_n[0]+j;
		    	G2(r2+i+j*atom_n[0], rcut, rcut2, G2_eta+k, F, f_index);
			}
			
			/*
			for (k = 0; k < nG2_eta[0]; k++) { # tercera especie
				for (k = 0; k < nG2_eta[0]; k++) { #gamma
					for (k = 0; k < nG2_eta[0]; k++) { #zetta
		    			f_index[0] = t[i]*(G2_n[0]+G4_n[0])+G2_n[0]+G4_zeta_n*nf* gama + n2*N-n2+n3;
		    			G4(r2+i+j*atom_n[0], rcut, rcut2, G2_eta+k, F, f_index);
			*/
		}

	}
}




