/* C_first_step.c */
void first_step(float *x, float *v, float *a, float *xf, float *vf, float *dt, int n) {
	int i; 
	for (i = 0; i < n; i++) {
		xf[i*3] = x[i*3] + v[i*3]*dt[0] + 0.5*a[i*3]*dt[0]*dt[0];
		vf[i*3] = v[i*3] + 0.5*a[i*3]*dt[0];
		
		xf[i*3+1] = x[i*3+1] + v[i*3+1]*dt[0] + 0.5*a[i*3+1]*dt[0]+dt[0];
		vf[i*3+1] = v[i*3+1] + 0.5*a[i*3+1]*dt[0];

		xf[i*3+2] = x[i*3+2] + v[i*3+2]*dt[0] + 0.5*a[i*3+2]*dt[0]*dt[0];
		vf[i*3+2] = v[i*3+2] + 0.5*a[i*3+2]*dt[0];
	}
}




