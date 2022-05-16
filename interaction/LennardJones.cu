/* VelVerlet.c */
#include <math.h>

void pair_force(float *ds, float *r2, float *f1, float *f2, int rcut) {
	float r8, r14;
	float r5, r11;
	float fg1, fg2, fg3;

	if (r2[0] < rcut*rcut) { 
	    r8 = r2[0]*r2[0]*r2[0]*r2[0];
	    r14 = r8*r2[0]*r2[0]*r2[0];
	        
	    fg1 = (48/r14 - 24/r8)*(ds[0]);
        fg2 = (48/r14 - 24/r8)*(ds[1]);
	    fg3 = (48/r14 - 24/r8)*(ds[2]);

	    f1[0] += fg1;
	    f1[1] += fg2;
       	f1[2] += fg3;
            
	    f2[0] -= fg1;
	    f2[1] -= fg2;
	    f2[2] -= fg3;
	}	
}

void pair_energ( float *r2, float *E, float rcut, float Eshift)  {
    float r6, r12;
    if ( r2[0] < rcut*rcut) {
	    r6 = r2[0]*r2[0]*r2[0];
	    r12 = r6*r6;
		E[0] += 4/r12 - 4/r6;
    }
    else { E[0] += 0; }
}

__global__ void forces_cuda(float *x, float *F, float *E, int n, float *rcut, float *box, float *bc, float *bc) {
    int ytid = blockIdx.x * blockDim.x + threadIdx.x;
    int xtid = blockIdx.x * blockDim.x + threadIdx.x;

    if (ytid < n){
    if (xtid < n && xtid > ytid+1){
	    if ( (x[3*ytid]-x[3*xtid]) > rcut[0]) { bc[0]=box[0]; }
	    else if ( (x[3*xtid]-x[3*ytid]) > rcut[0]) { bc[0]=-box[0]; }
	    else { bc[0]=0; }

	    if ( (x[3*ytid+1]-x[3*xtid+1]) > rcut[0]) { bc[1]=box[1]; }
	    else if ( (x[3*xtid+1]-x[3*ytid+1]) > rcut[0]) { bc[1]=-box[1]; }
	    else { bc[1]=0; }

	    if ( (x[3*ytid+2]-x[3*xtid+2]) > rcut[0]) { bc[2]=box[2]; }
	    else if ( (x[3*xtid+2]-x[3*ytid+2]) > rcut[0]) { bc[2]=-box[2]; }
	    else { bc[2]=0; }

		ds[0] = x[3*ytid]   - x[3*xtid]   - bc[0];
		ds[1] = x[3*ytid+1] - x[3*xtid+1] - bc[1];
		ds[2] = x[3*ytid+2] - x[3*xtid+2] - bc[2];

	    r2[0] = (ds[0])*(ds[0]) + (ds[1])*(ds[1]) + (ds[2])*(ds[2]);
		
	    pair_force( ds, r2, F+3*ytid, F+3*xtid, rcut[0]);
   	    pair_energ( r2, E, rcut[0], Eshift); 
	}
    }

}

void forces(float *x, float *v, float *F, float *E, int n, float *rcut, float *box) {
    int i, j;
    float Eshift, rcut2, rcut6;
    float bc[3], ds[3], r2[1];
    float *d_x, *d_F, *d_E, *d_rcut, *d_box; 

    rcut2 = rcut[0]*rcut[0];
    rcut6 = rcut2*rcut2*rcut2;
    Eshift = 4/rcut6*rcut6 - 4/rcut6;

    // Allocate device memory 
    cudaMalloc((void**)&d_x, sizeof(float) * n*3);
    cudaMalloc((void**)&d_F, sizeof(float) * n*3);
    cudaMalloc((void**)&d_E, sizeof(float) * n);
    cudaMalloc((void**)&d_rcut, sizeof(float) * 1);
    cudaMalloc((void**)&d_box, sizeof(float) * 3*3);
    cudaMalloc((void**)&d_bc, sizeof(float) * 3*3);

    // Transfer data from host to device memory
    cudaMemcpy(d_x, x, sizeof(float) * n*3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_F, F, sizeof(float) * n*3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_E, a, sizeof(float) * n, cudaMemcpyHostToDevice);
    cudaMemcpy(d_rcut, rcut, sizeof(float) * 1, cudaMemcpyHostToDevice);
    cudaMemcpy(d_box, box, sizeof(float) * 3, cudaMemcpyHostToDevice);

    // Executing kernel 
    int block_size = 256;
    int grid_size = ((n*n + block_size) / block_size);
    forces_cuda<<<grid_size,block_size>>>(d_x, d_F, d_E, n, d_rcut, d_box, d_bc);
    
    // Transfer data back to host memory
    cudaMemcpy(out, d_out, sizeof(float) * N, cudaMemcpyDeviceToHost);

    // Deallocate device memory
    cudaFree(d_x);
    cudaFree(d_F);
    cudaFree(d_E);
    cudaFree(d_n);
    cudaFree(d_rcut);
    cudaFree(d_box);

}





