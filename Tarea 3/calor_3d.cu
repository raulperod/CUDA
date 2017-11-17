#include <stdio.h>
#include <stdlib.h>
#include <string.h>

__global__ void llenarMatriz(float* __restrict__ dev_u, int nx, int ny, int nz){ 
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

	if( x >= nx || y >= ny || z >= nz){ return; }
	int g = z + ( blockDim.z * gridDim.z ) * ( y + x * ( blockDim.y * gridDim.y ) );
    
	dev_u[g]=0.0;
	if( (x*x+y*y+z*z) < 30000 ){ dev_u[g]=3.0; }
}

__device__ void fronteras( float* __restrict__ dev_u, int i, int j, int k, int nx, int ny, int nz){
	// Bandas izquierda y derecha
	dev_u[k+nz*(1+i*ny)] = dev_u[k+nz*(3+i*ny)];
    dev_u[k+nz*((ny-1)+i*ny)] = dev_u[k+nz*((ny-2)+i*ny)];
	// Bandas superior e inferior
	dev_u[k+nz*(j+1*ny)] = dev_u[k+nx*(j+3*ny)];
    dev_u[k+nz*(j+(nx-1)*ny)] = dev_u[k+nz*(j+(nx-1)*ny)];
}

__global__ void euler( float* __restrict__ dev_u, float* __restrict__ dev_ut, int nx, int ny, int nz ){
	float gama = 0.001f;
	float dt = 0.001;
	float h = (float)(1.0f/nx); 
	float h2 = h * h;
	int i = threadIdx.x + blockIdx.x*blockDim.x; // 0 - 255
   	int j = threadIdx.y + blockIdx.y*blockDim.y; // 0 - 255
    int k = blockIdx.z * blockDim.z + threadIdx.z;

	fronteras( dev_u, i, j, k, nx, ny, nz)
	__syncthreads();
	/* Euler Scheme */
	if( i*(nx-i-1)*j*(ny-j-1)*k*(nz-k-1) != 0){
		dev_ut[k+nz*(j+i*ny)] =  dev_u[k+nz*(j+i*ny)] + ( dev_u[k+nz*(j+(i+1)*ny)] 
                        + dev_u[k+nz*(j+(i-1)*ny)] + dev_u[k+nz*((j-1)+i*ny)] + dev_u[k+nz*((j+1)+i*ny)]
                        + dev_u[(k+1)+nz*(j+i*ny)] + dev_u[(k-1)+nz*(j+i*ny)] - 6.*dev_u[k+nz*(j+i*ny)] ) * dt * gama / h2;
	}
}

void main(){
    int nx = 256 , ny = 256 , nz = 256;
    float *u;
   
    char nombre[20];
    FILE *fl;
   

    /* Main Iteration loop */
    for( ntime = 0; ntime < 500; ntime++){
        
        
        if(ntime%100 == 0){
            printf("%d iteraciones \n", ntime);
            sprintf(nombre,"calor_3d_%i.csv",ntime);
            puts(nombre);
        
            fl=fopen(nombre,"w");
            for(int i = 1; i < nx; i++ ){
                for(int j = 1; j < ny; j++ ){
                    for(int k = 1; k < nz; k++){
                        fprintf(fl, "%d , %d , %d , %f \n", i, j, k, u[k+nz*(j+i*ny)]); 
                    }
                }
            }
            close(fl);
        }  
    }
   
    return;
}
