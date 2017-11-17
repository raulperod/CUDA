#include <stdio.h>
#include <stdlib.h>
#include <string.h>

__global__ void llenarMatriz(float* __restrict__ dev_u, int nx, int ny){ 
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;

	if( x >= nx || y >= ny){ return; }
	int g = x + y * ( blockDim.x * gridDim.x );

	dev_u[g]=0.0;
	if( (x*x+y*y) < 15000 ){ dev_u[g]=3.0; }
}

__device__ void fronteras( float* __restrict__ dev_u, int i, int j, int nx, int ny){
	// Bandas izquierda y derecha
	dev_u[1+ny*j]=d_u[2+ny*j];
	dev_u[(nx)+ny*j]=d_u[(nx-1)+ny*j];
	// Bandas superior e inferior
	dev_u[i+ny]=d_u[i+ny*2];
	dev_u[i+ny*(ny-1)]=d_u[i+ny*ny]; 
}

__global__ void euler( float* __restrict__ dev_u, float* __restrict__ dev_ut, int nx, int ny ){
	float gama = 0.001f;
	float dt = 0.001;
	float h = (float)(1.0f/nx); 
	float h2 = h * h;
	int i = threadIdx.x + blockIdx.x*blockDim.x; // 0 - 256
   	int j = threadIdx.y + blockIdx.y*blockDim.y; // 0 - 255

	fronteras( dev_u, i, j, nx, ny)
	__syncthreads();
	/* Euler Scheme */
	if( i*(nx-i-1)*j*(ny-j-1) != 0){
		dev_ut[j+i*ny] = dev_u[i+ny*j] + ( dev_u[i-1+ny*j] + dev_u[i+1+ny*j] + dev_u[i+(j-1)*ny] 
						+ dev_u[i+(j+1)*ny] - 4.0 * dev_u[i+ny*j] ) * dt * gama / h2;
	}
}

int main(){
   	int nx = 256, ny = 256;
	int tnx = 32, tny = 32;   
	int nbx = (nx+(tnx-1))/tnx;
	int nby = (ny+(tny-1))/tny;
	float *u;
	float *dev_u, *dev_ut;
	char nombre[20];
	FILE *fl;
	
	cudaMalloc( &dev_u, nx * ny * sizeof(float) );
	cudaMalloc( &dev_ut, nx * ny * sizeof(float) );

	u = (float*)malloc( nx * ny * sizeof(float) );

	dim3 NB = dim3( nbx, nby, 1);
	dim3 TB = dim3( tnx, tny, 1); 
	
	llenarMatriz<<< NB , TB>>>( dev_u, nx, ny);
	/* Main Iteration loop */
	for(int ntime = 0; ntime < 5000 ; ntime++){

		euler<<< NB , TB>>>( dev_u, dev_ut, nx, ny );
		
		cudaMemcpy(u, dev_ut, nx * ny * sizeof(float) , cudaMemcpyDeviceToHost);
		cudaMemcpy(dev_u, u, nx * ny * sizeof(float) , cudaMemcpyHostToDevice);

		if(ntime%500 == 0){
			
			printf("%d iteraciones \n", ntime);
			sprintf(nombre,"calor_2d_%i.csv",ntime);
			puts(nombre);
			fl=fopen(nombre,"w");

			for(int i = 0; i < nx; i++ ){
				for(int j = 0; j < ny; j++ ){
					fprintf(fl, "%d, %d, %f\n", i,j, u[j+i*nx]);
				}
			}
			fclose(fl);
		}   
	}
	//------------------------------
	printf("Ya terminamos!!!!!! \n");
	free( u );
	cudaFree( dev_u ); cudaFree( dev_ut );
	return 0;
}