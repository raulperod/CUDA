#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N 256
#define T 32

__global__ void calcularCalor2D( float *u, float *ut, int nx, int ny ){

	float gama, lap,h,h2,dt;
	/* Parameters to be used in the model */
	gama=0.001f;
	/* Numerical constants for the euler method */
	h=(float)(1.0f/nx); //0.025;
	h2=h*h;
	dt=0.001;
	//int i = threadIdx.x + blockIdx.x*blockDim.x; // 0 - 256
   	//int j = threadIdx.y + blockIdx.y*blockDim.y; // 0 - 255
	int i = threadIdx.x; // 0 - 256
	/* Main Iteration loop */
	// ---------------------------------
	/* Boundary Conditions: No-flux*/
	u[1+i*ny] = u[3+i*ny];
	u[(ny-1)+i*ny] = u[(ny-2)+i*ny];
	__syncthreads();
	// ---------------------------------
	u[i+ny] = u[i+3*ny];
	u[i+(nx-1)*nx] = u[i+(nx-2)*nx];
	__syncthreads();
	// ---------------------------------
	/* Euler Scheme */
	for(int j=0 ; j < nx ; j++){
		lap = u[j+(i+1)*ny] + u[j+(i-1)*ny] + u[(j-1)+i*ny] + u[(j+1)+i*ny] - 4.*u[j+i*ny];
		ut[j+i*ny] = u[j+i*ny] + lap * dt * gama / h2;
	}
	__syncthreads();
	// ---------------------------------
	/* Update of the mesh */
	for(int j=0 ; j < nx ; j++){
		u[j+i*ny]=ut[j+i*ny];
	}
	// ---------------------------------
}

void llenarMatriz(float *a, int n, int m);

int main(){

   	int i,j,ntime,nsteps,nx=8,ny=8;
	float *u, *ut;
	float *dev_u, *dev_ut;
	char nombre[20];
	FILE *fl;
	
	cudaMalloc( &dev_u, (nx+1) * (ny+1) * sizeof(float) );
	cudaMalloc( &dev_ut, nx * ny * sizeof(float) );

	u = (float*)malloc( (nx+1) * (ny+1) * sizeof(float) );
	ut = (float*)malloc( nx * ny * sizeof(float) );

	llenarMatriz(u, nx, ny);

	//dim3 numeroHilos(T, T);
	//dim3 numeroBloques(N/T, N/T);
	
	cudaMemcpy( dev_ut, ut, nx * ny * sizeof(float), cudaMemcpyHostToDevice );
    // ---------------------------------
	/* Main Iteration loop */
	for( ntime = 0; ntime < 5000 ; ntime++){

		cudaMemcpy( dev_u, u, (nx+1) * (ny+1) * sizeof(float), cudaMemcpyHostToDevice );

		//calcularCalor2D<<< numeroBloques, numeroHilos>>>( dev_u, dev_ut, nx, ny );
		calcularCalor2D<<< 1, nx >>>( dev_u, dev_ut, nx, ny );

		cudaMemcpy( u, dev_u, (nx+1) * (ny+1) * sizeof(float), cudaMemcpyDeviceToHost );
		
		if(ntime%500 == 0){
			
			printf("%d iteraciones \n", ntime);
			sprintf(nombre,"calor_2d_%i.csv",ntime);
			puts(nombre);
			fl=fopen(nombre,"w");

			for( i = 0; i < nx; i++ ){
				for( j = 0; j < ny; j++ ){
					fprintf(fl, "%d, %d, %f\n", i,j, u[j+i*nx]);
				}
			}
			fclose(fl);
		}   
	}
	//------------------------------
	
	printf("Ya terminamos!!!!!! \n");
	free( u ); free( ut );
	cudaFree( dev_u ); cudaFree( dev_ut );
	return 0;
}

void llenarMatriz(float *a, int n, int m){
   float *v1;
   v1 = a;

   for( int i=1; i < n; i++ ){
      for( int j=1; j < m; j++ ){
         if ( i*i+j*j > 20 && i*i+j*j < 50 ){
            v1[j+i*m] = 8.0f;
         }else{
            v1[j+i*m] = 0.0f;
         }
      }
   }
}