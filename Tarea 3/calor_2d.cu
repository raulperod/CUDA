#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#define N 256
#define T 32

__global__ void calcularCalor2D( float *x, float *y, int nx, int ny ) {
   int i = threadIdx.x + blockIdx.x*blockDim.x; // 0 - 255
   int j = threadIdx.y + blockIdx.y*blockDim.y; // 0 - 255

   // ----------------------------------------------------------------------------
   /* Main Iteration loop */
   for( ntime = 0; ntime < 100000; ntime++){
      // ------------------------------------
      /* Boundary Conditions: No-flux*/
      u[1+i*ny] = u[3+i*ny];
      u[(ny-1)+j*ny] = u[(ny-2)+j*ny];
      __syncthreads();
      // ------------------------------------
      u[i+ny] = u[i+3*ny];
      u[j+(nx-1)] = u[j+(nx-2)];
      __syncthreads();
      // ------------------------------------
      /* Euler Scheme */
      lap = u[j+(i+1)*ny] + u[j+(i-1)*ny] + u[(j-1)+i*ny] + u[(j+1)+i*ny] - 4.*u[j+i*ny];
      ut[j+i*ny] = u[j+i*ny] + lap * dt * gama / h2;
      __syncthreads();
      // ------------------------------------
      /* Update of the mesh */
      u[j+i*ny]=ut[j+i*ny];
      __syncthreads();
      // ------------------------------------
   }
   // ----------------------------------------------------------------------------
}

void llenarMatriz(float *a, int n, int m);

void main(){

   int i,j,ntime,nsteps,nx=256,ny=256;
   float *u, *ut;
   float gama,PI;
   float lap,h,h2,dt;
   /* Parameters to be used in the model */
   PI=3.1416f;
   gama=0.001f;
   /* Numerical constants for the euler method */
   h=(float)(1.0f/nx); //0.025;
   h2=h*h;
   dt=0.001;
   nsteps=3000000;

   cudaEvent_t start, end;
   cudaEventCreate(&start); 
   cudaEventCreate(&end);

   cudaMalloc( &dev_u, (nx+1) * (ny+1) * sizeof(float) );
   cudaMalloc( &dev_ut, nx * ny * sizeof(float) );

   u = (float*)malloc( (nx+1) * (ny+1) * sizeof(float) );
   ut = (float*)malloc( nx * ny * sizeof(float) );

   llenarMatriz(u, nx, ny);

   cudaMemcpy( dev_u, u, (nx+1) * (ny+1) * sizeof(float), cudaMemcpyHostToDevice );
   cudaMemcpy( dev_ut, ut, nx * ny * sizeof(float), cudaMemcpyHostToDevice );

   dim3 numeroHilos(T, T);
   dim3 numeroBloques(N/T, N/T);

   calcularCalor2D<<< numeroBloques, numeroHilos>>>( dev_u, dev_ut, nx, ny );

   printf("Ya terminamos!!!!!! \n");
   free( u ); free( ut );
   cudaFree( dev_u ); cudaFree( dev_ut );
   return;
}


void llenarMatriz(float *a, int n, int m){
   float *v1;
   v1 = a;

   for( i=1; i < n; i++ ){
      for( j=1; j < m; j++ ){
         if ( i*i+j*j > 1500 && i*i+j*j < 2500 ){
            v1[j+i*m] = 3.0f;
         }else{
            v1[j+i*m] = 0.0f;
         }
      }
   }
}