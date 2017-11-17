// line to compile with dislin: nvcc calor2d.cu -o c2 -I/usr/local/dislin -ldislin

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "dislin.h"

// Globals

#define cNX 256
#define cNY 256
#define Nsteps0 35e3 //5e3
#define GAMA 0.001

const int cnthx= 32;
const int cnthy= 32;

const float dt=0.001;
const float h=(float)(1.0f/cNX); 
const float h2=h*h;

float *d_u;
float *d_ut;
float *ptr;
float *h_u;
float h2u[cNX][cNY];

#define SWAP(ptr,x,y) {ptr=&x[0]; x=&y[0]; y=ptr;}

dim3 TPB,NB;

// Device functions

__global__ void initField(float *d_u);
__global__ void Euler(float *d_u,float *d_ut);

__device__ void frontera(float *d_u,int i,int j);

// Function to print
void print_to_file(float* f);


// MAIN
int main(int argc, char *argv[]){

  int istep,k,l;
  float milliseconds;
  cudaEvent_t startEvent, stopEvent;

  // Allocate THREADS/BLK and BLKS/GRID
  int Nbx,Nby;
  TPB = dim3(cnthx,cnthy,1); 
  Nbx = (cNX+(cnthx-1))/cnthx;
  Nby = (cNY+(cnthy-1))/cnthy;
  NB  = dim3(Nbx,Nby,1);

  // Allocate device memory
  int memsize = cNX * cNY * sizeof(float);
  cudaMalloc( (void **)&d_u , memsize );
  cudaMalloc( (void **)&d_ut, memsize );

  //allocate host memory (we will use this to print into files)
  h_u = (float *)malloc( memsize );

      // Start timer
  cudaEventCreate(&startEvent);
  cudaEventCreate(&stopEvent);
  cudaEventRecord(startEvent, 0);

  // Initialize the field on the gpu
  initField<<<NB,TPB>>>(d_u);

  // Create initial condition and integrate up to Nsteps0 iterations
      for(istep=0;istep<Nsteps0;istep++){
        Euler<<<NB,TPB>>>(d_u,d_ut);
//        cudaMemcpy(d_u, d_ut, memsize ,cudaMemcpyDeviceToDevice);
        SWAP(ptr,d_u,d_ut);
        if (istep%250==0){
          printf("iteration= %d \n",istep);
          cudaMemcpy(h_u, d_u, memsize ,cudaMemcpyDeviceToHost);
          for ( k = 0; k < cNY; k++){
            for ( l = 0; l < cNX; l++ ){
              h2u[l][k] = h_u[l+cNX*k];
            }
          }
          //print_to_file(h_u);

// Dislin Plotting 

        metafl ("GL");
        disini ();
        pagera ();

        titlin("2-D Karma model",2);

        axspos (450, 1800);
        axslen (2200, 1200);

        name   ("X-axis", "x");
        name   ("Y-axis", "y");

	intax();
	autres(cNX,cNY);
	axspos(600,1850);
	ax3len(1500,1500,1500);
	
        graf3(0.0,cNX,0.0,100.0,1.0,cNY,0.0,100.0,
                 0.0,3.5,0.0,0.25);

	crvmat((float *) h2u,cNX,cNY,2,2);
	
	height(50);
	title();
	endgrf();
	erase();

        }
      }

    // Check runtime
    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&milliseconds, startEvent, stopEvent);
    
    cudaEventDestroy(startEvent);
    cudaEventDestroy(stopEvent);

    printf(" \n Grid Size: %d, %d \n",cNX,cNY);
    printf(" Total Iterations: %d", (int)Nsteps0);
    printf(" \nRuntime (s): %f\n\n", milliseconds/1000);

    cudaFree(d_u);
    cudaFree(d_ut);
    free(h_u);

    return 0;
}

__global__ void initField(float *d_u){
  
  const int xtid = blockIdx.x * blockDim.x + threadIdx.x;
  const int ytid = blockIdx.y * blockDim.y + threadIdx.y;

  if(xtid>=cNX || ytid>=cNY){return;}
  const int gtid = xtid+ytid*(blockDim.x*gridDim.x);
  
  d_u[gtid]=0.0;
  if(  (xtid*xtid+ytid*ytid)<15000) d_u[gtid]=3.0;
}


void print_to_file(float* f){
      int i;

      FILE *fp;
      fp=fopen("OUTPUT.dat","w");

  //    for ( j = 0; j < cNY; j++)
    //    if( j == cNY/2 ){ 
  	  for ( i = 0; i < cNX; i++)
	    fprintf(fp,"%d %d %f\n",i,int(cNY/2),f[i+(int(cNY/2))*cNX]);
      //  }
      fclose(fp);
}

__global__ void Euler(float *__restrict__ d_u,float *__restrict__ d_ut){
  
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

// Condición de frontera
    frontera(d_u,i,j);
    __syncthreads();

// Actualización
    if(i*(cNX-i-1)*j*(cNY-j-1)!=0){
        d_ut[i+cNY*j] = d_u[i+cNY*j] + (d_u[i-1+cNY*j]+d_u[i+1+cNY*j]+d_u[i+(j-1)*cNY]+d_u[i+(j+1)*cNY]
                        -4.0*d_u[i+cNY*j])*dt*GAMA/h2;
    }
}

__device__ void frontera(float *__restrict__ d_u,int i,int j){

   // Bandas izquierda y derecha
      d_u[1+cNY*j]=d_u[2+cNY*j];
      d_u[(cNX)+cNY*j]=d_u[(cNX-1)+cNY*j];

   // Bandas superior e inferior
      d_u[i+cNY]=d_u[i+cNY*2];
      d_u[i+cNY*(cNY-1)]=d_u[i+cNY*cNY];  

}
