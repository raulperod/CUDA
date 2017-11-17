#include <stdio.h>

#define M 32  
#define N 32 

float a[M][N], b[M][N], c[M][N]; 

__global__ void kernelSumaMatrices(float *a, float *b,int m, int n) {
   int i = threadIdx.x + blockIdx.x*blockDim.x; 
   int j = threadIdx.y + blockIdx.y*blockDim.y; 
     
    while(i<m){
      j = threadIdx.y + blockIdx.y*blockDim.y;
      while(j<n){
         a[i*n+j]+=b[i*n+j];
         j+= blockDim.y*gridDim.y;
      }
      i+=blockDim.x*gridDim.x;
   } 
} 


void sumaMatrices(float *a, float *b, float *c, int m, int n) {
    float *ad, *bd, *cd; 
    int size = sizeof(float)*m*n; 
    dim3 nb(2,2); 
    dim3 nh(8,16); 

    cudaMalloc(&ad, size); 
    cudaMalloc(&bd, size); 
    cudaMalloc(&cd, size); 

    cudaMemcpy(ad, a, size, cudaMemcpyHostToDevice); 
    cudaMemcpy(bd, b, size, cudaMemcpyHostToDevice); 
    kernelSumaMatrices<<<nb , nh>>>(ad, bd, m, n); 
    cudaMemcpy(c, ad, size, cudaMemcpyDeviceToHost); 
    
    cudaFree(ad); cudaFree(bd); 
}

int main() {
    int i, j; 

    
    for (i=0; i<M; i++) {
       for (j=0; j<N; j++) {
          a[i][j]=b[i][j]=i+j; 
       }
    }
   
    sumaMatrices((float *)a, (float *)b, (float *)c,M, N); 

    for (i=0; i<M; i++) {
       for (j=0; j<N; j++) {
         printf("%4.2f  ", c[i][j]);
       }
       printf("\n"); 
    }
     
}


