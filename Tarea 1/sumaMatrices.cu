#include <stdio.h>
#include <stdlib.h>
#define N 2048
#define THREADS_PER_BLOCK 512

__global__ void suma( int *a, int *b, int *c ) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    c[index] = a[index] + b[index]
}

void llenarMatriz(int*);

int main( void ) {
    int *a, *b, *c;
    int *dev_a, *dev_b, *dev_c; 
    int size = N * N * sizeof( int ); 
    
    cudaMalloc( (void**)&dev_a, size );
    cudaMalloc( (void**)&dev_b, size );
    cudaMalloc( (void**)&dev_c, size );
    a = (int*)malloc( size );
    b = (int*)malloc( size );
    c = (int*)malloc( size );
    
    llenarMatriz( a );
    llenarMatriz( b );
    
    cudaMemcpy( dev_a, a, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_b, b, size, cudaMemcpyHostToDevice );
    
    suma<<< (N*N+1)/THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>( dev_a, dev_b, dev_c );
    
    cudaMemcpy( c, dev_c, size , cudaMemcpyDeviceToHost );
    
    free( a ); 
    free( b ); 
    free( c );
    cudaFree( dev_a );
    cudaFree( dev_b );
    cudaFree( dev_c );
    return 0;
}

void llenarMatriz(int *m){
    for(int i=0 ; i < N ; i++){ 
        for(int j=0 ; j < N ; j++){
            m[j+i*N] = 1; 
        }
    }
}
