#include <stdio.h>
#include <stdlib.h>
#define N 2048

__global__ void suma( int *a, int *b, int *c ) {
    int index = blockIdx.x + blockIdx.y * blockDim.x;
    c[index] = a[index] + b[index]
}

void llenarMatriz(int*);

int main( void ) {
    int *a, *b, *c;
    int *dev_a, *dev_b, *dev_c; 
    int size = N * N * sizeof( int ); 
    
    cudaMalloc( &dev_a, size );
    cudaMalloc( &dev_b, size );
    cudaMalloc( &dev_c, size );
    a = (int*)malloc( size );
    b = (int*)malloc( size );
    c = (int*)malloc( size );
    
    llenarMatriz( a );
    llenarMatriz( b );
    
    cudaMemcpy( dev_a, a, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_b, b, size, cudaMemcpyHostToDevice );
    
    dim3 numeroHilos(1,1);
    dim3 numeroBloques(2048,2048);

    suma<<< numeroBloques, numeroHilos  >>>( dev_a, dev_b, dev_c );
    
    cudaMemcpy( c, dev_c, size , cudaMemcpyDeviceToHost );
    
    free( a ); free( b ); free( c );
    cudaFree( dev_a ); cudaFree( dev_b ); cudaFree( dev_c );
    return 0;
}

void llenarMatriz(int *m){
    for(int i=0 ; i < N ; i++){ 
        for(int j=0 ; j < N ; j++){
            m[j+i*N] = 1; 
        }
    }
}
