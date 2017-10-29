#include <stdio.h>
#include <stdlib.h>
#define N (2048*2048)
#define THREADS_PER_BLOCK 512

__global__ void dot( int *a, int *b, int *c ) {
    __shared__ int temp[THREADS_PER_BLOCK];
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    temp[threadIdx.x] = a[index] * b[index];
    __syncthreads();
    if( 0 == threadIdx.x ) {
        int sum = 0;
        for( int i = 0; i < THREADS_PER_BLOCK; i++ ){
            sum += temp[i];
        }
        atomicAdd( c , sum );
    }
}

void llenarVector(int*);

int main( void ) {
    int *a, *b, *c; // host copies of a, b, c
    int *dev_a, *dev_b, *dev_c; // device copies of a, b, c
    int size = N * sizeof( int ); // we need space for N ints
    // allocate device copies of a, b, c
    cudaMalloc( (void**)&dev_a, size );
    cudaMalloc( (void**)&dev_b, size );
    cudaMalloc( (void**)&dev_c, sizeof( int ) );
    a = (int *)malloc( size );
    b = (int *)malloc( size );
    c = (int *)malloc( sizeof( int ) );
    // inicializo los vectores a y b
    llenarVector( a );
    llenarVector( b );
    // mandas a y b al GPU
    cudaMemcpy( dev_a, a, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_b, b, size, cudaMemcpyHostToDevice );
    // realizas el calculo
    dot<<< (N+1)/THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>( dev_a, dev_b, dev_c );
    // obtienes el producto punto
    cudaMemcpy( c, dev_c, sizeof( int ) , cudaMemcpyDeviceToHost );
    printf("El producto punto es: %d", *c)
    free( a ); 
    free( b ); 
    free( c );
    cudaFree( dev_a );
    cudaFree( dev_b );
    cudaFree( dev_c );
    return 0;
}

void llenarVector(int *a){
    for(int i=0 ; i < N ; i++){ a[i] = 1; }
}
