#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// curand
#define N_TOTAL 1048576
#define N 1024
#define T 512

__global__ void aproximarPi( int *x, int *y, int *z ) {
    int i = threadIdx.x + blockIdx.x*blockDim.x; // 0 - 1023
    int j = threadIdx.y + blockIdx.y*blockDim.y; // 0 - 1023
    int index = j + i*N; // 0 - 1048576     
    
    if( (x[index] * x[index] + y[index] * y[index]) <= 1.0f){
        atomicAdd( z , 1 );
    }

}

void llenarRandom(float*);

int main( void ) {
    float *x, *y, *z;
    float *dev_x, *dev_y, *dev_z; 
    int size = N_TOTAL * sizeof( float );
    
    cudaEvent_t start, end;
    cudaEventCreate(&start); 
    cudaEventCreate(&end);
    
    cudaMalloc( &dev_x, size );
    cudaMalloc( &dev_y, size );
    cudaMalloc( &dev_z, sizeof( int ) );

    x = (float*)malloc( size );
    y = (float*)malloc( size );
    z = (int*)malloc( sizeof( int ) )
    
    llenarRandom( x );
    llenarRandom( y );
    
    cudaMemcpy( dev_x, x, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_y, y, size, cudaMemcpyHostToDevice );
    
    dim3 numeroHilos(T, T);
    dim3 numeroBloques(N/T, N/T);
    // --------------------
    cudaEventRecord(start);
    aproximarPi<<< numeroHilos, numeroBloques  >>>( dev_x, dev_y, dev_z );
    cudaEventRecord(end);
    //---------------------
    cudaMemcpy( z, dev_z, sizeof( int ) , cudaMemcpyDeviceToHost );
    // sincronizar
    cudaEventSynchronize(end);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, end);
    printf("Pi: %f\n", 4.0 * z / N_TOTAL );
    printf("Tiempo: %f\n", milliseconds);

    free( a ); free( b ); free( c );
    cudaFree( dev_a ); cudaFree( dev_b ); cudaFree( dev_c );
    return 0;
}

void llenarRandom(float *v){
    srand(time(NULL));
    for(int i = 0 ; i < N ; i++ ) {
        v[i] = rand() / (RAND_MAX + 1.0f);
    }
}