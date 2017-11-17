#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// curand
#define N_TOTAL 4194304
#define N 2048
#define T 32

__global__ void aproximarPi( float *x, float *y, int *z ) {
    int i = threadIdx.x + blockIdx.x*blockDim.x; // 0 - 2047
    int j = threadIdx.y + blockIdx.y*blockDim.y; // 0 - 2047
    int index = j + i*N; // 0 - 4194303     
    
    if( (x[index] * x[index] + y[index] * y[index]) <= 1.0f){
        atomicAdd(z, 1);
    }

}

void llenarRandom(float *a, float *b);

int main( void ) {
    float *x, *y, pi; 
    int *z;
    float *dev_x, *dev_y;
    int *dev_z; 
    int size = N_TOTAL * sizeof( float );
    
    cudaEvent_t start, end;
    cudaEventCreate(&start); 
    cudaEventCreate(&end);
    
    cudaMalloc( &dev_x, size );
    cudaMalloc( &dev_y, size );
    cudaMalloc( &dev_z, sizeof( int ) );

    x = (float*)malloc( size );
    y = (float*)malloc( size );
    z = (int*)malloc( sizeof( int ) );
    
    llenarRandom( x, y );
    *z = 0;
    
    cudaMemcpy( dev_x, x, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_y, y, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_z, z, sizeof( int ) , cudaMemcpyHostToDevice );
    
    dim3 numeroHilos(T, T);
    dim3 numeroBloques(N/T, N/T);
    // --------------------
    cudaEventRecord(start);
    aproximarPi<<< numeroBloques, numeroHilos  >>>( dev_x, dev_y, dev_z );
    cudaEventRecord(end);
    //---------------------
    cudaMemcpy( z, dev_z, sizeof( int ) , cudaMemcpyDeviceToHost );
    // sincronizar
    cudaEventSynchronize(end);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, end);

    pi = ((4.0 * (*z)) / N_TOTAL);
    printf("z: %d\n", *z );
    printf("Pi: %f\n", pi );
    printf("Tiempo: %f\n", milliseconds);

    free( x ); free( y ); free( z );
    cudaFree( dev_x ); cudaFree( dev_y ); cudaFree( dev_z );
    return 0;
}

void llenarRandom(float *a, float *b){
    float *v1, *v2;
    v1 = a;
    v2 = b;
    srand(time(NULL));
    for(int i = 0 ; i < N_TOTAL ; i++ ) {
        v1[i] = rand() / (RAND_MAX + 1.0f);
        v2[i] = rand() / (RAND_MAX + 1.0f);
    }
}