#include <stdio.h>
#include <stdlib.h>
#define N 512

__global__ void dot( int *a, int *b, int *c ){

    __shared__ int temp[N];
    temp[threadIdx.x] = a[threadIdx.x] * b[threadIdx.x];
    __syncthreads(); 
    if( 0 == threadIdx.x ){
        int sum = 0;
        for( int i = 0; i < N; i++ ){
            sum += temp[i];
        }
        *c = sum;
    }
}

void llenarVector(int*);

int main(){
    int *a, *b, *c; 
    int *dev_a, *dev_b, *dev_c; 
    int size = N * sizeof( int ); 
    
    cudaMalloc( (void**)&dev_a, size );
    cudaMalloc( (void**)&dev_b, size );
    cudaMalloc( (void**)&dev_c, sizeof( int ) );

    a = (int*)malloc( size );
    b = (int*)malloc( size );
    c = (int*)malloc( sizeof( int ) );
    // inicializo los vectores a y b
    llenarVector( a );
    llenarVector( b );
    // mandas a y b al GPU
    cudaMemcpy( dev_a, a, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_b, b, size, cudaMemcpyHostToDevice );
    // realizo el calculo
    dot<<< 1, N >>>( dev_a, dev_b, dev_c );
    // obtengo el valor del producto de la GPU
    cudaMemcpy( c, dev_c, sizeof( int ) , cudaMemcpyDeviceToHost );

    printf("El producto punto es: %d\n", *c);
    // libera la memoria
    free( a ); 
    free( b ); 
    free( c );
    cudaFree( dev_a );
    cudaFree( dev_b );
    cudaFree( dev_c );
    return 0;
}

void llenarVector(int *a){
    int *v;
    v = a;
    for(int i=0 ; i < N ; i++){ v[i] = 1; }
}
