/*
Primero hay que cargar bien las librerías

  export LD_LIBRARY_PATH=/usr/local/cuda-8.0/lib64

luego, compilamos

  nvcc ejemplo.cu -o ejemplo -lcurand

*/


/* This program uses the host CURAND API to generate 100
 * pseudorandom floats.
 */
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

int main(int argc, char *argv[])
{
  size_t n = 100;
  size_t i;
  curandGenerator_t gen;
  float *devData, *hostData;
  /* Allocate n floats on host */
  hostData = (float *)malloc(n*sizeof(float));

  /* Allocate n floats on device */
     cudaMalloc((void **)&devData, n * sizeof(float));

  /* Create pseudo-random number generator */
     curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);

  /* Set seed */
     curandSetPseudoRandomGeneratorSeed(gen, (unsigned long long)clock());
  /* Generate n floats on device */
     curandGenerateUniform(gen, devData, n);
  /* Copy device memory to host */
     cudaMemcpy(hostData, devData, n * sizeof(float),
                       cudaMemcpyDeviceToHost);
  /* Set seed */
  //   curandSetPseudoRandomGeneratorSeed(gen, 1234ULL);

  /* Generate n floats on device */
  //   curandGenerateUniform(gen, devData, n);

  /* Copy device memory to host */
  //   cudaMemcpy(hostData, devData, n * sizeof(float), cudaMemcpyDeviceToHost);

  /* Show result */
  for(i = 0; i < n; i++) {
    printf("%1.4f \n", hostData[i]);
  }
  printf("\n");
  /* Cleanup */
  curandDestroyGenerator(gen);
  cudaFree(devData);
  free(hostData);
  return 0;
}
