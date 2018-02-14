#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
// Para PI, Producto punto
#define N 20000
// Para las matrices
#define M_N 100
#define M_M 100
#define M_L 100
// Para los calor 2d
#define NX_2D 100
#define NY_2D 100  
const float gama = 0.001f;
const float dt = 0.001;
const float h = (float)(1.0f/NX_2D); 
const float h2 = h * h;
float *ptr;
#define SWAP(ptr,x,y) {ptr=&x[0]; x=&y[0]; y=ptr;}
// Para los Bloques
#define B_NC 10
#define B_M 10
#define B_NUM_BLOQUES 100
// Repeticiones y numero de n's
#define NUM_N 10
#define REP 100
// curand
#define THREADS_PER_BLOCK 10


__global__ void aproximarPi( float *x, float *y, int *z, int tam) {
    int i = threadIdx.x + blockIdx.x*blockDim.x; // 0 - 2047
    int j = threadIdx.y + blockIdx.y*blockDim.y; // 0 - 2047
    int index = j + i*tam; // 0 - 4194303     
    
    if( (x[index] * x[index] + y[index] * y[index]) <= 1.0f){
        atomicAdd(z, 1);
    }
}

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

__global__ void suma( int *a, int *b, int *c, int n, int m) {
    int index = blockIdx.x + blockIdx.y * blockDim.y;
    if(index < n*m){
        c[index] = a[index] + b[index];
    }  
}

__global__ void multiplicacion( int *a, int *b, int *c, int n, int m, int l ) {
    int i = threadIdx.x + blockIdx.x*blockDim.x; 
    int j = threadIdx.y + blockIdx.y*blockDim.y; 
        
    c[j+i*l] = 0;

    for(int k=0 ; k < m ; k++ ){
        c[j+i*l] += a[k+i*m] * b[j+k*l];
    }
}

__global__ void llenarMatriz(float* __restrict__ dev_u, int nx, int ny){ 
	int x = threadIdx.x * blockDim.x + threadIdx.x;
	int y = threadIdx.y * blockDim.y + threadIdx.y;

	if( x >= nx || y >= ny){ return; }
	int g = x + y * ( blockDim.x * gridDim.x );

	dev_u[g]=0.0;
	if( (x*x+y*y) < nx*ny ){ dev_u[g]=3.0; }
}

__device__ void fronteras( float* __restrict__ dev_u, int i, int j, int nx, int ny){
	// Bandas izquierda y derecha
	dev_u[1+ny*j]=dev_u[2+ny*j];
	dev_u[(nx)+ny*j]=dev_u[(nx-1)+ny*j];
	// Bandas superior e inferior
	dev_u[i+ny]=dev_u[i+ny*2];
	dev_u[i+ny*(ny-1)]=dev_u[i+ny*ny]; 
}

__global__ void euler( float* __restrict__ dev_u, float* __restrict__ dev_ut, int nx, int ny ){
	int i = threadIdx.x + blockIdx.x*blockDim.x; // 0 - 256
   	int j = threadIdx.y + blockIdx.y*blockDim.y; // 0 - 255

	fronteras( dev_u, i, j, nx, ny);
	__syncthreads();
	/* Euler Scheme */
	if( i*(nx-i-1)*j*(ny-j-1) != 0){
		dev_ut[j+i*ny] = dev_u[i+ny*j] + ( dev_u[i-1+ny*j] + dev_u[i+1+ny*j] + dev_u[i+(j-1)*ny] 
						+ dev_u[i+(j+1)*ny] - 4.0 * dev_u[i+ny*j] ) * dt * gama / h2;
	}
}

__global__ void calcularCRS(int *val, int *col_ind, int *row_ptr, int *u, int *resultado, int l ){
    int i = threadIdx.x + blockIdx.x*blockDim.x; // 0 - 9 
    int j = threadIdx.y + blockIdx.y*blockDim.y; // 0 - 9
    int suma = 0;

    for(int k = row_ptr[i]-1; k < row_ptr[i+1]-1; k++){
        suma += val[k] * u[j + ( (col_ind[k]-1) * l) ];
    }
    resultado[j+i*l] = suma;
}

__global__ void calcularBloques(int *matriz, int *u, int *resultado, int num_bloques, int nc, int m ){
    int index1 = threadIdx.x + blockIdx.x*blockDim.x; // 0 - 1 
    int index2 = threadIdx.y + blockIdx.y*blockDim.y; // 0 - 1
    int suma = 0;

    for(int i=0 ; i < num_bloques ; i++){ 
        suma = 0;   
        for(int l=0 ; l < nc ; l++){
            suma += matriz[l+index1*nc] * u[index2+m*(l+i*nc)];
        }
        resultado[index2 + m*(index1+i*nc)] = suma;
    }
}

void calcularTiempoPI(int);
void calcularTimepoProductoPunto(int);
void calcularTiempoSumaMatrices(int, int);
void calcularTiempoMultiplicacionMatrices(int, int, int);
void calcularTiempoCalor2D(int, int);
void calcularTiempoCRS(int, int, int);
void calcularTiempoBloques(int, int, int);

int main( void ) {

    for(int i=1 ; i < NUM_N+1 ; i++){
		calcularTiempoPI(N*i);
	}

    for(int i=1 ; i < NUM_N+1 ; i++){
		calcularTimepoProductoPunto(N*i);
	}

    for(int i=1 ; i < NUM_N+1 ; i++){
		calcularTiempoSumaMatrices(M_N*i, M_M*i);
	}

	for(int i=1 ; i < NUM_N+1 ; i++){
		calcularTiempoMultiplicacionMatrices(M_N*i, M_M*i, M_L*i);
	}

	for(int i=1 ; i < NUM_N+1 ; i++){
		calcularTiempoCalor2D(NX_2D, NY_2D);
	}

	for(int i=1 ; i < NUM_N+1 ; i++){
		calcularTiempoCRS(M_N*i, M_M*i, M_L*i);
	}

	for(int i=1 ; i < NUM_N+1 ; i++){
		calcularTiempoBloques(B_NC*i, B_M, B_NUM_BLOQUES*i);
	}

    return 0;
}

void calcularTiempoPI(int n){
    float *x, *y, pi; 
    int *z;
    float *dev_x, *dev_y;
    int *dev_z;
    // Para la generalizacion >:(
    int tam = ( (int)sqrt(n) ) / THREADS_PER_BLOCK ;
    tam *= THREADS_PER_BLOCK; 
    // ------------------------------
    int size = tam * sizeof( float );
    float milliseconds = 0;
    float runtime = 0;
    // Tiempo
    cudaEvent_t start, end;
    cudaEventCreate(&start); 
    cudaEventCreate(&end);
    // memoria device
    cudaMalloc( &dev_x, size );
    cudaMalloc( &dev_y, size );
    cudaMalloc( &dev_z, sizeof( int ) );
    // memoria host
    x = (float*)malloc( size );
    y = (float*)malloc( size );
    z = (int*)malloc( sizeof( int ) );
    // llenar vectores
    srand(time(NULL));
    for(int i = 0 ; i < tam*tam ; i++ ) {
        x[i] = rand() / (RAND_MAX + 1.0f);
        y[i] = rand() / (RAND_MAX + 1.0f);
    }
    *z = 0;
    // enviar a de host to divice
    cudaMemcpy( dev_x, x, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_y, y, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_z, z, sizeof( int ) , cudaMemcpyHostToDevice );
    
    dim3 numeroHilos(THREADS_PER_BLOCK, THREADS_PER_BLOCK);
    dim3 numeroBloques(tam/THREADS_PER_BLOCK, tam/THREADS_PER_BLOCK);
    // --------------------
    for(int r=0 ; r<REP ; r++){
        cudaEventRecord(start);
        aproximarPi<<< numeroBloques, numeroHilos  >>>( dev_x, dev_y, dev_z, tam );
        cudaEventRecord(end);
        //---------------------
        cudaMemcpy( z, dev_z, sizeof( int ) , cudaMemcpyDeviceToHost );
        // sincronizar
        cudaEventSynchronize(end);   
        milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, end);
        //pi = ((4.0 * (*z)) / (tam*tam) );
        //printf("Pi: %f\n", pi );
        runtime += milliseconds;
    }
    runtime /= REP;
    printf("PI ) Tiempo Total[N = %d]: %f\n", n, runtime);

    free( x ); free( y ); free( z );
    cudaFree( dev_x ); cudaFree( dev_y ); cudaFree( dev_z );
}

void calcularTimepoProductoPunto(int n){
    int *a, *b, *c; 
    int *dev_a, *dev_b, *dev_c; 
    int size = n * sizeof( int ); 
    float milliseconds = 0;
    float runtime = 0;
    // Tiempo
    cudaEvent_t start, end;
    cudaEventCreate(&start); 
    cudaEventCreate(&end);

    cudaMalloc( &dev_a, size );
    cudaMalloc( &dev_b, size );
    cudaMalloc( &dev_c, sizeof( int ) );

    a = (int*)malloc( size );
    b = (int*)malloc( size );
    c = (int*)malloc( sizeof( int ) );
    // inicializo los vectores a y b
    for(int i=0 ; i < n ; i++){ a[i] = 1; b[i] = 1; }
    // mandas a y b al GPU
    cudaMemcpy( dev_a, a, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_b, b, size, cudaMemcpyHostToDevice );
    for(int r=0 ; r<REP ; r++){
        // realizo el calculo
        cudaEventRecord(start);
        dot<<< n/100, 100 >>>( dev_a, dev_b, dev_c );
        cudaEventRecord(end);
        // obtengo el valor del producto de la GPU
        cudaMemcpy( c, dev_c, sizeof( int ) , cudaMemcpyDeviceToHost );
        // sincronizar
        cudaEventSynchronize(end);   
        milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, end);
        runtime += milliseconds;
        //printf("El producto punto es: %d\n", *c);
    }
    runtime /= REP;
    printf("PP ) Tiempo Total[N = %d]: %f\n", n, runtime);
    // libera la memoria
    free( a ); free( b ); free( c );
    cudaFree( dev_a ); cudaFree( dev_b ); cudaFree( dev_c );
}

void calcularTiempoSumaMatrices(int n, int m){
    int *a, *b, *c;
    int *dev_a, *dev_b, *dev_c; 
    int size = n * m * sizeof( int ); 
    float milliseconds = 0;
    float runtime = 0;
    // Tiempo
    cudaEvent_t start, end;
    cudaEventCreate(&start); 
    cudaEventCreate(&end);

    cudaMalloc( &dev_a, size );
    cudaMalloc( &dev_b, size );
    cudaMalloc( &dev_c, size );

    a = (int*)malloc( size );
    b = (int*)malloc( size );
    c = (int*)malloc( size );
    
    for(int i=0 ; i < n ; i++){ 
        for(int j=0 ; j < m ; j++){
            a[j+i*m] = 1; 
        }
    }

    for(int i=0 ; i < n ; i++){ 
        for(int j=0 ; j < m ; j++){
            b[j+i*m] = 1; 
        }
    }

    cudaMemcpy( dev_a, a, size, cudaMemcpyHostToDevice );
    cudaMemcpy( dev_b, b, size, cudaMemcpyHostToDevice );
    
    dim3 block(THREADS_PER_BLOCK, THREADS_PER_BLOCK);
    dim3 grid(n/THREADS_PER_BLOCK, m/THREADS_PER_BLOCK);

    for(int r=0 ; r<REP ; r++){
        cudaEventRecord(start);
        suma<<<grid,block>>>( dev_a, dev_b, dev_c, n, m );
        cudaEventRecord(end);
        // sincronizar
        cudaEventSynchronize(end);   
        milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, end);
        runtime += milliseconds;
        cudaMemcpy( c, dev_c, size , cudaMemcpyDeviceToHost );
        //printf("valores: C[0,0]= %d , C[2047,2047]= %d\n", c[0], c[N*N-1]);
    }
    runtime /= REP;
    printf("SM ) Tiempo Total[N = %d, M = %d]: %f\n", n, m, runtime);
    
    free( a ); free( b ); free( c );
    cudaFree( dev_a ); cudaFree( dev_b ); cudaFree( dev_c );
}

void calcularTiempoMultiplicacionMatrices(int n, int m, int l){
    int *a, *b, *c;
    int *dev_a, *dev_b, *dev_c; 
    float milliseconds = 0;
    float runtime = 0;
    // Tiempo
    cudaEvent_t start, end;
    cudaEventCreate(&start); 
    cudaEventCreate(&end);

    cudaMalloc( &dev_a, n * m * sizeof( int ) );
    cudaMalloc( &dev_b, m * l * sizeof( int ) );
    cudaMalloc( &dev_c, n * l * sizeof( int ) );

    a = (int*)malloc( n * m * sizeof( int ) );
    b = (int*)malloc( m * l * sizeof( int ) );
    c = (int*)malloc( n * l * sizeof( int ) );
    
    for(int i=0 ; i < n ; i++){ 
        for(int j=0 ; j < m ; j++){
            a[j+i*m] = 1; 
        }
    }

    for(int i=0 ; i < m ; i++){ 
        for(int j=0 ; j < l ; j++){
            m[j+i*l] = 1; 
        }
    }
    
    cudaMemcpy( dev_a, a, n * m * sizeof( int ) , cudaMemcpyHostToDevice );
    cudaMemcpy( dev_b, b, m * l * sizeof( int ) , cudaMemcpyHostToDevice );
    
    dim3 numeroHilos(THREADS_PER_BLOCK, THREADS_PER_BLOCK);
    dim3 numeroBloques(n/THREADS_PER_BLOCK, l/THREADS_PER_BLOCK);

    // --------------------
    for(int r=0 ; r<REP ; r++){
        cudaEventRecord(start);
        multiplicacion<<< numeroHilos, numeroBloques  >>>( dev_a, dev_b, dev_c, n, m, l);
        cudaEventRecord(end);
        //---------------------
        cudaMemcpy( c, dev_c, n * l * sizeof( int ) , cudaMemcpyDeviceToHost );
        // sincronizar
        cudaEventSynchronize(end);
        milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, end);
        runtime += milliseconds;
    }
    runtime /= REP;
    printf("MM ) Tiempo Total[N = %d, M = %d]: %f\n", n, m, runtime);
    //printf("c[0][0] = %d , c[99][99] = %d\n", c[0], c[N*N-1] );

    free( a ); free( b ); free( c );
    cudaFree( dev_a ); cudaFree( dev_b ); cudaFree( dev_c );
}

void calcularTiempoCalor2D(int nx, int ny){
    int nbx = (nx+(THREADS_PER_BLOCK-1))/THREADS_PER_BLOCK;
	int nby = (ny+(THREADS_PER_BLOCK-1))/THREADS_PER_BLOCK;
	float *u;
	float *dev_u, *dev_ut;
    float milliseconds = 0;
    float runtime = 0;
    // Tiempo
    cudaEvent_t start, end;
    cudaEventCreate(&start); 
    cudaEventCreate(&end);
	
	cudaMalloc( &dev_u, nx * ny * sizeof(float) );
	cudaMalloc( &dev_ut, nx * ny * sizeof(float) );

	u = (float*)malloc( nx * ny * sizeof(float) );

	dim3 NB = dim3( nbx, nby, 1);
	dim3 TB = dim3( THREADS_PER_BLOCK, THREADS_PER_BLOCK, 1); 
	
	llenarMatriz<<< NB , TB>>>(dev_u, nx, ny);
    
    for(int r=0 ; r<REP ; r++){
        /* Main Iteration loop */
        cudaEventRecord(start);
        for(int ntime = 0; ntime < 100000 ; ntime++){
            euler<<< NB , TB>>>( dev_u, dev_ut, nx, ny );
            SWAP(ptr, dev_u, dev_ut);
        }
        cudaEventRecord(end);
        // sincronizar
        cudaEventSynchronize(end);
        milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, end);
        runtime += milliseconds;
    }
	runtime /= REP;
    printf("Calor2D ) Tiempo Total[Nx = %d, Ny = %d]: %f\n", nx, ny, runtime);
    //------------------------------
	free( u );
	cudaFree( dev_u ); cudaFree( dev_ut );
}

void calcularTiempoCRS(int n, int m, int l){
    int nnz = 0, contador = 0, indice = 0, fila = -1;
    int *matriz, *val, *col_ind, *row_ptr, *u, *resultado;
    int *dev_val, *dev_col_ind, *dev_row_ptr, *dev_u, *dev_resultado;
    float milliseconds = 0;
    float runtime = 0;
    // Tiempo
    cudaEvent_t start, end;
    cudaEventCreate(&start); 
    cudaEventCreate(&end);
    // llenar la matriz
    matriz = (int*)malloc( n * m * sizeof(int));

    for(int i=0 ; i < n ; i++){
        for(int j=0; j < m; j++){
            if( j >= i ){
                matriz[j+i*m] = 1;
            }else{
                matriz[j+i*m] = 0;
            }
        }
    }
    // GENERAR CRS ----------------------------------------------------
    // contar los no ceros (nnz)
    for(int i=0 ; i < n ; i++){
        for(int j=0; j < m; j++){
            if( matriz[j+i*m] != 0){
                nnz++;
            }
        }
    }
    // llenar val y col_ptr
    val = (int*)malloc( nnz * sizeof(int));
    col_ind = (int*)malloc( nnz * sizeof(int) );
    for(int i=0 ; i < n ; i++){
        for(int j=0; j < m; j++){
            if( matriz[j+i*m] != 0){
                val[contador] = matriz[j+i*m];
                col_ind[contador] = j+1;
                contador++;
            }
        }
    }

    // llenar row_ptr
    row_ptr = (int*)malloc( (n+1) * sizeof(int) );
    contador = 0; // contador = 0, indice = 0, fila = -1
    for(int i=0 ; i < n ; i++){
        for(int j=0; j < m; j++){
            if( matriz[j+i*m] != 0){
                contador++;
                if( i != fila){
                    row_ptr[indice] = contador;
                    fila = i;
                    indice++;
                }
            }
        }
    }
    row_ptr[n] = nnz+1;

    // lleno la matriz
    u = (int*)malloc( m * l * sizeof(int));
    
    for(int i = 0 ; i < m ; i++){
        for(int j = 0 ; j < l ; j++){
            u[j+i*l] = 1;
        }
    }

    resultado = (int*)malloc( n * l * sizeof(int));
    
    dim3 numeroHilos(THREADS_PER_BLOCK, THREADS_PER_BLOCK);
    dim3 numeroBloques(n/THREADS_PER_BLOCK, l/THREADS_PER_BLOCK);

    cudaMalloc( &dev_val,  nnz * sizeof(int) );
    cudaMalloc( &dev_col_ind, nnz * sizeof(int) );
    cudaMalloc( &dev_row_ptr, (n+1) * sizeof(int) );
    cudaMalloc( &dev_u, m * l * sizeof(int) );
    cudaMalloc( &dev_resultado, n * l * sizeof(int) );

    cudaMemcpy( dev_val, val, nnz * sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_col_ind, col_ind, nnz * sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_row_ptr, row_ptr, (n+1) * sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_u, u, m * l * sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_resultado, resultado, n * l * sizeof(int), cudaMemcpyHostToDevice );
    for(int r=0 ; r<REP ; r++){
        cudaEventRecord(start);
        calcularCRS<<< numeroBloques , numeroHilos >>>(dev_val, dev_col_ind, dev_row_ptr, dev_u, dev_resultado, l );
        cudaEventRecord(end);
        // sincronizar
        cudaEventSynchronize(end);
        milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, end);
        runtime += milliseconds;

        cudaMemcpy( resultado, dev_resultado, n * l * sizeof(int), cudaMemcpyDeviceToHost );
        /*// imprimir resultado
        for(int i=0; i < n; i++){
            for(int j=0 ; j < l ; j++){
                printf("resultado[%d][%d] = %d\n", i, j, resultado[j+i*l]);
            }
        }
        */
    }
    runtime /= REP;
    printf("CRS ) Tiempo Total[N = %d, M = %d, L = %d]: %f\n", n, m, l, runtime);

    free(matriz); free(val);
    free(col_ind); free(row_ptr);
    free(u); free(resultado);
    cudaFree(dev_val);
    cudaFree(dev_col_ind); cudaFree(dev_row_ptr);
    cudaFree(dev_u); cudaFree(dev_resultado);
}

void calcularTiempoBloques(int nc, int m, int num_bloques){
    int *matriz, *u, *resultado;
    int *dev_matriz, *dev_u, *dev_resultado;
    float milliseconds = 0;
    float runtime = 0;
    // Tiempo
    cudaEvent_t start, end;
    cudaEventCreate(&start); 
    cudaEventCreate(&end);

    matriz = (int*)malloc( nc * nc * sizeof(int));
    u = (int*)malloc( nc * num_bloques * m * sizeof(int));
    resultado = (int*)malloc( nc * num_bloques * m * sizeof(int));

    cudaMalloc( &dev_matriz, nc * nc * sizeof(int) );
    cudaMalloc( &dev_u, nc * num_bloques * m * sizeof(int) );
    cudaMalloc( &dev_resultado, nc * num_bloques * m * sizeof(int) );

    for(int i=0 ; i < nc ; i++){
        for(int j=0; j < nc; j++){
            matriz[j+i*nc] = j+i*nc+1;
        }
    }
    
    for(int i=0 ; i < nc*num_bloques ; i++){
        for(int j=0; j < m; j++){
            u[j+i*m] = j+i*m+1;
        }
    }

    cudaMemcpy( dev_matriz, matriz, nc * nc * sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_u, u, nc * num_bloques * m * sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy( dev_resultado, resultado, nc * num_bloques * m * sizeof(int), cudaMemcpyHostToDevice );

    dim3 numeroHilos(THREADS_PER_BLOCK, THREADS_PER_BLOCK);
    dim3 numeroBloques(nc/THREADS_PER_BLOCK, m/THREADS_PER_BLOCK);
    for(int r=0 ; r<REP ; r++){
        cudaEventRecord(start);
        calcularBloques<<< numeroBloques, numeroHilos>>>(dev_matriz, dev_u, dev_resultado, num_bloques, nc, m );
        cudaMemcpy( resultado, dev_resultado, nc * num_bloques * m * sizeof(int), cudaMemcpyDeviceToHost );
        cudaEventRecord(end);
        // sincronizar
        cudaEventSynchronize(end);
        milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, end);
        runtime += milliseconds;

        /*
        for(int i=0 ; i < nc*num_bloques ; i++){
            for(int j=0 ; j < m ; j++){
                printf("R[%d][%d] = %d\n", i, j, resultado[j+i*m]);
            }
        }
        */
    }
    runtime /= REP;
    printf("BLOQUES ) Tiempo Total[NC = %d, M = %d, NUM_BLOQUES = %d]: %f\n", nc, m, num_bloques, runtime);

    free(matriz); free(u); free(resultado);
    cudaFree(dev_matriz); cudaFree(dev_u); cudaFree(dev_resultado);
}
