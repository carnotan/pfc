#include <thrust/host_vector.h>
#include "data_structures.h"
#include <stdio.h>
#include "auxiliares.h"

#ifndef MIN
#define MIN(x,y) ((x < y) ? x : y)
#endif

/**
 * Función auxiliar que determina se un enteiro é potencia de 2 . 
 * @param x Enteiro a evaluar.
 * @return 1 se é potencia de 2, 0 noutro caso.
 */
extern "C"
bool isPow2(unsigned int x) {

    return ((x & (x - 1)) == 0);

}

/**
 * Función auxiliar para atopar a potencia de 2 superior máis cercana.
 * @param x Enteiro a evaluar.
 * @return  O enteiro potencia de 2 máis próximo a x.
 */
unsigned int nextPow2(unsigned int x) {

    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;

}

/**
 * Kernel xenérico para a redución de arrays baseado na sexta implementación de 
 * NVIDIA cuda samples. Reduce múltiples elementos por thread para aumentar a 
 * carga computacional de cada un . Realiza a suma en memoria compartida.
 * 
 * @g_idata Punteiro a un array de datos de entrada en memoria global.
 * @g_odata Punteiro a un array de datos de saída en memoria global.
 * @n Tamaño do array que se vai a reducir.
 */
template <class T, unsigned int blockSize, bool nIsPow2>
__global__ void
reduce6(T *g_idata, T *g_odata, unsigned int n) {

    T *sdata = SharedMemory<T>();
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockSize * 2 + threadIdx.x;
    unsigned int gridSize = blockSize * 2 * gridDim.x;
    T mySum = 0;

    //reducimos varios elementos por thread (determinado polo número de bloques)
    while (i < n) {
        mySum += g_idata[i];
        //Evitamos ler fora de límites 
        if (nIsPow2 || i + blockSize < n)
            mySum += g_idata[i + blockSize];
        i += gridSize;
    }
    // cada thread escribe a súa suma a memoria compartida.
    sdata[tid] = mySum;
    __syncthreads();
    // Redución en memoria compartida. Só a primeira metade dos threads realizan
    // escrituras para evitar conflictos en bancos de memoria.
    if ((blockSize >= 512) && (tid < 256)) {
        sdata[tid] = mySum = mySum + sdata[tid + 256];
    }
    __syncthreads();
    if ((blockSize >= 256) &&(tid < 128)) {
        sdata[tid] = mySum = mySum + sdata[tid + 128];
    }
    __syncthreads();
    if ((blockSize >= 128) && (tid < 64)) {
        sdata[tid] = mySum = mySum + sdata[tid + 64];
    }
    __syncthreads();

#if (__CUDA_ARCH__ >= 300 )
    if (tid < 32) {
        //Esta instrucción só está dispoñibe en arquitecturas kepler ou 
        superiores.
        //Emprega a función __shfl_down para realizar a reducción a nivel de 
                //warp
        if (blockSize >= 64) mySum += sdata[tid + 32];
        for (int offset = warpSize / 2; offset > 0; offset /= 2) {
            mySum += __shfl_down(mySum, offset);
        }
    }
#else

    //En arquitecturas previas a kepler, facemos desplegue completo do último 
    //warp (loop unroll).Aumenta a carga computacional sacrificando un
    //executable maior (space-time tradeoff)
    if ((blockSize >= 64) && (tid < 32)) {
        sdata[tid] = mySum = mySum + sdata[tid + 32];
    }
    __syncthreads();
    if ((blockSize >= 32) && (tid < 16)) {
        sdata[tid] = mySum = mySum + sdata[tid + 16];
    }
    __syncthreads();
    if ((blockSize >= 16) && (tid < 8)) {
        sdata[tid] = mySum = mySum + sdata[tid + 8];
    }
    __syncthreads();
    if ((blockSize >= 8) && (tid < 4)) {
        sdata[tid] = mySum = mySum + sdata[tid + 4];
    }
    __syncthreads();
    if ((blockSize >= 4) && (tid < 2)) {
        sdata[tid] = mySum = mySum + sdata[tid + 2];
    }
    __syncthreads();
    if ((blockSize >= 2) && (tid < 1)) {
        sdata[tid] = mySum = mySum + sdata[tid + 1];
    }

    __syncthreads();
#endif

    // O thread 0 escribe o resultado de volta a memoria global.
    if (tid == 0) g_odata[blockIdx.x] = mySum;
}

/**
 * Función auxiliar que obtén o número de bloques e threads a empregar na
 *  redución.
 * 
 * @param n Tamaño do vector a reducir
 * @param maxBlocks Número máximo de bloques.
 * @param maxThreads Número máximo de threads por bloque.
 * @param blocks Número de bloques.
 * @param threads Número de threads.
 * @return 0 se se obtiveron os parámetros correctamente, -1 se non é posible
 * reducir.
 */
int getNumBlocksAndThreads(int n, int maxBlocks, int maxThreads, int &blocks, 
        int &threads, int device, cudaDeviceProp prop) {

    threads = (n < maxThreads * 2) ? nextPow2((n + 1) / 2) : maxThreads;
    blocks = (n + (threads * 2 - 1)) / (threads * 2);
    if ((float) threads * blocks > (float) prop.maxGridSize[0] 
            * prop.maxThreadsPerBlock) {
        printf("Tamaño do array demasiado grande para facer a redución\n");
        return -1;
    }
    blocks = MIN(maxBlocks, blocks);
    return 0;
}

/**
 * Función wrapper para escoller a instancia correcta do kernel de reducción.
 * @param size Tamaño do vector a reducir.
 * @param threads Número de threads por bloque.
 * @param blocks Número de bloques. 
 * @param d_idata Punteiro ao array de datos de entrada en memoria global.
 * @param d_odata Punteiro ao array de datos de saída en memoria global.
 */
template <class T>
void
reduce(int size, int threads, int blocks, T *d_idata, T *d_odata) {

    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    int smemSize = (threads <= 32) ? 2 * threads * sizeof (T) : threads *
    sizeof (T);

    if (isPow2(size)) {
        switch (threads) {
            case 512:
                reduce6<T, 512, true> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 256:
                reduce6<T, 256, true> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 128:
                reduce6<T, 128, true> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 64:
                reduce6<T, 64, true> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 32:
                reduce6<T, 32, true> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 16:
                reduce6<T, 16, true> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 8:
                reduce6<T, 8, true> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 4:
                reduce6<T, 4, true> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 2:
                reduce6<T, 2, true> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 1:
                reduce6<T, 1, true> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;
        }
    } else {
        switch (threads) {
            case 512:
                reduce6<T, 512, false> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 256:
                reduce6<T, 256, false> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 128:
                reduce6<T, 128, false> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 64:
                reduce6<T, 64, false> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 32:
                reduce6<T, 32, false> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 16:
                reduce6<T, 16, false> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 8:
                reduce6<T, 8, false> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 4:
                reduce6<T, 4, false> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 2:
                reduce6<T, 2, false> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;

            case 1:
                reduce6<T, 1, false> << < dimGrid, dimBlock, smemSize >>>
                        (d_idata, d_odata, size);
                break;
        }
    }

}

/**
 * Función auxiliar que chama repetidas veces ao kernel de reducción a ta que 
 * completa a lonxitude do array a reducir.
 * 
 * @param n Tamaño do array a reducir.
 * @param numThreads Número de threads por bloque.
 * @param numBlocks Número de bloques.
 * @param maxThreads Número máximo de threads.
 * @param maxBlocks Número máximo de bloques.
 * @param d_salida Punteiro ao array de datos en memoria global que almacena as 
 * sucesivas reducións (faise unha soa copia ao final).
 * @param d_idata Punteiro ao array de datos de entrada en memoria global.
 * @param d_odata Punteiro ao array de datos de saída en memoria global.
 * @param pos Posición na que escribimos os datos en d_salida
 */
template <class T>
void performReduction(int n, int numThreads, int numBlocks, int maxThreads,
        int maxBlocks, T *d_salida, T *d_idata, T *d_odata, int pos, int device,
        cudaDeviceProp prop) {

    int s = numBlocks;

    reduce<T>(n, numThreads, numBlocks, d_idata, d_odata);
    gpuErrchk(cudaPeekAtLastError());
    while (s > 1) {
        int threads = 0, blocks = 0;
        getNumBlocksAndThreads(s, maxBlocks, maxThreads, blocks, threads, device
                ,prop);
        reduce<T>(s, threads, blocks, d_odata, d_odata);
        s = (s + (threads * 2 - 1)) / (threads * 2);

    }
    // Facemos as copias en memoria global para , unha vez analizada toda a
    // poboación, traernos a estrutura completa (máis eficiente)
    cudaMemcpy(d_salida + pos, d_odata, sizeof (T), cudaMemcpyDeviceToDevice);

}

//Instanciamos as funcións de redución para os dous tipos necesarios. Esto
//indícalle ao compilador que ten que crear o código para estes tipos .
//(Instanciación explícita) Ver: http://www.cplusplus.com/forum/articles/14272/

template void
reduce<int>(int size, int threads, int blocks,
        int *d_idata, int *d_odata);

template void
reduce<float>(int size, int threads, int blocks,
        float *d_idata, float *d_odata);

template void
performReduction(int n, int numThreads, int numBlocks, int maxThreads,
        int maxBlocks, int *d_salida, int *d_idata, int *d_odata, int pos,
        int device, cudaDeviceProp prop);

template void
performReduction(int n, int numThreads, int numBlocks, int maxThreads,
        int maxBlocks, float *d_salida, float *d_idata, float *d_odata, 
        int pos,int device, cudaDeviceProp prop);
