#ifndef REDUCE_H
#define REDUCE_H


unsigned int nextPow2(unsigned int x);
bool isPow2(unsigned int x) ;

template <class T, unsigned int blockSize, bool nIsPow2>
__global__ void
reduce6(T *g_idata, T *g_odata, unsigned int n);

void getNumBlocksAndThreads(int n, int maxBlocks, int maxThreads, int &blocks,
        int &threads,int device, cudaDeviceProp prop);

template <class T>
void
reduce(int size, int threads, int blocks,
        T *d_idata, T *d_odata);


template <class T>
void performReduction(int n, int numThreads, int numBlocks, int maxThreads, 
        int maxBlocks, T *d_salida, T *d_idata, T *d_odata, int pos, int device, cudaDeviceProp prop) ;


#endif
