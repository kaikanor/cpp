// #include "stdafx.h"	
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <ctime>
#include <C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.0\include\cuda_runtime.h>
//#define N 65599
#define s1 13
int kol = 0;
int N = 1024;
//using namespace std;
__constant__ int const_flag[1];

__global__ void FUN_KERNEL(float *Ad, float *Bd, float *Cd)
{
    int idx_thread = blockIdx.x * blockDim.x + threadIdx.x;
    Ad[idx_thread] = Ad[idx_thread] * (Bd[idx_thread]*Cd[idx_thread] + 34*Ad[idx_thread]*Ad[idx_thread]*Ad[idx_thread]) / (s1*s1*s1);
}
__global__ void TIME_CUDA(int *flagd)
{
    int kkk = 5;
 //   printf("==%d==\n", const_flag[0]);
    __threadfence_system();
}

//int main(int argc, char **argv)
float go_cuda();
float go_cpu();
int main()
{
    int now = clock();
    int koln = 1;
    for (N = 16384; N < 10*1000*1000; N*=2) {
        printf("\nN = %d\n", N);
        printf("cpu: "); kol = 0;
        while (kol < koln) printf("%f  ", go_cpu());
        printf("\ncuda: "); kol = 0;
        while (kol < koln) printf("%f  ", go_cuda());
    }
    printf("\nThe end!   %d\n", (clock() - now));
    return 0;
}

float go_cpu()
{
    kol++;
    float *A = new float[N];
    float *B = new float[N];
    float *C = new float[N];
    float elapsedTime;
    cudaEvent_t start, stop;
    int er1 = cudaEventCreate(&start);
    int er2 = cudaEventCreate(&stop);
    for (int i = 0; i < N; i++) { A[i] = i; B[i] = i; C[i] = i; }
//    int time = clock();
//    system("pause");
    int *flag = new int; 
    int *flago;
    flag[0] = 0;
    cudaMalloc((void**) &flago, sizeof(int));
    int blocks = 1, blocksize = 1;
    cudaMemcpy(flago, flag, sizeof(int), cudaMemcpyHostToDevice);
    TIME_CUDA<<<blocks, blocksize>>>(flago);
    int host_flag[1];
    host_flag[0] = 111;
    cudaMemcpyToSymbol("const_flag", host_flag, sizeof(int), 0, cudaMemcpyHostToDevice);

//    cudaMemcpyFromSymbol(&kkk, "const_flag"const char * symbol, size_t count, size_t offset, enum cudaMemcpyKind kind);
//    cudaMemcpyToSymbol(const char * symbol, const void * src, size_t count, size_t offset, enum cudaMemcpyKind kind);
//    cudaMemcpyFromSymbol(void * dst, const char * symbol, size_t count, size_t offset, enum cudaMemcpyKind kind);

//    cudaStream_t stream;
//    cudaStreamCreate(&stream);
    int er3 = cudaEventRecord(start, 0);
    int time = clock();
    for (int jj = 0; jj < 10; jj++) {
        for (int i = 0; i < N; i++) A[i] = A[i] * (B[i]*C[i] + 34*A[i]*A[i]*A[i]) / (s1*s1*s1);
    }
    time = clock() - time;

//    cudaStreamSynchronize(0);
//    cudaStreamDestroy(stream);
    int er4 = cudaEventRecord(stop, 0);

//    flag[0] = 1;
//    cudaMemcpy(flago, flag, sizeof(int), cudaMemcpyHostToDevice);
//    TIME_CUDA<<<blocks, blocksize>>>(flago);

    int er5 = cudaEventSynchronize(stop);
    int er6 = cudaEventElapsedTime(&elapsedTime, start, stop);
//    for (int i = 0; i <10; i++) { if ((i % 1) == 0) printf("M[%d]=%.0f  ", i, A[i]); if (((i + 1) % 10000) == 0) printf("\n");  }
    return (float)time / 10;
//    return elapsedTime;
}

float go_cuda()
{
    kol++;
    float *A = new float[N];
    float *B = new float[N];
    float *C = new float[N];
    float elapsedTime;
    cudaEvent_t start, stop;
    int er7 = cudaEventCreate(&start);
    int er8 = cudaEventCreate(&stop);
    for (int i = 0; i < N; i++) { A[i] = i; B[i] = i; C[i] = i; }
	int blocks, blocksize;
    float *devA, *devB, *devC;
	cudaMalloc((void**) &devA, N*sizeof(float));
    cudaMalloc((void**) &devB, N*sizeof(float));
    cudaMalloc((void**) &devC, N*sizeof(float));
    blocksize = 65536;
    blocks = (int)N / blocksize; 
//    int time = clock();
    int er9 = cudaEventRecord(start, 0);

    cudaMemcpy(devA, A, N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devB, B, N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(devC, C, N*sizeof(float), cudaMemcpyHostToDevice);
//    system("pause");
    for (int jj = 0; jj < 1; jj++) {
        FUN_KERNEL<<<blocks, blocksize>>>(devA, devB, devC);
    }

    
    cudaMemcpy(A, devA, N*sizeof(float), cudaMemcpyDeviceToHost);
//    time = clock() - time;
    cudaFree(devA);    
    int era = cudaEventRecord(stop, 0);

    int erb = cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
//    for (int i = 0; i <10; i++) { if ((i % 1) == 0) printf("M[%d]=%.0f  ", i, A[i]); if (((i + 1) % 10000) == 0) printf("\n");  }
//    return time;
    return elapsedTime;
}