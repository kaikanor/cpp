// ConsoleApplication8.cpp: определяет точку входа для консольного приложения.
//
#include "stdafx.h"
#include <iostream>
#include <ctime>
#include "mpi.h"
using namespace std;

#define N1 200000
int main(int argc, char **argv)
{
	int rank, size, now, i;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	//float m1[N1], m2[N1], m3[N1];
	float *A = new float[N1];
	float *B = new float[N1];
	float *C = new float[N1];
	float *D = new float[N1];
	float *E = new float[N1];
	float *G = new float[N1];
	float *Y1 = new float[N1];
	float *Y2 = new float[N1];
	float *Y3 = new float[N1];
	//	float *T1 = new float[N1];
	//	float *T2 = new float[N1];
	//	float *T3 = new float[N1];
	//	float *T4 = new float[N1];
	//	float *T5 = new float[N1];
	float T1, T2, T3, T4, T5;
	float x[12];

	if (size == 2) {
		float sr = 0;
		int kol = 100;
		now = clock();
		for (int k = 0; k < kol; k++) {
			for (int i = 0; i < N1; i++) { A[i] = B[i] = C[i] = D[i] = E[i] = G[i] = 2;	Y1[i] = Y2[i] = Y3[i] = 1; }
			for (i = 0; i < N1; i++) {
				if (rank == 0) {
					T1 = A[i] - B[i];
					MPI_Send(&T1, 1, MPI_FLOAT, 1, 1, MPI_COMM_WORLD);
					T2 = T1 * G[i];
					MPI_Send(&T2, 1, MPI_FLOAT, 1, 2, MPI_COMM_WORLD);
					x[3] = T1 * T1;
					MPI_Recv(&T5, 1, MPI_FLOAT, 1, 3, MPI_COMM_WORLD, &status);
					x[4] = 2 * T2;
					x[5] = x[4] * C[i];
					x[6] = C[i] * C[i];
					x[7] = x[3] + x[5];
					x[8] = x[7] / T5;
					Y1[i] = x[8] + x[6];
					MPI_Recv(&T2, 1, MPI_FLOAT, 1, 4, MPI_COMM_WORLD, &status);
					Y2[i] = Y1[i] - T2;
					MPI_Recv(&Y3[i], 1, MPI_FLOAT, 1, 5, MPI_COMM_WORLD, &status);
				}
				if (rank == 1) {
					T3 = E[i] * D[i];
					MPI_Recv(&T1, 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
					T4 = T3 - B[i];
					MPI_Recv(&T2, 1, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &status);
					T5 = G[i] * G[i];
					MPI_Send(&T5, 1, MPI_FLOAT, 0, 3, MPI_COMM_WORLD);
					x[4] = T2 + C[i];
					x[5] = T3 + T1;
					x[6] = T4 * T4;
					x[7] = x[6] * T4;
					x[8] = x[4] * x[5];
					x[9] = x[8] * x[7];
					Y3[i] = x[9] - T5;
					T1 = C[i] + T2;
					MPI_Send(&T1, 1, MPI_FLOAT, 0, 4, MPI_COMM_WORLD);
					MPI_Send(&Y3[i], 1, MPI_FLOAT, 0, 5, MPI_COMM_WORLD);
				}
			}
		}
		sr = (clock() - now) / kol;
		cout << sr << "\n";
		//		if (rank == 0) {
		//			for (i = 1; i <= 1000000; i *= 100) { cout << "Y1[" << i - 1 << "] = " << Y1[i - 1] << "\n"; }
		//			for (i = 1; i <= 1000000; i *= 100) { cout << "Y2[" << i - 1 << "] = " << Y2[i - 1] << "\n"; }
		//			for (i = 1; i <= 1000000; i *= 100) { cout << "Y3[" << i - 1 << "] = " << Y3[i - 1] << "\n"; }
		//		}
	}
	else {
		cout << "size must == 4 \n exit \n";
	}

	MPI_Finalize();
	return 0;
}

