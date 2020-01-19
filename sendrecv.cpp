// ConsoleApplication8.cpp: определяет точку входа для консольного приложения.
//
#include "stdafx.h"
#include <iostream>
#include <ctime>
#include "mpi.h"
using namespace std;

#define N1 2000000
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
	float *T0 = new float[N1];
	float *T1 = new float[N1];

	if (size == 4) {
		float sr = 0;
		int kol = 1;
		now = clock();
		for (int k = 0; k < kol; k++) {
			for (int i = 0; i<N1; i++) { A[i] = B[i] = C[i] = D[i] = E[i] = G[i] = 2;	Y1[i] = Y2[i] = Y3[i] = T0[i] = T1[i] = 1; }
			if (rank == 0) {
				for (i = 0; i < N1; i++) { D[i] = A[i] - B[i]; }
				MPI_Send(D, N1, MPI_FLOAT, 2, 1, MPI_COMM_WORLD);
				for (i = 0; i < N1; i++) { E[i] = G[i] * D[i]; }
				MPI_Send(E, N1, MPI_FLOAT, 1, 1, MPI_COMM_WORLD);
				for (i = 0; i < N1; i++) { Y1[i] = E[i] + C[i]; }
				MPI_Send(Y1, N1, MPI_FLOAT, 3, 1, MPI_COMM_WORLD);
			}
			if (rank == 1) {
				for (i = 0; i < N1; i++) { T0[i] = 2 * C[i]; }
				for (i = 0; i < N1; i++) { A[i] = C[i] * C[i]; }
				MPI_Recv(T1, N1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
				for (i = 0; i < N1; i++) { B[i] = T1[i] * T0[i]; }
				MPI_Recv(T0, N1, MPI_FLOAT, 2, 1, MPI_COMM_WORLD, &status);
				for (i = 0; i < N1; i++) { D[i] = B[i] + T0[i]; }
				MPI_Recv(T0, N1, MPI_FLOAT, 2, 2, MPI_COMM_WORLD, &status);
				for (i = 0; i < N1; i++) { E[i] = D[i] / T0[i]; }
				for (i = 0; i < N1; i++) { Y1[i] = E[i] + A[i]; }
				for (i = 0; i < N1; i++) { G[i] = Y1[i] - C[i]; }
				for (i = 0; i < N1; i++) { Y2[i] = G[i] - T1[i]; }
				MPI_Recv(Y3, N1, MPI_FLOAT, 3, 1, MPI_COMM_WORLD, &status);
				for (i = 0; i < N1; i++) { Y3[i] = Y3[i] - T0[i]; }
			}
			if (rank == 2) {
				for (i = 0; i < N1; i++) { A[i] = G[i] * G[i]; }
				MPI_Recv(D, N1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
				for (i = 0; i < N1; i++) { T0[i] = D[i] * D[i]; }
				MPI_Recv(E, N1, MPI_FLOAT, 3, 1, MPI_COMM_WORLD, &status);
				for (i = 0; i < N1; i++) { C[i] = E[i] + D[i]; }
				MPI_Send(T0, N1, MPI_FLOAT, 1, 1, MPI_COMM_WORLD);
				MPI_Send(C, N1, MPI_FLOAT, 3, 1, MPI_COMM_WORLD);
				MPI_Send(A, N1, MPI_FLOAT, 1, 2, MPI_COMM_WORLD);
			}
			if (rank == 3) {
				for (i = 0; i < N1; i++) { A[i] = E[i] * D[i]; }
				for (i = 0; i < N1; i++) { C[i] = A[i] - B[i]; }
				MPI_Send(A, N1, MPI_FLOAT, 2, 1, MPI_COMM_WORLD);
				for (i = 0; i < N1; i++) { G[i] = C[i] * C[i]; }
				MPI_Recv(E, N1, MPI_FLOAT, 2, 1, MPI_COMM_WORLD, &status);
				for (i = 0; i < N1; i++) { Y1[i] = E[i] * G[i]; }
				MPI_Recv(D, N1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
				for (i = 0; i < N1; i++) { Y2[i] = D[i] * Y1[i]; }
				for (i = 0; i < N1; i++) { T0[i] = Y2[i] * C[i]; }
				MPI_Send(T0, N1, MPI_FLOAT, 1, 1, MPI_COMM_WORLD);
			}
			//cout << clock() - now << "\n";
		}
		sr = (clock() - now) / kol;
		cout << sr << "\n";
		if (rank == 1) {
			for (i = 1; i <= 1000000; i *= 100) { cout << "Y1[" << i - 1 << "] = " << Y1[i - 1] << "\n"; }
			for (i = 1; i <= 1000000; i *= 100) { cout << "Y2[" << i - 1 << "] = " << Y2[i - 1] << "\n"; }
			for (i = 1; i <= 1000000; i *= 100) { cout << "Y3[" << i - 1 << "] = " << Y3[i - 1] << "\n"; }
		}
	}
	else {
		cout << "size must == 4 \n exit \n";
	}

	MPI_Finalize();
	return 0;
}

