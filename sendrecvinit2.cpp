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
	MPI_Status statuses[5];
	MPI_Request request[5];
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	float *A = new float[N1];
	float *B = new float[N1];
	float *C = new float[N1];
	float *D = new float[N1];
	float *E = new float[N1];
	float *G = new float[N1];
	float *Y1 = new float[N1];
	float *Y2 = new float[N1];
	float *Y3 = new float[N1];
	float T[6];
	float x[12];
	
	for (i = 0; i < 6; i++) { T[i] = 1; }
	if (size == 2) {
		float sr = 0;
		int kol = 1;
		now = clock();
		for (int k = 0; k < kol; k++) {
			for (i = 0; i < N1; i++) { A[i] = B[i] = C[i] = D[i] = E[i] = G[i] = 2;	Y1[i] = Y2[i] = Y3[i] = 1; }
			if (rank == 0) {
				MPI_Send_init(&T[1], 1, MPI_FLOAT, 1, 1, MPI_COMM_WORLD, &request[0]);
				MPI_Send_init(&T[2], 1, MPI_FLOAT, 1, 2, MPI_COMM_WORLD, &request[1]);
				MPI_Recv_init(&T[0], 1, MPI_FLOAT, 1, 5, MPI_COMM_WORLD, &request[4]);
				MPI_Recv_init(&T[5], 1, MPI_FLOAT, 1, 3, MPI_COMM_WORLD, &request[2]);
				MPI_Recv_init(&T[2], 1, MPI_FLOAT, 1, 4, MPI_COMM_WORLD, &request[3]);
			}
			if (rank == 1) {
				MPI_Recv_init(&T[1], 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &request[0]);
				MPI_Recv_init(&T[2], 1, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &request[1]);
				MPI_Send_init(&T[5], 1, MPI_FLOAT, 0, 3, MPI_COMM_WORLD, &request[2]);
				MPI_Send_init(&T[0], 1, MPI_FLOAT, 0, 5, MPI_COMM_WORLD, &request[4]);
				MPI_Send_init(&T[1], 1, MPI_FLOAT, 0, 4, MPI_COMM_WORLD, &request[3]);
			}
			for (i = 0; i < N1; i++) {
				if (rank == 0) {
					MPI_Start(&request[2]);
					MPI_Start(&request[4]);
					T[1] = A[i] - B[i];
					MPI_Start(&request[0]);
					T[2] = T[1] * G[i];
					MPI_Start(&request[1]);
					x[3] = T[1] * T[1];
					x[4] = 2 * T[2];
					MPI_Start(&request[3]);
					x[5] = x[4] * C[i];
					x[6] = C[i] * C[i];
					x[7] = x[3] + x[5];
					MPI_Wait(&request[2], &statuses[0]);
					x[8] = x[7] / T[5];
					Y1[i] = x[8] + x[6];
					MPI_Wait(&request[3], &statuses[1]);
					Y2[i] = Y1[i] - T[2];
					MPI_Wait(&request[4], &statuses[2]);
					Y3[i] = T[0];
					MPI_Waitall(2, &request[0], &statuses[3]);
				}
				if (rank == 1) {
					MPI_Start(&request[0]);
					MPI_Start(&request[1]);
					T[3] = E[i] * D[i];
					T[4] = T[3] - B[i];
					T[5] = G[i] * G[i];
					MPI_Start(&request[2]);
					MPI_Wait(&request[1], &statuses[0]);
					x[4] = T[2] + C[i];
					MPI_Wait(&request[0], &statuses[1]);
					x[5] = T[3] + T[1];
					x[6] = T[4] * T[4];
					x[7] = x[6] * T[4];
					x[8] = x[4] * x[5];
					x[9] = x[8] * x[7];
					Y3[i] = x[9] - T[5];
					T[1] = C[i] + T[2];
					MPI_Start(&request[3]);
					T[0] = Y3[i];
					MPI_Start(&request[4]);
					MPI_Waitall(3, &request[2], &statuses[2]);
				}
			}
		}
		sr = (clock() - now) / kol;
		cout << sr << "\n";
//				if (rank == 0) {
//					for (i = 1; i <= 1000000; i *= 100) { cout << "Y1[" << i - 1 << "] = " << Y1[i - 1] << "\n"; }
//					for (i = 1; i <= 1000000; i *= 100) { cout << "Y2[" << i - 1 << "] = " << Y2[i - 1] << "\n"; }
//					for (i = 1; i <= 1000000; i *= 100) { cout << "Y3[" << i - 1 << "] = " << Y3[i - 1] << "\n"; }
//				}
	}
	else {
		cout << "size must == 4 \n exit \n";
	}

	MPI_Finalize();
	return 0;
}

