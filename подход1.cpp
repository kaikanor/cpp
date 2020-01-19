// ConsoleApplication15.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <ctime>
#include <cmath>
#include "omp.h"
//#define N 2000
using namespace std;
int N;
//// Y = (A-B)*C + A/C

int inv_matr(float **matr, float **matri)
{
	float koeff;
	for (int k = 0; k < N; k++) {
		koeff = matr[k][k];
		for (int j = 0; j < N; j++) {
			matr[k][j] = matr[k][j] / koeff;
			matri[k][j] = matri[k][j] / koeff;
		}
		for (int i = 0; i < N; i++) {
			if (i != k) {
				koeff = matr[i][k];
				for (int j = 0; j < N; j++) {
					matr[i][j] = matr[i][j] - matr[k][j] * koeff;
					matri[i][j] = matri[i][j] - matri[k][j] * koeff;
				}
			}
		}
	}
	return 0;
}




int main(int argc, char **argv)
{
	N = 1024; // cin >> N;
	int kol = 1;
	cin >> kol;
	if ((kol < 1) || (kol > 12)) { return 0; }
	omp_set_num_threads(kol);
	int i, j, k, now = clock();
	int rank, size;
	float **A = new float*[N];
	float **B = new float*[N];
	float **C = new float*[N];
	float **Y = new float*[N];
	float **T1 = new float*[N];
	float **T2 = new float*[N];
	float **T3 = new float*[N];
	for (i = 0; i < N; i++) {
		A[i] = new float[N];
		B[i] = new float[N];
		C[i] = new float[N];
		Y[i] = new float[N];
		T1[i] = new float[N];
		T2[i] = new float[N];
		T3[i] = new float[N];
	}
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			A[i][j] = 1;
			B[i][j] = 2;
			C[i][j] = float(rand() % 999999) + 1;
		}
	}

#pragma omp parallel private(rank, size, i, j)
	{
		size = omp_get_num_threads();
		rank = omp_get_thread_num();
		for (i = rank * N / size; i < (rank + 1)*N / size; i++) {
			for (j = 0; j < N; j++) {
				T1[i][j] = A[i][j] - B[i][j];
			}
		}
	}

#pragma omp parallel private(rank, size, i, j, k)
	{
		size = omp_get_num_threads();
		rank = omp_get_thread_num();
		for (i = rank * N / size; i < (rank + 1)*N / size; i++) {
			for (j = 0; j < N; j++) {
				T2[i][j] = 0;
				for (k = 0; k < N; k++) {
					T2[i][j] += T1[i][k] * C[k][j];
				}
			}
		}
	}

#pragma omp parallel private(rank, size, i, j)
	{
		size = omp_get_num_threads();
		rank = omp_get_thread_num();
		for (i = rank * N / size; i < (rank + 1)*N / size; i++) {
			for (j = 0; j < N; j++) {
				Y[i][j] = C[i][j];
				if (i == j) { T3[i][j] = 1; }
				else { T3[i][j] = 0; }
			}
		}
	}
	int ttt = inv_matr(Y, T3);

#pragma omp parallel private(rank, size, i, j, k)
	{
		size = omp_get_num_threads();
		rank = omp_get_thread_num();
		for (i = rank * N / size; i < (rank + 1)*N / size; i++) {
			for (j = 0; j < N; j++) {
				T1[i][j] = 0;
				for (k = 0; k < N; k++) {
					T1[i][j] += A[i][k] * T3[k][j];
				}
			}
		}

	}

#pragma omp parallel private(rank, size, i, j) 
	{
		size = omp_get_num_threads();
		rank = omp_get_thread_num();
		for (i = rank * N / size; i < (rank + 1)*N / size; i++) {
			for (j = 0; j < N; j++) {
				Y[i][j] = T2[i][j] + T1[i][j];
			}
		}
	}


	cout << "N = " << N << "\n";
	float sum = 0;
	for (i = 0; i < 10; i++) { sum += abs(Y[i][i]); }
	cout << sum << "\n";
	cout << "clock = " << clock() - now << "\n";
//	system("pause");
	return 0;
}