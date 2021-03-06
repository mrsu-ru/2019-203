﻿#include "kiselevis.h"

/**
 * Введение в дисциплину
 */
void kiselevis::lab1()
{
	cout << "Hello, world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kiselevis::lab2()
{
double epsilon = 1.e-10;
	int *BeforeSwaped = new int[N];
	for (int i = 0; i < N; i++)
	{
		BeforeSwaped[i] = i;
	}
	for (int i = 0; i < N; i++)
	{
		int MaxIndex = i;
		for (int j = i+1; j < N; j++)
		{
			if (abs(A[j][i]) > abs(A[MaxIndex][i]))
			{
				MaxIndex = j;
			}
		}
		if (MaxIndex != i)
		{
			swap(A[MaxIndex], A[i]);
			swap(b[MaxIndex], b[i]);
			swap(BeforeSwaped[MaxIndex], BeforeSwaped[i]);
		}
		if (abs(A[i][i]) < epsilon)
		{
			cout << "NULL" << endl;
			return;
		}
		b[i] /= A[i][i];
		for (int j = N - 1; j >= i; j--)
		{
			A[i][j] /= A[i][i];
		}
		
		for (int j = 0; j < N; j++)
		{
			if (j != i)
			{
				b[j] -= A[j][i] * b[i];
				for (int k = N-1; k >= i; k--)
				{
					A[j][k] -= A[j][i] * A[i][k];
				}
				
			}

		}
	}
	for (int i = 0; i < N; i++) {
		x[i] = b[i];
	}
	for (int i = 0; i < N; i++) {
		if (BeforeSwaped[i] != i) {
			swap(A[BeforeSwaped[i]], A[i]);
			swap(b[BeforeSwaped[i]], b[i]);
		}
	}
	
	delete[] BeforeSwaped;
}



/**
 * Метод прогонки
 */
void kiselevis::lab3()
{
	double *p = new double[N];
	double *q = new double[N];

	// Прямой ход
	p[0] = -A[0][1] / A[0][0];
	q[0] = b[0] / A[0][0];
	for (int i = 1; i < N; i++) {
		p[i] = -A[i][i + 1] / (A[i][i] + A[i][i - 1] * p[i - 1]);
		q[i] = (b[i] - A[i][i - 1] * q[i - 1]) / (A[i][i] + A[i][i - 1] * p[i - 1]);
	}

	// Обратный ход
	x[N - 1] = q[N - 1];
	for (int i = N - 2; i >= 0; i--) {
		x[i] = p[i] * x[i + 1] + q[i];
	}

	delete[] p;
	delete[] q;
}



/**
 * Метод простых итераций
 */
void kiselevis::lab4()
{
	double eps = 1e-20;
	double tau = 1e-5;

	double* prevX = new double[N];
	
	while (true) {

		for (int i = 0; i < N; i++)
			prevX[i] = x[i];

		for (int i = 0; i < N; i++) {
			double sum = 0;
			for (int j = 0; j < N; j++) 
				sum += A[i][j] * prevX[j];
			x[i] = prevX[i] - tau * (sum - b[i]);
		}

		double maxErr = abs(x[0] - prevX[0]);
		for (int i = 1; i < N; i++)
			if (abs(x[i] - prevX[i]) > maxErr)
				maxErr = abs(x[i] - prevX[i]);

		if (maxErr < eps)
			break;

	}

	delete[] prevX;
}



/**
 * Метод Якоби или Зейделя
 */
void kiselevis::lab5()
{
	double eps = 1e-20;
	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}

	double *prev_x = new double[N];

	double norma = 0;
	do {
		for (int i = 0; i < N; i++) {
			prev_x[i] = x[i];
		}

		for (int i = 0; i < N; i++) {
			double result = b[i];
			for (int j = 0; j < N; j++) {
				if (i != j) {
					result -= (A[i][j] * prev_x[j]);
				}
			}

			x[i] = result / A[i][i];
		}

		norma = 0;
		for (int i = 0; i < N; i++) {
			if (abs(prev_x[i] - x[i]) > norma) {
				norma = abs(prev_x[i] - x[i]);
			}
		}
	} while (norma > eps);

	delete[] prev_x;
}




/**
 * Метод минимальных невязок
 */
void kiselevis::lab6()
{
	double eps = 1e-20;
	double *px = new double[N];
	for (int i = 0; i < N; i++) {
		px[i] = 0;
	}
	
	double *r = new double[N];

	while (true) {
		for (int i = 0; i < N; i++) {
			r[i] = b[i];

			for (int j = 0; j < N; j++) {
				r[i] -= (A[i][j] * px[j]);
			}
		}
		double tau = 0;
		double tmp = 0;
		for (int i = 0; i < N; i++) {
			double Ar = 0;

			for (int j = 0; j < N; j++) {
				Ar += (A[i][j] * r[j]);
			}
			tau += (Ar * r[i]);
			tmp += (Ar * Ar);
		}
		tau /= tmp;

		for (int i = 0; i < N; i++) {
			x[i] = px[i] + tau * r[i];
		}

		double error = 0;
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - px[i]) > error) {
				error = abs(x[i] - px[i]);
			}
		}

		if (error < eps) {
			break;
		}

		for (int i = 0; i < N; i++) {
			px[i] = x[i];
		}
	}


	delete[] px;
	delete[] r;
}



/**
 * Метод сопряженных градиентов
 */
void kiselevis::lab7()
{
	double eps = 1e-20;

	double* prevX = new double[N];
	double* prevR = new double[N];
	double* r = new double[N];
	double* z = new double[N];
	for (int i = 0; i < N; i++) {
		r[i] = b[i];
		z[i] = r[i];
	}

	while (true) {

		for (int i = 0; i < N; i++) {
			prevR[i] = r[i];
			prevX[i] = x[i];
		}

		double alpha = 0, denAlpha = 0;

		for (int i = 0; i < N; i++) {
			double Az = 0;
			for (int j = 0; j < N; j++) {
				Az += A[i][j] * z[j];
			}
			alpha += prevR[i] * prevR[i];
			denAlpha += Az * z[i];
		}
		alpha /= denAlpha;

		for (int i = 0; i < N; i++) {
			x[i] = prevX[i] + alpha * z[i];
		}

		double maxErr = abs(x[0] - prevX[0]);
		for (int i = 1; i < N; i++)
			if (abs(x[i] - prevX[i]) > maxErr)
				maxErr = abs(x[i] - prevX[i]);

		if (maxErr < eps)
			break;

		for (int i = 0; i < N; i++) {
			double Az = 0;

			for (int j = 0; j < N; j++) {
				Az += A[i][j] * z[j];
			}

			r[i] = prevR[i] - alpha * Az;
		}

		double beta = 0, denBeta = 0;
		for (int i = 0; i < N; i++) {
			beta += r[i] * r[i];
			denBeta += prevR[i] * prevR[i];
		}
		beta /= denBeta;

		for (int i = 0; i < N; i++) {
			z[i] = r[i] + beta * z[i];
		}
	}

	delete[] prevX;
	delete[] r;
	delete[] prevR;
	delete[] z;
}


void kiselevis::lab8()
{
	double eps = 1e-20;
	double **H = new double*[N];
	for (int i = 0; i < N; i++) 
		H[i] = new double[N];

	
	while(true){
		double n = 0;
		int i_max = 0, j_max = 1;
		for (int i = 0; i < N; i++)
			for (int j = i + 1; j < N; j++)
			{
				if (abs(A[i][j]) > abs(A[i_max][j_max]))
				{
					i_max = i;
					j_max = j;
				}

				n += A[i][j] * A[i][j];
			}

		if (sqrt(n) < eps) break;

		double fi = 0.5 * atan(2 * A[i_max][j_max] / (A[i_max][i_max] - A[j_max][j_max]));
		for (int i = 0; i < N; i++)
		{
			H[i][i_max] = A[i][i_max] * cos(fi) + A[i][j_max] * sin(fi);
			H[i][j_max] = A[i][j_max] * cos(fi) - A[i][i_max] * sin(fi);
		}

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				if (j != i_max && j != j_max) H[i][j] = A[i][j];

		for (int j = 0; j < N; j++)
		{
			A[i_max][j] = H[i_max][j] * cos(fi) + H[j_max][j] * sin(fi);
			A[j_max][j] = H[j_max][j] * cos(fi) - H[i_max][j] * sin(fi);
		}

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				if (i != i_max && i != j_max) A[i][j] = H[i][j];

	}

	for (int i = 0; i < N; i++) 
		x[i] = A[i][i];

	for (int i = 0; i < N; i++) 
		delete[] H[i];
	delete[] H;
}


void kiselevis::lab9()
{
	double eps = 1e-20;
	double lambda = 0;
	double plambda = 0;

	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	x[0] = 1;

	double *y = new double[N];

	while (true) {

		for (int i = 0; i < N; i++) {
			y[i] = 0;
			for (int j = 0; j < N; j++) {
				y[i] += (A[i][j] * x[j]);
			}
		}

		lambda = 0;
		for (int i = 0; i < N; i++) {
			lambda += (x[i] * y[i]);
		}

		if (abs(plambda - lambda) < eps) {
			break;
		}

		plambda = lambda;


		double norma_y = 0;
		for (int i = 0; i < N; i++) {
			norma_y += (y[i] * y[i]);
		}
		norma_y = sqrt(norma_y);

		for (int i = 0; i < N; i++) {
			x[i] = y[i] / norma_y;
		}
	}


	cout << lambda << endl;

	delete[] y;
}


std::string kiselevis::get_name()
{
  return "Kiselev Ilya Sergeevich";
}
