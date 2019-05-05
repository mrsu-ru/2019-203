#include "kiselevis.h"

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

}



/**
 * Метод минимальных невязок
 */
void kiselevis::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void kiselevis::lab7()
{

}


void kiselevis::lab8()
{

}


void kiselevis::lab9()
{

}


std::string kiselevis::get_name()
{
  return "Kiselev Ilya Sergeevich";
}
