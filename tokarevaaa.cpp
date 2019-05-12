#include "tokarevaaa.h"

/**
 * Введение в дисциплину
 */
void tokarevaaa::lab1()
{
 cout << "Hello, World!";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void tokarevaaa::lab2()
{
    double eps = 1.e-10;
	int *f = new int[N];
	for (int i = 0; i < N; i++) {
		f[i] = i;
	}

	for (int k = 0; k < N; k++) {
		int s = k;
		for (int i = k + 1; i < N; i++) {
			if (abs(A[s][k]) < abs(A[i][k])) {
				s = i;
			}
		}
		if (s != k) {
			swap(b[s], b[k]);
			swap(f[s], f[k]);
			for (int j = k; j < N; j++) {
				swap(A[s][j], A[k][j]);
			}
		}
		if (abs(A[k][k]) < eps) {
			cout << "Система не определена." << endl;
			return;
        }
		double t = A[k][k];
		b[k] /= t;
		for (int j = k; j < N; j++) {
			A[k][j] /= t;
		}
		for (int i = k + 1; i < N; i++) {
			t = A[i][k];

			b[i] -= b[k] * t;
			for (int j = k; j < N; j++) {
				A[i][j] -= A[k][j] * t;
			}
		}
	}
	// Обратный ход
	for (int k = N - 1; k > 0; k--) {
		for (int i = 0; i < k; i++) {
			b[i] -= A[i][k] * b[k];
		}
	}

	for (int i = 0; i < N; i++) {
		x[i] = b[f[i]];
	}

	delete[] f;

}



/**
 * Метод прогонки
 */
void tokarevaaa::lab3()
{
    double *alpha = new double[N];
	double *beta = new double[N];
	// Прямой ход
	alpha[0] = -A[0][1] / A[0][0];
	beta[0] = b[0] / A[0][0];
	for (int i = 1; i < N; i++) {
		alpha[i] = -A[i][i + 1] / (A[i][i] + A[i][i - 1] * alpha[i - 1]);
		beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / (A[i][i] + A[i][i - 1] * alpha[i - 1]);
	}
	// Обратный ход
	x[N - 1] = beta[N - 1];
	for (int i = N - 2; i >= 0; i--) {
		x[i] = alpha[i] * x[i + 1] + beta[i];
	}
	delete[] alpha;
	delete[] beta;
}



/**
 * Метод простых итераций
 */
void tokarevaaa::lab4()
{
	double eps = 1.e-14;
	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	const double t = 0.01;
	double *p_x = new double[N];
	double err = 0;
	do {
		for (int i = 0; i < N; i++) {
			p_x[i] = x[i];
		}
		for (int i = 0; i < N; i++) {
			double sum = 0;
			for (int j = 0; j < N; j++) {
				sum += (A[i][j] * p_x[j]);
			}
			x[i] += (t * (b[i] - sum));
		}
		err = 0;
		for (int i = 0; i < N; i++) {
			if (abs(p_x[i] - x[i]) > err) {
				err = abs(p_x[i] - x[i]);
			}
		}
	} while (err > eps);
	delete[] p_x;
}



/**
 * Метод Якоби 
 */
void tokarevaaa::lab5()
{
	double eps = 1.e-14;
	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	double *p_x = new double[N];
	double err = 0;
	do {
		for (int i = 0; i < N; i++) {
			p_x[i] = x[i];
		}
		for (int i = 0; i < N; i++) {
			double sum = 0;
			for (int j = 0; j < N; j++) {
				if (i != j) {
					sum += (A[i][j] * p_x[j]);
				}
			}

			x[i] = (b[i] - sum) / A[i][i];
		}
		err = 0;
		for (int i = 0; i < N; i++) {
			if ((p_x[i] - x[i]) > err) {
				err = abs(p_x[i] - x[i]);
			}
		}
	} while (err > eps);

	delete[] p_x;
}



/**
 * Метод минимальных невязок
 */
void tokarevaaa::lab6()
{
    double eps = 1.e-14;
    for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	double *p_x = new double[N];
	double *r = new double[N];
	double norma = 0;
	do {
		for (int i = 0; i < N; i++) {
			p_x[i] = x[i];
		}
		for (int i = 0; i < N; i++) {
			r[i] = b[i];
			for (int j = 0; j < N; j++) {
				r[i] -= (A[i][j] * x[j]);
			}
		}
		double tau = 0;
		double div = 0;
		for (int i = 0; i < N; i++) {
			double Ar = 0;
			for (int j = 0; j < N; j++) {
				Ar += (A[i][j] * r[j]);
			}
			tau += (Ar * r[i]);
			div += (Ar * Ar);
		}
		tau /= div;
		for (int i = 0; i < N; i++) {
			x[i] += (tau * r[i]);
		}
		norma = 0;
		for (int i = 0; i < N; i++) {
			if (abs(p_x[i] - x[i]) > norma) {
				norma = abs(p_x[i] - x[i]);
			}
		}
	} while (norma > eps);
	delete[] p_x;
	delete[] r;
}



/**
 * Метод сопряженных градиентов
 */
void tokarevaaa::lab7()
{

}


void tokarevaaa::lab8()
{

}


void tokarevaaa::lab9()
{

}


std::string tokarevaaa::get_name()
{
  return "Tokareva Alina Aleksandrovna";
}
