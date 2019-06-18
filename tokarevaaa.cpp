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
	double eps = 1.e-14;
	double *p_x = new double[N];
	double *p_r = new double[N];
	double *r = new double[N];
	double *z = new double[N];
    for (int i = 0; i < N; i++) {
		x[i] = 0.;
	}
	for (int i = 0; i < N; i++) {
		r[i] = b[i];
		for (int j = 0; j < N; j++) {
			r[i] -= (A[i][j] * x[j]);
		}
		z[i] = r[i];
	}

	double t = 0.;
	while (true) {
		t = 0.;
		double alph = 0.;
		for (int i = 0; i < N; i++) {
			double Az = 0.;
			for (int j = 0; j < N; j++) {
				Az += (A[i][j] * z[j]);
			}
			alph += (r[i] * r[i]);
			t += (Az * z[i]);
		}
		alph /= t;

		for (int i = 0; i < N; i++) {
			p_x[i] = x[i];
		}

		for (int i = 0; i < N; i++) {
			x[i] += (alph * z[i]);
		}

		// Находим погрешность
		double norma = 0.;
		for (int i = 0; i < N; i++) {
			if (abs(p_x[i] - x[i]) > norma) {
				norma = abs(p_x[i] - x[i]);
			}
		}
		if (norma < eps) {
			break;
		}

		for (int i = 0; i < N; i++) {
			p_r[i] = r[i];
		}

		for (int i = 0; i < N; i++) {
			double Az = 0.;
			for (int j = 0; j < N; j++) {
				Az += (A[i][j] * z[j]);
			}
			r[i] -= (alph * Az);
		}

		t = 0.;
		double beta = 0.;
		for (int i = 0; i < N; i++) {
			beta += (r[i] * r[i]) ;
			t += (p_r[i] * p_r[i]);
		}
		beta /= t;

		for (int i = 0; i < N; i++) {
			z[i] = r[i] + beta * z[i];
		}
	}
	delete[] p_x;
 	delete[] p_r;
	delete[] r;
    delete[] z;

}

/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void tokarevaaa::lab8()
{
	double eps = 1.e-14;
	double **T = new double*[N];
	for (int i = 0; i < N; i++) {
		T[i] = new double[N];
	}

	while (true) {
		int im = 0;
		int jm = 0;
		double norma = 0.;
		for (int i = 0; i < N - 1; i++) {
			for (int j = i + 1; j < N; j++) {
				norma += (A[i][j] * A[i][j]);
				if (abs(A[i][j]) > abs(A[im][jm])) {
					im = i;
					jm = j;
				}
			}
		}

		if (sqrt(norma) < eps) {
			break;
		}

		double f = .5 * atan(2. * A[im][jm] / (A[im][im] - A[jm][jm]));
		double c = cos(f);
		double s = sin(f);

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				T[i][j] = A[i][j];
			}
		}

		for (int k = 0; k < N; k++) {
			T[k][im] = A[k][im] * c + A[k][jm] * s;
			T[k][jm] = A[k][jm] * c - A[k][im] * s;
		}

		for (int k = 0; k < N; k++) {
			A[im][k] = T[im][k] * c + T[jm][k] * s;
			A[jm][k] = T[jm][k] * c - T[im][k] * s;
		}

		for (int i = 0; i < N; i++) {
			if ((i != im) && (i != jm)) {
				for (int j = 0; j < N; j++) {
					A[i][j] = T[i][j];
				}
			}
		}
	}

	for (int i = 0; i < N; i++) {
		x[i] = A[i][i];
	}

	for (int i = 0; i < N; i++) {
		delete[] T[i];
	}
    delete[] T;

}

/**
 * Нахождение наибольшего по модулю собственного значения матрицы
 */
void tokarevaaa::lab9()
{
	double eps = 1.e-14;
	for (int i = 0; i < N; i++) {
		x[i] = 0.;
	}
	x[0] = 1.;

	double *z = new double[N];
	double res = 0.;

	while (true) {
		for (int i = 0; i < N; i++) {
			z[i] = 0.;
			for (int j = 0; j < N; j++) {
				z[i] += (A[i][j] * x[j]);
			}
		}

		double num = 0.;
		double den = 0.;
		for (int i = 0; i < N; i++) {
			num += (z[i] * x[i]);
			den += (x[i] * x[i]);
		}
		double t = res;
		res = num / den;


		if (abs(t - res) < eps) {
			break;
		}

		double norma = sqrt(den);
		for (int i = 0; i < N; i++) {
			x[i] = z[i] / norma;
		}
	}
	cout << res << endl;
    delete[] z;

}


std::string tokarevaaa::get_name()
{
  return "Tokareva Alina Aleksandrovna";
}
