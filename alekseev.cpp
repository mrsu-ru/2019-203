#include "alekseev.h"

/**
 * Введение в дисциплину
 */
void alekseev::lab1()
{
	cout << "It's done" << endl;
	system("pause");
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void alekseev::lab2()
{
	double* temp = b;
	for (int i = 0, j = 0; (i < N - 1) && (j < N - 1); i++, j++)
	{
		int cur = i;
		double max = A[j][j];
		for (int num = j; num < N; num++)
			if (abs(max) < abs(A[num][j])) max = A[num][j];

		while (A[i][j] != max)
			i++;

		if (cur != i)
		{
			swap(temp[cur], temp[i]);
			swap(A[cur], A[i]);
			i = cur;
		}

		double t = 1 / A[i][j];
		temp[i] *= t;
		for (int cj = j; cj < N; cj++)
			A[i][cj] *= t;

		for (int ci = i + 1; ci < N; ci++)
		{
			temp[ci] += temp[i] * (-A[ci][j]);
			t = -A[ci][j];
			for (int cj = j; cj < N; cj++)
			{
				A[ci][cj] += A[i][cj] * t;
			}
		}
	}

	for (int i = N - 1, j = N - 1; (i > 0) && (j > 0); i--, j--)
	{
		int cur = i;
		temp[i] *= 1 / A[i][j];
		A[i][j] *= 1 / A[i][j];
		for (int ci = i - 1; ci >= 0; ci--)
		{
			temp[ci] += temp[i] * (-A[ci][j]);
			A[ci][j] = 0;
		}
	}

	x = temp;
}



/**
 * Метод прогонки
 */
void alekseev::lab3()
{
	double* alpha = new double[N - 1];
	double* gamma = new double[N];
	for (int i = 0; i < N; i++)
	{
		gamma[i] = A[i][i];
		if (i != 0) gamma[i] += A[i][i - 1] * alpha[i - 1];
		if (i != N - 1) alpha[i] = -A[i][i + 1] / gamma[i];
		x[i] = b[i] / gamma[i];
		if (i != 0) x[i] -= (A[i][i - 1] * x[i - 1]) / gamma[i];
	}

	for (int i = N - 2; i >= 0; i--)
		x[i] += alpha[i] * x[i + 1];

	delete[] alpha;
	delete[] gamma;
}



/**
 * Метод простых итераций
 */
void igaykiniv::lab4()
{
	double *new_x = new double[N], tau = 0.001, eps = 0.0000001;
	for (int i = 0; i < N; i++)
		x[i] = 0;

	do
	{
		for (int i = 0; i < N; i++)
		{
			double temp = 0;
			for (int j = 0; j < N; j++)
				temp += A[i][j] * x[j];

			new_x[i] = x[i] + tau * (b[i] - temp);
		}

		double maxdif = 0;
		for (int i = 0; i < N; i++)
		{
			if (fabs(x[i] - new_x[i]) > maxdif) maxdif = fabs(x[i] - new_x[i]);
			x[i] = new_x[i];
		}

		if (maxdif < eps) break;
	} while (true);

	delete[] new_x;
}



/**
 * Метод Якоби или Зейделя
 */
void alekseev::lab5()
{
	double *new_x = new double[N], eps = 0.0000001;
	bool condition;
	for (int i = 0; i < N; i++)
		x[i] = 1;

	do
	{
		condition = false;
		for (int i = 0; i < N; i++)
		{
			new_x[i] = b[i];
			for (int j = 0; j < N; j++)
			{
				if (i == j) continue;
				new_x[i] -= A[i][j] * x[j];
			}

			new_x[i] /= A[i][i];
			if (!condition) condition = (fabs(new_x[i] - x[i]) > eps);
			x[i] = new_x[i];
		}

	} while (condition);

	delete[] new_x;
}



/**
 * Метод минимальных невязок
 */
void alekseev::lab6()
{
	double *new_x = new double[N], *r = new double[N], eps = 0.0000001;
	for (int i = 0; i < N; i++)
		x[i] = 0;

	do
	{
		for (int i = 0; i < N; i++)
		{
			r[i] = b[i];
			for (int j = 0; j < N; j++)
				r[i] -= A[i][j] * x[j];
		}

		double tau, P = 0, Q = 0, t;
		for (int i = 0; i < N; i++)
		{
			t = 0;
			for (int j = 0; j < N; j++)
				t += A[i][j] * r[j];


			P += r[i] * t;
			Q += t * t;
		}

		tau = P / Q;
		for (int i = 0; i < N; i++)
			new_x[i] = x[i] + tau * r[i];

		double maxdif = 0;
		for (int i = 0; i < N; i++)
		{
			if (fabs(x[i] - new_x[i]) > maxdif) maxdif = fabs(x[i] - new_x[i]);
			x[i] = new_x[i];
		}

		if (maxdif < eps) break;
	} while (true);

	delete[] new_x;
	delete[] r;
}



/**
 * Метод сопряженных градиентов
 */
void alekseev::lab7()
{
	double *new_x = new double[N], *r = b, *new_r = new double[N], eps = 0.0000001;
	for (int i = 0; i < N; i++)
		x[i] = 0;

	double *z = new double[N];
	for (int i = 0; i < N; i++)
		z[i] = r[i];

	do
	{
		double tau1, tau2, P = 0, Q = 0, t;
		for (int i = 0; i < N; i++)
		{
			t = 0;
			for (int j = 0; j < N; j++)
				t += A[i][j] * r[j];


			P += r[i] * r[i];
			Q += t * z[i];
		}

		tau1 = P / Q;
		for (int i = 0; i < N; i++)
			new_x[i] = x[i] + tau1 * z[i];

		for (int i = 0; i < N; i++)
		{
			double temp = 0;
			for (int j = 0; j < N; j++)
				temp += A[i][j] * z[j];

			new_r[i] = r[i] - tau1 * temp;
		}

		double maxdif = 0;
		for (int i = 0; i < N; i++)
		{
			if (fabs(x[i] - new_x[i]) > maxdif) maxdif = fabs(x[i] - new_x[i]);
			x[i] = new_x[i];
		}

		if (maxdif < eps) break;

		P = 0; Q = 0;
		for (int i = 0; i < N; i++)
		{
			P += new_r[i] * new_r[i];
			Q += r[i] * r[i];
		}

		tau2 = P / Q;
		for (int i = 0; i < N; i++)
			z[i] = new_r[i] + tau2 * z[i];

		for (int i = 0; i < N; i++)
			r[i] = new_r[i];

	} while (true);

	delete[] new_x;
	delete[] new_r;
	delete[] r;
	delete[] z;
}


void alekseev::lab8()
{
	double **H = new double*[N], eps = 1.e-10;
	for (int i = 0; i < N; i++) H[i] = new double[N];

	do
	{
		double n = 0;
		int i_max = 0, j_max = 1;
		for (int i = 0; i < N; i++)
			for (int j = i + 1; j < N; j++)
			{
				if (fabs(A[i][j]) > fabs(A[i_max][j_max]))
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

	} while (true);

	for (int i = 0; i < N; i++) x[i] = A[i][i];

	for (int i = 0; i < N; i++) delete[] H[i];
	delete[] H;
}


void alekseev::lab9()
{
	double *new_x = new double[N], lam = 0, eps = 1.e-10;
	for (int i = 1; i < N; i++) x[i] = 0;

	x[0] = 1;
	do
	{
		double new_lam = 0;
		for (int i = 0; i < N; i++)
		{
			new_x[i] = 0;

			for (int j = 0; j < N; j++)
				new_x[i] += A[i][j] * x[j];

			new_lam += x[i] * new_x[i];
		}

		if (fabs(new_lam - lam) < eps) break;

		lam = new_lam;
		double n = 0;
		for (int i = 0; i < N; i++) n += new_x[i] * new_x[i];

		n = sqrt(n);
		for (int i = 0; i < N; i++) x[i] = new_x[i] / n;

	} while (true);

	x[0] = lam;
	delete[] new_x;
}


std::string alekseev::get_name()
{
  return "AlekseevAN";
}
