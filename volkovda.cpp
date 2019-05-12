#include "volkovda.h"

const double EPSILON = 1.e-15;



/**
 * Введение в дисциплину
 */
void volkovda::lab1()
{
	std::cout << "Hello, World!!!\n";
}



/**
 * Метод Гаусса с выбором главного элемента
 */
void volkovda::lab2()
{
	// Отслеживаем позиции вектора x
	int *p = new int[N];
	for (int i = 0; i < N; i++) {
		p[i] = i;
	}

	for (int k = 0; k < N; k++) {

		// Выбор ведущего элемента
		int d = k;
		for (int i = k + 1; i < N; i++) {
			if (std::abs(A[d][k]) < std::abs(A[i][k])) {
				d = i;
			}
		}

		// Меняем строки
		if (d != k) {
			for (int j = k; j < N; j++) {
				std::swap(A[d][j], A[k][j]);
			}
			std::swap(p[d], p[k]);
			std::swap(b[d], b[k]);
		}

		// Система не определена
		if (std::abs(A[k][k]) < EPSILON) {
			std::cout << "System of linear equations has no solution.\n";
			return;
		}

		// Нормируем строку k
		b[k] /= A[k][k];
		for (int j = N - 1; j >= k; j--) {
			A[k][j] /= A[k][k];
		}

		// Приводим матрицу к диагональному виду на позиции k
		for (int i = 0; i < N; i++) {
			if (i != k) {
				b[i] -= (A[i][k] * b[k]);
				for (int j = N - 1; j >= k; j--) {
					A[i][j] -= (A[i][k] * A[k][j]);
				}
			}
		}
	}

	// Записываем результат в вектор x с учетом измененных позиций
	for (int i = 0; i < N; i++) {
		x[i] = b[p[i]];
	}

	delete[] p;
}



/**
 * Метод прогонки
 */
void volkovda::lab3()
{
	double *v = new double[N];
	double *w = new double[N];

	// Прямой ход
	v[0] = -A[0][1] / A[0][0];
	w[0] = b[0] / A[0][0];
	for (int i = 1; i < N; i++) {
		v[i] = -A[i][i + 1] / (A[i][i] + A[i][i - 1] * v[i - 1]);
		w[i] = (b[i] - A[i][i - 1] * w[i - 1]) / (A[i][i] + A[i][i - 1] * v[i - 1]);
	}

	// Обратный ход
	x[N - 1] = w[N - 1];
	for (int i = N - 2; i >= 0; i--) {
		x[i] = v[i] * x[i + 1] + w[i];
	}

	delete[] v;
	delete[] w;
}



/**
 * Метод простых итераций
 */
void volkovda::lab4()
{
	// Вектор значений x на предыдущий итерации
	double *px = new double[N];
	for (int i = 0; i < N; i++) {
		px[i] = 0.0;
	}

	// Итерационный параметр
	const double tau = 0.01;

	int iteration = 0;
	while (true) {
		iteration++;

		// Посчитаем значения неизвестных на текущей итерации
		for (int i = 0; i < N; i++) {
			double var = 0.0;

			// Вычислим сумму
			for (int j = 0; j < N; j++) {
				var += (A[i][j] * px[j]);
			}

			x[i] = px[i] + tau * (b[i] - var);
		}

		// Посчитаем максимальную по модулю погрешность
		double error = 0.0;
		for (int i = 0; i < N; i++) {
			if (std::abs(x[i] - px[i]) > error) {
				error = std::abs(x[i] - px[i]);
			}
		}

		// При достижении необходимой точности завершаем процесс
		if (error < EPSILON) {
			break;
		}

		// Текущее зачение итерации представим как предыдущее
		for (int i = 0; i < N; i++) {
			px[i] = x[i];
		}
	}

	std::cout << "Number of iterations : " << iteration << '\n';

	delete[] px;
}



/**
 * Метод Якоби или Зейделя
 */
void volkovda::lab5()
{
	// Метод Зейделя (модификация метода простых итераций)

	// Вектор значений x на предыдущий итерации
	double *px = new double[N];
	for (int i = 0; i < N; i++) {
		px[i] = 0.0;
	}

	int iteration = 0;
	while (true) {
		iteration++;

		// Посчитаем значения неизвестных на текущей итерации
		for (int i = 0; i < N; i++) {
			double var = 0.0;

			// Вычисляем сумму k+1 -го шага
			for (int j = 0; j < i; j++) {
				var += (A[i][j] * x[j]);
			}
			// Вычисялем сумму k -го шага
			for (int j = i + 1; j < N; j++) {
				var += (A[i][j] * px[j]);
			}

			x[i] = (b[i] - var) / A[i][i];
		}

		// Посчитаем максимальную по модулю погрешность
		double error = 0.0;
		for (int i = 0; i < N; i++) {
			if (std::abs(x[i] - px[i]) > error) {
				error = std::abs(x[i] - px[i]);
			}
		}

		// При достижении необходимой точности завершаем процесс
		if (error < EPSILON) {
			break;
		}

		// Текущее зачение итерации представим как предыдущее
		for (int i = 0; i < N; i++) {
			px[i] = x[i];
		}
	}

	std::cout << "Number of iterations : " << iteration << '\n';

	delete[] px;
}



/**
 * Метод минимальных невязок
 */
void volkovda::lab6()
{
	// Вектор значений x на предыдущий итерации
	double *px = new double[N];
	for (int i = 0; i < N; i++) {
		px[i] = 0.0;
	}
	
	// Вектор невязок
	double *r = new double[N];

	int iteration = 0;
	while (true) {
		iteration++;

		// Рассчитываем вектор невязок
		for (int i = 0; i < N; i++) {
			r[i] = b[i];

			for (int j = 0; j < N; j++) {
				r[i] -= (A[i][j] * px[j]);
			}
		}

		// Рассчитываем итерационный параметр
		double tau = 0.0;
		double tmp = 0.0;
		for (int i = 0; i < N; i++) {
			double Ar = 0.0;

			for (int j = 0; j < N; j++) {
				Ar += (A[i][j] * r[j]);
			}

			tau += (Ar * r[i]);
			tmp += (Ar * Ar);
		}
		tau /= tmp;

		// Рассчитывается новое приближение к вектору x
		for (int i = 0; i < N; i++) {
			x[i] = px[i] + tau * r[i];
		}

		// Посчитаем максимальную по модулю погрешность
		double error = 0.0;
		for (int i = 0; i < N; i++) {
			if (std::abs(x[i] - px[i]) > error) {
				error = std::abs(x[i] - px[i]);
			}
		}

		// При достижении необходимой точности завершаем процесс
		if (error < EPSILON) {
			break;
		}

		// Текущее зачение итерации представим как предыдущее
		for (int i = 0; i < N; i++) {
			px[i] = x[i];
		}
	}

	std::cout << "Number of iterations : " << iteration << '\n';

	delete[] px;
	delete[] r;
}



/**
 * Метод сопряженных градиентов
 */
void volkovda::lab7()
{
	// Вектор значений x на предыдущий итерации
	double *px = new double[N];
	for (int i = 0; i < N; i++) {
		px[i] = 0.0;
	}

	// Вектор невязок на предыдущий итерации
	double *pr = new double[N];

	// Вектор невязок
	double *r = new double[N];

	// Вектор приближения
	double *d = new double[N];

	// Расчитаем невязку и приближение на 0 -вой итерации
	for (int i = 0; i < N; i++) {
		pr[i] = b[i];

		for (int j = 0; j < N; j++) {
			pr[i] -= (A[i][j] * px[j]);
		}

		d[i] = pr[i];
	}

	int iteration = 0;
	while (true) {
		iteration++;

		double temp;

		// Вычисли итерационный параметр alph
		temp = 0.0;
		double alph = 0.0;
		for (int i = 0; i < N; i++) {
			double Ad = 0.0;

			for (int j = 0; j < N; j++) {
				Ad += (A[i][j] * d[j]);
			}

			alph += (pr[i] * pr[i]);
			temp += (Ad * d[i]);
		}
		alph /= temp;

		// Рассчитывается новое приближение к вектору x
		for (int i = 0; i < N; i++) {
			x[i] = px[i] + alph * d[i];
		}

		// Посчитаем максимальную по модулю погрешность
		double error = 0.0;
		for (int i = 0; i < N; i++) {
			if (std::abs(x[i] - px[i]) > error) {
				error = std::abs(x[i] - px[i]);
			}
		}

		// При достижении необходимой точности завершаем процесс
		if (error < EPSILON) {
			break;
		}

		// Расчитаем новый вектор невязок для последующей итерации
		for (int i = 0; i < N; i++) {
			double Ad = 0.0;

			for (int j = 0; j < N; j++) {
				Ad += (A[i][j] * d[j]);
			}

			r[i] = pr[i] - alph * Ad;
		}

		// Вычисли итерационный параметр beta для последующего приближения
		temp = 0.0;
		double beta = 0.0;
		for (int i = 0; i < N; i++) {
			beta += (r[i] * r[i]);
			temp += (pr[i] * pr[i]);
		}
		beta /= temp;

		// Приближение i+1 -ой итерации
		for (int i = 0; i < N; i++) {
			d[i] = r[i] + beta * d[i];
		}

		// Текущее зачение итерации представим как предыдущее
		for (int i = 0; i < N; i++) {
			px[i] = x[i];
			pr[i] = r[i];
		}
	}

	std::cout << "Number of iterations : " << iteration << '\n';

	delete[] px;
	delete[] pr;
	delete[] r;
	delete[] d;
}



/*
* Метод вращения для нахождения собственных значений матрицы
*/
void volkovda::lab8()
{
	// Промежуточная матрица для вычислений
	double **B = new double*[N];
	for (int i = 0; i < N; i++) {
		B[i] = new double[N];
	}

	int iteration = 0;
	while (true) {
		iteration++;

		// Вычисляем норму и находим наибольший не диагональный элемент
		int iMax = 0;
		int jMax = 0;
		double norma = 0.0;
		for (int i = 0; i < N - 1; i++) {
			for (int j = i + 1; j < N; j++) {
				if (std::abs(A[i][j]) > std::abs(A[iMax][jMax])) {
					iMax = i;
					jMax = j;
				}
				norma += (A[i][j] * A[i][j]);
			}
		}

		// Матрица являеться диагональной
		if (std::sqrt(norma) < EPSILON) {
			break;
		}

		// Вычисляем значения переменных матрицы вражения
		double fi = 0.5 * std::atan(2.0 * A[iMax][jMax] / (A[iMax][iMax] - A[jMax][jMax]));
		double cos = std::cos(fi);
		double sin = std::sin(fi);

		// Находим новое приближение как A = B' * A * B

		for (int k = 0; k < N; k++) {
			B[k][iMax] = A[k][iMax] * cos + A[k][jMax] * sin;
			B[k][jMax] = A[k][jMax] * cos - A[k][iMax] * sin;
		}

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if ((j != iMax) && (j != jMax)) {
					B[i][j] = A[i][j];
				}
			}
		}

		for (int k = 0; k < N; k++) {
			A[iMax][k] = B[iMax][k] * cos + B[jMax][k] * sin;
			A[jMax][k] = B[jMax][k] * cos - B[iMax][k] * sin;
		}

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if ((i != iMax) && (i != jMax)) {
					A[i][j] = B[i][j];
				}
			}
		}
	}

	std::cout << "Number of iterations : " << iteration << '\n';

	// Запишем в вектор решений вектор собственных значений матрицы
	std::cout << "The vector X is the vector of eigenvalues of the matrix.\n";
	for (int i = 0; i < N; i++) {
		x[i] = A[i][i];
	}

	for (int i = 0; i < N; i++) {
		delete[] B[i];
	}
	delete[] B;

	//***********************************************************************************

	// Нахождение наибольшего по модулю собственного значения матрицы
	// см. lab9
	
	double max = 0.0;
	for (int i = 0; i < N; i++) {
		if (max < std::abs(x[i])) {
			max = x[i];
		}
	}
	std::cout << "The largest eigenvalue of the matrix: " << max << '\n';
}



/*
* Нахождение наибольшего по модулю собственного значения матрицы
*/
void volkovda::lab9()
{
	// Степенной метод без сдвига

	// Искомое наибольшее по модулю собственное значение матрицы
	double lambda = 0.0;
	double plambda = 0.0;

	// Изначельное приближение, полагаем ||x|| = 1.0
	for (int i = 0; i < N; i++) {
		x[i] = 0.0;
	}
	x[0] = 1.0;

	// Вектор последующего приближения
	double *y = new double[N];

	int iteration = 0;
	while (true) {
		iteration++;

		// Вычисляем i -ое приближение
		for (int i = 0; i < N; i++) {
			y[i] = 0.0;
			for (int j = 0; j < N; j++) {
				y[i] += (A[i][j] * x[j]);
			}
		}

		// Вычисляем новое значение
		lambda = 0.0;
		for (int i = 0; i < N; i++) {
			lambda += (x[i] * y[i]);
		}

		// Достаточно близки к ответу
		if (std::abs(plambda - lambda) < EPSILON) {
			break;
		}

		plambda = lambda;

		// "Увеличение" степени

		double norma_y = 0.0;
		for (int i = 0; i < N; i++) {
			norma_y += (y[i] * y[i]);
		}
		norma_y = std::sqrt(norma_y);

		for (int i = 0; i < N; i++) {
			x[i] = y[i] / norma_y;
		}
	}

	std::cout << "Number of iterations : " << iteration << '\n';

	std::cout << "The largest eigenvalue of the matrix: " << lambda << '\n';

	delete[] y;
}


std::string volkovda::get_name()
{
  return "Volkov Dmitry Alekseevich";
}
