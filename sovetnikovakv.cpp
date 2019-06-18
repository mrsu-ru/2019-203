#include "sovetnikovakv.h"

const double EPSILON = 1.e-15;

/**
 * Введение в дисциплину
 */
void sovetnikovakv::lab1()
{
	std::cout << "Hello world!" << std::endl;
}



/**
 * Метод Гаусса с выбором главного элемента
 */
void sovetnikovakv::lab2()
{
	// Запоминаем позиции x при перестановке строк
	int *p = new int[N];
	for (int i = 0; i < N; i++) {
		p[i] = i;
	}

	// Прямой ход
	for (int k = 0; k < N; k++) {
		// Поиск ведущего элемента
		int s = k;
		for (int i = k + 1; i < N; i++) {
			if (std::abs(A[s][k]) < std::abs(A[i][k])) {
				s = i;
			}
		}

		// Меняем строки
		if (s != k) {
			std::swap(b[s], b[k]);
			std::swap(p[s], p[k]);
			for (int j = k; j < N; j++) {
				std::swap(A[s][j], A[k][j]);
			}
		}

		// Система не определена
		if (std::abs(A[k][k]) < EPSILON) {
			std::cout << "Система не определена." << std::endl;
			return;
		}

		// Нормируем k-ю строчку
		double tmp = A[k][k];

		b[k] /= tmp;
		for (int j = k; j < N; j++) {
			A[k][j] /= tmp;
		}

		// Приводим матрицу к верхне треугольному виду
		for (int i = k + 1; i < N; i++) {
			tmp = A[i][k];

			b[i] -= (b[k] * tmp);
			for (int j = k; j < N; j++) {
				A[i][j] -= (A[k][j] * tmp);
			}
		}
	}

	// Обратный ход
	for (int k = N - 1; k > 0; k--) {
		for (int i = 0; i < k; i++) {
			b[i] -= (A[i][k] * b[k]);
		}
	}

	// Записываем результат в вектор x
	for (int i = 0; i < N; i++) {
		x[i] = b[p[i]];
	}

	delete[] p;
}



/**
 * Метод прогонки
 */
void sovetnikovakv::lab3()
{
	double *alph = new double[N];
	double *beta = new double[N];

	// Прямой ход
	alph[0] = -A[0][1] / A[0][0];
	beta[0] = b[0] / A[0][0];
	for (int i = 1; i < N; i++) {
		alph[i] = -A[i][i + 1] / (A[i][i - 1] * alph[i - 1] + A[i][i]);
		beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / (A[i][i - 1] * alph[i - 1] + A[i][i]);
	}

	// Обратный ход
	x[N - 1] = beta[N - 1];
	for (int i = N - 2; i >= 0; i--) {
		x[i] = alph[i] * x[i + 1] + beta[i];
	}

	delete[] alph;
	delete[] beta;
}



/**
 * Метод простых итераций
 */
void sovetnikovakv::lab4()
{
	// Начальное приближение
	for (int i = 0; i < N; i++) {
		x[i] = 0.;
	}

	// Итерационный параметр
	const double tau = 0.01;

	// Вектор решений на предыдущей итерации
	double *prev_x = new double[N];

	double norma = 0.;
	do {
		// Запишем предыдущий вектор
		for (int i = 0; i < N; i++) {
			prev_x[i] = x[i];
		}

		// Найдем новый вектор решений на текущей итерации
		for (int i = 0; i < N; i++) {
			// Подставим текущее решение
			double result = b[i];
			for (int j = 0; j < N; j++) {
				result -= (A[i][j] * prev_x[j]);
			}

			// Приближаем вектор решение к ответу
			x[i] += (tau * result);
		}

		// Находим погрешность
		norma = 0.;
		for (int i = 0; i < N; i++) {
			if (std::abs(prev_x[i] - x[i]) > norma) {
				norma = std::abs(prev_x[i] - x[i]);
			}
		}
	} while (norma > EPSILON);

	delete[] prev_x;
}



/**
 * Метод Якоби или Зейделя
 */
void sovetnikovakv::lab5()
{
	// Метод Якоби

	// Начальное приближение
	for (int i = 0; i < N; i++) {
		x[i] = 0.;
	}

	// Вектор решений на предыдущей итерации
	double *prev_x = new double[N];

	double norma = 0.;
	do {
		// Запишем передыдущий вектор
		for (int i = 0; i < N; i++) {
			prev_x[i] = x[i];
		}

		// Найдем новый вектор решений на текущей итерации
		for (int i = 0; i < N; i++) {
			// Подставим текущее решение
			double result = b[i];
			for (int j = 0; j < N; j++) {
				if (i != j) {
					result -= (A[i][j] * prev_x[j]);
				}
			}

			x[i] = result / A[i][i];
		}

		// Находим погрешность
		norma = 0.;
		for (int i = 0; i < N; i++) {
			if (std::abs(prev_x[i] - x[i]) > norma) {
				norma = std::abs(prev_x[i] - x[i]);
			}
		}
	} while (norma > EPSILON);

	delete[] prev_x;
}



/**
 * Метод минимальных невязок
 */
void sovetnikovakv::lab6()
{
	// Начальное приближение
	for (int i = 0; i < N; i++) {
		x[i] = 0.;
	}

	// Вектор решений на предыдущей итерации
	double *prev_x = new double[N];

	// Вектор невязок
	double *r = new double[N];

	double norma = 0.;
	do {
		// Запишем передыдущий вектор
		for (int i = 0; i < N; i++) {
			prev_x[i] = x[i];
		}

		// Найдем новый вектор невязок
		for (int i = 0; i < N; i++) {
			r[i] = b[i];
			for (int j = 0; j < N; j++) {
				r[i] -= (A[i][j] * x[j]);
			}
		}

		// Найдем итерационный параметр
		double tau = 0.;
		double div = 0.;
		for (int i = 0; i < N; i++) {
			double Ar = 0.;
			for (int j = 0; j < N; j++) {
				Ar += (A[i][j] * r[j]);
			}
			tau += (Ar * r[i]);
			div += (Ar * Ar);
		}
		tau /= div;

		// Найдем новый вектор решений на текущей итерации
		for (int i = 0; i < N; i++) {
			x[i] += (tau * r[i]);
		}

		// Находим погрешность
		norma = 0.;
		for (int i = 0; i < N; i++) {
			if (std::abs(prev_x[i] - x[i]) > norma) {
				norma = std::abs(prev_x[i] - x[i]);
			}
		}
	} while (norma > EPSILON);

	delete[] prev_x;
	delete[] r;
}



/**
 * Метод сопряженных градиентов
 */
void sovetnikovakv::lab7()
{
	// Начальное приближение
	for (int i = 0; i < N; i++) {
		x[i] = 0.;
	}

	// Вектор решений на предыдущей итерации
	double *prev_x = new double[N];

	// Вектор невязок на предыдущей итерации
	double *prev_r = new double[N];

	// Вектор невязок
	double *r = new double[N];

	// Вектор приближения
	double *z = new double[N];
	
	// Найдем начальную невязку и приближение
	for (int i = 0; i < N; i++) {
		r[i] = b[i];
		for (int j = 0; j < N; j++) {
			r[i] -= (A[i][j] * x[j]);
		}
		z[i] = r[i];
	}

	double temp = 0.;
	while (true) {
		// Найдем 1-й итерационный параметр
		temp = 0.;
		double alph = 0.;
		for (int i = 0; i < N; i++) {
			double Az = 0.;
			for (int j = 0; j < N; j++) {
				Az += (A[i][j] * z[j]);
			}
			alph += (r[i] * r[i]);
			temp += (Az * z[i]);
		}
		alph /= temp;

		// Запишем предыдущее решение
		for (int i = 0; i < N; i++) {
			prev_x[i] = x[i];
		}

		// Найдем новый вектор решений
		for (int i = 0; i < N; i++) {
			x[i] += (alph * z[i]);
		}

		// Находим погрешность
		double norma = 0.;
		for (int i = 0; i < N; i++) {
			if (std::abs(prev_x[i] - x[i]) > norma) {
				norma = std::abs(prev_x[i] - x[i]);
			}
		}

		// Достижение необходимой точности
		if (norma < EPSILON) {
			break;
		}

		// Запишем предыдущую невязку
		for (int i = 0; i < N; i++) {
			prev_r[i] = r[i];
		}

		// Найдем новую невязку для последующей итерации
		for (int i = 0; i < N; i++) {
			double Az = 0.;
			for (int j = 0; j < N; j++) {
				Az += (A[i][j] * z[j]);
			}
			r[i] -= (alph * Az);
		}

		// Найдем 2-й итерационный параметр
		temp = 0.;
		double beta = 0.;
		for (int i = 0; i < N; i++) {
			beta += (r[i] * r[i]) ;
			temp += (prev_r[i] * prev_r[i]);
		}
		beta /= temp;

		// Расчитаем приближение для последующей итерации
		for (int i = 0; i < N; i++) {
			z[i] = r[i] + beta * z[i];
		}
	}

	delete[] prev_x;
 	delete[] prev_r;
	delete[] r;
	delete[] z;
}



/*
* Метод вращения для нахождения собственных значений матрицы
*/
void sovetnikovakv::lab8()
{
	// Матрица для промежуточных вычислений
	double **T = new double*[N];
	for (int i = 0; i < N; i++) {
		T[i] = new double[N];
	}

	while (true) {
		// Находим наибольший внедиагональный элемент
		// Вычисляем норму внедиагональных элементов
		int im = 0;
		int jm = 0;
		double norma = 0.;
		for (int i = 0; i < N - 1; i++) {
			for (int j = i + 1; j < N; j++) {
				norma += (A[i][j] * A[i][j]);
				if (std::abs(A[i][j]) > std::abs(A[im][jm])) {
					im = i;
					jm = j;
				}
			}
		}

		// Матрица является диагональной
		if (std::sqrt(norma) < EPSILON) {
			break;
		}

		// Находим значения переменных матрицы поворота
		double f = .5 * std::atan(2. * A[im][jm] / (A[im][im] - A[jm][jm]));
		double c = std::cos(f);
		double s = std::sin(f);

		// Копируем матрицу A
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				T[i][j] = A[i][j];
			}
		}

		// Осуществляем поворот A * T
		for (int k = 0; k < N; k++) {
			T[k][im] = A[k][im] * c + A[k][jm] * s;
			T[k][jm] = A[k][jm] * c - A[k][im] * s;
		}

		// Осуществялем поворот T' * (A * T)
		for (int k = 0; k < N; k++) {
			A[im][k] = T[im][k] * c + T[jm][k] * s;
			A[jm][k] = T[jm][k] * c - T[im][k] * s;
		}

		// Записываем все расчеты от первого поворота в матрицу A
		for (int i = 0; i < N; i++) {
			if ((i != im) && (i != jm)) {
				for (int j = 0; j < N; j++) {
					A[i][j] = T[i][j];
				}
			}
		}
	}

	// Вектор собственных значений матрицы
	for (int i = 0; i < N; i++) {
		x[i] = A[i][i];
	}

	for (int i = 0; i < N; i++) {
		delete[] T[i];
	}
	delete[] T;
}


/*
* Нахождение наибольшего по модулю собственного значения матрицы
*/
void sovetnikovakv::lab9()
{
	// Начальное приближение
	for (int i = 0; i < N; i++) {
		x[i] = 0.;
	}
	x[0] = 1.;

	// Вектор последующего приближения
	double *z = new double[N];

	// Наибольшее собственнное значение матрицы
	double result = 0.;

	while (true) {
		// Находим новое приближение, "увеличиваем степень"
		for (int i = 0; i < N; i++) {
			z[i] = 0.;
			for (int j = 0; j < N; j++) {
				z[i] += (A[i][j] * x[j]);
			}
		}

		// Находим собственное значение
		double num = 0.;
		double den = 0.;
		for (int i = 0; i < N; i++) {
			num += (z[i] * x[i]);
			den += (x[i] * x[i]);
		}
		double temp = result;
		result = num / den;

		// Достижение необходимой точности
		if (std::abs(temp - result) < EPSILON) {
			break;
		}

		// Запишем предыдущее приближение
		double norma = std::sqrt(den);
		for (int i = 0; i < N; i++) {
			x[i] = z[i] / norma;
		}
	}

	std::cout << "Наибольшее по модулю собственное значение матрицы: " << result << std::endl;

	delete[] z;
}


std::string sovetnikovakv::get_name()
{
  return "Советникова Карина Владиславовна";
}
