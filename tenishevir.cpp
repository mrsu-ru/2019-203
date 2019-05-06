#include "tenishevir.h"

/**
 * Введение в дисциплину
 */
void tenishevir::lab1()
{
  cout << "Hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void tenishevir::lab2()
{
    double p;
	int maxn;

    for (int k=0; k<N-1; k++){
        maxn = k;
        for (int i=k+1; i<N; i++)
			if(abs(A[i][k]) > abs(A[maxn][k])) maxn = i; ///Выбор главного элемента
        std::swap(A[maxn], A[k]); ///Меняем строки местами
        std::swap(b[maxn], b[k]);

        for (int i=k+1; i<N; i++){
            p = A[i][k]/A[k][k];
            for (int j=k; j<N; j++)
                A[i][j] -= p*A[k][j];
            b[i] -= p*b[k];
        }
    }

    for(int i = 0; i<N; i++){
        x[i]=b[i];
    }

    for (int i=N-1; i>=0; i--){
        for (int j=i+1;j<N;j++)
            x[i]-=A[i][j]*x[j];
        x[i] /= A[i][i];
    }
}



/**
 * Метод прогонки
 */
void tenishevir::lab3()
{
    double *P = new double [N]; ///Коэффициенты "альфа"
	double *Q = new double [N]; ///Коэффициенты "бетта"

	P[0] = -A[0][1]/A[0][0];
	Q[0] = b[0]/A[0][0];

	for(int i=1; i<N; i++) ///Определяем прогоночные коэффициенты
	{
		P[i] = A[i][i+1]/(-A[i][i] - A[i][i-1]*P[i-1]);
		Q[i] = (-b[i] + A[i][i-1]*Q[i-1])/(-A[i][i] - A[i][i-1]*P[i-1]);
	}

	x[N-1] = Q[N-1];
	for(int i=N-2; i>=0; i--) ///Определяем решение
		x[i] = P[i]*x[i+1] + Q[i];

	delete [] P;
	delete [] Q;
}



/**
 * Метод простых итераций
 */
void tenishevir::lab4()
{
    double eps = 1e-15;
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
void tenishevir::lab5()
{
  double eps = 1e-15;

	double* prevX = new double[N];

	while (true) {

		for (int i = 0; i < N; i++)
			prevX[i] = x[i];

		for (int i = 0; i < N; i++) {
			double sum = 0;
			for (int j = 0; j < i; j++)
				sum += A[i][j] * x[j];
			for (int j = i + 1; j < N; j++)
				sum += A[i][j] * prevX[j];
			x[i] = (b[i] - sum) / A[i][i];
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
 * Метод минимальных невязок
 */
void tenishevir::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void tenishevir::lab7()
{

}


void tenishevir::lab8()
{

}


void tenishevir::lab9()
{

}


std::string tenishevir::get_name()
{
  return "Tenishev Ilya";
}
