#include "borisovaem.h"

/**
 * Введение в дисциплину
 */
void borisovaem::lab1()
{
cout << "Hello, world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void borisovaem::lab2()
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
void borisovaem::lab3()
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
void borisovaem::lab4()
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
void borisovaem::lab5()
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
void borisovaem::lab6()
{
    double eps = 1e-15;

	double* prevX = new double[N];
	double* r = new double[N];

	while (true) {

		for (int i = 0; i < N; i++)
			prevX[i] = x[i];

		for (int i = 0; i < N; i++) {
			r[i] = b[i];

			for (int j = 0; j < N; j++) {
				r[i] -= A[i][j] * x[j];
			}
		}

		double tau = 0;
		double denomTau = 0;

		for (int i = 0; i < N; i++) {
			double Ar = 0;

			for (int j = 0; j < N; j++) {
				Ar += A[i][j] * r[j];
			}

			tau += Ar * r[i];
			denomTau += Ar * Ar;
		}

		tau /= denomTau;

		for (int i = 0; i < N; i++) {
			x[i] = prevX[i] + tau * r[i];
		}

		double maxErr = abs(x[0] - prevX[0]);
		for (int i = 1; i < N; i++)
			if (abs(x[i] - prevX[i]) > maxErr)
				maxErr = abs(x[i] - prevX[i]);

		if (maxErr < eps)
			break;

	}

	delete[] prevX;
	delete[] r;
}



/**
 * Метод сопряженных градиентов
 */
void borisovaem::lab7()
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

    }while(true);

    delete[] new_x;
    delete[] new_r;
    delete[] r;
    delete[] z;
}


void borisovaem::lab8()
{

}


void borisovaem::lab9()
{
  double eps = 1e-15;
	double* y = new double[N];
	double lambda = 0;
	x[0] = 1;

	while (true)
	{
		double newLambda = 0;
		for (int i = 0; i < N; i++) {
			y[i] = 0;

			for (int j = 0; j < N; j++) {
				y[i] += A[i][j] * x[j];
			}

			newLambda += y[i] * x[i];
		}

		if (abs(newLambda - lambda) < eps) break;

		lambda = newLambda;

		double n = 0;
		for (int i = 0; i < N; i++) {
			n += y[i] * y[i];
		}
		n = sqrt(n);

		for (int i = 0; i < N; i++) {
			x[i] = y[i] / n;
		}
	}
	x[0] = lambda;

	delete[] y;
}


std::string borisovaem::get_name()
{
  return "Borisova E.M.";
}
