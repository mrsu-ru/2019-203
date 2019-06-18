#include "eremkinnv.h"

/**
 * Введение в дисциплину
 */
void eremkinnv::lab1()
{
    cout<<"Hallo world" <<endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void eremkinnv::lab2()
{
double p;
    int maxn;

    for (int k=0; k<N-1; k++)
    {

        maxn = k;
        for (int i=k+1; i<N; i++)
            if(abs(A[i][k]) > abs(A[maxn][k])) maxn = i;
        std::swap(A[maxn], A[k]);
        std::swap(b[maxn], b[k]);

        for (int i=k+1; i<N; i++)
        {
            p = A[i][k]/A[k][k];
            for (int j=k; j<N; j++)
                A[i][j] -= p*A[k][j];
            b[i] -= p*b[k];
        }
    }

    for(int i = 0; i<N; i++)
    {
        x[i]=b[i];
    }

    for (int i=N-1; i>=0; i--)
    {
        for (int j=i+1;j<N;j++)
            x[i]-=A[i][j]*x[j];
        x[i] /= A[i][i];
    }
}



/**
 * Метод прогонки
 */
void eremkinnv::lab3()
{
	double* alpha = new double[N - 1];
	double* gamma = new double[N];

	for (int i = 0; i < N; i++) {

		gamma[i] = A[i][i];
		if (i != 0) gamma[i] += A[i][i - 1] * alpha[i - 1];

		if (i != N - 1) alpha[i] = -A[i][i + 1] / gamma[i];

		x[i] = b[i] / gamma[i];
		if (i != 0) x[i] -= A[i][i - 1] * x[i - 1] / gamma[i];

	}

	for (int i = N - 2; i >= 0; i--)
		x[i] += alpha[i] * x[i + 1];

	delete[] alpha;
	delete[] gamma;
}



/**
 * Метод простых итераций
 */
void eremkinnv::lab4()
{
double Eps=1e-15;
double Err;
double *nx = new double[N];

for (int i=0;i<N;i++)
	x[i]=b[i];
int step=0;

do{
step++;
  for(int i=0;i < N;i++)
  {
   nx[i]=-b[i];

   for(int j=0;j < N;j++)
   {
    if(i!=j)
     nx[i]+=A[i][j]*x[j];
   }

   nx[i]/=-A[i][i];
  }
  Err=0;
for(int i=0; i<N; i++) {
if(std::abs(x[i]-nx[i]) > Err)
Err = std::abs(x[i]-nx[i]);
}
for(int i=0; i<N; i++)
	x[i]=nx[i];
std::cout<<step<<"    "<<Err<<endl;
}while (Err>Eps);

delete[] nx;
}



/**
 * Метод Якоби или Зейделя
 */
void eremkinnv::lab5()
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
                new_x[i] -= A[i][j]*x[j];
            }

            new_x[i] /= A[i][i];
            if (!condition) condition = (fabs(new_x[i] - x[i]) > eps);
            x[i] = new_x[i];
        }

    }while(condition);

    delete[] new_x;
}



/**
 * Метод минимальных невязок
 */
void eremkinnv::lab6()
{
double Eps = 1e-18;
	double Del, Res, Abs;

	double *K = new double[N];
	double *L = new double[N];
	double *xrez = new double[N];



	for (int i = 0; i<N; i++)
		xrez[i] = 0;


	do{

		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * xrez[j];
		}


		for (int i = 0; i < N; i++) {
			L[i] = K[i] - b[i];
		}


		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * L[j];
		}

		Res = 0;
		Abs = 0;


		for (int i = 0; i < N; i++) {
			Res += K[i] * L[i];
			Abs += K[i] * K[i];
		}

		if (Res==Abs) Res=1;
		else {
		Res = Res / Abs;
		}

		for (int i = 0; i < N; i++)
			x[i] = xrez[i] - Res*L[i];


		Del = abs(x[0] - xrez[0]);
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - xrez[i])>Del)
				Del = abs(x[i] - xrez[i]);
			xrez[i] = x[i];
		}
	} while (Eps < Del);
}



/*
 * Метод сопряженных градиентов
 */
void eremkinnv::lab7()
{
double eps = 1e-15;

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


/*
* Метод вращения для нахождения собственных значений матрицы
*/
void eremkinnv::lab8()
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


/*
* Нахождение максимального по модулю собственного значения матрицы
*/
void eremkinnv::lab9()
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


std::string eremkinnv::get_name()
{
  return "Eremkin";
}
