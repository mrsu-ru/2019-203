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

}



/**
 * Метод минимальных невязок
 */
void eremkinnv::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void eremkinnv::lab7()
{

}


void eremkinnv::lab8()
{

}


void eremkinnv::lab9()
{

}


std::string eremkinnv::get_name()
{
  return "Eremkin";
}
