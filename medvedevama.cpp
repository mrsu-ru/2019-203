#include "medvedevama.h"
#include <iomanip>
#include <cmath>
#include <locale>



/**
 * Введение в дисциплину
 */
void medvedevama::lab1()
{
cout<<"Hello, World!"<<endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void medvedevama::lab2()
{
double temp;

for (int k = 0; k < N; k++) { 
	int max=k;
		for(int i=k+1;i<N;i++)
			if(abs(A[i][k]) > abs(A[max][k]))
				max=i;
	for(int i=0;i<N;i++)
	std::swap(A[k][i],A[max][i]);
	std::swap(b[k],b[max]);

//Прямой ход
temp = A[k][k];
for (int j = 0; j < N; j++)
	A[k][j] = A[k][j] / temp;
b[k] =b[k]/temp;

for (int i = k + 1; i < N; i++) {
	temp = A[i][k];
	for (int j = 0; j < N; j++) {
		A[i][j] =A[i][j]- A[k][j] * temp;
	}
b[i] =b[i]- b[k] * temp;
}
}

//Обратный ход
for (int k = N - 1; k > 0; k--)
{
  for (int i = k - 1; i >= 0; i--)
    {
       temp = A[i][k];
       for (int j = 0; j < N; j++)
           A[i][j] =A[i][j]- A[k][j] * temp;
       b[i] =b[i] - b[k] * temp;
    }
}

for(int i=0; i<N; i++)
x[i]=b[i];
}



/**
 * Метод прогонки
 */
void medvedevama::lab3()
{
  	double *upper, *middle, *lower; // верхняя, главная и нижняя диагонали
	upper = new double[N];
	middle = new double[N];
	lower = new double[N];
	double k;

	lower[0] = 0;
	upper[N-1] = 0;

	for (int i = 0; i < N; i++)//Заполнение"диагональных" массивов
	{
		if (i - 1 >= 0 && i - 1 < N)
		upper[i] = A[i-1][i];
		middle[i] = A[i][i];
		if (i + 1 >= 0 && i + 1 < N)
		lower[i] = A[i+1][i]; //нижняя
	}
	//Прямой ход
	for (int i = 1; i < N; i++) //Вычисляем коэффициенты прогонки
	{
		k = lower[i]/middle[i-1];
		middle[i] = middle[i] - k*upper[i-1];
		b[i] = b[i] - k*b[i-1];
	}

	//Обратный ход
	x[N-1] = b[N-1]/middle[N-1]; //Вычисляется решение

	for (int i = N - 2; i >= 0; i--) //рекуррентная формула для вычисления остальных неизвестных
		x[i]=(b[i]-upper[i]*x[i+1])/middle[i];

	delete[] upper, middle, lower;
}



/**
 * Метод простых итераций
 */
void medvedevama::lab4()
{
  	double mis=1e-15;//порядок ошибки
	double error;
	double *nx = new double[N];//для хранения промежуточных значений

	for (int i=0;i<N;i++)//для первичного приближения возьмём столбец свободных членов
		x[i]=b[i];
	int step=0;

	do
	{
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
	error=0;
	for(int i=0; i<N; i++)
	{
		if(std::abs(x[i]-nx[i]) > error)//Максимальная разница между элементами решения
			error = std::abs(x[i]-nx[i]);
	}
	for(int i=0; i<N; i++)
		x[i]=nx[i];
	std::cout<<step<<"    "<<error<<endl;
	}
	while (error>mis);
	delete[] nx;
}



/**
 * Метод Якоби или Зейделя
 */
void medvedevama::lab5()
{
//Метод Якоби
double *oldx = new double[N];

for (int i=0; i<N; i++) {
    x[i]=0; // первоначальное решение
    }

double Err=0.0;
double eps=1e-20;
int k=0;

do {
  k++;
  Err=0.0;
  for(int i=0; i<N; i++)
      oldx[i]=x[i]; //предыдущее решение сюда
  for(int i=0; i<N; i++)
  {
    double s=0; //вычисляем s, но не берём диагональные элементы
    for(int j=0; j<i; j++)
        s += A[i][j] * oldx[j];
    for(int j=i+1; j<N; j++)
        s += A[i][j] * oldx[j];
     x[i]=(b[i] - s)/A[i][i]; // вычисляется новое решение
    }
Err= std::abs(oldx[0]-x[0]);
for(int i=0; i<N; i++)
{
  if(std::abs(oldx[i]-x[i]) > Err)
     Err = std::abs(oldx[i]-x[i]);//максимальная допустимая разница между предыдущим решением и текущим.
}
} while(Err >= eps);
delete [] oldx;
}



/**
 * Метод минимальных невязок
 */
void medvedevama::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void medvedevama::lab7()
{

}


void medvedevama::lab8()
{

}


void medvedevama::lab9()
{

}


std::string medvedevama::get_name()
{
  return "Medvedeva M.A.";
}
