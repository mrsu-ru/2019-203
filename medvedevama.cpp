#include "medvedevama.h"
#include <iomanip>
#include <cmath>
#include <locale>

double * Gauss(double **A, double *b, int N)
{
  double *x, max;
  int k, index;
  const double eps = 0.00001;
  x = new double[N]; double change[N]; int L;
  k = 0;
  while (k < N)
  {
    max = abs(A[k][k]);
    index = k;
    for (int i = k + 1; i < N; i++)
    {
      if (abs(A[i][k]) > max)
      {
        max = abs(A[i][k]);
        index = i;
      }
    }

    if (max < eps)
    {
      cout << "Решение получить невозможно из-за нулевого столбца ";
      cout << index << " матрицы A" << endl;
      return 0;
    }
    for (int i = 0; i < N; i++)
    {
    change[i]=index;
    if(index!=i) {
		swap(A[i],A[index]);
        swap(b[i],b[index]);
		gity}

    for (int j = 0; j < N; j++)
    {
      double temp = A[k][j];

      A[k][j] = A[index][j];
      A[index][j] = temp;
    }
    double temp = b[k];
    b[k] = b[index];
    b[index] = temp;

    for (int i = k; i < N; i++)
    {
      double temp = A[i][k];
      if (abs(temp) < eps) continue;
      for (int j = 0; j < N; j++)
        A[i][j] = A[i][j] / temp;
      b[i] = b[i] / temp;
      if (i == k)  continue;
      for (int j = 0; j < N; j++)
        A[i][j] = A[i][j] - A[k][j];
      b[i] = b[i] - b[k];
    }
    k++;
  }
  for (k = N - 1; k >= 0; k--)
  {
    x[k] = b[k];
    for (int i = 0; i < k; i++)
      b[i] = b[i] - A[i][k] * x[k];
  }
   for (int i = N-1; i >= 0; i--)
    {
       if (index!=i)
	      { L = change[i];
            swap(A[L], A[i]);
            swap(b[L], b[i]);
			}
	}
  return x;
}
}


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
 system("chcp 1251");
  system("cls");
  x = Gauss(A, b, N);
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
