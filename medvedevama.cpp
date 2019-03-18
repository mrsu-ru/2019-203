#include "medvedevama.h"
#include <iomanip>
#include <cmath>
#include <locale>

double * gauss(double **A, double *b, int N)
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
    if(index!=i) {swap(A[i],A[index]);
                  swap(b[i],b[index]); }

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
   for (int i = N-1; i >= 0; i--) { if (index!=i) {L = change[i];
                                              swap(A[L], A[i]);
                                              swap(b[L], b[i]); } }
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
  x = gauss(A, b, N);
  }



/**
 * Метод прогонки
 */
void medvedevama::lab3()
{

}



/**
 * Метод простых итераций
 */
void medvedevama::lab4()
{

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
