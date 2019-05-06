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

}



/**
 * Метод простых итераций
 */
void eremkinnv::lab4()
{

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
