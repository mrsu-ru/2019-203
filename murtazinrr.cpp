#include "murtazinrr.h"

/**
 * Введение в дисциплину
 */
void murtazinrr::lab1()
{
cout <<"Hello world" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void murtazinrr::lab2()
{
 double p;
	int maxn;

    for (int k=0; k<N-1; k++)
    {
        maxn = k;
        for (int i=k+1; i<N; i++)
			if(abs(A[i][k]) > abs(A[maxn][k])) maxn = i; ///Выбор главного элемента
        std::swap(A[maxn], A[k]); ///Меняем строки местами
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
void murtazinrr::lab3()
{
double *up, *mid, *low; 
up = new double[N]; 
mid = new double[N]; 
low = new double[N]; 
double k; 

low[0] = 0; 
up[N-1] = 0; 

for (int i = 0; i < N; i++)
{ 
	if (i - 1 >= 0 && i - 1 < N) 
	up[i] = A[i-1][i]; 
	mid[i] = A[i][i]; 
	if (i + 1 >= 0 && i + 1 < N) 
	low[i] = A[i+1][i]; 
} 
//Прямая прогонка 
for (int i = 1; i < N; i++)
{ 
	k = low[i]/mid[i-1]; 
	mid[i] = mid[i] - k*up[i-1]; 
	b[i] = b[i] - k*b[i-1]; 
} 


x[N-1] = b[N-1]/mid[N-1]; 

for (int i = N - 2; i >= 0; i--) 
	x[i]=(b[i]-up[i]*x[i+1])/mid[i]; 

delete[] up;
delete[] mid; 
delete[] low; 
}



/**
 * Метод простых итераций
 */
void murtazinrr::lab4()
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
void murtazinrr::lab5()
{
//Метод Якоби
double *oldx = new double[N]; 

for (int i=0; i<N; i++) { 
x[i]=0; 
} 

double Err=0.0; 
double eps=1e-20; 
int k=0; 

do { 
k++; 
Err=0.0; 
for(int i=0; i<N; i++) 
oldx[i]=x[i]; 
for(int i=0; i<N; i++) 
{ 
double s=0; 
for(int j=0; j<i; j++) 
s += A[i][j] * oldx[j]; 
for(int j=i+1; j<N; j++) 
s += A[i][j] * oldx[j]; 
x[i]=(b[i] - s)/A[i][i]; 
} 
Err= std::abs(oldx[0]-x[0]); 
for(int i=0; i<N; i++) 
{ 
if(std::abs(oldx[i]-x[i]) > Err) 
Err = std::abs(oldx[i]-x[i]); 
} 
} while(Err >= eps); 
delete [] oldx; 
}



/**
 * Метод минимальных невязок
 */
void murtazinrr::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void murtazinrr::lab7()
{

}


void murtazinrr::lab8()
{

}


void murtazinrr::lab9()
{

}


std::string murtazinrr::get_name()
{
  return "Murtazin R.R.";
}
