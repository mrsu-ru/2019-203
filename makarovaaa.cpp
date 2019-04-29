#include "makarovaaa.h"

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
void makarovaaa::lab1()
{
cout<<"Hello, World!"<<endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void makarovaaa::lab2()
{
  system("chcp 1251");
  system("cls");
  x = gauss(A, b, N);

}



/**
 * Метод прогонки
 */
void makarovaaa::lab3()
{
//Это частный случай для метода Гаусса, используется, когда матрица трёхдиагональная
double *up, *mid, *low;  
up = new double[N]; 
mid = new double[N]; 
low = new double[N]; 
double k; 

low[0] = 0; 
up[N-1] = 0; 

for (int i = 0; i < N; i++)//Заполняем "диагональные" массивы 
{ 
	if (i - 1 >= 0 && i - 1 < N) 
	up[i] = A[i-1][i];  
	mid[i] = A[i][i];  
	if (i + 1 >= 0 && i + 1 < N) 
	low[i] = A[i+1][i];  
} 
//Прямая прогонка 
for (int i = 1; i < N; i++) //Вычисляем коэффициенты прогонки 
{ 
	k = low[i]/mid[i-1]; 
	mid[i] = mid[i] - k*up[i-1]; 
	b[i] = b[i] - k*b[i-1]; 
} 

//Обратная прогонка 
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
void makarovaaa::lab4()
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
if(abs(x[i]-nx[i]) > Err)//Максимальная разница между элементами решения 
Err = abs(x[i]-nx[i]);
}
for(int i=0; i<N; i++) 
	x[i]=nx[i];
cout<<step<<"    "<<Err<<endl;
}while (Err>Eps);

delete[] nx;
}



/**
 * Метод Якоби или Зейделя
 */
void makarovaaa::lab5()
{
//Метод Зейделя
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
Err= abs(oldx[0]-x[0]); 
for(int i=0; i<N; i++) 
{ 
if(abs(oldx[i]-x[i]) > Err) 
Err = abs(oldx[i]-x[i]);//максимальная разница между предыдущим решением и текущим. 
} 
} while(Err >= eps); 
delete [] oldx;
}



/**
 * Метод минимальных невязок
 */
void makarovaaa::lab6()
{

}



/**
 * Метод сопряженных градиентов
 */
void makarovaaa::lab7()
{

}


void makarovaaa::lab8()
{

}


void makarovaaa::lab9()
{

}


std::string makarovaaa::get_name()
{
  return "Makarova Anna Alekseevna";
}
