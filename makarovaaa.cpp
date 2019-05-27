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
double *new_x = new double[N], tau = 0.001, eps = 1.e-9;
    for (int i = 0; i < N; i++)
        x[i] = 0;

    do
    {
        for (int i = 0; i < N; i++)
        {
            double temp = 0;
            for (int j = 0; j < N; j++)
                temp += A[i][j] * x[j];

            new_x[i] = x[i] + tau * (b[i] - temp);
        }

        double maxdif = 0;
        for (int i = 0; i < N; i++)
        {
            if (fabs(x[i] - new_x[i]) > maxdif) maxdif = fabs(x[i] - new_x[i]);
            x[i] = new_x[i];
        }

        if (maxdif < eps) break;
    }while(true);

    delete[] new_x;
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
double Eps = 1e-18;
	double Del, Res, Abs;//погрешность, невязка, модуль

	double *K = new double[N];
	double *L = new double[N];
	double *xrez = new double[N];
	
	
	//задаём первоначальное приближение
	for (int i = 0; i<N; i++)
		xrez[i] = 0;

	
	do{
		//находим редуцированную систему
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * xrez[j];
		}

		
		for (int i = 0; i < N; i++) {
			L[i] = K[i] - b[i]; //нахождение вектора невязки
		}

		
		//нахождение скалярного произведения матрицы системы и вектора невязки
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * L[j];
		}

		Res = 0;
		Abs = 0;
		
		//нахождение значения итерационного параметра
		for (int i = 0; i < N; i++) {
			Res += K[i] * L[i];
			Abs += K[i] * K[i];
		}
		
		if (Res==Abs) Res=1;
		else {
		Res = Res / Abs;
		}
		//получение приближения решения
		for (int i = 0; i < N; i++)
			x[i] = xrez[i] - Res*L[i];
		
		//Проверка на уменьшение погрешности
		Del = abs(x[0] - xrez[0]);
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - xrez[i])>Del)
				Del = abs(x[i] - xrez[i]);
			xrez[i] = x[i];
		}
	} while (Eps < Del);

}



/**
 * Метод сопряженных градиентов
 */
void makarovaaa::lab7()
{
double Eps = 1e-20;
	double Del, s, sAbs;//погрешность итерации, скалярный шаг, модуль шага


	double *K = new double[N];
	double *L = new double[N];
	double *M = new double[N];
	double *xrez = new double[N];//итерационные решения
	
	
	//задание начального приближения 
	for (int i = 0; i<N; i++){
		xrez[i] = 0;
	}
	
	
	do {
		//нахождение скалярного произведения матрицы системы и вектора приближенного решения
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * xrez[j];
		}
		
		//нахождение градиента
		for (int i = 0; i < N; i++) {
			L[i] = K[i] - b[i];
		}
		
		//скалярное произведение матрицы системы и градиента
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * L[j];
		}
		
		
		for (int i = 0; i < N; i++) {
			M[i] = 0;
			for (int j = 0; j < N; j++) {
				M[i] += A[i][j] * K[j];
			}
		}
		
		s = 0;
		sAbs = 0;
		
		//нахождение величины смещения по направлению градиента
		for (int i = 0; i < N; i++) {
			s += K[i] * L[i];
			sAbs += M[i] * K[i];
		}
		if (s == sAbs)
			s = 1;
		else 
			s = s / sAbs;
		// новое приближенное решение
		for (int i = 0; i < N; i++)
			x[i] = xrez[i] - s*L[i];
		
		//проверка на уменьшение погрешности
		Del = abs(x[0] - xrez[0]);
		
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - xrez[i])>Del)
				Del = abs(x[i] - xrez[i]);
				xrez[i] = x[i];
		}
	} while (Eps < Del);
}


void makarovaaa::lab8()
{

}


void makarovaaa::lab9()
{
double * Y = new double[N];//первый вектор приближения
	double * M = new double[N];//второй вектор приближения
	double maxL, L, sum;
	double EPS = 1e-15;
	
	
	for (int i = 0; i < N; i++)
		Y[i] = 0;
	Y[0] = 1;
	
	do{
		sum = 0;
		 
		for (int i = 0; i < N; i++)
			sum += Y[i] * Y[i];
		
		L = sqrt(sum);
		
		//построение последовательности векторов
		for (int i = 0; i < N; i++)
		{
			M[i] = 0;
			for (int j = 0; j < N; j++)
				M[i] += A[i][j] * Y[j] / L;
		}
		sum = 0;
		
		//сравнение нормы полученного вектора с заданной погрешностью
		for (int i = 0; i < N; i++)
			sum += M[i] * M[i];
		maxL = sqrt(sum);
		
		for (int i = 0; i<N; i++)
			Y[i] = M[i];
	} while (abs(maxL - L)>EPS);

	//cout << maxL << endl;
}

double static f(double x)
{
    double f=x*x*x-2*x*x-5*x+6;
	return f;
}
//производная этой функции 
double static df(double x)
{
    double df=3*x*x-4*x-5;
	return df;
}

//сжимающее отображение 
double static g(double x)
{
    return x - f(x)/df(x);
}


std::string makarovaaa::get_name()
{
  return "Makarova Anna Alekseevna";
}
