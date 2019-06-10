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



/**
 * Метод сопряженных градиентов
 */
void murtazinrr::lab7()
{
double Eps = 1e-20;
	double Del, s, sAbs;


	double *K = new double[N];
	double *L = new double[N];
	double *M = new double[N];
	double *xrez = new double[N];
	
	
	
	for (int i = 0; i<N; i++){
		xrez[i] = 0;
	}
	
	
	do {
		
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
		
		
		for (int i = 0; i < N; i++) {
			M[i] = 0;
			for (int j = 0; j < N; j++) {
				M[i] += A[i][j] * K[j];
			}
		}
		
		s = 0;
		sAbs = 0;
		
		
		for (int i = 0; i < N; i++) {
			s += K[i] * L[i];
			sAbs += M[i] * K[i];
		}
		if (s == sAbs)
			s = 1;
		else 
			s = s / sAbs;
		
		for (int i = 0; i < N; i++)
			x[i] = xrez[i] - s*L[i];
		
		
		Del = abs(x[0] - xrez[0]);
		
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - xrez[i])>Del)
				Del = abs(x[i] - xrez[i]);
				xrez[i] = x[i];
		}
	} while (Eps < Del);
}


void murtazinrr::lab8()
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

    }while(true);

    for (int i = 0; i < N; i++) x[i] = A[i][i];

    for (int i = 0; i < N; i++) delete[] H[i];
    delete[] H;
}


void murtazinrr::lab9()
{
double eps=pow(10,-10), lambda=0, newlambda=0, y[N], mody;
    x[0]=1;
    do{
        for (int i=0; i<N; i++){
            y[i]=0;
            for (int j=0; j<N; j++){
                y[i]+=A[i][j]*x[j];
            }
        }
        lambda=newlambda;
        newlambda=0;
        for (int i=0; i<N; i++){
            newlambda+=y[i]*x[i];
        }
        mody=0;
        for (int i=0; i<N; i++){
            mody+=pow(y[i], 2);
        }
        mody=sqrt(mody);
        for (int i=0; i<N; i++){
            x[i]=y[i]/mody;
        }
    }while (abs(lambda-newlambda)>eps);
    x[0]=newlambda;
}


std::string murtazinrr::get_name()
{
  return "Murtazin R.R.";
}
