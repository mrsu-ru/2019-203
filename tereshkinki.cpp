#include "tereshkinki.h"

/**
 * Введение в дисциплину
 */
void tereshkinki::lab1()
{
cout<<"Hello World!!";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void tereshkinki::lab2()
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
void tereshkinki::lab3()
{
	double *P = new double [N]; ///Коэффициенты "альфа"
	double *Q = new double [N]; ///Коэффициенты "бетта"

	P[0] = -A[0][1]/A[0][0];
	Q[0] = b[0]/A[0][0];

	for(int i=1; i<N; i++) ///Определяем прогоночные коэффициенты
	{
		P[i] = A[i][i+1]/(-A[i][i] - A[i][i-1]*P[i-1]);
		Q[i] = (-b[i] + A[i][i-1]*Q[i-1])/(-A[i][i] - A[i][i-1]*P[i-1]);
	}

	x[N-1] = Q[N-1];
	for(int i=N-2; i>=0; i--) ///Определяем решение
		x[i] = P[i]*x[i+1] + Q[i];

	delete [] P;
	delete [] Q;
}



/**
 * Метод простых итераций
 */
void tereshkinki::lab4()
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
void tereshkinki::lab5()
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
void tereshkinki::lab6()
{
	
	double Eps = 1e-18;//погрешность
	double Del, Res, Abs;//погрешность, невязка, модуль

	double *K = new double[N];//w
	double *L = new double[N];//v
	double *xrez = new double[N];
	
	
	//задаём первоначальное приближение
	for (int i = 0; i<N; i++)
		xrez[i] = 0;

	//цикл для нахождения корней
	do{
		//находим редуцированную систему(одна часть)
		for (int i = 0; i < N; i++) {
			K[i] = 0;
			for (int j = 0; j < N; j++)
				K[i] += A[i][j] * xrez[j];
		}

		//находим редуцированную систему(вторая часть)
		for (int i = 0; i < N; i++) {
			L[i] = K[i] - b[i];//нахождение вектора невязки
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
void tereshkinki::lab7()
{
	double Eps = 1e-20;//заданная погрешность
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
		
		//нахождение величины смещения по направлению градиента(скалярного шага)
		for (int i = 0; i < N; i++) {
			s += K[i] * L[i];
			sAbs += M[i] * K[i];
		}
		if (s == sAbs)
			s = 1;
		else 
			s = s / sAbs;
		//записываем новое приближенное решение
		for (int i = 0; i < N; i++)
			x[i] = xrez[i] - s*L[i];
		
		//проверка на уменьшение погрешности
		Del = abs(x[0] - xrez[0]);
		
		for (int i = 0; i < N; i++) {
			if (abs(x[i] - xrez[i])>Del) {
				Del = abs(x[i] - xrez[i]);
            }
			xrez[i] = x[i];
		}
	} while (Eps < Del);
}


void tereshkinki::lab8()
{

}


void tereshkinki::lab9()
{

}


std::string tereshkinki::get_name()
{
  return "Tereshkin Kirill Igorevich";
}
