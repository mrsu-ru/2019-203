#include "malovava.h"

/**
 * Введение в дисциплину
 */
void malovava::lab1()
{
cout<<"Hello world!!!"<<endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void malovava::lab2()
{
double Q = 0;

    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            Q = A[j][i] / A[i][i];
            for (int k = i; k < N; k++)
            {
                A[j][k] -= Q * A[i][k];
            }
            b[j] -= Q * b[i];
        }
    }

    for(int i = 0; i < N; i++)
	{
        x[i] = b[i];
	}

    for (int i = N - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < N; j++)
		{
			x[i] -= A[i][j] * x[j];
		}

     x[i] /= A[i][i];
	}

}



/**
 * Метод прогонки
 */
void malovava::lab3()
{

	double* P = new double[N]; // "альфа"
    double* Q = new double[N]; // "бетта"

    for (int i = 0; i < N; i++)
    {
        P[i] = 0;
        Q[i] = 0;
    }

    P[0] = -A[0][1] / A[0][0];
    Q[0] = b[0] / A[0][0];

    for(int i = 1; i < N; i++)
    {
        P[i] = A[i][i + 1] / (-A[i][i - 1] * P[i - 1] - A[i][i]);
        Q[i] = (-b[i] + A[i][i - 1] * Q[i - 1]) / (-A[i][i - 1] * P[i - 1] - A[i][i]);
    }

    x[N - 1] = Q[N - 1];

    for(int i = N - 2; i >= 0; i--)
	{
        x[i] = P[i] * x[i + 1] + Q[i];
	}

	delete[] P;
	delete[] Q;
}



/**
 * Метод простых итераций
 */
void malovava::lab4()
{
double eps = 1e-15;
	double tau = 1e-5;

	double* prevX = new double[N];

	while (true) {

		for (int i = 0; i < N; i++)
			prevX[i] = x[i];

		for (int i = 0; i < N; i++) {
			double sum = 0;
			for (int j = 0; j < N; j++)
				sum += A[i][j] * prevX[j];
			x[i] = prevX[i] - tau * (sum - b[i]);
		}

		double maxErr = abs(x[0] - prevX[0]);
		for (int i = 1; i < N; i++)
			if (abs(x[i] - prevX[i]) > maxErr)
				maxErr = abs(x[i] - prevX[i]);

		if (maxErr < eps)
			break;

	}

	delete[] prevX;


}



/**
 * Метод Якоби или Зейделя
 */
void malovava::lab5()
{   double eps=1e-20;
	for (int i = 0; i < N; i++) {
		x[i] = 0;
	}
	double *prev_x = new double[N];

	double norma = 0.0;
	do {
		for (int i = 0; i < N; i++) {
			prev_x[i] = x[i];           //записываем предыдущее решение
		}
		for (int i = 0; i < N; i++) {   // найдем новое решение на текущей итерации
			double result = b[i];
			for (int j = 0; j < N; j++) {
				if (i != j) {
					result -= (A[i][j] * prev_x[j]);
				}
			}
			x[i] = result / A[i][i];
		}
		norma = 0.0;    //найдем погрешность
		for (int i = 0; i < N; i++) {
			if (std::abs(prev_x[i] - x[i]) > norma) {
				norma = std::abs(prev_x[i] - x[i]);   //максимальная разница между предыдущим решением и текущим
			}
		}
	} while (norma > eps);

	delete[] prev_x;
}



/**
 * Метод минимальных невязок
 */
void malovava::lab6()
{
   double *r=new double[N];
    for (int i=0; i<N; i++)
              r[i]=0;

    double eps=10e-16;
    double nx=0;
    double tay=0;

       for(;;){
            double differ=0, sum1=0, sum2=0;
            for (int i=0; i<N; i++){
                r[i]=b[i];
                for (int j=0; j<N; j++)
                    r[i]-=A[i][j]*x[j];    // найдем вектор невязок
            }
            for (int i=0; i<N; i++){
                double vec=0;
                for (int k=0; k<N; k++)
                     vec+=A[i][k]*r[k];
                sum1+=r[i]*vec;
                sum2+=vec*vec;
            }
            tay=sum1/sum2;               //итерационный параметр
            for (int i=0; i<N; i++){
                nx=x[i];
                x[i]+=r[i]*tay;       // новый вектор решений на новой итерации

                if(abs(x[i]-nx)>differ)
                    differ=abs(x[i]-nx);
            }
            if(differ<eps) break;
            }
            delete[] r;
}



/**
 * Метод сопряженных градиентов
 */
void malovava::lab7()
{
double eps = 1e-15;

	double* prevX = new double[N];
	double* prevR = new double[N];
	double* r = new double[N];
	double* z = new double[N];
	for (int i = 0; i < N; i++) {
		r[i] = b[i];
		z[i] = r[i];
	}

	while (true) {

		for (int i = 0; i < N; i++) {
			prevR[i] = r[i];
			prevX[i] = x[i];
		}

		double alpha = 0, Alp = 0;

		for (int i = 0; i < N; i++) {
			double Az = 0;
			for (int j = 0; j < N; j++) {
				Az += A[i][j] * z[j];
			}
			alpha += prevR[i] * prevR[i];
			Alp += Az * z[i];
		}
		alpha /= Alp;

		for (int i = 0; i < N; i++) {
			x[i] = prevX[i] + alpha * z[i];
		}

		double maxErr = abs(x[0] - prevX[0]);
		for (int i = 1; i < N; i++)
			if (abs(x[i] - prevX[i]) > maxErr)
				maxErr = abs(x[i] - prevX[i]);

		if (maxErr < eps)
			break;

		for (int i = 0; i < N; i++) {
			double Az = 0;

			for (int j = 0; j < N; j++) {
				Az += A[i][j] * z[j];
			}

			r[i] = prevR[i] - alpha * Az;
		}

		double beta = 0, Bet = 0;
		for (int i = 0; i < N; i++) {
			beta += r[i] * r[i];
			Bet += prevR[i] * prevR[i];
		}
		beta /= Bet;

		for (int i = 0; i < N; i++) {
			z[i] = r[i] + beta * z[i];
		}
	}

	delete[] prevX;
	delete[] r;
	delete[] prevR;
	delete[] z;
}


void malovava::lab8()
{

}


void malovava::lab9()
{
 double eps=pow(10,-10);
 double lambda=0, newlambda=0, y[N], md;
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
        md=0;
        for (int i=0; i<N; i++){
            md+=pow(y[i], 2);
        }
        md=sqrt(md);
        for (int i=0; i<N; i++){
            x[i]=y[i]/md;
        }
    } while (abs(lambda-newlambda)>eps);

    x[0]=newlambda;
}


std::string malovava::get_name()
{
  return "Malova V.A.";
}
