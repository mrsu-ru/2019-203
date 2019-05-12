#include "maslovma.h"

/**
 * Введение в дисциплину
 */
void maslovma::lab1()
{
  cout << "Hello world!" << endl;
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void maslovma::lab2()
{

  double temp;
    int maxInd;

      for(int i=0;i<N-1;i++)
        {
         maxInd = i;

          for(int k=i+1;k<N;k++)
           {
            if(abs(A[k][i]) > abs(A[maxInd][i]))
              {
                maxInd = k;
              }
           }

           if (maxInd != i)
            {
              for(int k=i;k<N;k++)
               {
                std::swap(A[i][k],A[maxInd][k]);
               }
               std::swap(b[i],b[maxInd]);
            }

         for(int j=i+1;j<N;j++)
          {
             temp = A[j][i]/A[i][i];

             b[j] -= temp*b[i];

              for(int q=i;q<N;q++)
                {
                  A[j][q] -= temp*A[i][q];
                }
           }
        }

        for(int i=N-1;i>=0;i--)
         {
            b[i]/=A[i][i];

             for(int j=i-1;j>=0;j--)
              {
                b[j] -= b[i]*A[j][i];
              }

         }
    for(int i=0;i<N;i++)
         x[i] = b[i];
}



/**
 * Метод прогонки
 */
void maslovma::lab3()
{
   double* P = new double[N-1];
    double* Q = new double[N];

     P[0] = -A[0][1]/A[0][0];
     Q[0] = b[0]/A[0][0];

      for(int i=1;i<N-1;i++)
       {
         P[i] = -A[i][i+1]/(A[i][i]+A[i-1][i]*P[i-1]);
         Q[i] = (b[i] - A[i][i-1]*Q[i-1])/(A[i][i] + A[i][i-1]*P[i-1]);
       }
        Q[N-1] = (b[N-1] - A[N-1][N-2]*Q[N-2])/(A[N-1][N-1] + A[N-1][N-2]*P[N-2]);


       x[N-1] = Q[N-1];
        for(int i=N-2;i>=0;i--)
          {
            x[i] = P[i]*x[i+1]+Q[i];
          }


     delete [] P;
     delete [] Q;

}



/**
 * Метод простых итераций
 */
void maslovma::lab4()
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
void maslovma::lab5()
{
  const double eps = 10E-20;

    double* D = new double[N];
    double** LU = new double*[N];

     for(int i=0;i<N;i++)
      LU[i] = new double[N];

    for(int i=0;i<N;i++)
     for(int j=0;j<N;j++)
       {
         if(i!=j)
          {
            LU[i][j] = A[i][j];
          }
         else
          {
            D[i] = A[i][j];
            LU[i][j] = 0;
          }
       }

  double* new_x = new double[N];
  double S;
  double temp;
  double k=0;

  for(int i=0;i<N;i++)
        new_x[i] = D[i];

   do{
      for(int i=0;i<N;i++)
       {
         S = 0;
        for(int j=0;j<N;j++)
         {
           if(j!=i)
            S +=  LU[i][j]*new_x[j];
         }
         x[i] = (b[i] - S)/D[i];
       }

        temp = fabs(new_x[0] - x[0]);

       for(int i=1;i<N;i++)
        if( fabs(new_x[i] - x[i]) > temp )
         temp = fabs(new_x[i] - x[i]);

        for(int i=0;i<N;i++)
            new_x[i] = x[i];

     k++;
   }while(temp >= eps);

  delete [] new_x;

  for(int i=0;i<N;i++)
     delete LU[i];

    delete [] LU;
    delete [] D;

}


/**
 * Метод минимальных невязок
 */
void MatrVect(int N, double **M, double *V, double *R)
{
    for(int i=0; i<N; i++)
    {
        R[i]=0;
        for(int j=0; j<N; j++)
              R[i]+= M[i][j]*V[j];
    }
}

double ScalarVect(int N, double* v1, double* v2)
{
	double result=0;
	for (int i=0; i<N; i++)
        result+=(v1[i]*v2[i]);
	return result;
}


void maslovma::lab6()
{
   const double eps = 1.e-6;

    int count = 0; ///  количество итераций
    double *U = new double [N];
    double *r = new double [N];
    double *TempX = new double[N];
    double *p = new double[N];
    double Tau = 0.0;

    for (int i=0; i<N; i++)
        TempX[i]=0; /// первое приближение задаём нулевым

    do
    {
        MatrVect(N, A, TempX, U);

        for(int i=0; i<N; i++)
        {
            r[i] = U[i]-b[i]; /// Вектор невязок
        }

        MatrVect(N, A, r, U);

        double TempTau1 = ScalarVect(N, U, r);
        double TempTau2 = ScalarVect(N, U, U);
        if (TempTau2 == 0) break;

        Tau = TempTau1/TempTau2; /// Итерационный параметр

        for(int i=0; i<N; i++)
            x[i] = TempX[i] - Tau*r[i];

        for(int i=0; i<N; i++)
            p[i] = x[i]-TempX[i];

        count++;
    } while ((sqrt(ScalarVect(N, p, p)) >= eps) && (count < 500000));

    delete[] U;
    delete[] r;
    delete[] p;
    delete[] TempX;

}



/**
 * Метод сопряженных градиентов
 */
void maslovma::lab7()
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

		double alpha = 0, denAlpha = 0;

		for (int i = 0; i < N; i++) {
			double Az = 0;
			for (int j = 0; j < N; j++) {
				Az += A[i][j] * z[j];
			}
			alpha += prevR[i] * prevR[i];
			denAlpha += Az * z[i];
		}
		alpha /= denAlpha;

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

		double beta = 0, denBeta = 0;
		for (int i = 0; i < N; i++) {
			beta += r[i] * r[i];
			denBeta += prevR[i] * prevR[i];
		}
		beta /= denBeta;

		for (int i = 0; i < N; i++) {
			z[i] = r[i] + beta * z[i];
		}
	}

	delete[] prevX;
	delete[] r;
	delete[] prevR;
	delete[] z;
}


void maslovma::lab8()
{

}

/*Степенной метод нахождения 
максимального по модулю собственного значения
 и соответствующего ему собственного вектора*/

void maslovma::lab9()
{
   double *new_x = new double[N], lam = 0, eps = 1.e-10;
    for (int i = 1; i < N; i++) x[i] = 0;

    x[0] = 1;
    do
    {
        double new_lam = 0;
        for (int i = 0; i < N; i++)
        {
            new_x[i] = 0;

            for (int j = 0; j < N; j++)
                new_x[i] += A[i][j] * x[j];

            new_lam += x[i] * new_x[i];
        }

        if (fabs(new_lam - lam) < eps) break;

        lam = new_lam;
        double n = 0;
        for (int i = 0; i < N; i++) n += new_x[i] * new_x[i];

        n = sqrt(n);
        for (int i = 0; i < N; i++) x[i] = new_x[i] / n;

    } while (true);

    x[0] = lam;
    delete[] new_x;
}


std::string maslovma::get_name()
{
  return "MaslovMA";
}
