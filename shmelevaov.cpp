#include "shmelevaov.h"
const double eps = 1.e-15;
/**
 * Введение в дисциплину
 */
void shmelevaov::lab1()
{
	cout << "Hello World!";
}


/**
 * Метод Гаусса с выбором главного элемента
 */

 void change_lines (double **A, double *b, int first_line, int second_line, int N, int *help)
{
    for (int i = first_line; i < N; i++)
    {
        swap(A[first_line][i], A[second_line][i]);
    }

    swap(b[first_line], b[second_line]);

    swap(help[second_line], help[first_line]);
}

void shmelevaov::lab2()
{
    int *help = new int[N];
	for (int i = 0; i < N; i++)
    {
		help[i] = i;
    }

    for (int i = 0; i < N; i++)
    {
        int maximum = i;
        for (int j = i; j < N; j++)
        {
            if (fabs(A[j][i]) > fabs(A[maximum][i]))
            {
                maximum = j;
            }
        }

        if (fabs(A[maximum][i]) > eps)
        {
           change_lines(A, b, i, maximum, N, help);
        }

        double diagonal = A[i][i];

        A[i][i] = 1;

        for (int j = i + 1; j < N; j++)
        {
            A[i][j] /= diagonal;
        }

        b[i] /= diagonal;

        for (int j = i + 1; j < N; j++)
        {
            double temp = A[j][i];
            A[j][i] = 0;
            for(int k = i + 1; k < N; k++)
            {
                A[j][k] -= A[i][k] * temp;
            }

            b[j] -= b[i] * temp;
        }
    }

    for (int i = N - 1; i > 0; i--)
	{
		for (int j = i - 1; j >= 0; j--)
		{
			double temp = A[j][i];
			A[j][i] -= A[i][i] * temp;
			b[j] -= b[i] * temp;
		}
	}

    for (int i = 0; i < N; i++)
    {
		x[help[i]] = b[i];
	}
}



/**
 * Метод прогонки
 */
void shmelevaov::lab3()
{
	double alpha[N - 1], beta[N];
    alpha[0] = - A[0][1] / A[0][0];
    beta[0] = b[0] / A[0][0];

    for (int i = 1; i < N; i++)
	{
		double y = A[i][i] + A[i][i - 1] * alpha[i - 1];
        alpha[i] = - A[i][i + 1] / y;
        beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
    }

    x[N - 1] = beta[N - 1];
    for (int i = N - 2; i >= 0; i--)
	{
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
}



/**
 * Метод простых итераций
 */
void shmelevaov::lab4()
{
	const double lambda = 0.01;
	double prev[N];
	for (int i = 0; i < N; i++)
	{
		prev[i] = 0;
	}

	double difference;
	do
	{
		for (int i = 0; i < N; i++)
		{
			double sum = 0;
			for (int j = 0; j < N; j++)
			{
				sum += A[i][j] * prev[j];
			}

			x[i] = prev[i] - lambda * (sum - b[i]);
		}

		difference = 0;
		for (int i = 0; i < N; i++)
		{
			if (fabs(x[i] - prev[i]) > difference)
			{
				difference = fabs(x[i] - prev[i]);
			}
		}

		for (int i = 0; i < N; i++)
		{
			prev[i] = x[i];
		}

	} while (difference > eps);
}



/**
 * Метод Якоби или Зейделя
 */
void shmelevaov::lab5()
{
	double prev[N];
	for (int i = 0; i < N; i++)
	{
		prev[i] = 0;
	}

	double difference;
	do
	{
		for (int i = 0; i < N; i++)
		{
			double sum = 0;
			for (int j = 0; j < i; j++)
			{
				sum += A[i][j] * prev[j];
			}

			for (int j = i + 1; j < N; j++)
			{
				sum += A[i][j] * prev[j];
			}

			x[i] = (b[i] - sum) / A[i][i];
		}

		difference = 0;
		for (int i = 0; i < N; i++)
		{
			if (fabs(x[i] - prev[i]) > difference)
			{
				difference = fabs(x[i] - prev[i]);
			}
		}

		for (int i = 0; i < N; i++)
		{
			prev[i] = x[i];
		}

	} while (difference > eps);
}



/**
 * Метод минимальных невязок
 */
void shmelevaov::lab6()
{
    double prev[N];
	for (int i = 0; i < N; i++)
	{
		prev[i] = 0;
	}

	double r[N];
	double difference;
	do
	{
		for (int i = 0; i < N; i++)
		{
			r[i] = b[i];

			for (int j = 0; j < N; j++)
			{
				r[i] -= (A[i][j] * prev[j]);
			}
		}

		double numerator_tau = 0;
		double denominator_tau = 0;
		for (int i = 0; i < N; i++)
		{
			double Ar = 0.0;

			for (int j = 0; j < N; j++)
			{
				Ar += (A[i][j] * r[j]);
			}

			numerator_tau += (Ar * r[i]);
			denominator_tau += (Ar * Ar);
		}
		double tau = numerator_tau / denominator_tau;

		for (int i = 0; i < N; i++)
		{
			x[i] = prev[i] + tau * r[i];
		}

		difference = 0;
		for (int i = 0; i < N; i++)
		{
			if (fabs(x[i] - prev[i]) > difference)
			{
				difference = fabs(x[i] - prev[i]);
			}
		}

		for (int i = 0; i < N; i++)
		{
			prev[i] = x[i];
		}
	} while (difference > eps);
}



/**
 * Метод сопряженных градиентов
 */
void shmelevaov::lab7()
{
	double prev_x[N];
	for (int i = 0; i < N; i++)
	{
		prev_x[i] = 0;
	}

	double r[N];
	double prev_r[N];
	double p[N];

	for (int i = 0; i < N; i++)
	{
		prev_r[i] = b[i];

		for (int j = 0; j < N; j++)
		{
			prev_r[i] -= (A[i][j] * prev_x[j]);
		}

		p[i] = prev_r[i];
	}

	double difference;
	do
	{
		double alpha = 0;
		double alpha_denominator = 0;
		for (int i = 0; i < N; i++)
		{
			double Ap = 0;
			for (int j = 0; j < N; j++)
			{
				Ap += (A[i][j] * p[j]);
			}

			alpha += (prev_r[i] * prev_r[i]);
			alpha_denominator += (Ap * p[i]);
		}

		alpha /= alpha_denominator;

		for (int i = 0; i < N; i++)
		{
			x[i] = prev_x[i] + alpha * p[i];
		}

		difference = 0;
		for (int i = 0; i < N; i++)
		{
			if (fabs(x[i] - prev_x[i]) > difference)
			{
				difference = fabs(x[i] - prev_x[i]);
			}
		}

		for (int i = 0; i < N; i++)
		{
			double Ap = 0;
			for (int j = 0; j < N; j++)
			{
				Ap += (A[i][j] * p[j]);
			}

			r[i] = prev_r[i] - alpha * Ap;
		}

		double beta = 0;
		double beta_denominator = 0;
		for (int i = 0; i < N; i++)
		{
			beta += (r[i] * r[i]);
			beta_denominator += (prev_r[i] * prev_r[i]);
		}

		beta /= beta_denominator;

		for (int i = 0; i < N; i++)
		{
			p[i] = r[i] + beta * p[i];
		}

		for (int i = 0; i < N; i++)
		{
			prev_x[i] = x[i];
			prev_r[i] = r[i];
		}
	}while (difference > eps);
}


void shmelevaov::lab8()
{
	double Ja[N][N];
	while (true) 
	{
		int i_max = 0;
		int j_max = 0;
		for (int i = 0; i < N - 1; i++) 
		{
			for (int j = i + 1; j < N; j++) 
			{
				if (fabs(A[i][j]) > fabs(A[i_max][j_max])) 
				{
					i_max = i;
					j_max = j;
				}
			}
		}
		
		if (fabs(A[i_max][j_max]) < eps)
		{
			break;
		}

		double fi = 0.5 * atan(2 * A[i_max][j_max] / (A[i_max][i_max] - A[j_max][j_max]));
		double cos_fi = cos(fi);
		double sin_fi = sin(fi);

		for (int k = 0; k < N; k++) 
		{
			Ja[k][i_max] = A[k][i_max] * cos_fi + A[k][j_max] * sin_fi;
			Ja[k][j_max] = A[k][j_max] * cos_fi - A[k][i_max] * sin_fi;
		}

		for (int i = 0; i < N; i++) 
		{
			for (int j = 0; j < N; j++) 
			{
				if ((j != i_max) && (j != j_max)) 
				{
					Ja[i][j] = A[i][j];
				}
			}
		}

		for (int k = 0; k < N; k++) 
		{
			A[i_max][k] = Ja[i_max][k] * cos_fi + Ja[j_max][k] * sin_fi;
			A[j_max][k] = Ja[j_max][k] * cos_fi - Ja[i_max][k] * sin_fi;
		}

		for (int i = 0; i < N; i++) 
		{
			for (int j = 0; j < N; j++) 
			{
				if ((i != i_max) && (i != j_max)) 
				{
					A[i][j] = Ja[i][j];
				}
			}
		}
	}
	
	for (int i = 0; i < N; i++) 
	{
		x[i] = A[i][i];
	}
}


void shmelevaov::lab9()
{
    double y[N], lambda = 0;
    for (int i = 1; i < N; i++) 
	{
		x[i] = 0;
	}

    x[0] = 1;
    while (true)
    {
        double new_lambda = 0;
        for (int i = 0; i < N; i++)
        {
            y[i] = 0;
            for (int j = 0; j < N; j++)
                y[i] += A[i][j] * x[j];
        }
		
		new_lambda = y[0] / x[0];
        if (fabs(new_lambda - lambda) < eps) 
		{
			break;
		}

        lambda = new_lambda;
        double norma = 0;
        for (int i = 0; i < N; i++) 
		{
			norma += y[i] * y[i];
		}

        norma = sqrt(norma);
        for (int i = 0; i < N; i++) 
		{
			x[i] = y[i] / norma;
		}
    };

    cout << "answer : " << lambda;
}


std::string shmelevaov::get_name()
{
  return "Shmeleva O. V.";
}
