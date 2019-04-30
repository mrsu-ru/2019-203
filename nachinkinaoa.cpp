#include "nachinkinaoa.h"

/**
 * Введение в дисциплину
 */
void nachinkinaoa::lab1()
{
   cout <<"Hello word";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void nachinkinaoa::lab2()
{
    double eps=10e-7;
    int num;
    int massiv[N];
    double ved;
    for (int i=0; i<N; i++){
                num=i;
        for (int j=i+1; j<N; j++){
            if(abs(A[j][i])>abs(A[num][i]))
               num=j;}
        if(abs(A[num][i])<eps) {cout << "No solution"; break;}
        if(num!=i){
            swap(A[i],A[num]);
            swap(b[i],b[num]);
        }
        massiv[i]=num;
        ved=A[i][i];//ведущий элемент находится на диагонали
        b[i]/=ved;
        for (int j=i; j<N; A[i][j++]/=ved);//делим по строке на ведущий элемент
        for (int j=i+1; j<N; j++){
                ved = A[j][i];
                for (int k=i; k<N; k++)
                    A[j][k]-=A[i][k]*ved;
                 b[j]-=b[i]*ved;
             }
        for (int j=0; j<i; j++){
                ved = A[j][i];
                for (int k=i; k<N; k++)
                    A[j][k]-=A[i][k]*ved;
                 b[j]-=b[i]*ved;
             }
        }
       for (int i=0; i<N; i++)
          x[i]=b[i];
       for (int i=N-1; i>=0; i--){
          if(massiv[i]!=i){
             swap(A[i], A[massiv[i]]);
             swap(b[i], b[massiv[i]]);}
        }

}



/**
 * Метод прогонки
 */
void nachinkinaoa::lab3()
{
     double *alfa = new double[N-1];
       double *beta = new double[N];
       alfa[0]=-A[0][1]/A[0][0];
       beta[0]=b[0]/A[0][0];
       double gamma;
       for (int i=1; i<N; i++){
        gamma=A[i][i]+alfa[i-1]*A[i][i-1];
        if(i!=N-1)
            alfa[i]=-A[i][i+1]/gamma;
        beta[i]=(b[i]-A[i][i-1]*beta[i-1])/gamma;
       }
       x[N-1]=beta[N-1];
       for (int i=N-2; i>=0; i--)
        x[i]=alfa[i]*x[i+1]+beta[i];
        delete[] alfa;
	    delete[] beta;

}



/**
 * Метод простых итераций
 */
void nachinkinaoa::lab4()
{
 for (int i=0; i<N; i++)
              x[i]=0;
       double eps=10e-15;
       double tay=0.01;
       double var=0;
       for(;;){
             double differ=0;
         for (int i=0; i<N; i++){
            double ay=0;
            var=x[i];
             for (int j=0; j<N; j++)
              ay+=A[i][j]*x[j];
            x[i]-=tay*(ay-b[i]);// новые значения на новой итерации
          if(abs(x[i]-var)>differ)
                differ=abs(x[i]-var);
         }
        if(differ<eps) break;
       }
}



/**
 * Метод Якоби или Зейделя
 */
void nachinkinaoa::lab5()
{
    for (int i=0; i<N; i++)
              x[i]=0;
       double eps=10e-15;
       double var=0;
       for(;;){
          double differ=0;
         for (int i=0; i<N; i++){
            double ay=0;
            var=x[i];
            for (int j=0; j<i; j++)
              ay+=A[i][j]*x[j];
            for (int j=i+1; j<N; j++)
              ay+=A[i][j]*x[j];
            x[i]=(b[i]-ay)/A[i][i];// новые значения на новой итерации
          if(abs(x[i]-var)>differ)
                differ=abs(x[i]-var);
         }
        if(differ<eps) break;
       }

}



/**
 * Метод минимальных невязок
 */
void nachinkinaoa::lab6()
{
    double *r=new double[N];
    for (int i=0; i<N; i++)
              r[i]=0;
    double eps=10e-16, var=0, tay=0;
       for(;;){
            double differ=0, sum1=0, sum2=0;
            for (int i=0; i<N; i++){
                r[i]=b[i];
                for (int j=0; j<N; j++)
                    r[i]-=A[i][j]*x[j];//вектор невязки
            }
            for (int i=0; i<N; i++){
                double vec=0;
                for (int k=0; k<N; k++)
                     vec+=A[i][k]*r[k];
                sum1+=r[i]*vec;
                sum2+=vec*vec;
            }
            tay=sum1/sum2;
            for (int i=0; i<N; i++){
                var=x[i];
                x[i]+=r[i]*tay;// новые значения на новой итерации
                if(abs(x[i]-var)>differ)
                    differ=abs(x[i]-var);
            }
            if(differ<eps) break;
            }
            delete[] r;
}



/**
 * Метод сопряженных градиентов
 */
void nachinkinaoa::lab7()
{
       for (int i=0; i<N; i++)
              x[i]=0;
       double *r=new double[N];
       for (int i=0; i<N; i++){
            r[i]=b[i];
         for (int j=0; j<N; j++)
            r[i]-=A[i][j]*x[j];
       }
       double *z=new double[N];
         for (int i=0; i<N; i++)
            z[i]=r[i];
       double eps=10e-16, var=0, alfa=0;
       for(;;){
            double differ=0, sum1=0, sum2=0, vec=0;
             for (int i=0; i<N; i++){
                vec=0;
                for (int k=0; k<N; k++)
                     vec+=A[i][k]*z[k];
                sum1+=r[i]*r[i];
                sum2+=vec*z[i];
            }
            alfa=sum1/sum2;
             for (int i=0; i<N; i++){
                    var=x[i];
            x[i]+=alfa*z[i];
                if(abs(x[i]-var)>differ)
                    differ=abs(x[i]-var);
                            }
            if(differ<eps) break;
            sum2=0;
            sum1=0;
             for (int i=0; i<N; i++){
                  vec=0;
                sum2+=r[i]*r[i];
                for (int j=0; j<N; j++)
                 vec+=A[i][j]*z[j];
                 r[i]-=alfa*vec;
                 sum1+=r[i]*r[i];

              }
            for (int i=0; i<N; i++){
                sum1+=r[i]*r[i];
                sum2+=vec*z[i];
            }
            for (int i=0; i<N; i++)
            z[i]=r[i]+sum1*z[i]/sum2;
      }
      delete[] r;
      delete[] z;

}


void nachinkinaoa::lab8()
{
       double **B=new double*[N];
       for (int i=0; i<N; i++)
              B[i]=new double [N];
       double eps=1.e-10;
       int i_var, j_var;
       for(;;){
              i_var=0;
              j_var=1;
              double n=0;
              for (int i=0; i<N-1; i++)
                  for (int j=i+1; j<N; j++){
                    if(abs(A[i][j])>abs(A[i_var][j_var])){
                       i_var=i;
                       j_var=j;
                    }
                    n+=A[i][j]*A[i][j];
                    }
              if(sqrt(n)<eps) break;
            double fi=0.5*atan(2*A[i_var][j_var]/(A[i_var][i_var]-A[j_var][j_var]));
            for (int i=0; i<N; i++){
                 B[i][i_var]= A[i][i_var]*cos(fi)+A[i][j_var]*sin(fi);
                 B[i][j_var]=-A[i][i_var]*sin(fi)+A[i][j_var]*cos(fi);
            }
            for (int i=0; i<N; i++)
                  for (int j=0; j<N; j++)
                    if(j!=j_var&&j!=i_var) B[i][j]=A[i][j];
            for (int i=0; i<N; i++){
                 A[i_var][i]= B[i_var][i]*cos(fi)+B[j_var][i]*sin(fi);
                 A[j_var][i]=-B[i_var][i]*sin(fi)+B[j_var][i]*cos(fi);
            }
            for (int i=0; i<N; i++)
                  for (int j=0; j<N; j++)
                    if(i!=j_var&&i!=i_var) A[i][j]=B[i][j];
            }
       for (int i=0; i<N; i++)
          x[i]=A[i][i];

        for (int i=0; i<N; i++)
              delete[] B[i];
      delete[] B;
}


void nachinkinaoa::lab9()
{
       for (int i=0; i<N; i++)
              x[i]=0;
              x[0]=1;
       double *y=new double[N];
       double eps=10e-16, pr_lyamba, lyamba=0, sum;
       for(;;){
            sum=0;
            pr_lyamba=lyamba;
            lyamba=0;
            for (int i=0; i<N; i++){
                y[i]=0;
                for (int k=0; k<N; k++)
                     y[i]+=A[i][k]*x[k];
                lyamba+=y[i]*x[i];
                sum+=y[i]*y[i];
            }
           sum=pow(sum,0.5);
             for (int i=0; i<N; i++)
                x[i]=y[i]/sum;
             if(abs(lyamba-pr_lyamba)<eps) break;
            }
            x[0] = lyamba;
       delete[] y;
}


std::string nachinkinaoa::get_name()
{
  return "Начинкина Олеся Алексеевна";
}
