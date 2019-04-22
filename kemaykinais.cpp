#include "kemaykinais.h"

/**
 * Введение в дисциплину
 */
void kemaykinais::lab1()
{
    cout << "Hello world!";
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void kemaykinais::lab2()
{
    double eps=0.00001;
    int lead;
    double maxElem, keep;
    int change[N];
    for (int i=0; i<N; i++){
        lead=i;
        for (int j=i+1; j<N; j++){
            if (abs(A[j][i])>abs(A[lead][i])){
                lead=j;
            }
        }
        if (abs(A[lead][i])<eps){
            cout<<"Нет решений";
            break;
        }
        change[i]=lead;
        if (lead!=i){
            swap(A[i],A[lead]);
            swap(b[i],b[lead]);
        }
        maxElem=A[i][i];
        A[i][i]=1;
        for (int j=i+1; j<N; j++){
            A[i][j]/=maxElem;
        }
        b[i]/=maxElem;
        for (int j=i+1; j<N; j++){
            keep=A[j][i];
            A[j][i]=0;
            for (int k=i+1; k<N; k++){
                A[j][k]-=A[i][k]*keep;
            }
            b[j]-=b[i]*keep;
        }
    }
    for (int i=N-1; i>=0; i--){
        for (int j=i-1; j>=0; j--){
            keep=A[j][i];
            A[j][i]=0;
            for (int k=i-1; k>=0; k--){
                A[j][k]-=A[i][k]*keep;
            }
            b[j]-=b[i]*keep;
        }
    }
    for (int i=0; i<N; i++){
        x[i]=b[i];
    }
    for (int i=N-1; i>=0; i--){
        if (change[i]!=i){
            swap(A[change[i]], A[i]);
            swap(b[change[i]], b[i]);
        }
    }
}



/**
 * Метод прогонки
 */
void kemaykinais::lab3()
{
    double coeff1[N-1], coeff2[N];
    coeff1[0]=-A[0][1]/A[0][0];
    coeff2[0]=b[0]/A[0][0];
    for (int i=1; i<N-1; i++){
        coeff1[i]=-A[i][i+1]/(A[i][i]+A[i][i-1]*coeff1[i-1]);
        coeff2[i]=(b[i]-A[i][i-1]*coeff2[i-1])/(A[i][i]+A[i][i-1]*coeff1[i-1]);
    }
    coeff2[N-1]=(b[N-1]-A[N-1][N-2]*coeff2[N-2])/(A[N-1][N-1]+A[N-1][N-2]*coeff1[N-2]);

    x[N-1]=coeff2[N-1];
    for (int i=N-2; i>=0; i--){
        x[i]=coeff1[i]*x[i+1]+coeff2[i];
    }
}



/**
 * Метод простых итераций
 */
void kemaykinais::lab4()
{
    double eps=1e-20;
    double t=0.0001;
    double y[N];
    double sum;
    double maxRes=1;
    for (;maxRes>eps;){
        for (int i=0; i<N; i++){
            y[i]=x[i];
        }
        for (int i=0; i<N; i++){
            sum=0;
            for (int j=0; j<N; j++){
                sum+=A[i][j]*y[j];
            }
            x[i]=y[i]-t*(sum-b[i]);
        }
        maxRes=abs(x[0]-y[0]);
        for (int i=1; i<N; i++){
            if (abs(x[i]-y[i])>maxRes){
                maxRes=abs(x[i]-y[i]);
            }
        }
    }
}



/**
 * Метод Зейделя
 */
void kemaykinais::lab5()
{
    double eps=1e-20;
    double y[N];
    double sum;
    double maxRes=1;
    for (;maxRes>eps;){
        for (int i=0; i<N; i++){
            y[i]=x[i];
        }
        for (int i=0; i<N; i++){
            sum=0;
            for (int j=0; j<i; j++){
                sum+=A[i][j]*x[j];
            }
            for (int j=i+1; j<N; j++){
                sum+=A[i][j]*y[j];
            }
            x[i]=(1/A[i][i])*(b[i]-sum);
        }
        maxRes=abs(x[0]-y[0]);
        for (int i=1; i<N; i++){
            if (abs(x[i]-y[i])>maxRes){
                maxRes=abs(x[i]-y[i]);
            }
        }
    }
}



/**
 * Метод минимальных невязок
 */
void kemaykinais::lab6()
{
    double eps=pow(10,-16), err;
    double Ax, Ar[N];
    double r0[N];
    double tau, xprev;
    double sum1, sum2;
    do{
        for (int i=0; i<N; i++){
            Ax=0;
            for (int j=0; j<N; j++){
                Ax+=A[i][j]*x[j];
            }
            r0[i]=b[i]-Ax;
        }
        for (int i=0; i<N; i++){
            Ar[i]=0;
            for (int j=0; j<N; j++){
                Ar[i]+=A[i][j]*r0[j];
            }
        }
        sum1=0; sum2=0;
        for (int i=0; i<N; i++){
            sum1+=r0[i]*Ar[i];
            sum2+=Ar[i]*Ar[i];
        }
        tau=sum1/sum2;
        err=0;
        for(int i=0; i<N; i++){
            xprev=x[i];
            x[i]+=tau*r0[i];
            if (abs(x[i]-xprev)>err){
                err=abs(x[i]-xprev);
            }
        }
    }while (err>eps);
}



/**
 * Метод сопряженных градиентов
 */
void kemaykinais::lab7()
{
    double alpha, beta, temp, err, eps=pow(10, -10);
    double sum1, sum2, Ax;
    double r0[N], z0[N], Az[N];
    for (int i=0; i<N; i++){
        Ax=0;
        for (int j=0; j<N; j++){
            Ax+=A[i][j]*x[j];
        }
        r0[i]=b[i]-Ax;
        z0[i]=r0[i];
    }
    do{
        for (int i=0; i<N; i++){
            Az[i]=0;
            for (int j=0; j<N; j++){
                Az[i]+=A[i][j]*z0[j];
            }
        }
        sum1=0; sum2=0;
        for (int i=0; i<N; i++){
            sum1+=r0[i]*r0[i];
            sum2+=Az[i]*z0[i];
        }
        alpha=sum1/sum2;
        err=0;
        for (int i=0; i<N; i++){
            temp = x[i];
            x[i]=x[i]+alpha*z0[i];
            if (abs(temp-x[i]) > err) err = abs(temp-x[i]);
        }
        sum1=0; sum2=0;
        for (int i=0; i<N; i++){
            sum2+=r0[i]*r0[i];
            r0[i]=r0[i]-alpha*Az[i];
            sum1+=r0[i]*r0[i];
        }
        beta=sum1/sum2;
        for (int i=0; i<N; i++){
            z0[i]=r0[i]+beta*z0[i];
        }
    }while(err>eps);
}



/**
* Метод вращения для нахождения собственных значений матрицы
*/
void kemaykinais::lab8()
{

}



/**
* Нахождение наибольшего по модулю собственного значения матрицы
*/
void kemaykinais::lab9()
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


std::string kemaykinais::get_name()
{
  return "Kemaykina Irina";
}
