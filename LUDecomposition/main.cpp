/*Asim Amjad
 LU Decomposition Project
 Assignment 2

*/

#include<iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <math.h>
#include <random>
#include <ctime>

using namespace std;

int n =484;

int main(int argc, char **argv)
{
    // function to make LU Decomposition
    void lu(double[][100], double[][100], double[][100], int n);
    //  function to print out for two dimensional matrix
    void output(double[][100], int);
    //  function to print out for one dimensional matrix
    void output2(double[100], int);
    // function to compute LY=B
    void LYCompute(double [][100], double [100], double [100], int n);
    // function to compute UY=X
    void UXCompute (double [][100], double [100], double [100], int n);
    // function to compute error
    void ActualX(double [100], double [100],double [100][100],double [100], int n);

    // function to calculate error between actual and imaginary
    void ErrorTwo(double [100],double [100], double [100]);

    // function to calculate error
    double CalError(double [100]);

    // declaration of arrays a, l and u
    double a[100][100], l[100][100], u[100][100];  // a=original matrix, l =L matrix, u = U matrix
    double b[100]; // for random solution of B
    double YMatrix[100]; // for Y values
    double XMatrix[100];  // for X values
    double AM[100]; // for actual matrix
    double EM[100];// error matrix

    int  i = 0, j = 0; // integers

    // for random numbers
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
//
//    cout<<"Please Enter Size of Matrix? ";
//    cin>>n;
//    cout<<endl;
    myclock::duration d = myclock::now() - beginning;

    unsigned seed2 = d.count();

   minstd_rand0 generator (seed2);
   uniform_int_distribution<int> distribution(-9999999,9999999);
// to print out for original matrix
  cout<<"Original Random Matrix "<<endl;
  for(int r=0; r<n; r++)
  {
    for(int c=0; c<n; c++)
    {
        a[r][c]=double(distribution(generator))/10000.0;
        cout<<a[r][c]<<"  ";
    }
    cout<<endl;
  }

     cout<<"Random Matrix for B "<<endl;
  for(int r=0; r<n; r++)
  {
    for(int c=0; c<1; c++)
    {
        b[r]=double(distribution(generator))/100000.0;
        cout<<b[r]<<"  ";
    }
    cout<<endl;
  }
    // function call to make LU Decomposition
    lu(a, l, u, n);

    // print for L Decomposition Matrix with array l
    cout << "\nL Decomposition\n\n";
    output(l, n);

    // print for U Decomposition Matrix with array u
    cout << "\nU Decomposition\n\n";
    output(u, n);
    // call function LYCompute
    LYCompute(l, b,YMatrix, n);
     //print out Y matrix
//     cout<<endl;
//     cout<<"Y Matrix "<<endl;
//     output2(YMatrix,n);

    // function to call UXCompute
     UXCompute(u,YMatrix,XMatrix,n);

    // print out X Matrix
     cout<<"X Matrix "<<endl;
     output2(XMatrix, n);
     cout<<endl;
     // call the function Actual matrix
     ActualX(XMatrix,b,a,AM,n);
     // print out Error Matrix
     cout<<"Actual Matrix "<<endl;
     output2(AM, n);
     cout<<endl;


    // call the error matrix
    ErrorTwo(AM,b, EM);
    // print error matrix
    cout<<"Error Matrix"<<endl;
    output2(EM,n);

    // call callerror function
    double TotalError =CalError(EM);
    cout<<"Error Vector "<<TotalError<<endl;
    // Two Norm Form
    double TwoNorm = TotalError/n;
    cout<<fixed << showpoint;
    cout << setprecision(20);
    cout << "Two Norm Form:" << TwoNorm;
    cout<< endl;

    return 0;
}

// function method for making LU Decomposition
void lu(double a[][100], double l[][100], double u[][100], int n)
{
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (j < i)
                l[j][i] = 0;
            else
            {
                l[j][i] = a[j][i];
                for (k = 0; k < i; k++)
                {
                    l[j][i] = l[j][i] - l[j][k] * u[k][i];
                }
            }
        }
        for (j = 0; j < n; j++)
        {
            if (j < i)
                u[i][j] = 0;
            else if (j == i)
                u[i][j] = 1;
            else
            {
                u[i][j] = a[i][j] / l[i][i];
                for (k = 0; k < i; k++)
                {
                    u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
                }
            }
        }
    }
}

// print out for L and U Decomposition Matrices
void output(double x[][100], int n)
{
    int i = 0, j = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            cout<<x[i][j]<<" ";
        }
        cout << "\n";
    }
}

void output2(double x[100], int n)
{
    int i = 0;
    for (i = 0; i < n; i++)
    {
        cout<<x[i];

        cout << "\n";
    }
}


// forward substitution for LY=B
void LYCompute(double L[][100], double B[100], double Y[100], int n)
{
 for (int i=0; i<n; i++)
 {
     double sum=0;
     double value=L[i][i];
     for(int j=0; j<n; j++)
     {
         double temp = L[i][j]*Y[j];
         sum +=temp;

     }

     Y[i] = (B[i] - sum)/value;
 }

}


// Backward substitution for Y=UX


void UXCompute (double U[][100], double Y[100], double X[100], int n)

{


    for(int i = n-1; i >=0; i--){
        if(i == (n-1))
            X[n-1] = Y[n-1];
        else{
            double value = Y[i];
            double temp = 0;
            for(int j = n-1; j > 0; j--){
                temp -= U[i][j] * X[j];
            }
            X[i] = temp + value;
        }
    }

}

// Make error vector


void ActualX(double X[100], double B[100],double a[100][100],double AMatrix[100], int n)
{

    for(int i = 0; i < n; i++){
        double sum = 0;
        for(int j = 0; j < n; j++){
            double temp = a[i][j]*X[j];
            sum += temp;
        }
        AMatrix[i] = sum;
    }


}

// calculate error by function
void ErrorTwo(double A[100],double X[100], double E[100])
{
 for(int i = 0; i < n; i++)


        E[i] = A[i] - X[i];
}
//calculate error matrix
double CalError(double E[100])
{

    double sum = 0;
    for(int i = 0; i < n; i++ ){
        sum += sqrt( pow(E[i],2) );
    }

return sum;
}

