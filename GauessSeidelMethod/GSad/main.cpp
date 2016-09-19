/*Asim Amjad
Gauss Seidel Iteration Project
 Assignment 3 (Part II)

*/

#include<iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <math.h>
#include <random>
#include <ctime>

using namespace std;


int n =676;



static double a[676][676];
static  double temp[676];
static double  FinalMatrix[676];
static double b[676]; // for random solution of B

int main(int argc, char **argv)
{



    // for random numbers
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();

//    cout<<"Please Enter Size of Matrix? ";
//    cin>>n;
//    cout<<endl;
    myclock::duration d = myclock::now() - beginning;

    unsigned seed2 = d.count();

   minstd_rand0 generator (seed2);

   uniform_int_distribution<int> distribution(10,100);

/// to print out for original matrix
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
cout<<endl;

  /// Make Diagonal Dominant Matrix


    for(int i = 0; i < n; i++){
        double dValue = 0;
        for(int k = 0; k < n; k++){
            if(a[i][k] < 0)
                dValue += abs(a[i][k]);
            else
                dValue += a[i][k];
        }
        a[i][i] = dValue;
    }
  cout<<"Diagonal Dominant Matrix "<<endl;
  for(int r=0; r<n; r++)
  {
    for(int c=0; c<n; c++)
    {

        cout<<a[r][c]<<"  ";
    }
    cout<<endl;
  }
cout<<endl;

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

cout<<endl;



/// We need to make a imaginary matrix which has all values 1
cout<<"Imaginary Matrix by giving 1 Value to all columns "<<endl;
for(int r=0; r<n; r++)
{

    temp[r]=1;
    cout<<temp[r]<<endl;

}

cout<<endl;

/// Now we need to use this formula to apply Gauss Seidel Method
/// X(new) = D^-1 ( b - L X (new) - U X(old) )
/// We also have to tell how many iteration we want but we put only one iteration
int iterations=1;

while (iterations > 0)
    {
        for (int i = 0; i < n; i++)
        {
            FinalMatrix[i] = (b[i] / a[i][i]);
            for (int j = 0; j < n; j++)
            {
                if (j == i)
                    continue;
                FinalMatrix[i] = FinalMatrix[i] - ((a[i][j] / a[i][i]) * temp[j]);
                temp[i] = FinalMatrix[i];
            }

            cout<<"X"<<i+1<<" = "<<FinalMatrix[i]<<endl;

        }
        cout << "\n";
        iterations--;
    }

return 0;
}  /// end of main


