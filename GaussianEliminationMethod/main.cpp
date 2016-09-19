/*
Asim Amjad
Cs 417
Assignment 1
Make a program which will calculate the Gaussian Elimination Method with Diagonal Dominance
Then Calculate the error with original and estimated value and lastly calculate Two Norm Method.
*/

#include <chrono>
#include <iomanip>
#include <vector>
#include <math.h>
#include <iostream>
#include <random>
#include <ctime>


using namespace std;


// Create a matrix
void MatrixMod(vector <vector<double> >&V2, vector<double> V4){
    int n = V2.size();
    for(int i = 0; i < n; i++){
        V2[i][n] = V4[i];
    }

}
// Make a Diagonal Dominance
void DomDia(vector <vector<double> > &V2){
    int n = V2.size();

    for(int i = 0; i < n; i++){
        double dValue = 0;
        for(int k = 0; k < n; k++){
            if(V2[i][k] < 0)
                dValue += abs(V2[i][k]);
            else
                dValue += V2[i][k];
        }
        V2[i][i] = dValue;
    }
}

// Random Solution
vector<double> V1(vector <vector<double> >V2, vector< vector<double> >V3){
    int n = V2.size();
    vector<double> V4(n);

    for(int i = 0; i < n; i++){
        double sum = 0;
        for(int j = 0; j < n; j++){
            double dValue = V2[i][j]*V3[j][0];
            sum += dValue;
        }
        V4[i] = sum;

    }
    return V4;
}


// function for echelon form
void ECF(vector <vector<double> > &V2)
{
    int n = V2.size();

    for(int i = 0; i < n; i++){
        double Value = abs(V2[i][i]);
        int max_row = i;
        for(int j = i+1; j < n; j++){
            if((abs(V2[j][i]) > Value)){
                Value = abs(V2[j][i]);
                max_row = j;
            }
        }

        for(int j = i; j < n+1; j++){
            double Max = V2[max_row][j];
            V2[max_row][j] = V2[i][j];
            V2[i][j] = Max;
        }
        for(int k = i+1; k < n; k++){
            double dValue = -V2[k][i]/V2[i][i];
            for(int j = i; j < n+1; j++){
                if(i == j){
                    V2[k][j] = 0;
                }
                else{
                    V2[k][j] += dValue * V2[i][j];
                }
            }
        }
        for(int k = i; k < n; k++){
            double Max = V2[i][k];
                for(int j = i; j < n+1; j++){
                    if(k == i){
                        V2[i][j] /= Max;
                    }
                }
        }

    }

}
 // for X solution
vector<double> XSolu(vector< vector<double> > &V2){

    int n = V2.size();
    vector<double> X(n);
    for(int i = n-1; i >= 0; i--){
        double dValue = V2[i][n]/V2[i][i];
        X[i] = dValue;
        for(int j = i-1; j >= 0; j--){
            V2[j][n] -= V2[j][i] * X[i] ;
        }
    }
    return X;


}

// for B solution

vector<double> BSol(vector< vector<double> > C, vector<double> X){
    int n = C.size();
    vector<double> M(n);

        for(int i = 0; i < n; i++){
        double sum = 0;
        for(int j = 0; j < n; j++){
            double dValue = C[i][j]*X[j];

            sum += dValue;

        }
        M[i] = sum;
    }

    return M;
}


// for error
vector<double> Error(vector <double>V4, vector <double> M, int n){

    vector<double> E(n);
    for(int i = 0; i < n; i ++){
        double dValue;
        dValue += V4[i] - M[i];
        E[i]= dValue;

    }


    return E;
}

// for two norm
double Norm(vector<double> E){
    int n = E.size();
    double sum = 0;
    for(int i = 0; i < n; i++ ){
        sum += (sqrt( pow(E[i],2))/n );

    }

return sum;
}




void Print1(vector< vector<double> > V2){
    cout << left;
    int n = V2.size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            cout << setw(12) << V2[i][j];
        }
        cout<< endl;
    }
}


void Print2(vector< vector<double> > V3){
    int n = V3.size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < 1; j++){
            cout << V3[i][j] << "  ";
        }
        cout<< endl;
    }
}

void Print3(vector<double> V4){
    int n = V4.size();
    for(int i = 0; i < n; i++){
        cout << V4[i] << "\n";
    }
    cout<< endl;
}


void Print4(vector< vector<double> > &V2){
    cout<<fixed << showpoint;
    cout << setprecision(6);
    int n = V2.size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n + 1; j++){
            cout << V2[i][j] << "\t";
            if(j == n - 1)
                cout << "|  ";
        }
        cout<< endl;
    }
    cout << endl;
}

void Print5(vector<double> X){
    cout<<fixed << showpoint;
    cout << setprecision(4);
    int n = X.size();
    for(int i = 0; i < n; i++){
        cout << X[i] << "\n";
    }
    cout<< endl;
}


void Print6(vector< vector<double> > C){
    cout << left;
    int n = C.size();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            cout << setw(10) << C[i][j];
        }
        cout << endl;
    }

}


void print_M(vector<double> M){
    int n = M.size();
    for(int i = 0; i < n; i++){
        cout << M[i] << "\n";
    }
    cout<< endl;

}

void printE(vector< double> E){
    int n = E.size();
    for(int i = 0; i < n; i++){
        cout << E[i] << "\n";
    }
    cout<< endl;

}

/*           End of Functions   */

int main()
{

    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();

    int n;
    cout<<"Please enter size of matrix ";
    cin>>n;


    cout<<endl;

    myclock::duration d = myclock::now() - beginning;
    unsigned seed2 = d.count();
    minstd_rand0 generator(seed2);
    uniform_int_distribution<int> distribution(-9999999,9999999);

    vector< vector<double> > V2(n);
    for(int r = 0; r < n; r++){
        for(int c = 0; c < n; c++){
            V2[r].push_back(r*c);
            V2[r][c]=double(distribution(generator))/10000.0;

        }
    }

    uniform_int_distribution<int> distrib(-20,20);
    vector< vector<double> > V3(n);

    for(int r=0; r < n; r++){
        for(int c = 0; c < 1; c++){
            V3[r].push_back(r*c);
            V3[r][c]= double(distrib(generator))/100.000;
        }
   }


    cout<<"Random Matrix: " << endl;
    Print1(V2);
    cout << endl;

    DomDia(V2);
    cout<<"Diagonal Dominance matrix:" << endl;
    Print1(V2);
    cout<< endl;
    vector< vector<double> >C(V2);
    cout<<"Random Vector Y:" << endl;
    Print2(V3);
    cout<<endl;

    vector<double>V4 = V1(V2, V3);
    cout << "Random Solution for B:" << endl;
    Print3(V4);
    cout << endl;

    MatrixMod(V2, V4);


    ECF(V2);
    cout<<"Row Echelon Form:" << endl;
    Print4(V2);

    // XSolu is for back solving
    vector<double> X = XSolu(V2);
    cout<<"Back Solving:" << endl;
    Print5(X);
    cout << endl;

    cout<<"Matrix Values "<<endl;
    Print4(V2);
    cout << endl;

    cout<<"Original Matrix X:" << endl;
    vector<double>M = BSol(C,X);


    print_M(M);
    cout << endl;

    vector<double> E = Error(V4, M, n);
    cout << "Error Vector: " << endl;
    printE(E);
    cout << endl;

    //Two Norm
    double norm= Norm(E);
    cout << "Two Norm Form:"<<endl;
    cout << norm;
    cout << endl;

    return 0;
}
