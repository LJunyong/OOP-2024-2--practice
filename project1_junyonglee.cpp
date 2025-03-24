#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

/*
************************************************************************************************

윈도우 사이즈 적절히 reasonable 하게 키워서 구현해도 된다고 하셔서 setw(11)을 사용했습니다.
output case의 solution만 맞으면 된다고 하셔서 Gaussian Elimination result의 형태가 다소 다릅니다. 
하지만 결론적으로 output case와 같은 행렬입니다.

*************************************************************************************************
*/


// print the matrix, use a window size of 3 and right align
void printMatrix(double** a, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(a[i][j]) < 1e-10) {
                a[i][j] = 0;
            }
            if (j == 0)
                cout << ' ' << a[i][j];
            else
                cout << right << setw(11) << a[i][j];
        }
        cout << right << setw(11) << a[i][n] << endl;
    }
}

// Swap two rows of a matrix
void swapRows(double** A, int a, int b, int n) {
    for (int j = 0; j <= n; j++) {
        double temp = A[a][j];
        A[a][j] = A[b][j];
        A[b][j] = temp;
    }
}

void performGaussianElimination(double** A, int n) {
    for (int i = 0; i < n; i++) {
        if (A[i][i] == 0) {
            for (int j = i + 1; j < n; j++) {
                if (A[j][i] != 0) {
                    swapRows(A, i, j, n);
                    break;
                }
            }
        }
        double m = A[i][i];
        for (int j = 0; j <= n; j++) {
            A[i][j] = (A[i][j] / m);
        }
        for (int h = i + 1; h < n; h++) {
            double c = A[h][i];
            for (int j = i; j <= n; j++) {
                A[h][j] -= c * A[i][j];
            }
        }
    }
    cout << "Gaussian Elimination result is:" << endl;
    printMatrix(A, n);
}
void backSubstitution(double** M, int n, double* s) {
    double x, y, z;
    z = M[2][3] / M[2][2];
    y = (M[1][3] - (M[1][2]*z)) / M[1][1];
    x = (M[0][3] - (M[0][2] * z) - (M[0][1] * y)) / M[0][0];
    s[0] = x;
    s[1] = y;
    s[2] = z;
    for (int i = 0; i < n; i++) {
        if (fabs(s[i]) < 1e-10) {
            s[i] = 0;
        }
    }
}

void solveSystem(double** matrix, int n) {
    performGaussianElimination(matrix, n);
    double* solution = new double[n];
    backSubstitution(matrix, n, solution);

    std::cout << "The solution to the system of linear equations is:" << std::endl;
    std::cout << std::setprecision(6);
    for (int i = 0; i < n; i++) {
        std::cout << "x[" << i << "] = " << solution[i] << std::endl;
    }
    std::cout << std::endl;
    delete[] solution;
}



int main() {
    // create the augmented matrix
    int n = 3; // size of the matrix
    double** A = new double* [n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n + 1];
    }

    // initialize the matrix with values (does not require row exchange)
    A[0][0] = 2; A[0][1] = 3; A[0][2] = -1; A[0][3] = 1;
    A[1][0] = 4; A[1][1] = 4; A[1][2] = 3; A[1][3] = 3;
    A[2][0] = 2; A[2][1] = -3; A[2][2] = 1; A[2][3] = -1;

    // print original matrix
    std::cout << "The original matrix (the last column is augmented) is:" << std::endl;
    printMatrix(A, n);

    // solve the system of linear equations
    solveSystem(A, n);

    // initialize the matrix with values (requires row exchange)
    A[0][0] = 0; A[0][1] = 3; A[0][2] = -1; A[0][3] = 1;
    A[1][0] = 4; A[1][1] = 4; A[1][2] = 3; A[1][3] = 3;
    A[2][0] = 2; A[2][1] = -3; A[2][2] = 1; A[2][3] = -1;

    // print original matrix
    std::cout << "The original matrix (the last column is augmented) is:" << std::endl;
    printMatrix(A, n);

    // solve the system of linear equations
    solveSystem(A, n);

    // free the memory used by the matrix
    for (int i = 0; i < n; i++) {
        delete[] A[i];
    }
    delete[] A;

    return 0;
}
