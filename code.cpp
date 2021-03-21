#include <stdio.h>
#include <iostream>
#include <random>
#include <string.h>
#include <vector>
#include <algorithm>
#include <chrono>
#include <limits>
#include <fstream>
#include <string>
#include <array>

using namespace std;


vector<int> naive_multiply (vector<int>& A, vector<int>& B, int n) {
    vector<int> ans(n*n,0);
    int entry = 0;
    for (int i=0; i<n; i++) {
        for (int k=0; k<n; k++) {
            entry = 0;
            for (int j=0; j<n; j++) {
                entry += (A[j+i*n]) * (B[k+j*n]);
            }
            ans[k+i*n] = entry;
        }
    }
    return ans;
}


void strassen (vector<int>* matrix_list, int depth) {
    // base case 
    if (depth == 0) {
        matrix_list[2][0] = matrix_list[0][0]*matrix_list[1][0];
        return;
    }

    // otherwise

    // on this level the matrices are 2^depth size
    // factors are matrix_list[3*depth], matrix_list[3*depth+1], answer is matrix_list[3*depth+2]
    int n = (int)pow(2, depth);

    // factors on the level below: matrix_list[3*depth-3], matrix_list[3*depth-2]

    // compute p1
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*n/2] = matrix_list[3*depth][j+i*n]; // A
            matrix_list[3*depth-2][j+i*n/2] = matrix_list[3*depth+1][(j+n/2) + i*n] - matrix_list[3*depth+1][(j+n/2)+ (i+n/2)*n]; // F-H
        }
    }
    strassen(matrix_list, depth-1);
    // result of multiplication is matrix_list[3*depth-1]
    // add matrix_list[3*depth-1] to second and fourth quadrant of matrix_list[3*depth+2] 
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][(j+n/2)+i*n] += matrix_list[3*depth-1][j+i*n/2];
            matrix_list[3*depth+2][(j+n/2)+(i+n/2)*n] += matrix_list[3*depth-1][j+i*n/2];
        }
    }


    // compute p2
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*n/2] = matrix_list[3*depth][j+i*n] + matrix_list[3*depth][(j+n/2)+i*n]; // A+B
            matrix_list[3*depth-2][j+i*n/2] = matrix_list[3*depth+1][(j+n/2)+ (i+n/2)*n]; // H
        }
    }
    strassen(matrix_list, depth-1);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][j+i*n] -= matrix_list[3*depth-1][j+i*n/2]; // - first quadrant
            matrix_list[3*depth+2][(j+n/2)+i*n] += matrix_list[3*depth-1][j+i*n/2]; // + second quadrant
        }
    }

    // compute p3
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*n/2] = matrix_list[3*depth][j+(i+n/2)*n] + matrix_list[3*depth][(j+n/2)+(i+n/2)*n]; // C+D
            matrix_list[3*depth-2][j+i*n/2] = matrix_list[3*depth+1][j+i*n]; // E
        }
    }
    strassen(matrix_list, depth-1);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][j+(i+n/2)*n] += matrix_list[3*depth-1][j+i*n/2]; // + third quadrant
            matrix_list[3*depth+2][(j+n/2)+(i+n/2)*n] -= matrix_list[3*depth-1][j+i*n/2]; // - fourth quadrant
        }
    }

    // compute p4
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*n/2] = matrix_list[3*depth][(j+n/2)+(i+n/2)*n]; // D
            matrix_list[3*depth-2][j+i*n/2] = matrix_list[3*depth+1][j+(i+n/2)*n]-matrix_list[3*depth+1][j+i*n]; // G-E
        }
    }
    strassen(matrix_list, depth-1);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][j+i*n] += matrix_list[3*depth-1][j+i*n/2]; // + first quadrant
            matrix_list[3*depth+2][j+(i+n/2)*n] += matrix_list[3*depth-1][j+i*n/2]; // + third quadrant
        }
    }


    // compute p5
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*n/2] = matrix_list[3*depth][j+i*n]+matrix_list[3*depth][(j+n/2)+(i+n/2)*n]; // A+D
            matrix_list[3*depth-2][j+i*n/2] = matrix_list[3*depth+1][j+i*n]+matrix_list[3*depth+1][(j+n/2)+(i+n/2)*n]; // E+H
        }
    }
    strassen(matrix_list, depth-1);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][j+i*n] += matrix_list[3*depth-1][j+i*n/2]; // + first quadrant
            matrix_list[3*depth+2][(j+n/2)+(i+n/2)*n] += matrix_list[3*depth-1][j+i*n/2]; // + fourth quadrant
        }
    }


    // compute p6
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*n/2] = matrix_list[3*depth][(j+n/2)+i*n] - matrix_list[3*depth][(j+n/2)+(i+n/2)*n]; // B-D
            matrix_list[3*depth-2][j+i*n/2] = matrix_list[3*depth+1][j+(i+n/2)*n] + matrix_list[3*depth+1][(j+n/2)+(i+n/2)*n]; // G+H
        }
    }
    strassen(matrix_list, depth-1);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][j+i*n] += matrix_list[3*depth-1][j+i*n/2]; // + first quadrant
        }
    }

    // compute p7
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*n/2] = matrix_list[3*depth][j+i*n] - matrix_list[3*depth][j+(i+n/2)*n]; // A-C
            matrix_list[3*depth-2][j+i*n/2] = matrix_list[3*depth+1][j+i*n] + matrix_list[3*depth+1][(j+n/2)+i*n]; // E+F
        }
    }
    strassen(matrix_list, depth-1);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][(j+n/2)+(i+n/2)*n] -= matrix_list[3*depth-1][j+i*n/2]; // - fourth quadrant
        }
    }
}




//int argc, char** argv
int main() {
    
    int dimension;
    cin >> dimension;
    const int count = ceil(log2(dimension))+1;
    // cout << count;

    // create array of matrices


    vector<int> matrix_list[3*count];
    for (int d=0; d<count; d++) {
        for (int j=0; j<3; j++) {
            matrix_list[3*d+j]=vector<int>((int)pow(2,2*d), 0);
        }
    }

    // cout << matrix_list[0][0] << endl;
    // for(int i=0; i<count; i++) cout << matrix_list[3*i].size() << " ";

    // assign matrix_list[3*count - 3] = M1, matrix_list[3*count - 2] = M2. 
    // Answer will be matrix_list[3*count-1] ..

    // stress test
    srand((unsigned)time(NULL));
    vector<int> bm1(dimension*dimension, 0), bm2(dimension*dimension, 0);
    for (int i=0; i<dimension; i++) {
        for (int j=0; j<dimension; j++) {
            bm1[j+i*dimension] = rand()%2;
            bm2[j+i*dimension] = rand()%2;
        }
    }

    matrix_list[3*count - 3] = bm1;
    matrix_list[3*count - 2] = bm2;

    auto t1 = chrono::high_resolution_clock::now();
    strassen(matrix_list, count-1);

    // for (int i=0; i<dimension; i++) cout << matrix_list[3*count-1][i] << " ";
    // cout << endl;

    auto t2 = chrono::high_resolution_clock::now();
    auto d1 = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
    cout << "strassen took " << d1.count() << " ms | ";

    vector<int> naive = naive_multiply(matrix_list[3*count-3], matrix_list[3*count-2], dimension);
    // for (int i=0; i<dimension; i++) cout << naive[i] << " ";
    // cout << endl;

    auto t3 = chrono::high_resolution_clock::now();
    auto d2 = chrono::duration_cast<chrono::milliseconds>(t3 - t2);
    cout << "naive took " << d2.count() << " ms | ";


    
    // strassen(matrix_list)

    // vector<int> A = {1,2,3,4}, B={0,1,1,0};
    // vector<int> ans = naive_multiply(A,B, 2);
    // for (int i:ans ) cout << i << " " ;


    // for (int i:bm1[0]) cout << i << " ";

    // auto t1 = chrono::high_resolution_clock::now();
    
    // vector<vector<int>> ans = strassen_pure(bm1, bm2, n);

    // auto t2 = chrono::high_resolution_clock::now();
    // auto d1 = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
    // cout << "strassen took " << d1.count() << " ms | ";

    // ans =  naive_multiply(bm1, bm2, n);

    // auto t3 = chrono::high_resolution_clock::now();
    // auto d2 = chrono::duration_cast<chrono::milliseconds>(t3 - t2);
    // cout << "naive took " << d2.count() << " ms | ";
}