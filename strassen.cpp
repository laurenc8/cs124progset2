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
    for (int i=0; i<n; i++) {
        for (int k=0; k<n; k++) {
            for (int j=0; j<n; j++) {
                ans[k+i*n] += (A[j+i*n]) * (B[k+j*n]);
            }
        }
    }
    return ans;
}

void naive_multiply_recursion (vector<int>* matrix_list, int depth, vector<int>& sizes) {
    int n = sizes[depth]; 
    for (int &i: matrix_list[3*depth+2]) i=0;
    for (int i=0; i<n; i++) {
        for (int k=0; k<n; k++) {
            for (int j=0; j<n; j++) {
                matrix_list[3*depth+2][k+i*n] += (matrix_list[3*depth][j+i*n]) * (matrix_list[3*depth+1][k+j*n]);
            }
        }
    }
}

const int n0 = 100;

void strassen (vector<int>* matrix_list, int depth, vector<int>& sizes) {
    // base case 
    if (sizes[depth] <n0) {
        naive_multiply_recursion(matrix_list, depth, sizes); 
        return;
    }

    // factors are matrix_list[3*depth], matrix_list[3*depth+1], answer is matrix_list[3*depth+2]
    int n = sizes[depth], m =sizes[depth-1];

    // factors on the level below: matrix_list[3*depth-3], matrix_list[3*depth-2]

    // compute p1
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*m] = matrix_list[3*depth][j+i*n]; // A
            matrix_list[3*depth-2][j+i*m] = matrix_list[3*depth+1][(j+n/2) + i*n] - matrix_list[3*depth+1][(j+n/2)+ (i+n/2)*n]; // F-H
        }
    }
    strassen(matrix_list, depth-1, sizes);
    // result of multiplication is matrix_list[3*depth-1]
    // add matrix_list[3*depth-1] to second and fourth quadrant of matrix_list[3*depth+2] 
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][(j+n/2)+i*n] = matrix_list[3*depth-1][j+i*m]; // + second quadrant
            matrix_list[3*depth+2][(j+n/2)+(i+n/2)*n] = matrix_list[3*depth-1][j+i*m]; // + fourth quadrant
        }
    }


    // compute p2
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*m] = matrix_list[3*depth][j+i*n] + matrix_list[3*depth][(j+n/2)+i*n]; // A+B
            matrix_list[3*depth-2][j+i*m] = matrix_list[3*depth+1][(j+n/2)+ (i+n/2)*n]; // H
        }
    }
    strassen(matrix_list, depth-1, sizes);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][j+i*n] = - matrix_list[3*depth-1][j+i*m]; // - first quadrant
            matrix_list[3*depth+2][(j+n/2)+i*n] += matrix_list[3*depth-1][j+i*m]; // + second quadrant
        }
    }

    // compute p3
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*m] = matrix_list[3*depth][j+(i+n/2)*n] + matrix_list[3*depth][(j+n/2)+(i+n/2)*n]; // C+D
            matrix_list[3*depth-2][j+i*m] = matrix_list[3*depth+1][j+i*n]; // E
        }
    }
    strassen(matrix_list, depth-1, sizes);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][j+(i+n/2)*n] = matrix_list[3*depth-1][j+i*m]; // + third quadrant
            matrix_list[3*depth+2][(j+n/2)+(i+n/2)*n] -= matrix_list[3*depth-1][j+i*m]; // - fourth quadrant
        }
    }

    // compute p4
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*m] = matrix_list[3*depth][(j+n/2)+(i+n/2)*n]; // D
            matrix_list[3*depth-2][j+i*m] = matrix_list[3*depth+1][j+(i+n/2)*n]-matrix_list[3*depth+1][j+i*n]; // G-E
        }
    }
    strassen(matrix_list, depth-1, sizes);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][j+i*n] += matrix_list[3*depth-1][j+i*m]; // + first quadrant
            matrix_list[3*depth+2][j+(i+n/2)*n] += matrix_list[3*depth-1][j+i*m]; // + third quadrant
        }
    }


    // compute p5
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*m] = matrix_list[3*depth][j+i*n]+matrix_list[3*depth][(j+n/2)+(i+n/2)*n]; // A+D
            matrix_list[3*depth-2][j+i*m] = matrix_list[3*depth+1][j+i*n]+matrix_list[3*depth+1][(j+n/2)+(i+n/2)*n]; // E+H
        }
    }
    strassen(matrix_list, depth-1, sizes);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][j+i*n] += matrix_list[3*depth-1][j+i*m]; // + first quadrant
            matrix_list[3*depth+2][(j+n/2)+(i+n/2)*n] += matrix_list[3*depth-1][j+i*m]; // + fourth quadrant
        }
    }


    // compute p6
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*m] = matrix_list[3*depth][(j+n/2)+i*n] - matrix_list[3*depth][(j+n/2)+(i+n/2)*n]; // B-D
            matrix_list[3*depth-2][j+i*m] = matrix_list[3*depth+1][j+(i+n/2)*n] + matrix_list[3*depth+1][(j+n/2)+(i+n/2)*n]; // G+H
        }
    }
    strassen(matrix_list, depth-1, sizes);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][j+i*n] += matrix_list[3*depth-1][j+i*m]; // + first quadrant
        }
    }

    // compute p7
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*m] = matrix_list[3*depth][j+i*n] - matrix_list[3*depth][j+(i+n/2)*n]; // A-C
            matrix_list[3*depth-2][j+i*m] = matrix_list[3*depth+1][j+i*n] + matrix_list[3*depth+1][(j+n/2)+i*n]; // E+F
        }
    }
    strassen(matrix_list, depth-1, sizes);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][(j+n/2)+(i+n/2)*n] -= matrix_list[3*depth-1][j+i*m]; // - fourth quadrant
        }
    }
}

void strassen_time (vector<int>* matrix_list, int depth, vector<int> &sizes) {
    // on this level the matrices are 2^depth size
    // factors are matrix_list[3*depth], matrix_list[3*depth+1], answer is matrix_list[3*depth+2]
    int n = sizes[depth];

    // factors on the level below: matrix_list[3*depth-3], matrix_list[3*depth-2]

    // compute p1
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*n/2] = matrix_list[3*depth][j+i*n]; // A
            matrix_list[3*depth-2][j+i*n/2] = matrix_list[3*depth+1][(j+n/2) + i*n] - matrix_list[3*depth+1][(j+n/2)+ (i+n/2)*n]; // F-H
        }
    }
    naive_multiply_recursion(matrix_list, depth-1, sizes);
    // result of multiplication is matrix_list[3*depth-1]
    // add matrix_list[3*depth-1] to second and fourth quadrant of matrix_list[3*depth+2] 
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][(j+n/2)+i*n] = matrix_list[3*depth-1][j+i*n/2]; // + second quadrant
            matrix_list[3*depth+2][(j+n/2)+(i+n/2)*n] = matrix_list[3*depth-1][j+i*n/2]; // + fourth quadrant
        }
    }


    // compute p2
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth-3][j+i*n/2] = matrix_list[3*depth][j+i*n] + matrix_list[3*depth][(j+n/2)+i*n]; // A+B
            matrix_list[3*depth-2][j+i*n/2] = matrix_list[3*depth+1][(j+n/2)+ (i+n/2)*n]; // H
        }
    }
    naive_multiply_recursion(matrix_list, depth-1, sizes);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][j+i*n] = - matrix_list[3*depth-1][j+i*n/2]; // - first quadrant
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
    naive_multiply_recursion(matrix_list, depth-1, sizes);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][j+(i+n/2)*n] = matrix_list[3*depth-1][j+i*n/2]; // + third quadrant
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
    naive_multiply_recursion(matrix_list, depth-1, sizes);
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
    naive_multiply_recursion(matrix_list, depth-1, sizes);
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
    naive_multiply_recursion(matrix_list, depth-1, sizes);
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
    naive_multiply_recursion(matrix_list, depth-1, sizes);
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            matrix_list[3*depth+2][(j+n/2)+(i+n/2)*n] -= matrix_list[3*depth-1][j+i*n/2]; // - fourth quadrant
        }
    }
}



//
int main(int argc, char** argv) {
    if (argc < 2 || argc > 4) {
        cout << "usage : ./strassen flag dimension inputfile" << endl;
        return 0;
    }

    // non-negative flag will use inputfile, negative flag will use RNG
    int dimension = stoi(argv[2]), flag = stoi(argv[1]);

    // timing
    if (flag == -10) {
        for (int n=1; n<=dimension; n++) {
            // generate 2 random matrices
            srand((unsigned)time(NULL));

            vector<int> sizes(2);
            sizes[1] = n + n%2;
            sizes[0] = sizes[1]/2;
            sizes[0] += sizes[0]%2;

            vector<int> matrix_list[6];
            for (int d=0; d<2; d++) {
                for (int j=0; j<3; j++) {
                    matrix_list[3*d+j] = vector<int>(sizes[d]*sizes[d], 0);
                }
            }

            for (int i=0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    matrix_list[3][j+i*n] = rand() % 2;
                    matrix_list[4][j+i*n] = rand() % 2;
                }
            }
            
            // do the timing and compare
            auto t1 = chrono::high_resolution_clock::now();
            for (int i=0; i<10; i++) {
                strassen_time(matrix_list, 1, sizes);
            }            
            auto t2 = chrono::high_resolution_clock::now();
            auto d1 = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
            
            for (int i=0; i<10; i++) {
                naive_multiply_recursion(matrix_list, 1, sizes);
            }   
            
            auto t3 = chrono::high_resolution_clock::now();
            auto d2 = chrono::duration_cast<chrono::milliseconds>(t3 - t2);

            cout << "case n = " << n <<": strassen took " << d1.count() << " ms | naive took " << d2.count() << " ms | " << endl;
        }

        return 0;
    } 
    else {
        const int count = ceil(log2(dimension))+1;

        // create array of matrices
        vector<int> sizes(count, 1);
        sizes[count-1] = dimension + dimension%2;
        for (int i=count-2; i> 0; i--) {
            sizes[i] = sizes[i+1]/2;
            sizes[i] += sizes[i]%2;
        }

        // for (int s:sizes) cout << s << " ";
        // cout << endl;

        vector<int> matrix_list[3*count];
        for (int d=0; d<count; d++) {
            for (int j=0; j<3; j++) {
                matrix_list[3*d+j] = vector<int>(sizes[d]*sizes[d], 0);
            }
        }


        if (flag >= 0) {
            ifstream INFILE;
            INFILE.open(argv[3]);
            if(!INFILE) {
                cout << "File cannot be opened" << endl;
                return -1;
            }
            for (int num=0; num<2; num++) {
                for (int i=0; i<dimension; i++) {
                    for (int j=0; j<dimension; j++) {
                        INFILE >> matrix_list[3*count - 3 + num][j+sizes[count-1]*i];
                    }
                }
            }
            INFILE.close();
            
            if (flag > 0) {
                for (int i: matrix_list[3*count-3]) cout << i << " ";
                cout << endl;
                for (int i: matrix_list[3*count-2]) cout << i << " ";
                cout << endl;
            }
        }
        if (flag == -1 || flag == -2) { // stress test
            srand((unsigned)time(NULL));
            
            for (int num=0; num<2; num++) {
                for (int i=0; i<dimension; i++) {
                    for (int j=0; j<dimension; j++) {
                        matrix_list[3*count - 3 + num][j+sizes[count-1]*i] = rand()%2;
                    }
                }
            }
        }
        

        auto t1 = chrono::high_resolution_clock::now();
        strassen(matrix_list, count-1, sizes);

        if (flag == 1 || flag == -1) {
            for (int i: matrix_list[3*count-1]) cout << i << " ";
            cout << endl;
        }
        

        auto t2 = chrono::high_resolution_clock::now();
        auto d1 = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
        

        if (flag) {
            naive_multiply_recursion(matrix_list, count-1, sizes);
            if (flag == 1 || flag == -1) {    
                for (int i: matrix_list[3*count-1]) cout << i << " ";
                cout << endl;
            }
        }
        
        

        auto t3 = chrono::high_resolution_clock::now();
        auto d2 = chrono::duration_cast<chrono::milliseconds>(t3 - t2);
        if (flag) {
            cout << "strassen took " << d1.count() << " ms | " << endl;
            cout << "naive took " << d2.count() << " ms | " << endl;
        }
        
        // final answer
        for (int i=0; i<dimension; i++) {
            cout << matrix_list[3*count-1][i+i*sizes[count-1]] << endl;
        }
    }
}