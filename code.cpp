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

using namespace std;


vector<vector<int>> naive_multiply (vector<vector<int>>& A, vector<vector<int>>& B, int n) {
    vector<vector<int>> ans(n, vector<int>(0));
    int entry = 0;
    for (int i=0; i<n; i++) {
        for (int k=0; k<n; k++) {
            entry = 0;
            for (int j=0; j<n; j++) {
                entry += (A[i][j]) * (B[j][k]);
            }
            ans[i].emplace_back(entry);
        }
    }
    return ans;
}

vector<vector<int>> strassen_pure (vector<vector<int>>& M1, vector<vector<int>>& M2, int n) {
    vector<vector<int>> ans(n, vector<int>(0));
    // vector<vector<int>> p1(n/2, vector<int>(0)), p2(n/2, vector<int>(0)), p3(n/2, vector<int>(0)), p4(n/2, vector<int>(0)),
    // p5(n/2, vector<int>(0)), p6(n/2, vector<int>(0)), p7(n/2, vector<int>(0));

    if (n==1) {
        ans[0].emplace_back(M1[0][0]*M2[0][0]);
    } else {
        vector<vector<int>> p1,p2,p3,p4,p5,p6,p7;
        vector<vector<int>> m1(n/2, vector<int>(0)), m2(n/2, vector<int>(0));

        // this is for p1
        for (int i=0; i<n/2; i++) {
            for (int j=0; j<n/2; j++) {
                m1[i].emplace_back(M1[i][j]);
            }
        }
        for (int i=0; i<n/2; i++) {
            for (int j=0; j<n/2; j++) {
                m2[i].emplace_back(M2[i][j+n/2]);
                m2[i][j] -= M2[i+n/2][j+n/2];
            }
        }
        p1 = strassen_pure(m1, m2, n/2);

        // p2
        for (int i=0; i<n/2; i++) {
            for (int j=0; j<n/2; j++) {
                m1[i][j] = M1[i][j] + M1[i][j+n/2];
                m2[i][j] = M2[i+n/2][j+n/2];
            }
        }
        p2 = strassen_pure(m1, m2, n/2);

        // p3
        for (int i=0; i<n/2; i++) {
            for (int j=0; j<n/2; j++) {
                m1[i][j] = M1[i+n/2][j] + M1[i+n/2][j+n/2];
                m2[i][j] = M2[i][j];
            }
        }
        p3 = strassen_pure(m1, m2, n/2);

        // p4
        for (int i=0; i<n/2; i++) {
            for (int j=0; j<n/2; j++) {
                m1[i][j] = M1[i+n/2][j+n/2];
                m2[i][j] = M2[i+n/2][j] - M2[i][j];
            }
        }
        p4 = strassen_pure(m1, m2, n/2);

        // p5
        for (int i=0; i<n/2; i++) {
            for (int j=0; j<n/2; j++) {
                m1[i][j] = M1[i][j] + M1[i+n/2][j+n/2];
                m2[i][j] = M2[i][j] + M2[i+n/2][j+n/2];
            }
        }
        p5 = strassen_pure(m1, m2, n/2);

        // p6
        for (int i=0; i<n/2; i++) {
            for (int j=0; j<n/2; j++) {
                m1[i][j] = M1[i][j+n/2] - M1[i+n/2][j+n/2];
                m2[i][j] = M2[i+n/2][j] + M2[i+n/2][j+n/2];
            }
        }
        p6 = strassen_pure(m1, m2, n/2);

        // p7
        for (int i=0; i<n/2; i++) {
            for (int j=0; j<n/2; j++) {
                m1[i][j] = M1[i][j] - M1[i+n/2][j];
                m2[i][j] = M2[i][j] + M2[i][j+n/2];
            }
        }
        p7 = strassen_pure(m1, m2, n/2);

        // recombine
        for (int i=0; i<n/2; i++) {
            for (int j=0; j<n/2; j++) {
                ans[i].emplace_back(p4[i][j] + p5[i][j] + p6[i][j] - p2[i][j]);
            }
            for (int j=0; j<n/2; j++) {
                ans[i].emplace_back(p1[i][j] + p2[i][j]);
            }
            for (int j=0; j<n/2; j++) {
                ans[i+n/2].emplace_back(p3[i][j] + p4[i][j]);
            }
            for (int j=0; j<n/2; j++) {
                ans[i+n/2].emplace_back(p1[i][j] + p5[i][j] - p3[i][j] - p7[i][j]);
            }
        }
    }
    return ans;
}


void strassen_helper (vector<vector<int>>& M1, vector<vector<int>>& M2, int n, vector<vector<int>> &ans, vector<vector<int>> &ops) {
    vector<vector<int>> subans(n, vector<int>(n, 0));

    // calculate subans
    vector<vector<int>> m1(n/2, vector<int>(0)), m2(n/2, vector<int>(0));

    // this is for p1
    // before every step reset subans to 0
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++) {
            subans[i][j]=0;
        }
    }
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            m1[i].emplace_back(M1[i][j]);
        }
    }
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            m2[i].emplace_back(M2[i][j+n/2]);
            m2[i][j] -= M2[i+n/2][j+n/2];
        }
    }
    strassen_helper(m1, m2, n/2, subans, {{2,1},{4,1}});
    for (auto op:ops) {
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (op[0]==1) {
                    ans[i][j] += op[1]*subans[i][j];
                } else if (op[0]==2) {
                    ans[i][j+n] += op[1]*subans[i][j];
                } else if (op[0]==3) {
                    ans[i+n][j] += op[1]*subans[i][j];
                } else {
                    ans[i+n][j+n] += op[1]*subans[i][j];
                }
            }
        }
    }


    // this is for p2
    // before every step reset subans to 0
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++) {
            subans[i][j]=0;
        }
    }
    for (int i=0; i<n/2; i++) {
            for (int j=0; j<n/2; j++) {
                m1[i][j] = M1[i][j] + M1[i][j+n/2];
                m2[i][j] = M2[i+n/2][j+n/2];
            }
        }
    strassen_helper(m1, m2, n/2, subans, {{1,-1},{2,1}});
    for (auto op:ops) {
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (op[0]==1) {
                    ans[i][j] += op[1]*subans[i][j];
                } else if (op[0]==2) {
                    ans[i][j+n] += op[1]*subans[i][j];
                } else if (op[0]==3) {
                    ans[i+n][j] += op[1]*subans[i][j];
                } else {
                    ans[i+n][j+n] += op[1]*subans[i][j];
                }
            }
        }
    }


    // this is for p3
    // before every step reset subans to 0
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++) {
            subans[i][j]=0;
        }
    }
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            m1[i][j] = M1[i+n/2][j] + M1[i+n/2][j+n/2];
            m2[i][j] = M2[i][j];
        }
    }
    strassen_helper(m1, m2, n/2, subans, {{3,1},{4,-1}});
    for (auto op:ops) {
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (op[0]==1) {
                    ans[i][j] += op[1]*subans[i][j];
                } else if (op[0]==2) {
                    ans[i][j+n] += op[1]*subans[i][j];
                } else if (op[0]==3) {
                    ans[i+n][j] += op[1]*subans[i][j];
                } else {
                    ans[i+n][j+n] += op[1]*subans[i][j];
                }
            }
        }
    }


    // this is for p4
    // before every step reset subans to 0
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++) {
            subans[i][j]=0;
        }
    }
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            m1[i][j] = M1[i+n/2][j+n/2];
            m2[i][j] = M2[i+n/2][j] - M2[i][j];
        }
    }
    strassen_helper(m1, m2, n/2, subans, {{1,1},{3,1}});
    for (auto op:ops) {
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (op[0]==1) {
                    ans[i][j] += op[1]*subans[i][j];
                } else if (op[0]==2) {
                    ans[i][j+n] += op[1]*subans[i][j];
                } else if (op[0]==3) {
                    ans[i+n][j] += op[1]*subans[i][j];
                } else {
                    ans[i+n][j+n] += op[1]*subans[i][j];
                }
            }
        }
    }

    // this is for p5
    // before every step reset subans to 0
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++) {
            subans[i][j]=0;
        }
    }
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            m1[i][j] = M1[i][j] + M1[i+n/2][j+n/2];
            m2[i][j] = M2[i][j] + M2[i+n/2][j+n/2];
        }
    }
    strassen_helper(m1, m2, n/2, subans, {{1,1},{4,1}});
    for (auto op:ops) {
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (op[0]==1) {
                    ans[i][j] += op[1]*subans[i][j];
                } else if (op[0]==2) {
                    ans[i][j+n] += op[1]*subans[i][j];
                } else if (op[0]==3) {
                    ans[i+n][j] += op[1]*subans[i][j];
                } else {
                    ans[i+n][j+n] += op[1]*subans[i][j];
                }
            }
        }
    }


    // this is for p6
    // before every step reset subans to 0
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++) {
            subans[i][j]=0;
        }
    }
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            m1[i][j] = M1[i][j+n/2] - M1[i+n/2][j+n/2];
            m2[i][j] = M2[i+n/2][j] + M2[i+n/2][j+n/2];
        }
    }
    strassen_helper(m1, m2, n/2, subans, {{1,1}});
    for (auto op:ops) {
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (op[0]==1) {
                    ans[i][j] += op[1]*subans[i][j];
                } else if (op[0]==2) {
                    ans[i][j+n] += op[1]*subans[i][j];
                } else if (op[0]==3) {
                    ans[i+n][j] += op[1]*subans[i][j];
                } else {
                    ans[i+n][j+n] += op[1]*subans[i][j];
                }
            }
        }
    }


    // this is for p7
    // before every step reset subans to 0
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++) {
            subans[i][j]=0;
        }
    }
    for (int i=0; i<n/2; i++) {
        for (int j=0; j<n/2; j++) {
            m1[i][j] = M1[i][j] - M1[i+n/2][j];
            m2[i][j] = M2[i][j] + M2[i][j+n/2];
        }
    }
    strassen_helper(m1, m2, n/2, subans, {{4,-1}});
    for (auto op:ops) {
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (op[0]==1) {
                    ans[i][j] += op[1]*subans[i][j];
                } else if (op[0]==2) {
                    ans[i][j+n] += op[1]*subans[i][j];
                } else if (op[0]==3) {
                    ans[i+n][j] += op[1]*subans[i][j];
                } else {
                    ans[i+n][j+n] += op[1]*subans[i][j];
                }
            }
        }
    }



    // modify approrpriate values in ans?


}

int main() {
    // naive test
    // vector<vector<int>> A = {{1,2},{3,4}}, B = {{0,1},{1,0}};
    // vector<vector<int>> ans = naive_multiply(A,B, 2);
    // for (int i=0; i<2; i++) {
    //     for (int j=0; j<2; j++) {
    //         cout << ans[i][j] << " ";
    //     } cout << endl;
    // }

    // strassen test
    // int n=4;
    // vector<vector<int>> A = {{1,2,3,4},{5,6,7,8},{9,10,11,12},{13,14,15,16}}, B = {{1000, 0, 0, 1000}, {1,1,1,1}, {-100, 99, 99, 99},{0,1,2,3}};
    // vector<vector<int>> ans1 = naive_multiply(A,B,n);
    // for (int i=0; i<n; i++) {
    //     for (int j=0; j<n; j++) {
    //         cout << ans1[i][j] << " ";
    //     } cout << endl;
    // }
    // vector<vector<int>> ans2 = strassen_pure(A,B,n);
    // for (int i=0; i<n; i++) {
    //     for (int j=0; j<n; j++) {
    //         cout << ans2[i][j] << " ";
    //     } cout << endl;
    // }

    // stress test
    srand((unsigned)time(NULL));
    int n=256;
    vector<vector<int>> bm1(n, vector<int>(0)), bm2(n, vector<int>(0));
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            bm1[i].emplace_back(rand()%2);
            bm2[i].emplace_back(rand()%2);
        }
    }
    // for (int i:bm1[0]) cout << i << " ";

    auto t1 = chrono::high_resolution_clock::now();
    
    vector<vector<int>> ans = strassen_pure(bm1, bm2, n);

    auto t2 = chrono::high_resolution_clock::now();
    auto d1 = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
    cout << "strassen took " << d1.count() << " ms | ";

    ans =  naive_multiply(bm1, bm2, n);

    auto t3 = chrono::high_resolution_clock::now();
    auto d2 = chrono::duration_cast<chrono::milliseconds>(t3 - t2);
    cout << "naive took " << d2.count() << " ms | ";
}