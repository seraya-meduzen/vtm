#include <iomanip>
#include <iostream>
#include <cstdio>
#include <list>
#include <vector>
#include <map>
#include <string>
#include <set>
#include <sys/types.h>
#include <unistd.h>
#include <chrono>
#include <cmath>
#include <cstring>

extern "C" {
    extern int dgesdd_(char*,int*,int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*, int*);
}

struct SVD {
    std::vector<double> u;
    std::vector<double> s;
    std::vector<double> vt;
};

//m, n - shape
//lda - leading dimension
struct SVD SingularValueDecomposition(int m, int n, int lda, double *a) {
    int sigma_size = m < n ? m : n;

    std::vector<double> s(sigma_size);

    std::vector<double> u(m * m);
    std::vector<double> vt(n * n);

    double workSize;
    double *work = &workSize;
    int lwork = -1;
    std::vector<int> iwork(8 * sigma_size);
    int info = 0;

    char q[] = "A";

    dgesdd_(q, &m, &n, a, &lda, &s[0], &u[0], &m, &vt[0], &n, work, &lwork, &iwork[0], &info);

    lwork = workSize;
    work = (double*) malloc(lwork * sizeof(work[0]));

    dgesdd_(q, &m, &n, a, &lda, &s[0], &u[0], &m, &vt[0], &n, work, &lwork, &iwork[0], &info);

    return SVD{.u = u, .s = s, .vt = vt};
}

int main(int argc, char** argv){
    std::vector<double> A = {1, 5, -4, -8, -10, 4, 5, 6, 1, 3, 5, 8, 7, 9, -4, 2};

    int m = 4;
    int n = 4;
    int lda = 4;

    auto res = SingularValueDecomposition(m, n, lda, &A[0]);

    std::cout << "U = ";
    for (auto x : res.u) std::cout << x << " ";
    std::cout << std::endl;
    std::cout << "S = ";
    for (auto x : res.s) std::cout << x << " ";
    std::cout << std::endl;
    std::cout << "V.T = ";
    for (auto x : res.vt) std::cout << x << " ";
    std::cout << std::endl;

    return 0;
}
