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

    int m = std::stoi(argv[1]);
    int n = m;
    int lda = m;

    std::vector<double> A(m * n, 1);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    auto res = SingularValueDecomposition(m, n, lda, &A[0]);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << std::endl << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl << std::endl;
    return 0;
}
