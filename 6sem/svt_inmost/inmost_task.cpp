#include <omp.h>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <string>
#include <algorithm>
#include <bitset>
#include <vector>
#include <chrono>
#include "inmost.h"

using namespace INMOST;

using std::cin, std::cout, std::endl, std::stoll;
using std::string, std::vector, std::to_string;

int main(int argc, char **argv)
{
    Solver::Initialize(&argc, &argv);

    Sparse::Matrix A;

    A.Load("/home/meduzen/Documents/vtm codes/cpp/6sem/svt_inmost/0/A_g2_2.mtx");

    Sparse::Vector rhs;

    rhs.Load("/home/meduzen/Documents/vtm codes/cpp/6sem/svt_inmost/0/rhs_g2_2.mtx");

    Solver S("inner_ilu2");
    S.SetParameter("verbosity", "3");
    S.SetParameter("drop_tolerance", argv[1]);

    std::chrono::steady_clock::time_point begin_1 = std::chrono::steady_clock::now(), begin_2, end_1, end_2;

    S.SetMatrix(A);

    end_1 = std::chrono::steady_clock::now();

    Sparse::Vector sol = rhs;

    begin_2 = std::chrono::steady_clock::now();

    S.Solve(sol, rhs);

    end_2 = std::chrono::steady_clock::now();

    std::cout << endl
              << "Preconditioning = " << std::chrono::duration_cast<std::chrono::milliseconds>(end_1 - begin_1).count() << "[ms]" << endl
              << endl;
    std::cout << endl
              << "Solving = " << std::chrono::duration_cast<std::chrono::milliseconds>(end_2 - begin_2).count() << "[ms]" << endl
              << endl;

    std::cout << "iterations: " << S.Iterations() << std::endl;

    Solver::Finalize();

    return 0;
}
