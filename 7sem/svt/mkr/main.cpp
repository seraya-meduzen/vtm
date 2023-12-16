#include <cstdio>
#include <vector>
#include <stdexcept>
#include <iostream>

#include "inmost.h"
#include <cmath>

using namespace INMOST;

double f(double x, double y)
{
    return sin(4 * M_PI * x) * sin(M_PI * y);
}

double u(double x, double y, double dx, double dy)
{
    return sin(4 * M_PI * x) * sin(M_PI * y) / ((16 * dx + dy) * M_PI * M_PI);
}

double g_bond(double x, double y)
{
    return 0;
}

void save_solution(size_t n, const std::vector<double> &sol, const std::vector<double> &exact_sol, std::string &vtk_fname)
{
    size_t N = n * n;
    double h = 1. / n;
    FILE *fd = fopen(vtk_fname.c_str(), "w");
    if (!fd)
    {
        throw std::runtime_error("can't open file " + vtk_fname);
    }
    fprintf(fd, "# vtk DataFile Version 2.0\n");
    fprintf(fd, "solution by FDM\n");
    fprintf(fd, "ASCII\nDATASET UNSTRUCTURED_GRID\n");
    fprintf(fd, "POINTS %lu double\n", N);
    for (int i = 0; i * i < N; ++i)
    {
        for (int j = 0; j * j < N; ++j)
        {
            fprintf(fd, "%lf %lf %lf\n", i * h, j * h, 0.0);
        }
    }
    fprintf(fd, "POINT_DATA %lu\n", N);
    fprintf(fd, "SCALARS U double 1\n");
    fprintf(fd, "LOOKUP_TABLE default\n");
    for (int i = 0; i < N; ++i)
    {
        fprintf(fd, "%lf\n", sol[i]);
    }
    fputc('\n', fd);
    if (!exact_sol.empty())
    {
        fprintf(fd, "SCALARS U_true double 1\n");
        fprintf(fd, "LOOKUP_TABLE default\n");
        for (int i = 0; i < N; ++i)
        {
            fprintf(fd, "%lf\n", exact_sol[i]);
        }
        fputc('\n', fd);
    }
}

int main(int argc, char **argv)
{
    Solver::Initialize(&argc, &argv);

    int n = strtol(argv[1], nullptr, 10), N = n * n;
    double h = 1. / n;
    double dx = std::stod(argv[2]), dy = std::stod(argv[3]);

    std::cout << "h: " << h << " n: " << n << std::endl;
    std::cout << "dx: " << dx << " dy: " << dy << std::endl;

    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector sol;

    A.SetInterval(0, N);
    b.SetInterval(0, N);
    sol.SetInterval(0, N);

    for (int k = 0; k < N; ++k)
    {
        int i = k / n, j = k % n;
        double x = h * i, y = h * j;

        if (i == 0 || j == 0 || i == n - 1 || j == n - 1)
        {
            A[k][k] = 1;
            b[k] = g_bond(x, y);
        }
        else
        {
            A[k][k] = 2 * (dx + dy);
            A[k][k + n] = -dx;
            A[k][k - n] = -dx;
            A[k][k - 1] = -dy;
            A[k][k + 1] = -dy;
            b[k] = f(x, y) * h * h;
        }
    }

    Solver S(Solver::INNER_ILU2);

    S.SetParameter("absolute_tolerance", "1e-12");
    S.SetParameter("verbosity", "2");
    S.SetParameter("relative_tolerance", "1e-7");
    S.SetMatrix(A);
    S.Solve(b, sol);

    std::vector<double> full_sol(N);
    for (int i = 0; i < N; ++i)
    {
        full_sol[i] = sol[i];
    }

    std::vector<double> exact_sol(N);
    for (int k = 0; k < N; ++k)
    {
        int i = k / n, j = k % n;
        exact_sol[k] = u(i * h, j * h, dx, dy);
    }

    std::string res_fname = "../res_" + std::to_string(n) + ".vtk";

    double L_2 = 0., C_ = 0;
    for (int i = 0; i < N; ++i)
    {
        L_2 += (exact_sol[i] - sol[i]) * (exact_sol[i] - sol[i]);
        C_ = std::max(C_, abs(exact_sol[i] - sol[i]));
    }
    L_2 = sqrt(L_2) / N;

    std::cout << "L2: " << L_2 << " C-norm: " << C_ << std::endl;

    save_solution(n, full_sol, exact_sol, res_fname);

    Solver::Finalize();
}