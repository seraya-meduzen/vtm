#include <omp.h>
#include <iostream>
#include <random>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> floatDist(-1, 0);

using namespace std;

double inner(double *a, double *b, int n)
{
	double sum = 0;

#pragma omp parallel for reduction(+ \
								   : sum)
	for (int i = 0; i < n; i++)
	{
		double tmp = a[i] * b[i];
		sum += tmp;
	}
	return sum;
}

int main(int argc, char *argv[])
{

	int n = 1000000;
	double *a = new double[n];
	double *b = new double[n];
	for (int i = 0; i < n; i++)
	{
		a[i] = (-floatDist(gen));
		b[i] = (-floatDist(gen));
	}

	cout << inner(a, b, n);
}
