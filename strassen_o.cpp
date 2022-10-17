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

using std::cout;
using std::vector;
using std::setw, std::endl;
using std::string;

const int b = 32;

unsigned int mix(unsigned int a, unsigned int b, unsigned int c) {
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c;
}

void multiply_block(int *X, int* Y, int *res, int n) {
    for (int j = 0; j < b; ++j)
        for (int s = 0; s < b; ++s)
            for (int i = 0; i < b; ++i)
                res[i + j * b] += X[i + s * b] * Y[s + j * n];
}

void copy_a(int *block_a, int* A, int n) {
    for (int j = 0; j < b; ++j)
        for (int i = 0; i < b; ++i)
            block_a[i + j * b] = A[i + j * n];
}

void copy_c(int *block_c, int* C, int n) {
    for (int j = 0; j < b; ++j)
        for (int i = 0; i < b; ++i)
            C[i + j * n] = block_c[i + j * b];
}

vector<int> matmul(vector<int>& A, vector<int>& B) {
    const int n = sqrt(A.size());

    vector<int> C(n * n);

    static int block_c[b * b];
    static int block_a[b * b];

    for (int j = 0; j < n; j += b) {
        for (int i = 0; i < n; i += b) {
            memset(block_c, 0, b * b * sizeof(int));
            for (int k = 0; k < n; k += b) {
                copy_a(block_a, &(A[i + k * n]), n);
                multiply_block(block_a, &(B[k + j * n]), block_c, n);
            }
            copy_c(block_c, &(C[i + j * n]), n);
        }
    }

    return C;
}


vector<int> mtx_sum(vector<int> mtx_A, vector<int> mtx_B, int split_index, int multiplier = 1) {
    for (auto i = 0; i < split_index; ++i)
		for (auto j = 0; j < split_index; ++j)
			mtx_A[i * split_index + j] = mtx_A[i * split_index + j] + (multiplier * mtx_B[i * split_index + j]);

	return mtx_A;
}

vector<int> mtx_mult(vector<int> mtx_A, vector<int> mtx_B) {
	int col_1 = sqrt(mtx_A.size());
	int row_1 = sqrt(mtx_A.size());
	int col_2 = sqrt(mtx_B.size());
	int row_2 = sqrt(mtx_B.size());

    if (col_1 <= 1024) {
        return matmul(mtx_B, mtx_A);  //block multiplication somehow counts B*A, not A*B....(
    }

	vector<int> result_matrix(mtx_A.size(), 0);

    if (col_1 == 1) {
        result_matrix[0] = mtx_A[0] * mtx_B[0];
    } else {
		int split_index = col_1 / 2;

		vector<int> a_11(split_index * split_index);
		vector<int> a_12(split_index * split_index);
		vector<int> a_21(split_index * split_index);
		vector<int> a_22(split_index * split_index);
		vector<int> b_11(split_index * split_index);
		vector<int> b_12(split_index * split_index);
		vector<int> b_21(split_index * split_index);
    	vector<int> b_22(split_index * split_index);


		for (auto i = 0; i < split_index; i++) {
			for (auto j = 0; j < split_index; j++) {
				a_11[i * split_index + j] = mtx_A[i * col_1 + j];
				a_12[i * split_index + j] = mtx_A[i * col_1 + j + split_index];
				a_21[i * split_index + j] = mtx_A[(split_index + i) * col_1 + j];
				a_22[i * split_index + j] = mtx_A[(i + split_index) * col_1 + j + split_index];
				b_11[i * split_index + j] = mtx_B[i * col_1 + j];
				b_12[i * split_index + j] = mtx_B[i * col_1 + j + split_index];
				b_21[i * split_index + j] = mtx_B[(split_index + i)  * col_1 + j];
				b_22[i * split_index + j] = mtx_B[(i + split_index) * col_1 + j + split_index];
			}
        }

        vector<int> D(mtx_mult(mtx_sum(a_11, a_22, split_index), mtx_sum(b_11, b_22, split_index)));
        vector<int> D_1(mtx_mult(mtx_sum(a_12, a_22, split_index, -1), mtx_sum(b_21, b_22, split_index)));
        vector<int> D_2(mtx_mult(mtx_sum(a_21, a_11, split_index, -1), mtx_sum(b_11, b_12, split_index)));
        vector<int> H_2(mtx_mult(mtx_sum(a_21, a_22, split_index), b_11));
        vector<int> H_1(mtx_mult(mtx_sum(a_11, a_12, split_index), b_22));
        vector<int> V_1(mtx_mult(a_22, mtx_sum(b_21, b_11, split_index, -1)));
		vector<int> V_2(mtx_mult(a_11, mtx_sum(b_12, b_22, split_index, -1)));


		vector<int> result_matrix_00(mtx_sum(mtx_sum(mtx_sum(D, D_1, split_index), V_1, split_index), H_1, split_index, -1));
		vector<int> result_matrix_01(mtx_sum(V_2, H_1, split_index));
    	vector<int> result_matrix_10(mtx_sum(V_1, H_2, split_index));
		vector<int> result_matrix_11(mtx_sum(mtx_sum(mtx_sum(D, D_2, split_index), V_2, split_index), H_2, split_index, -1));


		for (auto i = 0; i < split_index; i++) {
			for (auto j = 0; j < split_index; j++) {
				result_matrix[i * col_1 + j] = result_matrix_00[i * split_index + j];
				result_matrix[i * col_1 + j + split_index] = result_matrix_01[i * split_index + j];
				result_matrix[(split_index + i) * col_1 + j] = result_matrix_10[i * split_index + j];
				result_matrix[(i + split_index) * col_1 + j + split_index] = result_matrix_11[i * split_index + j];
			}
        }
	}

	return result_matrix;
}


vector<int> Generator(int size) {
    srand(mix(clock(), time(NULL), getpid()));

    vector<int> res(size * size);

    for (int i = 0; i < size * size; ++i) {
        res[i] = rand();
    }

    return res;
}

int main(int argc, char** argv) {
    if (argc < 2) return 1;

    const int size = std::stoi(argv[1]);

    vector<int> mtx_A = Generator(size);
    vector<int> mtx_B = Generator(size);

    // vector<int> mtx_A = {1, 5, -4, -8, -10, 4, 5, 6, 1, 3, 5, 8, 7, 9, -4, 2};
    // vector<int> mtx_B = {4, 5, -4, -8, -10, 4, 6, 6, 1, 3, 2, 8, 7, 10, -4, 2};


    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    vector<int> result_matrix(mtx_mult(mtx_A, mtx_B));


    // for (auto x : result_matrix) cout << x << " ";

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << endl << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << endl << endl;


    return 0;
}
