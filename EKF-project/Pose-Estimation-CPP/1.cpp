#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

// #define max(a,b)		((a) > (b) ? (a) : (b))
// #define min(a,b)		((a) < (b) ? (a) : (b))

void print_matrix(const gsl_matrix *m) {
	size_t i, j;

	for (i = 0; i < m->size1; i++) {
		for (j = 0; j < m->size2; j++) {
			printf("%f\t", gsl_matrix_get(m, i, j));
		}
		printf("\n");
	}
}
//
// void print_vector (const gsl_vector * v) {
// 	size_t i;
//
// 	for (i = 0; i < v->size; i++) {
// 		printf("%f\t", gsl_vector_get (v, i));
// 	}
// }

gsl_matrix* moore_penrose_pinv(gsl_matrix *A, const double rcond) {

	gsl_matrix *V, *Sigma_pinv, *U, *A_pinv;
	gsl_matrix *_tmp_mat = NULL;
	gsl_vector *_tmp_vec;
	gsl_vector *u;
	double x, cutoff;
	size_t i, j;
	int n = A->size1;
	int m = A->size2;
	bool was_swapped = false;


	if (m > n) {
		was_swapped = true;
		_tmp_mat = gsl_matrix_alloc(m, n);
		gsl_matrix_transpose_memcpy(_tmp_mat, A);
		A = _tmp_mat;
		i = m;
		m = n;
		n = i;
	}

	V = gsl_matrix_alloc(m, m);
	u = gsl_vector_alloc(m);
	_tmp_vec = gsl_vector_alloc(m);
	gsl_linalg_SV_decomp(A, V, u, _tmp_vec);
	gsl_vector_free(_tmp_vec);

	Sigma_pinv = gsl_matrix_alloc(m, n);
	gsl_matrix_set_zero(Sigma_pinv);
	cutoff = rcond * gsl_vector_max(u);

	for (i = 0; i < m; ++i) {
		if (gsl_vector_get(u, i) > cutoff) {
			x = 1. / gsl_vector_get(u, i);
		}
		else {
			x = 0.;
		}
		gsl_matrix_set(Sigma_pinv, i, i, x);
	}

	U = gsl_matrix_alloc(n, n);
	gsl_matrix_set_zero(U);

	for (i = 0; i < n; ++i) {
		for (j = 0; j < m; ++j) {
			gsl_matrix_set(U, i, j, gsl_matrix_get(A, i, j));
		}
	}

	if (_tmp_mat != NULL) {
		gsl_matrix_free(_tmp_mat);
	}

	_tmp_mat = gsl_matrix_alloc(m, n);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., V, Sigma_pinv, 0., _tmp_mat);

	if (was_swapped) {
		A_pinv = gsl_matrix_alloc(n, m);
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., U, _tmp_mat, 0., A_pinv);
	}
	else {
		A_pinv = gsl_matrix_alloc(m, n);
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., _tmp_mat, U, 0., A_pinv);
	}

	gsl_matrix_free(_tmp_mat);
	gsl_matrix_free(U);
	gsl_matrix_free(Sigma_pinv);
	gsl_vector_free(u);
	gsl_matrix_free(V);

	return A_pinv;
}


int main() {
	const int N = 3;
	const int M = 3;
	const double rcond = 8.881784197001252e-16;


	gsl_matrix *A = gsl_matrix_alloc(N, M);
	gsl_matrix *A_pinv;

	gsl_matrix_set(A, 0, 0, 0.);
	gsl_matrix_set(A, 0, 1, 0.);
	gsl_matrix_set(A, 0, 2, 0.);
	gsl_matrix_set(A, 1, 0, 0.);
	gsl_matrix_set(A, 1, 1, 0.);
	gsl_matrix_set(A, 1, 2, 0.);
	gsl_matrix_set(A, 2, 0, 0.);
	gsl_matrix_set(A, 2, 1, 0.);
	gsl_matrix_set(A, 2, 2, 0.01);

	A_pinv = moore_penrose_pinv(A, rcond);


	print_matrix(A_pinv);

	gsl_matrix_free(A);
	gsl_matrix_free(A_pinv);

	return 0;
}
