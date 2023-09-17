#include "NumCpp.hpp"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <iomanip>
#include <iostream>
#include <cstdio>
#include <list>
#include <vector>
#include <map>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cfloat>
#include <limits>
#include <stdlib.h>

using std::vector;
using std::list;
using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::stringstream;
using std::map;
using std::pair;

#define EPS 8.881784197001252e-16

enum { OriCovSize = 3, PoseCovSize = 6, TwistCovSize = 6 };


map<string, vector<int>> _AXES2TUPLE = {
    {"sxyz", vector<int>{0, 0, 0, 0}}, {"sxyx", vector<int>{0, 0, 1, 0}}, {"sxzy", vector<int>{0, 1, 0, 0}},
    {"sxzx", vector<int>{0, 1, 1, 0}}, {"syzx", vector<int>{1, 0, 0, 0}}, {"syzy", vector<int>{1, 0, 1, 0}},
    {"syxz", vector<int>{1, 1, 0, 0}}, {"syxy", vector<int>{1, 1, 1, 0}}, {"szxy", vector<int>{2, 0, 0, 0}},
    {"szxz", vector<int>{2, 0, 1, 0}}, {"szyx", vector<int>{2, 1, 0, 0}}, {"szyz", vector<int>{2, 1, 1, 0}},
    {"rzyx", vector<int>{0, 0, 0, 1}}, {"rxyx", vector<int>{0, 0, 1, 1}}, {"ryzx", vector<int>{0, 1, 0, 1}},
    {"rxzx", vector<int>{0, 1, 1, 1}}, {"rxzy", vector<int>{1, 0, 0, 1}}, {"ryzy", vector<int>{1, 0, 1, 1}},
    {"rzxy", vector<int>{1, 1, 0, 1}}, {"ryxy", vector<int>{1, 1, 1, 1}}, {"ryxz", vector<int>{2, 0, 0, 1}},
    {"rzxz", vector<int>{2, 0, 1, 1}}, {"rxyz", vector<int>{2, 1, 0, 1}}, {"rzyz", vector<int>{2, 1, 1, 1}}
};

int _NEXT_AXIS[4] = {1, 2, 0, 1};


struct Gaussian {
    nc::NdArray<double> mean;
    nc::NdArray<double> covariance;

    Gaussian() {}

    Gaussian(nc::NdArray<double>& mean_, nc::NdArray<double>& covariance_) {
        mean = mean_;
        covariance = covariance_;
    }

    Gaussian(const Gaussian& temp) {
        mean = temp.mean;
        covariance = temp.covariance;
    }

    Gaussian& operator=(const Gaussian& x) {
        if (&x == this) {
            return *this;
        }

        mean = x.mean;
        covariance = x.covariance;

        return *this;
    }

    ~Gaussian() {
        // mean.~NdArray(); covariance.~NdArray();
    }

    nc::NdArray<double> get_mean() const { return mean; }
    nc::NdArray<double> get_covariance() const { return covariance; }
};

nc::NdArray<double> sum_mtx(nc::NdArray<double>& Mat1, nc::NdArray<double>& Mat2) {

    auto m = nc::shape(Mat1).rows;
    auto n = nc::shape(Mat2).cols;

    auto rslt = nc::zeros<double>(m, n);

    #pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            #pragma omp atomic
            rslt(i, j) = Mat1(i, j) + Mat2(i, j);
        }
    }

    return rslt;
}

nc::NdArray<double> sub_mtx(nc::NdArray<double>& Mat1, nc::NdArray<double>& Mat2) {

    auto m = nc::shape(Mat1).rows;
    auto n = nc::shape(Mat2).cols;

    auto rslt = nc::zeros<double>(m, n);

    #pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            #pragma omp atomic
            rslt(i, j) = Mat1(i, j) - Mat2(i, j);
        }
    }

    return rslt;
}

double dot(nc::NdArray<double>& vec_1, nc::NdArray<double>& vec_2) {
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2] + vec_1[3] * vec_2[3]);
}

nc::NdArray<double> quaternion_matrix(nc::NdArray<double>& quaternion) {
    auto q = quaternion;

    auto nq = dot(q, q);

    if (nq < EPS) {
        return nc::eye<double>(4);
    }

    q[0] *= sqrtf(2.0 / nq);
    q[1] *= sqrtf(2.0 / nq);
    q[2] *= sqrtf(2.0 / nq);
    q[3] *= sqrtf(2.0 / nq);

    auto tmp_q = nc::outer(q, q);

    return nc::NdArray<double>{
        {1.0 - tmp_q(1, 1) - tmp_q(2, 2), tmp_q(0, 1) - tmp_q(2, 3), tmp_q(0, 2) + tmp_q(1, 3), 0.0},
        {tmp_q(0, 1) + tmp_q(2, 3), 1.0 - tmp_q(0, 0) - tmp_q(2, 2), tmp_q(1, 2) - tmp_q(0, 3), 0.0},
        {tmp_q(0, 2) - tmp_q(1, 3), tmp_q(1, 2) + tmp_q(0, 3), 1.0 - tmp_q(0, 0) - tmp_q(1, 1), 0.0},
        {0.0, 0.0, 0.0, 1.0}
    };
}

nc::NdArray<double> euler_from_matrix(nc::NdArray<double> matrix, const string axes = "sxyz") {
    auto fprf = _AXES2TUPLE[axes];

    auto firstaxis = fprf[0], parity = fprf[1], repetition = fprf[2], frame = fprf[3];

    auto i = firstaxis;
    auto j = _NEXT_AXIS[i + parity];
    auto k = _NEXT_AXIS[i - parity + 1];

    auto M = matrix(nc::Slice(0, 3), nc::Slice(0, 3));

    double ax, ay, az, cy, sy;

    if (repetition) {
        sy = sqrtf(M(i, j) * M(i, j) + M(i, k) * M(i, k));

        if (sy > EPS) {
            ax = atan2(M(i, j), M(i, k));
            ay = atan2(sy, M(i, i));
            az = atan2(M(j, i), -M(k, i));
        } else {
            ax = atan2(-M(j, k), M(j, j));
            ay = atan2(sy, M(i, i));
            az = 0.0;
        }
    } else {
        cy = sqrtf(M(i, i) * M(i, i) + M(j, i) * M(j, i));
        if (cy > EPS) {
            ax = atan2(M(k, j), M(k, k));
            ay = atan2(-M(k, i), cy);
            az = atan2(M(j, i), M(i, i));
        } else {
            ax = atan2(-M(j, k), M(j, j));
            ay = atan2(-M(k, i), cy);
            az = 0.0;
        }
    }

    if (parity) {
        ax = -ax; ay = -ay; az = -az;
    }
    if (frame) {
        std::swap(ax, az);
    }

    return nc::NdArray<double>{ax, ay, az};
}

nc::NdArray<double> euler_matrix(double ai, double aj, double ak, string axes = "sxyz") {
    auto fprf = _AXES2TUPLE[axes];

    auto firstaxis = fprf[0], parity = fprf[1], repetition = fprf[2], frame = fprf[3];

    auto i = firstaxis;
    auto j = _NEXT_AXIS[i + parity];
    auto k = _NEXT_AXIS[i - parity + 1];

    if (frame) {
        std::swap(ai, ak);
    }
    if (parity) {
        ai = -ai; aj = -aj; ak = -ak;
    }

    auto si = sin(ai);
    auto sj = sin(aj);
    auto sk = sin(ak);
    auto ci = cos(ai);
    auto cj = cos(aj);
    auto ck = cos(ak);

    auto cc = ci * ck;
    auto cs = ci * sk;
    auto sc = si * ck;
    auto ss = si * sk;


    nc::NdArray<double> M = nc::eye<double>(4);

    if (repetition) {
        M(i, i) = cj;
        M(i, j) = sj * si;
        M(i, k) = sj * ci;
        M(j, i) = sj * sk;
        M(j, j) = -cj * ss + cc;
        M(j, k) = -cj * cs - sc;
        M(k, i) = -sj * ck;
        M(k, j) = cj * sc + cs;
        M(k, k) = cj * cc - ss;
    } else {
        M(i, i) = cj * ck;
        M(i, j) = sj * sc - cs;
        M(i, k) = sj * cc + ss;
        M(j, i) = cj * sk;
        M(j, j) = sj * ss + cc;
        M(j, k) = sj * cs - sc;
        M(k, i) = -sj;
        M(k, j) = cj * si;
        M(k, k) = cj * ci;
    }

    return M;
}

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


struct ExtendedKalmanFilter {
// private:
    Gaussian estimated_state;
// public:
    ExtendedKalmanFilter () {}
    ExtendedKalmanFilter(Gaussian& prior_state) {
        estimated_state = Gaussian(prior_state.mean, prior_state.covariance);
    }

    static nc::NdArray<double> tfm_to_state(nc::NdArray<double>& tfm) {
        auto euler = euler_from_matrix(tfm);

        return nc::NdArray<double>{tfm(0, 3), tfm(1, 3), tfm(2, 3), euler[0], euler[1], euler[2]};
    }

    static nc::NdArray<double> state_to_tfm(nc::NdArray<double>& state) {
        auto mat = euler_matrix(state[3], state[4], state[5]);

        mat(0, 3) = state[0];
        mat(1, 3) = state[1];
        mat(2, 3) = state[2];

        return mat;
    }

    static nc::NdArray<double> t_jacobian(nc::NdArray<double>& state, nc::NdArray<double>& control) {
        auto j = nc::eye<double>(6);

        auto v = control[0];
        auto r = state[3];
        auto p = state[4];
        auto y = state[5];

        j(0, 4) = -v * cos(y) * sin(p);
        j(0, 5) = -v * cos(p) * sin(y);
        j(1, 4) = -v * sin(p) * sin(y);
        j(1, 5) = v * cos(p) * cos(y);
        j(2, 4) = -v * cos(p);

        return j;
    }

    static nc::NdArray<double> w_jacobian(nc::NdArray<double>& state, nc::NdArray<double>& control) {
        auto j = nc::zeros<double>(6, 2);

        auto r = state[3];
        auto p = state[4];
        auto y = state[5];

        j(0, 0) = cos(p) * cos(y);
        j(1, 0) = cos(p) * sin(y);
        j(2, 0) = -sin(p);
        j(5, 1) = 1;

        return j;
    }


    nc::NdArray<double> transition(nc::NdArray<double>& state, nc::NdArray<double>& control) {
        auto tf_prev = state_to_tfm(state);
        auto mat = euler_matrix(0, 0, control[1], "rxyz");

        mat(0, 3) = control[0];

        auto tmp = nc::matmul(tf_prev, mat);

        return tfm_to_state(tmp);
    }

    Gaussian predict(const Gaussian& control) {
        auto control_m = control.mean;
        auto control_c = control.covariance;

        auto state_m = estimated_state.mean;
        auto state_c = estimated_state.covariance;

        auto new_state_m = transition(state_m, control_m);

        auto j_t = t_jacobian(state_m, control_m);
        auto j_w = w_jacobian(state_m, control_m);

        auto j_t_state_c = nc::matmul(j_t, state_c);
        auto j_w_control_c = nc::matmul(j_w, control_c);

        auto j_t_state_c_j_t = nc::matmul(j_t_state_c, nc::transpose(j_t));
        auto j_w_control_c_j_w = nc::matmul(j_w_control_c, nc::transpose(j_w));

        auto new_state_c = sum_mtx(j_t_state_c_j_t, j_w_control_c_j_w);

        estimated_state = Gaussian(new_state_m, new_state_c);

        return estimated_state;
    }

    Gaussian update(nc::NdArray<double>& h, nc::NdArray<double>& v, Gaussian& measurement) {
        auto measure_m = nc::transpose(measurement.mean);
        auto measure_c = measurement.covariance;

        auto state_m = nc::transpose(estimated_state.mean);
        auto state_c = estimated_state.covariance;
        auto size = nc::size(state_m);

        auto h_state_m = nc::matmul(h, state_m);

        auto innovation = sub_mtx(measure_m, h_state_m);

        auto h_state_c = nc::matmul(h, state_c);
        auto h_state_c_h = nc::matmul(h_state_c, nc::transpose(h));
        auto v_measure_c = nc::matmul(v, measure_c);
        auto v_measure_c_v = nc::matmul(v_measure_c, nc::transpose(v));

        auto innovation_cov = sum_mtx(h_state_c_h, v_measure_c_v);


        const int N = nc::shape(innovation_cov).rows;
        const int M = nc::shape(innovation_cov).cols;

        gsl_matrix *A = gsl_matrix_alloc(N, M);

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                gsl_matrix_set(A, i, j, innovation_cov(i, j));
            }
        }

    	gsl_matrix* A_pinv = moore_penrose_pinv(A, EPS);


        auto pinv_innovation_cov = nc::zeros<double>(A_pinv->size1, A_pinv->size2);

        for (int i = 0; i < A_pinv->size1; ++i) {
            for (int j = 0; j < A_pinv->size2; ++j) {
                pinv_innovation_cov(i, j) = gsl_matrix_get(A_pinv, i, j);
            }
        }

        gsl_matrix_free(A);
    	gsl_matrix_free(A_pinv);


        auto state_c_h = nc::matmul(state_c, nc::transpose(h));

        auto kalman_gain = nc::matmul(state_c_h, pinv_innovation_cov);

        auto kalman_gain_innovation = nc::matmul(kalman_gain, innovation);

        auto new_state_m = sum_mtx(state_m, kalman_gain_innovation);

        auto eye = nc::eye<double>(size);
        auto kalman_gain_h = nc::matmul(kalman_gain, h);

        auto new_state_c = nc::matmul(sub_mtx(eye, kalman_gain_h), state_c);

        nc::NdArray<double> temp_state = {new_state_m(0, 0), new_state_m(1, 0), new_state_m(2, 0),
                                          new_state_m(3, 0), new_state_m(4, 0), new_state_m(5, 0)};

        // cout << temp_state << endl;

        estimated_state = Gaussian(temp_state, new_state_c);

        return estimated_state;
    }
};

int main(int argc, char const** argv) {

    nc::Timer<std::chrono::milliseconds> TikTac;


    int index = 0;
    uint64_t start_time = 0, time, time_odo, time_imu, dt;
    string line;

    uint64_t current_time = 0;

    double lv_x, w, trash;

    double roll, pitch, yaw;

    nc::NdArray<double> quaternion = nc::NdArray<double>(1, 4);

    nc::NdArray<double> h_mat = nc::zeros<double>(3, 6);
    h_mat(0, 3) = 1;
    h_mat(1, 4) = 1;
    h_mat(2, 5) = 1;

    nc::NdArray<double> v_mat = nc::eye<double>(3);

    nc::NdArray<double> m_noise = nc::zeros<double>(3, 3);

    Gaussian prior, state;
    ExtendedKalmanFilter ekf_filter;


    std::ifstream fin("./telemetry_new.txt");
    std::ofstream pos_out("./pos_out.txt", std::fstream::out | std::fstream::trunc);

    TikTac.tic();

    while (std::getline(fin, line)) {
        stringstream streamObj(line);

        streamObj >> quaternion[3] >> quaternion[0] >> quaternion[1] >> quaternion[2];
        streamObj >> lv_x >> trash >> trash;
        streamObj >> trash >> trash >> w;

        auto euler = euler_from_matrix(quaternion_matrix(quaternion));

        if (index == 0) {
            for (int i = 0; i < PoseCovSize; ++i)
                for (int j = 0; j < PoseCovSize; ++j)
                    streamObj >> trash;

            for (int i = 0; i < TwistCovSize; ++i)
                for (int j = 0; j < TwistCovSize; ++j)
                    streamObj >> trash;

            for (int i = 0; i < OriCovSize; ++i)
                for (int j = 0; j < OriCovSize; ++j)
                    streamObj >> m_noise(i, j);
        } else {
            double trash;
            for (int i = 0; i < PoseCovSize * PoseCovSize + TwistCovSize * TwistCovSize + OriCovSize * OriCovSize; ++i) {
                streamObj >> trash;
            }
        }

        streamObj >> time;

        if (index == 0) start_time = time;

        time_odo = time - start_time;
        time_imu = time_odo;

        roll = euler[0];
        pitch = euler[1];
        yaw = euler[2];

        if (index == 0) {
            auto Gaussian_mean = nc::NdArray<double>{0, 0, 0, roll, pitch, yaw};
            auto Gaussian_cov = nc::eye<double>(6, 6);

            #pragma omp parallel for
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    #pragma omp atomic
                    Gaussian_cov(i, j) *= 0.01;
                }
            }

            prior = Gaussian(Gaussian_mean, Gaussian_cov);
            ekf_filter = ExtendedKalmanFilter(prior);
            state = prior;
        }

        dt = time_odo - current_time;
        current_time = time_odo;

        pos_out << std::setprecision(12) << std::fixed << state.mean[0] << " " << state.mean[1] << endl;

        auto control_mean = nc::NdArray<double>{lv_x * dt, w * dt};

        auto control_cov = nc::NdArray<double>{{0.01, 0}, {0, 0.01}};
        Gaussian control(control_mean, control_cov);

        ekf_filter.predict(control);

        auto measurment_mean = nc::NdArray<double>{roll, pitch, yaw};
        Gaussian measurement(measurment_mean, m_noise);

        state = ekf_filter.update(h_mat, v_mat, measurement);

        ++index;
    }

    TikTac.toc(true);

    pos_out.close();
    fin.close();


    return 0;
}
