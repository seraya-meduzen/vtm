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

#include <fcntl.h>
#include <linux/limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/wait.h>

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

const string ORIENTATION = "orientation";
const string LINEAR_TWIST = "linear_twist";
const string ANGULAR_TWIST = "angular_twist";
const string POSISITON_COV = "position_cov";
const string TWIST_COV = "twist_cov";
const string ORIENTATION_COV = "orientation_cov";
const string TIME = "time";


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
    double* mean;
    double** covariance = nullptr;
    int arr_size;
    int m;
    int n;

    Gaussian() {
        arr_size = 0;
        m = 0;
        n = 0;
    }

    Gaussian(double* mean_, int arr_size_, double** covariance_, int m_, int n_) {
        m = m_; n = n_; arr_size = arr_size_;

        mean = new double[arr_size_];
        covariance = new double*[m];

        for (int i = 0; i < m; ++i) covariance[i] = new double[n];


        std::memcpy(mean, mean_, sizeof(double) * arr_size_);

        for (int i = 0; i < m; ++i) {
            covariance[i] = new double[n];
            std::memcpy(covariance[i], covariance_[i], sizeof(double) * n);
        }

    }

    Gaussian(Gaussian&& temp) {
        arr_size = temp.arr_size;
        m = temp.m;
        n = temp.n;
        temp.mean = nullptr;
        temp.covariance = nullptr;
    }

    Gaussian(const Gaussian& temp) {
        arr_size = temp.arr_size;
        m = temp.m;
        n = temp.n;

        covariance = new double*[m];
        mean = new double[arr_size];

        for (int i = 0; i < arr_size; ++i) mean[i] = temp.mean[i];


        for (int i = 0; i < m; ++i) {
            covariance[i] = new double[n];
            for (int j = 0; j < n; ++j) {
                covariance[i][j] = temp.covariance[i][j];
            }
        }
    }

    Gaussian& operator=(const Gaussian& x) {
        if (&x == this) {
            return *this;
        }

        arr_size = x.arr_size;
        m = x.m;
        n = x.n;

        mean = new double[arr_size];

        for (int i = 0; i < arr_size; ++i) {
            mean[i] = x.mean[i];
        }

        covariance = new double*[m];

        for (int i = 0; i < m; ++i) {
            covariance[i] = new double[n];
            for (int j = 0; j < n; ++j) {
                covariance[i][j] = x.covariance[i][j];
            }
        }

        return *this;
    }

    ~Gaussian() {
        if (arr_size) delete [] mean;

        if (n && m) {
            for (int i = 0; i < m; ++i) delete [] covariance[i];
            delete [] covariance;
        }
    }

    double* get_mean() const { return mean; }
    double** get_covariance() const { return covariance; }
};




struct Odometry {
    double** orientation;
    double** linear_twist;
    double** angular_twist;
    double** position_cov;
    double** twist_cov;
    double** orientation_cov;
    uint64_t* time;

    int message_cnt;

    Odometry(const int msgcnt = 0) {
        message_cnt = msgcnt;

        orientation = new double*[message_cnt];
        for (int i = 0; i < message_cnt; ++i) {
            orientation[i] = new double[3];
        }

        linear_twist = new double*[message_cnt];
        for (int i = 0; i < message_cnt; ++i) {
            linear_twist[i] = new double[3];
        }

        angular_twist = new double*[message_cnt];
        for (int i = 0; i < message_cnt; ++i) {
            angular_twist[i] = new double[3];
        }

        position_cov = new double*[PoseCovSize];
        for (int i = 0; i < PoseCovSize; ++i) {
            position_cov[i] = new double[PoseCovSize];
        }

        twist_cov = new double*[TwistCovSize];
        for (int i = 0; i < TwistCovSize; ++i) {
            twist_cov[i] = new double[TwistCovSize];
        }

        orientation_cov = new double*[OriCovSize];
        for (int i = 0; i < OriCovSize; ++i) {
            orientation_cov[i] = new double[OriCovSize];
        }

        time = new uint64_t[message_cnt];
    }

    ~Odometry() {
        for (int i = 0; i < message_cnt; ++i) delete [] orientation[i];
        for (int i = 0; i < message_cnt; ++i) delete [] linear_twist[i];
        for (int i = 0; i < message_cnt; ++i) delete [] angular_twist[i];
        for (int i = 0; i < PoseCovSize; ++i) delete [] position_cov[i];
        for (int i = 0; i < TwistCovSize; ++i) delete [] twist_cov[i];
        for (int i = 0; i < OriCovSize; ++i) delete [] orientation_cov[i];

        delete [] orientation;
        delete [] linear_twist;
        delete [] angular_twist;
        delete [] position_cov;
        delete [] twist_cov;
        delete [] orientation_cov;
        delete [] time;

    }
};

/*==============================================================================

                                    Prototypes

==============================================================================*/

int GetMsgCnt(const string FileName);
Odometry GetValues(const string FileName);
double dot(double* vec_1, double* vec_2);
double** outer(double* u, double* v);
double** quaternion_matrix(double* quaternion);


/*==============================================================================

                                    Functions

==============================================================================*/

double dot(double* vec_1, double* vec_2) {
    return (vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2] + vec_1[3] * vec_2[3]);
}

double** outer(double* u, double* v) {
    double** res = new double*[4];
    for (int i = 0; i < 4; ++i) {
        res[i] = new double[4];
    }

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            res[i][j] = u[i] * v[j];
        }
    }

    return res;
}


double** quaternion_matrix(double* quaternion) {
    double** res = new double*[4];
    for (int i = 0; i < 4; ++i) {
        res[i] = new double[4];
    }


    double* q = new double[4];
    std::memcpy(q, quaternion, sizeof(double) * 4);


    double nq = dot(q, q);

    if (nq < EPS) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                res[i][j] = 0;
                if (i == j) res[i][i] = 1.0f;
            }
        }
    }

    for (int i = 0; i < 4; ++i) {
        q[i] *= sqrtf(2.0 / nq);
    }


    auto tmp = outer(q, q);

    res[0][0] = 1.0 - tmp[1][1] - tmp[2][2]; res[0][1] = tmp[0][1] - tmp[2][3]; res[0][2] = tmp[0][2] + tmp[1][3]; res[0][3] = 0.0;
    res[1][0] = tmp[0][1] + tmp[2][3]; res[1][1] = 1.0 - tmp[0][0] - tmp[2][2]; res[1][2] = tmp[1][2] - tmp[0][3]; res[1][3] = 0.0;
    res[2][0] = tmp[0][2] - tmp[1][3]; res[2][1] = tmp[1][2] + tmp[0][3]; res[2][2] = 1.0 - tmp[0][0] - tmp[1][1]; res[2][3] = 0.0;
    res[3][0] = 0.0; res[3][1] = 0.0; res[3][2] = 0.0; res[3][3] = 1.0;

    delete [] q;

    return res;
}


double* euler_from_matrix(double** matrix, const string axes = "sxyz") {
    auto fprf = _AXES2TUPLE[axes];

    int firstaxis = fprf[0], parity = fprf[1], repetition = fprf[2], frame = fprf[3];

    int i = firstaxis;
    int j = _NEXT_AXIS[i + parity];
    int k = _NEXT_AXIS[i - parity + 1];


    double** M = new double*[3];
    for (int i = 0; i < 3; ++i) {
        M[i] = new double[3];
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            M[i][j] = matrix[i][j];
        }
    }

    double ax, ay, az, cy, sy;

    if (repetition) {
        sy = sqrtf(M[i][j] * M[i][j] + M[i][k] * M[i][k]);

        if (sy > EPS) {
            ax = atan2(M[i][j], M[i][k]);
            ay = atan2(sy, M[i][i]);
            az = atan2(M[j][i], -M[k][i]);
        } else {
            ax = atan2(-M[j][k], M[j][j]);
            ay = atan2(sy, M[i][i]);
            az = 0.0;
        }
    } else {
        cy = sqrtf(M[i][i] * M[i][i] + M[j][i] * M[j][i]);
        if (cy > EPS) {
            ax = atan2(M[k][j], M[k][k]);
            ay = atan2(-M[k][i], cy);
            az = atan2(M[j][i], M[i][i]);
        } else {
            ax = atan2(-M[j][k], M[j][j]);
            ay = atan2(-M[k][i], cy);
            az = 0.0;
        }
    }

    if (parity) {
        ax = -ax; ay = -ay; az = -az;
    }
    if (frame) {
        std::swap(ax, az);
    }

    for (int i = 0; i < 3; ++i) delete [] M[i];
    delete [] M;


    double* res = new double[3];

    res[0] = ax; res[1] = ay; res[2] = az;

    return res;
}

double** euler_matrix(double ai, double aj, double ak, string axes = "sxyz") {
    auto fprf = _AXES2TUPLE[axes];

    int firstaxis = fprf[0], parity = fprf[1], repetition = fprf[2], frame = fprf[3];

    int i = firstaxis;
    int j = _NEXT_AXIS[i + parity];
    int k = _NEXT_AXIS[i - parity + 1];

    if (frame) {
        std::swap(ai, ak);
    }
    if (parity) {
        ai = -ai; aj = -aj; ak = -ak;
    }

    double si = sin(ai);
    double sj = sin(aj);
    double sk = sin(ak);
    double ci = cos(ai);
    double cj = cos(aj);
    double ck = cos(ak);

    double cc = ci * ck;
    double cs = ci * sk;
    double sc = si * ck;
    double ss = si * sk;

    double** M = new double*[4];
    for (int i = 0; i < 4; ++i) {
        M[i] = new double[4];
        for (int j = 0; j < 4; ++j) {
            M[i][j] = 0;
            if (i == j) M[i][i] = 1;
        }
    }

    if (repetition) {
        M[i][i] = cj;
        M[i][j] = sj * si;
        M[i][k] = sj * ci;
        M[j][i] = sj * sk;
        M[j][j] = -cj * ss + cc;
        M[j][k] = -cj * cs - sc;
        M[k][i] = -sj * ck;
        M[k][j] = cj * sc + cs;
        M[k][k] = cj * cc - ss;
    } else {
        M[i][i] = cj * ck;
        M[i][j] = sj * sc - cs;
        M[i][k] = sj * cc + ss;
        M[j][i] = cj * sk;
        M[j][j] = sj * ss + cc;
        M[j][k] = sj * cs - sc;
        M[k][i] = -sj;
        M[k][j] = cj * si;
        M[k][k] = cj * ci;
    }

    return M;
}



Odometry GetValues(const string FileName) {
    std::ifstream fin(FileName);

    auto size = GetMsgCnt(FileName);
    Odometry odometry(size);


    int index = 0;
    uint64_t start_time = 0, time;
    string line;

    double* quaternion = new double[4];

    while (std::getline(fin, line)) {
        stringstream streamObj(line);

        streamObj >> quaternion[3] >> quaternion[0] >> quaternion[1] >> quaternion[2];
        streamObj >> odometry.linear_twist[index][0] >> odometry.linear_twist[index][1] >> odometry.linear_twist[index][2];
        streamObj >> odometry.angular_twist[index][0] >> odometry.angular_twist[index][1] >> odometry.angular_twist[index][2];


        auto euler = euler_from_matrix(quaternion_matrix(quaternion));

        odometry.orientation[index][0] = euler[0];
        odometry.orientation[index][1] = euler[1];
        odometry.orientation[index][2] = euler[2];

        if (index == 0) {
            for (int i = 0; i < PoseCovSize; ++i)
                for (int j = 0; j < PoseCovSize; ++j)
                    streamObj >> odometry.position_cov[i][j];

            for (int i = 0; i < TwistCovSize; ++i)
                for (int j = 0; j < TwistCovSize; ++j)
                    streamObj >> odometry.twist_cov[i][j];

            for (int i = 0; i < OriCovSize; ++i)
                for (int j = 0; j < OriCovSize; ++j)
                    streamObj >> odometry.orientation_cov[i][j];

        } else {
            double trash;
            for (int i = 0; i < PoseCovSize * PoseCovSize + TwistCovSize * TwistCovSize + OriCovSize * OriCovSize; ++i) {
                streamObj >> trash;
            }
        }

        streamObj >> time;

        if (index == 0) start_time = time;

        odometry.time[index] = time - start_time;

        ++index;
    }



    delete [] quaternion;

    fin.close();
    return odometry;
}



int GetMsgCnt(const string FileName) {
    std::ifstream fin(FileName);

    int number_of_lines = 0;
    string line;

    while (std::getline(fin, line)) ++number_of_lines;


    fin.close();
    return number_of_lines;
}

double** transpose(double** mtx, int n, int m) {
    double** res = new double*[m];
    for (int i = 0; i < m; ++i) {
        res[i] = new double[n];
    }

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            res[i][j] = mtx[j][i];
        }
    }

    return res;
}

double** multiply(double** Mat1, double** Mat2, int l, int m, int n) {
    double** rslt = new double*[l];
    for (int i = 0; i < l; ++i) rslt[i] = new double[n];


    for (int i = 0; i < l; i++) {
        for (int j = 0; j < n; j++) {
            rslt[i][j] = 0;
            for (int k = 0; k < m; k++) {
                rslt[i][j] += Mat1[i][k] * Mat2[k][j];
            }
        }
    }

    return rslt;
}

double** sum_mtx(double** Mat1, double** Mat2, int m, int n) {
    double** rslt = new double*[m];
    for (int i = 0; i < m; ++i) rslt[i] = new double[n];

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            rslt[i][j] = Mat1[i][j] + Mat2[i][j];
        }
    }

    return rslt;
}

double** sub_mtx(double** Mat1, double** Mat2, int m, int n) {
    double** rslt = new double*[m];
    for (int i = 0; i < m; ++i) rslt[i] = new double[n];

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            rslt[i][j] = Mat1[i][j] - Mat2[i][j];
        }
    }

    return rslt;
}

struct ExtendedKalmanFilter {
private:
    Gaussian estimated_state;
public:
    ExtendedKalmanFilter(const Gaussian& prior_state) {
        estimated_state = Gaussian(prior_state.mean, prior_state.arr_size, prior_state.covariance, prior_state.m, prior_state.n);
    }

    static double* tfm_to_state(double** tfm) {
        auto euler = euler_from_matrix(tfm);

        double* res = new double[6];
        res[0] = tfm[0][3];
        res[1] = tfm[1][3];
        res[2] = tfm[2][3];
        res[3] = euler[0];
        res[4] = euler[1];
        res[5] = euler[2];

        return res;
    }

    static double** state_to_tfm(double* state) {
        auto mat = euler_matrix(state[3], state[4], state[5]);

        mat[0][3] = state[0];
        mat[1][3] = state[1];
        mat[2][3] = state[2];

        return mat;
    }

    static double** t_jacobian(double* state, double* control) {
        double** j = new double*[6];
        for (int i = 0; i < 6; ++i) {
            j[i] = new double[6];
            for (int k = 0; k < 6; ++k) {
                j[i][k] = 0;
                if (i == k) j[i][i] = 1;
            }
        }

        double v = control[0];
        double r = state[3];
        double p = state[4];
        double y = state[5];

        j[0][4] = -v * cos(y) * sin(p);
        j[0][5] = -v * cos(p) * sin(y);
        j[1][4] = -v * sin(p) * sin(y);
        j[1][5] = v * cos(p) * cos(y);
        j[2][4] = -v * cos(p);

        return j;
    }

    static double** w_jacobian(double* state, double* control) {
        double** j = new double*[6];
        for (int i = 0; i < 6; ++i) {
            j[i] = new double[2];

            j[i][0] = 0.0; j[i][1] = 0.0;
        }

        double r = state[3];
        double p = state[4];
        double y = state[5];

        j[0][0] = cos(p) * cos(y);
        j[1][0] = cos(p) * sin(y);
        j[2][0] = -sin(p);
        j[5][1] = 1;

        return j;
    }


    double* transition(double* state, double* control) {
        auto tf_prev = state_to_tfm(state);
        auto mat = euler_matrix(0, 0, control[1], "rxyz");

        mat[0][3] = control[0];

        return tfm_to_state(multiply(tf_prev, mat, 4, 4, 4));
    }

    Gaussian predict(const Gaussian& control) {
        auto control_m = control.get_mean();
        auto control_c = control.get_covariance();

        auto state_m = estimated_state.get_mean();
        auto state_c = estimated_state.get_covariance();

        auto new_state_m = transition(state_m, control_m);

        auto j_t = t_jacobian(state_m, control_m);
        auto j_w = w_jacobian(state_m, control_m);

        auto new_state_c = sum_mtx(multiply(multiply(j_t, state_c, 6, estimated_state.m, estimated_state.n), transpose(j_t, 6, 6), 6, estimated_state.n, 6), multiply(multiply(j_w, control_c, 6, control.m, control.n), transpose(j_w, 6, 2), 6, control.m, 6), 6, 6);

        estimated_state = Gaussian(new_state_m, 6, new_state_c, 6, 6);

        return estimated_state;
    }

    Gaussian update(double** h, double** v, Gaussian measurement) {
        double** measure_m = new double*[measurement.arr_size];

        for (int i = 0; i < measurement.arr_size; ++i) {
            measure_m[i] = new double[1];
            measure_m[i][0] = measurement.get_mean()[i];
        }

        double** state_m = new double*[estimated_state.arr_size];

        for (int i = 0; i < estimated_state.arr_size; ++i) {
            state_m[i] = new double[1];
            state_m[i][0] = estimated_state.get_mean()[i];
        }

        auto measure_c = measurement.get_covariance();
        auto state_c = estimated_state.get_covariance();
        auto size = estimated_state.arr_size;

        auto innovation = sub_mtx(measure_m, multiply(h, state_m, 3, 6, 1), 3, 1);

        auto innovation_cov = sum_mtx(multiply(multiply(h, state_c, 3, 6, estimated_state.n), transpose(h, 3, 6), 3, 6, 3), multiply(multiply(v, measure_c, 3, 3, measurement.n), transpose(v, 3, 3), 3,  measurement.n, 3), 3, 3);


        pid_t pinv;

        if ((pinv = fork()) == 0) {
            std::ofstream pinvout("./pinv.txt", std::fstream::out | std::fstream::trunc);

            for (int i = 0 ; i < 3; i++){
                for (int j = 0; j < 3; j++) {
                    pinvout << innovation_cov[i][j] << " ";
                }
            }

            pinvout.close();

            execlp("python3", "python3", "pinv.py", NULL);

            _exit(1);
        }

        wait(NULL);

        double** innovation_cov_pinv = new double*[3];
        for (int i = 0; i < 3; i++) innovation_cov_pinv[i] = new double[3];


        std::ifstream pinvin("./pinvnew.txt");

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                pinvin >> innovation_cov_pinv[i][j];
            }
        }

        pinvin.close();
        remove("./pinvnew.txt");
        remove("./pinv.txt");

        double** eye = new double*[6];
        for (int i = 0; i < 6; ++i) {
            eye[i] = new double[6];
            for (int j = 0; j < 6; ++j) {
                eye[i][j] = 0;
                if (i == j) eye[i][i] = 1;
            }
        }


        auto kalman_gain = multiply(multiply(state_c, transpose(h, 3, 6), 6, 6, 3), innovation_cov_pinv, 6, 3, 3);

        auto new_state_m = sum_mtx(state_m, multiply(kalman_gain, innovation, 6, 3, 1), 6, 1);

        auto new_state_c = multiply(sub_mtx(eye, multiply(kalman_gain, h, 6, 3, 6), 6, 6), state_c, 6, 6, 6);


        double* temp_state = new double[6];

        temp_state[0] = new_state_m[0][0];
        temp_state[1] = new_state_m[1][0];
        temp_state[2] = new_state_m[2][0];
        temp_state[3] = new_state_m[3][0];
        temp_state[4] = new_state_m[4][0];
        temp_state[5] = new_state_m[5][0];


        estimated_state = Gaussian(temp_state, 6, new_state_c, 6, 6);

        // delete [] temp_state;
        //
        // for (int i = 0; i < 6; ++i) delete [] eye[i];
        // delete [] eye;
        //
        // for (int i = 0; i < 3; ++i) delete [] innovation_cov_pinv[i];
        // delete [] innovation_cov_pinv;
        //
        // for (int i = 0; i < estimated_state.arr_size; ++i) delete [] state_m[i];
        // delete [] state_m;
        //
        // for (int i = 0; i < measurement.arr_size; ++i) delete [] measure_m[i];
        // delete [] state_m;

        return estimated_state;
    }
};


int main(int argc, char const** argv) {

    Odometry odometry = GetValues("./telemetry_new.txt");
    auto size = odometry.message_cnt;

    double* lv_x = new double[size];
    double* w = new double[size];

    double* roll_tmp = new double[size];
    double* pitch_tmp = new double[size];
    double* yaw_tmp = new double[size];

    for (int i = 0; i < size; ++i) {
        lv_x[i] = odometry.linear_twist[i][0];
        w[i] = odometry.angular_twist[i][2];

        roll_tmp[i] = odometry.orientation[i][0];
        pitch_tmp[i] = odometry.orientation[i][1];
        yaw_tmp[i] = odometry.orientation[i][2];
    }

    uint64_t* time_odo = odometry.time;
    uint64_t* time_imu = odometry.time;

    /*==========================================================================

                CREATING PYTHON SCRIPT FOR LINEAR INTERPOLATION

    ==========================================================================*/

    setbuf(stdin, NULL);

    char* buf = (char*) calloc (1500, sizeof(char));

    snprintf(buf, 1500, "import numpy as np\n"
                            "import sys\n"
                            "message_count = sum(1 for line in open(sys.argv[1]))\n"
                            "time_odo = np.zeros(message_count)\n"
                            "time_imu = np.zeros(message_count)\n"
                            "third = np.zeros(message_count)\n"
                            "index = 0\n"
                            "with open(sys.argv[1]) as f:\n"
                            "   for line in f:\n"
                            "       time_odo[index], time_imu[index], third[index] = line.split(' ')\n"
                            "       index += 1\n"
                            "towrite = np.interp(time_odo, time_imu, third)\n"
                            "np.savetxt(sys.argv[2], towrite)\n");

    int fd = open("./interp1.py", O_WRONLY | O_TRUNC | O_CREAT, 0660);
    dprintf(fd, "%s", buf);
    free(buf);
    close(fd);


    pid_t roll_pid, pitch_pid, yaw_pid;

    if ((roll_pid = fork()) == 0) {
        std::ofstream roll_out("./roll.txt", std::fstream::out | std::fstream::trunc);

        for (int i = 0; i < size; ++i) {
            roll_out << time_odo[i] << " " << time_imu[i] << " " << roll_tmp[i] << endl;
        }

        roll_out.close();

        execlp("python3", "python3", "interp1.py", "./roll.txt", "./rollinterp.txt", NULL);

        _exit(1);
    } else {
        if ((pitch_pid = fork()) == 0) {
            std::ofstream pitch_out("./pitch.txt", std::fstream::out | std::fstream::trunc);

            for (int i = 0; i < size; ++i) {
                pitch_out << time_odo[i] << " " << time_imu[i] << " " << pitch_tmp[i] << endl;
            }

            pitch_out.close();

            execlp("python3", "python3", "interp1.py", "./pitch.txt", "./pitchinterp.txt", NULL);

            _exit(1);
        } else {
            if ((yaw_pid = fork()) == 0) {
                std::ofstream yaw_out("./yaw.txt", std::fstream::out | std::fstream::trunc);

                for (int i = 0; i < size; ++i) {
                    yaw_out << time_odo[i] << " " << time_imu[i] << " " << yaw_tmp[i] << endl;
                }

                yaw_out.close();

                execlp("python3", "python3", "interp1.py", "./yaw.txt", "./yawinterp.txt", NULL);

                _exit(1);
            }
        }
    }

    while (wait(NULL) != -1) {};

    delete [] roll_tmp;
    delete [] pitch_tmp;
    delete [] yaw_tmp;

    double* roll = new double[size];
    double* pitch = new double[size];
    double* yaw = new double[size];

    std::ifstream roll_in("./rollinterp.txt");
    std::ifstream pitch_in("./pitchinterp.txt");
    std::ifstream yaw_in("./yawinterp.txt");

    for (int i = 0; i < size; ++i) {
        roll_in >> roll[i];
        pitch_in >> pitch[i];
        yaw_in >> yaw[i];
    }

    roll_in.close();
    pitch_in.close();
    yaw_in.close();

    remove("./interp1.py");
    remove("./roll.txt");
    remove("./pitch.txt");
    remove("./yaw.txt");
    remove("./rollinterp.txt");
    remove("./pitchinterp.txt");
    remove("./yawinterp.txt");

    /*==========================================================================

                                END INTERPOLATION

    ==========================================================================*/

    double** h_mat = new double*[3];
    double** v_mat = new double*[3];

    for (int i = 0 ; i < 3; ++i) {
        h_mat[i] = new double[6];
        v_mat[i] = new double[3];
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (j < 3) v_mat[i][j] = 0;
            h_mat[i][j] = 0;
        }
    }

    h_mat[0][3] = 1;
    h_mat[1][4] = 1;
    h_mat[2][5] = 1;
    v_mat[0][0] = 1;
    v_mat[1][1] = 1;
    v_mat[2][2] = 1;


    double** m_noise = odometry.orientation_cov;

    double** tmp_mtx;
    double* tmp_arr;

    tmp_mtx = new double*[6];
    for (int i = 0; i < 6; ++i) {
        tmp_mtx[i] = new double[6];

        for (int j = 0; j < 6; ++j) {
            tmp_mtx[i][j] = 0.0;
            if (i == j) tmp_mtx[i][i] = 0.01;
        }
    }

    tmp_arr = new double[6];
    tmp_arr[0] = 0.0; tmp_arr[1] = 0.0; tmp_arr[2] = 0.0;
    tmp_arr[3] = roll[0]; tmp_arr[4] = pitch[0]; tmp_arr[5] = yaw[0];

    Gaussian prior(tmp_arr, 6, tmp_mtx, 6, 6);

    ExtendedKalmanFilter ekf_filter(prior);


    delete [] tmp_arr;
    for (int i = 0; i < 6; ++i) delete [] tmp_mtx[i];
    delete [] tmp_mtx;


    Gaussian state = prior;


    size = odometry.message_cnt;
    uint64_t current_time = 0;

    double** pos = new double*[size];
    for (int i = 0; i < size; ++i) {
        pos[i] = new double[3];
        for (int j = 0; j < 3; ++j) pos[i][j] = 0;
    }

    /*==========================================================================

                                    MAIN LOOP

    ==========================================================================*/

    for (int i = 0; i < size; ++i) {
        uint64_t dt = time_odo[i] - current_time;
        current_time = time_odo[i];

        pos[i][0] = state.get_mean()[0]; pos[i][1] = state.get_mean()[1]; pos[i][2] = state.get_mean()[2];
        tmp_arr = new double[2];
        tmp_arr[0] = lv_x[i] * dt; tmp_arr[1] = w[i] * dt;



        tmp_mtx = new double*[2];
        for (int k = 0; k < 2; ++k) {
            tmp_mtx[k] = new double[2];
        }

        tmp_mtx[0][0] = 0.01; tmp_mtx[0][1] = 0.0;
        tmp_mtx[1][0] = 0.0; tmp_mtx[1][1] = 0.01;

        Gaussian control(tmp_arr, 2, tmp_mtx, 2, 2);

        ekf_filter.predict(control);

        delete [] tmp_arr;
        for (int k = 0; k < 2; ++k) delete [] tmp_mtx[k];
        delete [] tmp_mtx;


        tmp_arr = new double[3];
        tmp_arr[0] = roll[i]; tmp_arr[1] = pitch[i]; tmp_arr[2] = yaw[i];

        Gaussian measurement(tmp_arr, 3, m_noise, 3, 3);

        state.~Gaussian();

        state = ekf_filter.update(h_mat, v_mat, measurement);
    }

    std::ofstream pos_out("./pos_out.txt", std::fstream::out | std::fstream::trunc);

    for (int i = 0; i < size; ++i) {
        pos_out << std::setprecision(12) << std::fixed << pos[i][0] << " " << pos[i][1] << endl;
    }

    pos_out.close();

    return 0;
}
