#include <omp.h>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <string>
#include <algorithm>
#include <bitset>
#include <vector>

using std::string, std::vector, std::to_string;
using std::cin, std::cout, std::endl, std::stoll;

const int K = 3; //Toom-3

string subtract(string lhs, string rhs);
string add(string lhs, string rhs);

unsigned int digit_to_decimal(char digit) {
	return digit - '0';
}

char decimal_to_digit(unsigned int decimal) {
	return decimal + '0';
}

string trim_leading_zeros(string num) {
	bool wasNeg = false;
	string result = "";

	if (num.find("-") != string::npos) {
		wasNeg = true;
		num.erase(0,1);
	}

	size_t i = 0;
	while (num.at(i) == '0') {
		i++;
		if (i == num.length()) {
			return "0";
		}
	}

	while (i != num.length()) {
		result += num.at(i);
		i++;
	}

	if (wasNeg) {
		result.insert(0,"-");
	}

    // return result.erase(0, std::min(result.find_first_not_of('0'), result.size()-1));
    return result;
}

string add(string lhs, string rhs) {
	int carryOver = 0;
	int leftSide = 0;
	int rightSide = 0;
	int sumOfLeftAndRight = 0;
	string result = "";

	bool carryOverExists = true;
	bool negFound = false;

	int leftIterate = lhs.length()-1;
	int rightIterate = rhs.length()-1;

	if ((lhs.find('-') != string::npos && rhs.find('-') != string::npos) ||
		(rhs.find('-') != string::npos && lhs.find('0') == 0) ||
		(lhs.find('-') != string::npos && rhs.find('0') == 0)
	  )
	{
		lhs.erase(0,1);
		rhs.erase(0,1);
		negFound = true;
		leftIterate -= 1;
		rightIterate -= 1;
	} else if (lhs.find('-') != string::npos && rhs.find('-') == string::npos) {
		lhs.erase(0,1);
		result = subtract(lhs,rhs);

		if ((lhs.length() == rhs.length() && rhs > lhs) || lhs.length() < rhs.length())
			result.erase(0,1);
		else
			result.insert(0,"-");
		return trim_leading_zeros(result);
	}

	else if (lhs.find('-') == string::npos && rhs.find('-') != string::npos) {
		rhs.erase(0,1);
		result = subtract(lhs,rhs);
		if ((lhs.length() == rhs.length() && rhs < lhs) || lhs.length() > rhs.length()) {
			if (result.find('-') != string::npos)
				result.erase(0,1);
		}
		return result;
	}

	while (carryOverExists) {
		if (leftIterate >= 0) {
			leftSide = digit_to_decimal(lhs.at(leftIterate));
		}
		if (rightIterate >= 0) {
			rightSide = digit_to_decimal(rhs.at(rightIterate));
		}
		if (leftIterate < 0) {
			leftSide = 0;
		}
		if (rightIterate < 0) {
			rightSide = 0;
		}

		sumOfLeftAndRight = leftSide + rightSide + carryOver;
		carryOver = 0;

		if (sumOfLeftAndRight >= 10) {
			sumOfLeftAndRight = sumOfLeftAndRight % 10;
			result += decimal_to_digit(sumOfLeftAndRight);
			carryOver = 1;
		} else if (leftIterate < 0 && rightIterate < 0) {
			carryOverExists = false;
			result += decimal_to_digit(sumOfLeftAndRight);
		} else {
			result += decimal_to_digit(sumOfLeftAndRight);
		}
		leftIterate -= 1;
		rightIterate -= 1;
	}

    string finalResult = "";
	int i = result.length() - 1;

	while (i >= 0) {
		finalResult += result.at(i);
		i -= 1;
	}

	finalResult = trim_leading_zeros(finalResult);

	if (negFound) {
		finalResult.insert(0,"-");
	}

    // return finalResult.erase(0, std::min(finalResult.find_first_not_of('0'), finalResult.size()-1));
	return finalResult;
}

string subtract(string lhs, string rhs) {
	int negativeCarryOver = 0;
	int leftSide = 0;
	int rightSide = 0;
	int leftMinusRight = 0;
	string temp = "";
	string result = "";
	string addResult = "";

	bool carryOverExists = true;
	bool negFound = false;

	int leftIterate = lhs.length()-1;
	int rightIterate = rhs.length()-1;

	if (lhs.find('-') != string::npos && rhs.find('-') != string::npos) {
		lhs.erase(0,1);
		rhs.erase(0,1);
		leftIterate -= 1;
		rightIterate -= 1;

		if ((lhs.length() == rhs.length() && rhs >= lhs) || rhs.length() > lhs.length()) {
			lhs.swap(rhs);
			std::swap(leftIterate,rightIterate);
		}
        else
			negFound = true;
	} else if (lhs.find('-') == string::npos && rhs.find('-') == string::npos) {
		if ((lhs.length() == rhs.length() && rhs > lhs) || rhs.length() > lhs.length()) {
			negFound = true;
			lhs.swap(rhs);
			std::swap(leftIterate,rightIterate);
		}
	}

	else if (lhs.find('-') != string::npos && rhs.find('-') == string::npos) {
		lhs.erase(0,1);
		result = add(lhs,rhs);
		result.insert(0,"-");
		return result;
	} else if (lhs.find('-') == string::npos && rhs.find('-') != string::npos) {
		rhs.erase(0,1);
		result = add(lhs,rhs);
		return result;
	}

	while (carryOverExists) {
		if (leftIterate >= 0) {
			leftSide = digit_to_decimal(lhs.at(leftIterate));
		}
		if (rightIterate >= 0) {
			rightSide = digit_to_decimal(rhs.at(rightIterate));
		}
		if (leftIterate < 0) {
			leftSide = 0;
		}
		if (rightIterate < 0) {
			rightSide = 0;
		}

		leftMinusRight = leftSide - rightSide - negativeCarryOver;

		if (leftMinusRight < 0) {
			leftMinusRight = (leftSide + 10) - rightSide - negativeCarryOver;
			result += decimal_to_digit(leftMinusRight);
			negativeCarryOver = 1;
		}

		else if (leftIterate < 0 && rightIterate < 0) {
			carryOverExists = false;
			result += decimal_to_digit(leftMinusRight);
		}
		else
		{
			negativeCarryOver = 0;
			result += decimal_to_digit(leftMinusRight);
		}
		leftIterate -= 1;
		rightIterate -= 1;
	}

    string finalResult = "";
	int i = result.length() - 1;

	while (i >= 0) {
		finalResult += result.at(i);
		i -= 1;
	}

	finalResult = trim_leading_zeros(finalResult);

	if (negFound) {
		finalResult.insert(0,"-");
	}
	return finalResult;
}


vector<long long int> matmul(vector<double> matrix, vector<string> vector_) {
    const int matrix_size = std::sqrt(matrix.size());
    vector<long long int> results(matrix_size, 0);

    #pragma omp parallel num_threads(omp_get_num_threads())
    {
        int y, i;
        vector<double> results_private(matrix_size);

        for (y = 0; y < matrix_size ; y++) {
            #pragma omp for
            for (i = 0; i < matrix_size; i++) {
                results_private[y] += std::stoll(vector_[i]) * matrix[i + matrix_size * y];
            }
        }
        #pragma omp critical
        {
            for (y = 0; y < matrix_size; y++) results[y] += round(results_private[y]);
        }
    }

    return results;
}


int eval_i(string m, string n) {
    return std::max(int(int(m.length() / 4) / 3), int(n.length() / 4) / 3) + 1;
}

string toom_cook(string m, string n){
    if ((m.size() <= 9 && n.size() <= 9) || (m[0] == '-' && n[0] == '-' && m.size() <= 10 && n.size() <= 10)){
        return to_string(std::stoll(m) * std::stoll(n));
    }

    bool m_minus = false;
    bool n_minus = false;

    if (m[0] == '-') {
        m.erase(0, 1);
        m_minus = true;
    }
    if (n[0] == '-') {
        n.erase(0, 1);
        n_minus = true;
    }


    const int B = 4 * eval_i(m, n);

    vector<string> p(K, "0");
    vector<string> q(K, "0");

    int k = 0;
    string copy(m);
    std::reverse(copy.begin(), copy.end());

    for (unsigned i = 0; i < copy.size(); i += B) {
        p[k] = copy.substr(i, B);
        std::reverse(p[k].begin(), p[k].end());
        k++;
    }

    copy = n;
    std::reverse(copy.begin(), copy.end());

    k = 0;
    for (unsigned i = 0; i < copy.size(); i += B) {
        q[k] = copy.substr(i, B);
        std::reverse(q[k].begin(), q[k].end());
        k++;
    }

    vector<string> p_p = {p[0], add(add(p[0], p[1]), p[2]), add(subtract(p[0], p[1]), p[2]),  add(subtract(p[0], to_string(2 * stoll(p[1]))), to_string(4 * stoll(p[2]))), p[2]};
    vector<string> q_q = {q[0], add(add(q[0], q[1]), q[2]), add(subtract(q[0], q[1]), q[2]),  add(subtract(q[0], to_string(2 * stoll(q[1]))), to_string(4 * stoll(q[2]))), q[2]};
    //
    // vector<string> q_q = {q[0], to_string(std::stoll(q[0]) + std::stoll(q[1]) + std::stoll(q[2])), to_string(std::stoll(q[0]) - std::stoll(q[1]) + std::stoll(q[2])),
    // to_string(std::stoll(q[0]) - 2 * std::stoll(q[1]) + 4 * std::stoll(q[2])), q[2]};

    // vector<string> p_p = {p[0], to_string(std::stoll(p[0]) + std::stoll(p[1]) + std::stoll(p[2])), to_string(std::stoll(p[0]) - std::stoll(p[1]) + std::stoll(p[2])),           to_string(std::stoll(p[0]) - 2 * std::stoll(p[1]) + 4 * std::stoll(p[2])), p[2]};

    // vector<string> q_q = {q[0], to_string(std::stoll(q[0]) + std::stoll(q[1]) + std::stoll(q[2])), to_string(std::stoll(q[0]) - std::stoll(q[1]) + std::stoll(q[2])),
    // to_string(std::stoll(q[0]) - 2 * std::stoll(q[1]) + 4 * std::stoll(q[2])), q[2]};

    // vector<string> r_in_points = {to_string(std::stoll(p_p[0]) * std::stoll(q_q[0])), to_string(std::stoll(p_p[1]) * std::stoll(q_q[1])), to_string(std::stoll(p_p[2]) * std::stoll(q_q[2])), to_string(std::stoll(p_p[3]) * std::stoll(q_q[3])), to_string(std::stoll(p_p[4]) * std::stoll(q_q[4]))};

    // vector<string> r_in_points = {toom_cook(p_p[0], q_q[0]), toom_cook(p_p[1], q_q[1]), toom_cook(p_p[2], q_q[2]), to_string(stoll(p_p[3]) * stoll(q_q[3])), to_string(stoll(p_p[4]) * stoll(q_q[4]))};

    vector<string> r_in_points = {toom_cook(p_p[0], q_q[0]), toom_cook(p_p[1], q_q[1]), toom_cook(p_p[2], q_q[2]), toom_cook(p_p[3], q_q[3]), toom_cook(p_p[4], q_q[4])};

    vector<double> matrix = {1, 0, 0, 0, 0, 1.0 / 2, 1.0 / 3, -1, 1.0 / 6, -2, -1, 1.0 / 2, 1.0 / 2, 0, -1, -1.0 / 2, 1.0 / 6, 1.0 / 2, -1.0 / 6, 2, 0, 0, 0, 0, 1};

    auto res = matmul(matrix, r_in_points);

    vector<string> recompos = {to_string(res[0]), to_string(res[1]) + string(B, '0'), to_string(res[2]) + string(B * 2, '0'), to_string(res[3]) + string(B * 3, '0'), to_string(res[4]) + string(B * 4, '0')};

    if ((m_minus && n_minus) || (!m_minus && !n_minus)) {
        return add(add(add(add(recompos[4], recompos[3]), recompos[2]), recompos[1]), recompos[0]);
    }

    return "-" + add(add(add(add(recompos[4], recompos[3]), recompos[2]), recompos[1]), recompos[0]);
}


int main(int argc, char** argv) {
    // std::cout<<"Hello from "<< omp_get_thread_num() << " of " << omp_get_num_threads() << std::endl;
    //
    // #pragma omp parallel
    // {
    //     std::ostringstream oss;
    //     oss << "Hello from "<< omp_get_thread_num() << " of " << omp_get_num_threads() << std::endl;
    //     std::cout<<oss.str();
    // }

    cout << toom_cook("1234567890123456789012", "987654321987654321098") << endl;
    cout << toom_cook("1234567890123456789012", "-987654321987654321098") << endl;
    cout << toom_cook("-1234567890123456789012", "-987654321987654321098") << endl;
    cout << toom_cook("-1234567890123456789012", "987654321987654321098") << endl;
    cout << toom_cook("-1234567890123456789012123456789456", "2") << endl;
    cout << toom_cook("2", "-1234567890123456789012123456789456") << endl;
    cout << toom_cook("-118518517008", "2") << endl;
    cout << toom_cook("118518517008", "2") << endl;


    return 0;
}
