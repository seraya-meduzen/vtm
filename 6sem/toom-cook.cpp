#include <omp.h>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <string>
#include <algorithm>
#include <bitset>
#include <vector>
#include <utility>

using std::cin, std::cout, std::endl, std::stoll, std::make_pair, std::stoi;
using std::string, std::vector, std::to_string, std::pair;

const int K = 3; // Toom-3

string subtract(string lhs, string rhs);
string add(string lhs, string rhs);
string multiply(string lhs, string rhs);
string division(string number, int divisor);
string toom_cook(string lhs, string rhs);
char decimal_to_digit(unsigned int decimal);
unsigned int digit_to_decimal(char digit);

string multiply(string lhs, string rhs)
{
	return toom_cook(lhs, rhs);
}

unsigned int digit_to_decimal(char digit)
{
	return digit - '0';
}

char decimal_to_digit(unsigned int decimal)
{
	return decimal + '0';
}

string trim_leading_zeros(string num)
{
	bool wasNeg = false;

	if (num.find("-") != string::npos)
	{
		wasNeg = true;
		num.erase(0, 1);
	}

	if (wasNeg)
	{
		return "-" + num.erase(0, std::min(num.find_first_not_of('0'), num.size() - 1));
	}

	return num.erase(0, std::min(num.find_first_not_of('0'), num.size() - 1));
}

string add(string lhs, string rhs)
{
	int carryOver = 0;
	int leftSide = 0;
	int rightSide = 0;
	int sumOfLeftAndRight = 0;
	string result = "";

	bool carryOverExists = true;
	bool negFound = false;

	int leftIterate = lhs.length() - 1;
	int rightIterate = rhs.length() - 1;

	if ((lhs.find('-') != string::npos && rhs.find('-') != string::npos) ||
		(rhs.find('-') != string::npos && lhs.find('0') == 0) ||
		(lhs.find('-') != string::npos && rhs.find('0') == 0))
	{
		lhs.erase(0, 1);
		rhs.erase(0, 1);
		negFound = true;
		leftIterate -= 1;
		rightIterate -= 1;
	}
	else if (lhs.find('-') != string::npos && rhs.find('-') == string::npos)
	{
		lhs.erase(0, 1);
		result = subtract(lhs, rhs);

		if ((lhs.length() == rhs.length() && rhs > lhs) || lhs.length() < rhs.length())
			result.erase(0, 1);
		else
			result.insert(0, "-");
		return trim_leading_zeros(result);
	}

	else if (lhs.find('-') == string::npos && rhs.find('-') != string::npos)
	{
		rhs.erase(0, 1);
		result = subtract(lhs, rhs);
		if ((lhs.length() == rhs.length() && rhs < lhs) || lhs.length() > rhs.length())
		{
			if (result.find('-') != string::npos)
				result.erase(0, 1);
		}
		return result;
	}

	while (carryOverExists)
	{
		if (leftIterate >= 0)
		{
			leftSide = digit_to_decimal(lhs.at(leftIterate));
		}
		if (rightIterate >= 0)
		{
			rightSide = digit_to_decimal(rhs.at(rightIterate));
		}
		if (leftIterate < 0)
		{
			leftSide = 0;
		}
		if (rightIterate < 0)
		{
			rightSide = 0;
		}

		sumOfLeftAndRight = leftSide + rightSide + carryOver;
		carryOver = 0;

		if (sumOfLeftAndRight >= 10)
		{
			sumOfLeftAndRight = sumOfLeftAndRight % 10;
			result += decimal_to_digit(sumOfLeftAndRight);
			carryOver = 1;
		}
		else if (leftIterate < 0 && rightIterate < 0)
		{
			carryOverExists = false;
			result += decimal_to_digit(sumOfLeftAndRight);
		}
		else
		{
			result += decimal_to_digit(sumOfLeftAndRight);
		}
		leftIterate -= 1;
		rightIterate -= 1;
	}

	string finalResult = "";
	int i = result.length() - 1;

	while (i >= 0)
	{
		finalResult += result.at(i);
		i -= 1;
	}

	finalResult = trim_leading_zeros(finalResult);

	if (negFound)
	{
		finalResult.insert(0, "-");
	}

	// return finalResult.erase(0, std::min(finalResult.find_first_not_of('0'), finalResult.size()-1));
	return finalResult;
}

string subtract(string lhs, string rhs)
{
	int negativeCarryOver = 0;
	int leftSide = 0;
	int rightSide = 0;
	int leftMinusRight = 0;
	string temp = "";
	string result = "";
	string addResult = "";

	bool carryOverExists = true;
	bool negFound = false;

	int leftIterate = lhs.length() - 1;
	int rightIterate = rhs.length() - 1;

	if (lhs.find('-') != string::npos && rhs.find('-') != string::npos)
	{
		lhs.erase(0, 1);
		rhs.erase(0, 1);
		leftIterate -= 1;
		rightIterate -= 1;

		if ((lhs.length() == rhs.length() && rhs >= lhs) || rhs.length() > lhs.length())
		{
			lhs.swap(rhs);
			std::swap(leftIterate, rightIterate);
		}
		else
			negFound = true;
	}
	else if (lhs.find('-') == string::npos && rhs.find('-') == string::npos)
	{
		if ((lhs.length() == rhs.length() && rhs > lhs) || rhs.length() > lhs.length())
		{
			negFound = true;
			lhs.swap(rhs);
			std::swap(leftIterate, rightIterate);
		}
	}

	else if (lhs.find('-') != string::npos && rhs.find('-') == string::npos)
	{
		lhs.erase(0, 1);
		result = add(lhs, rhs);
		result.insert(0, "-");
		return result;
	}
	else if (lhs.find('-') == string::npos && rhs.find('-') != string::npos)
	{
		rhs.erase(0, 1);
		result = add(lhs, rhs);
		return result;
	}

	while (carryOverExists)
	{
		if (leftIterate >= 0)
		{
			leftSide = digit_to_decimal(lhs.at(leftIterate));
		}
		if (rightIterate >= 0)
		{
			rightSide = digit_to_decimal(rhs.at(rightIterate));
		}
		if (leftIterate < 0)
		{
			leftSide = 0;
		}
		if (rightIterate < 0)
		{
			rightSide = 0;
		}

		leftMinusRight = leftSide - rightSide - negativeCarryOver;

		if (leftMinusRight < 0)
		{
			leftMinusRight = (leftSide + 10) - rightSide - negativeCarryOver;
			result += decimal_to_digit(leftMinusRight);
			negativeCarryOver = 1;
		}

		else if (leftIterate < 0 && rightIterate < 0)
		{
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

	while (i >= 0)
	{
		finalResult += result.at(i);
		i -= 1;
	}

	finalResult = trim_leading_zeros(finalResult);

	if (negFound)
	{
		finalResult.insert(0, "-");
	}
	return finalResult;
}

string division(string number, int divisor)
{
	bool num_isNeg = false;
	bool div_isNeg = false;

	if (number[0] == '-')
	{
		num_isNeg = true;
		number.erase(0, 1);
	}

	if (divisor < 0)
	{
		div_isNeg = true;
		divisor = -divisor;
	}

	string ans;

	int idx = 0;
	int temp = number[idx] - '0';

	while (temp < divisor)
		temp = temp * 10 + (number[++idx] - '0');

	while (number.size() > idx)
	{
		ans += (temp / divisor) + '0';

		temp = (temp % divisor) * 10 + number[++idx] - '0';
	}

	if (ans.length() == 0)
		return "0";

	if (num_isNeg && div_isNeg || !num_isNeg && !div_isNeg)
		return ans;

	return "-" + ans;
}

vector<string> matmul(vector<pair<string, string>> matrix, vector<string> vector_)
{
	const int matrix_size = std::sqrt(matrix.size());
	vector<string> results(matrix_size, "0");

	int y, i;
	vector<pair<string, string>> results_private(matrix_size, make_pair("0", "1"));

#pragma omp parallel
	{
#pragma omp single
		for (y = 0; y < matrix_size; y++)
		{
#pragma omp task
			{
				for (i = 0; i < matrix_size; i++)
				{
					pair<string, string> tmp = make_pair(toom_cook(vector_[i], matrix[i + matrix_size * y].first), matrix[i + matrix_size * y].second);

					auto num = add(toom_cook(tmp.first, results_private[y].second), toom_cook(tmp.second, results_private[y].first));
					auto denom = toom_cook(tmp.second, results_private[y].second);

					results_private[y].first = num;
					results_private[y].second = denom;
				}
			}
#pragma omp taskwait
			results[y] = division(results_private[y].first, stoi(results_private[y].second));
		}
	}

	return results;
}

int eval_i(string m, string n)
{
	return std::max(int(int(m.length() / 4) / 3), int(int(n.length() / 4) / 3)) + 1;
}

string toom_cook(string m, string n)
{
	if ((m.size() <= 9 && n.size() <= 9) || (m[0] == '-' && n[0] == '-' && m.size() <= 10 && n.size() <= 10))
	{
		return to_string(std::stoll(m) * std::stoll(n));
	}

	bool m_minus = false;
	bool n_minus = false;

	if (m[0] == '-')
	{
		m.erase(0, 1);
		m_minus = true;
	}
	if (n[0] == '-')
	{
		n.erase(0, 1);
		n_minus = true;
	}

	const int B = 4 * eval_i(m, n);

	vector<string> p(K, "0");
	vector<string> q(K, "0");

	int k = 0;
	string copy(m);
	std::reverse(copy.begin(), copy.end());

#pragma omp parallel
	{
#pragma omp single
		{
			for (unsigned i = 0; i < copy.size(); i += B)
			{
#pragma omp task
				{
					p[k] = copy.substr(i, B);
					std::reverse(p[k].begin(), p[k].end());
					k++;
				}
			}
		}
	}

	copy = n;
	std::reverse(copy.begin(), copy.end());

	k = 0;

#pragma omp parallel
	{
#pragma omp single
		{
			for (unsigned i = 0; i < copy.size(); i += B)
			{
#pragma omp task
				{
					q[k] = copy.substr(i, B);
					std::reverse(q[k].begin(), q[k].end());
					k++;
				}
			}
		}
	}

	vector<string> p_p = {p[0], add(add(p[0], p[1]), p[2]), add(subtract(p[0], p[1]), p[2]), add(subtract(p[0], toom_cook("2", p[1])), toom_cook("4", p[2])), p[2]};
	vector<string> q_q = {q[0], add(add(q[0], q[1]), q[2]), add(subtract(q[0], q[1]), q[2]), add(subtract(q[0], toom_cook("2", q[1])), toom_cook("4", q[2])), q[2]};
	vector<string> r_in_points = {toom_cook(p_p[0], q_q[0]), toom_cook(p_p[1], q_q[1]), toom_cook(p_p[2], q_q[2]), toom_cook(p_p[3], q_q[3]), toom_cook(p_p[4], q_q[4])};

	vector<pair<string, string>> matrix = {make_pair("1", "1"), make_pair("0", "1"), make_pair("0", "1"), make_pair("0", "1"), make_pair("0", "1"),
										   make_pair("1", "2"), make_pair("1", "3"), make_pair("-1", "1"), make_pair("1", "6"), make_pair("-2", "1"),
										   make_pair("-1", "1"), make_pair("1", "2"), make_pair("1", "2"), make_pair("0", "1"), make_pair("-1", "1"),
										   make_pair("-1", "2"), make_pair("1", "6"), make_pair("1", "2"), make_pair("-1", "6"), make_pair("2", "1"),
										   make_pair("0", "1"), make_pair("0", "1"), make_pair("0", "1"), make_pair("0", "1"), make_pair("1", "1")};

	auto res = matmul(matrix, r_in_points);

	vector<string> recompos = {res[0], res[1] + string(B, '0'), res[2] + string(B * 2, '0'), res[3] + string(B * 3, '0'), res[4] + string(B * 4, '0')};

	if ((m_minus && n_minus) || (!m_minus && !n_minus))
	{
		return add(add(add(add(recompos[4], recompos[3]), recompos[2]), recompos[1]), recompos[0]);
	}

	return "-" + add(add(add(add(recompos[4], recompos[3]), recompos[2]), recompos[1]), recompos[0]);
}

int main(int argc, char **argv)
{

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
