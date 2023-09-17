#pragma execution_character_set("utf-8")

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

extern "C"
{
    const int K = 3; // Toom-3
    char *subtract(char *lhs_, char *rhs_);
    char *add(char *lhs_, char *rhs_);
    char *division(char *number_, int divisor);
    char decimal_to_digit(unsigned int decimal);
    unsigned int digit_to_decimal(char digit);

    unsigned int digit_to_decimal(char digit)
    {
        return digit - '0';
    }

    char decimal_to_digit(unsigned int decimal)
    {
        return decimal + '0';
    }

    char *trim_leading_zeros(char *num_)
    {
        string num(num_);

        bool wasNeg = false;

        if (num.find("-") != string::npos)
        {
            wasNeg = true;
            num.erase(0, 1);
        }

        if (wasNeg)
        {
            return string("-" + num.erase(0, std::min(num.find_first_not_of('0'), num.size() - 1))).data();
        }

        return num.erase(0, std::min(num.find_first_not_of('0'), num.size() - 1)).data();
    }

    char *add(char *lhs_, char *rhs_)
    {
        string lhs(lhs_), rhs(rhs_);

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
            result = string(subtract(lhs.data(), rhs.data()));

            if ((lhs.length() == rhs.length() && rhs > lhs) || lhs.length() < rhs.length())
                result.erase(0, 1);
            else
                result.insert(0, "-");
            return trim_leading_zeros(result.data());
        }

        else if (lhs.find('-') == string::npos && rhs.find('-') != string::npos)
        {
            rhs.erase(0, 1);
            result = string(subtract(lhs.data(), rhs.data()));
            if ((lhs.length() == rhs.length() && rhs < lhs) || lhs.length() > rhs.length())
            {
                if (result.find('-') != string::npos)
                    result.erase(0, 1);
            }
            return result.data();
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

        finalResult = string(trim_leading_zeros(finalResult.data()));

        if (negFound)
        {
            finalResult.insert(0, "-");
        }

        // return finalResult.erase(0, std::min(finalResult.find_first_not_of('0'), finalResult.size()-1));

        // cout << finalResult.data() << endl;

        return finalResult.data();
    }

    char *subtract(char *lhs_, char *rhs_)
    {
        string lhs(lhs_), rhs(rhs_);

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
            result = string(add(lhs.data(), rhs.data()));
            result.insert(0, "-");
            return result.data();
        }
        else if (lhs.find('-') == string::npos && rhs.find('-') != string::npos)
        {
            rhs.erase(0, 1);
            result = string(add(lhs.data(), rhs.data()));
            return result.data();
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

        finalResult = string(trim_leading_zeros(finalResult.data()));

        if (negFound)
        {
            finalResult.insert(0, "-");
        }
        return finalResult.data();
    }

    char *division(char *number_, int divisor)
    {
        string number(number_);

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
            return string("0").data();

        if (num_isNeg && div_isNeg || !num_isNeg && !div_isNeg)
            return ans.data();

        return string("-" + ans).data();
    }

    void pyadd(char *lhs_, char *rhs_, char *result)
    {
        auto tmp = string(add(lhs_, rhs_));
        memcpy(result, tmp.data(), tmp.size());
        result[tmp.size()] = '\0';
    }
    void pysubtract(char *lhs_, char *rhs_, char *result)
    {
        auto tmp = string(subtract(lhs_, rhs_));
        memcpy(result, tmp.data(), tmp.size());
        result[tmp.size()] = '\0';
    }
    void pydivision(char *lhs_, char *rhs_, char *result)
    {
        auto tmp = string(division(lhs_, std::stoi(rhs_)));
        memcpy(result, tmp.data(), tmp.size());
        result[tmp.size()] = '\0';
    }
}

int main(int argc, char const *argv[])
{
    string lhs = "132";
    string rhs = "132234";
    char *result = new char[46];

    pyadd(lhs.data(), rhs.data(), result);

    // cout << result << endl;

    return 0;
}
