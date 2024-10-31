#include <iostream>
#include <string>
#include <vector>
#include "central.h"
#include <limits>
#include <cmath>

using std::string;
using std::vector;

//long exp_solver(void) { return 2; }

int prec(string oper) {

	if (oper == "^") return 3;
	if (oper == "*") return 2;
	if (oper == "/") return 2;
	if (oper == "+") return 1;
	if (oper == "-") return 1;
	return 0;

}


vector<string> reverse(vector<string> list){
	string temp;
	int n = list.size();
	for (int i = 0; i < n / 2; i++) {
		temp = list[list.size() - 1 - i];
		list[list.size() - 1 - i] = list[i];
		list[i] = temp;
	}
	return list;
}

bool is_num(string num) {
	try {
		std::stod(num);
		return true;
	}
	catch (std::invalid_argument&) {
		return false;
	}
	catch (std::out_of_range&) {
		return false;
	}
}

bool in_arr_str(vector<string> lst, string c) {
	for (string i : lst) if (i == c) return true;
	return false;
}

bool is_alpha(string str) {
	for (char i : str) if (!std::isalpha(i)) return false;
	return true;
}



vector<string> post_fix_conv(string exp) {
	vector<string> stack_1;
	vector<string> OUT;
	vector<string> operators = { "^", "*", "/", "+", "-" };
	string c;
	int size_exp = exp.length();
	int i = 0;

	while (i < size_exp) {

		// The character currently being observed, it is converted to a string for compliance
		string c(1, exp[i]);

		// Handling paranthesis
		if (c == "(") stack_1.push_back(c);
		
		else if(c == ")") {

			while (!stack_1.empty() && stack_1.back() != "(") {
				OUT.push_back(stack_1.back());
				stack_1.pop_back();
			}
			stack_1.pop_back();
		}

		// Handling negative numbers
		else if (c == "-" && (i == 0 || in_arr_str(operators, ("" + exp[i - 1])) || exp[i - 1] == '(')) {
			OUT.push_back("-");
		}

		// Numeric values and special cases
		else if (is_num(c) || c == "." || is_alpha(c)) {

			// Handling scientific notation, incorporated from python
			if (c == "e") {
				OUT.back() += c + exp.substr(i + 1, 3);
				i += 3;
			}

			//else if () This will be for handling infinity TODO

			//else if () This will be for handling nan TODO

			else {
				if (OUT.size() &&
					(std::isdigit(exp[i - 1]) || exp[i - 1] == '.' ||
						(exp[i - 1] == '-' && (i == 1 || in_arr_str(operators, std::string(1, exp[i - 2])) || exp[i - 2] == '(')))) {

					OUT.back() += c;  // Append to current token
				}
				else OUT.push_back(c);  // Start a new token
			}
		}
		// Handling operators
		else if (in_arr_str(operators,c)){
			while (stack_1.size() && prec(stack_1.back()) >= prec(c)){

				OUT.push_back(stack_1.back());
				stack_1.pop_back();
			}
			stack_1.push_back(c);
		}

		i++;
	}
	stack_1 = reverse(stack_1);
	OUT.insert(OUT.end(), stack_1.begin(), stack_1.end());
	
	return OUT;

}


double post_fix_exp_solver(vector<string> stack_1) {
	
	vector<double> stack_2;
	vector<string> operators = { "^", "*", "/", "+", "-" };
	int n = stack_1.size();
	double f_num;
	double s_num;
	for (int i = 0; i < n; i++) {

		// Case where stack_1[i] is a number
		if (is_num(stack_1[i])) {
			stack_2.push_back(std::stod(stack_1[i]));
		}

		// Case where stack_1[i] is a operator
		else if (in_arr_str(operators, stack_1[i])) {

			s_num = stack_2.back();
			stack_2.pop_back();

			f_num = stack_2.back();
			stack_2.pop_back();

			// Operator cases
			if (stack_1[i] == "+") stack_2.push_back(f_num + s_num);
			if (stack_1[i] == "-") stack_2.push_back(f_num - s_num);
			if (stack_1[i] == "*") stack_2.push_back(f_num * s_num);

			// Operators with special cases
			if (stack_1[i] == "/") {
				if (s_num == 0) {
					stack_2.push_back(std::numeric_limits<float>::quiet_NaN());
				}
				else stack_2.push_back(f_num/s_num);	
			}

			if (stack_1[i] == "^") {

				if (f_num == 0 && s_num < 0) {
					stack_2.push_back(std::numeric_limits<float>::quiet_NaN());
				}
				else stack_2.push_back(std::pow(f_num, s_num));

			}

		}
		else std::cout << "ERROR :" << i;

	}

	return stack_2[0];
}