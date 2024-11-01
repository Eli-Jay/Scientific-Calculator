// Eli Enyart
// Central file used for calculation functions

#include <iostream>
#include <string>
#include <vector>
#include "central.h"
#include <limits>
#include <cmath>
#include <regex>

using std::string;
using std::vector;

// Returns the precidence of an operator
int prec(string oper) {

	if (oper == "^") return 3;
	if (oper == "*") return 2;
	if (oper == "/") return 2;
	if (oper == "+") return 1;
	if (oper == "-") return 1;
	return 0;

}

// Reverses a list
vector<string> reverse(vector<string> list){
	string temp;
	size_t n = list.size();
	for (int i = 0; i < n / 2; i++) {
		temp = list[list.size() - 1 - i];
		list[list.size() - 1 - i] = list[i];
		list[i] = temp;
	}
	return list;
}

// Checks if num is a valid double/float
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

// Checks if a string is in a vector of strings
bool in_arr_str(vector<string> lst, string c) {
	for (string i : lst) if (i == c) return true;
	return false;
}

bool is_alpha(string str) {
	for (char i : str) if (!std::isalpha(i)) return false;
	return true;
}


// Converts the expression to post fix
vector<string> post_fix_conv(string exp) {
	vector<string> stack_1;
	vector<string> OUT;
	vector<string> operators = { "^", "*", "/", "+", "-" };
	string c;
	size_t size_exp = exp.length();
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

		// Handling and detecting negatives from minus
		else if (c == "-" && (i == 0 || in_arr_str(operators, std::string(1, exp[i - 1])) || exp[i - 1] == '(')) {
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
			
			// Numeric case
			else {
				// Detecting if current number should be added to previous token (is a negative value, decimal, numbers greater than 1 digit)
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

// Solves the post fix expression
double post_fix_exp_solver(vector<string> stack_1) {
	
	vector<double> stack_2;
	vector<string> operators = { "^", "*", "/", "+", "-" };
	size_t n = stack_1.size();
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

// Solves string expressions
double exp_solver(string exp) {

	return post_fix_exp_solver(post_fix_conv(exp));
}

double f_x(string exp, char var, double x) {

	string point_str = std::to_string(x);
	string pattern = std::string(1, var);
	string temp_exp = std::regex_replace(exp, std::regex(pattern), point_str);
	return post_fix_exp_solver(post_fix_conv(temp_exp));

}

// Integration given limits (NO CAS)

double integration(double upper_bound, double lower_bound, string exp, char var) {

	string temp_exp;
	double x;
	string r_val;
	double res = 0;
	size_t search_i = exp.find(var);

	vector<double> nodes = {-0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
		-0.4333953941292472, -0.1488743389816312, 0.1488743389816312,
		0.4333953941292472, 0.6794095682990244, 0.8650633666889845,
		0.9739065285171717};

	vector<double> weights = {0.0666713443086881, 0.1494513491505806, 0.2190863625159820,
		0.2692667193099963, 0.2955242247147529, 0.2955242247147529,
		0.2692667193099963, 0.2190863625159820, 0.1494513491505806,
		0.0666713443086881};

	
	for (size_t i = 1; i <= 10; i++) {

		x = ((upper_bound - lower_bound) / 2) * nodes[i - 1] + (upper_bound + lower_bound) / 2;
		res += weights[i - 1] * f_x(exp, var, x);
	}
	
	return res * (upper_bound - lower_bound)/2;
}

// Calculates derivative at a point (NO CAS)
double derivative_point(string exp, double point, char var) {
	
	double h = 1e-6;
	double f_prime = -1 * f_x(exp, var, point + 2 * h) + 8 * f_x(exp, var, point + h)
		- 8 * f_x(exp, var, point - h) + f_x(exp, var, point - 2 * h);
	f_prime /= (12 * h);

	return f_prime;
}

// Expression solving (NO CAS)
double bisection_method(std::string equ, char var, double a, double b) {
	size_t small_iter = 0;
	size_t max_iter = 500;
	double precision = 1e-8;
	long double c = 0;

	long double f_a = f_x(equ, var, a);
	long double f_b = f_x(equ, var, b);

	// Ensure initial interval has a root
	while (f_a * f_b > 0) {
		a *= 2;
		b *= 2;
		f_a = f_x(equ, var, a);
		f_b = f_x(equ, var, b);
		small_iter++;
		if (small_iter >= 200) return 404.404; // No root found in a reasonable interval
	}

	// Main bisection loop
	for (size_t i = 0; i <= max_iter; i++) {
		c = (a + b) / 2;
		long double f_c = f_x(equ, var, c);

		// If f(c) is close enough to zero, return c
		if (std::abs(f_c) < precision) {
			return c;
		}

		// Choose the side to continue with
		if (f_a * f_c < 0) {
			b = c;
			f_b = f_c;
		}
		else {
			a = c;
			f_a = f_c;
		}

		// Stop if interval is smaller than desired precision
		if (std::abs(b - a) < precision) {
			return (a + b) / 2;
		}

		//std::cout << c << "\n";
	}

	return 34434.404;  // Return if max iterations reached without convergence
}