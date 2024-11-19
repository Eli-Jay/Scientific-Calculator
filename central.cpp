// Eli Enyart
// Central file used for calculation functions

#include <iostream>
#include <string>
#include <vector>
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
	if (oper == "%") return 2;
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

//Bubble sort (TODO: Use more efficient sorting method later)
vector<double> bubble_sort(vector<double> lst) {
	for (size_t i = 0; i < lst.size() - 1; i++) {
		for (size_t j = 0; j < lst.size() - 1; j++) {
			if (lst[j] >= lst[j + 1]) {
				std::swap(lst[j], lst[j+1]);
			}
		}
	}
	return lst;
}


// Checks if a string is in a vector of strings
bool in_arr_str(const vector<string>& lst, const string& c) {
    for (const auto& i : lst) if (i == c) return true;
    return false;
}

bool in_arr_flt(vector<double> lst, double c) {
	for (long double i : lst) {
		//printf("hi %lf %lf\n", i, c);
		if (std::abs(i-c) < 0.0001) {
			return true;
		}
	}
	return false;
}

bool condition_array(vector<bool> lst) {
	for (bool i : lst) if (i == true) return true;
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
	vector<string> operators = { "^", "*", "/", "+", "-", "%"};
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
	vector<string> operators = { "^", "*", "/", "+", "-", "%"};
	size_t n = stack_1.size();
	double f_num;
	double s_num;
	double result;
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
			if (stack_1[i] == "+") result = (f_num + s_num);
			if (stack_1[i] == "-") result = (f_num - s_num);
			if (stack_1[i] == "*") result = (f_num * s_num);
			//if (stack_1[i] == "%") stack_2.push_back(f_num % s_num); TODO FIX

			// Operators with special cases
			if (stack_1[i] == "/") {
				if (s_num == 0) {
					result = (std::numeric_limits<float>::quiet_NaN());
				}
				else result = (f_num/s_num);	
			}

			if (stack_1[i] == "^") {

				if (f_num == 0 && s_num < 0) {
					result = (std::numeric_limits<float>::quiet_NaN());
				}
				else result = (std::pow(f_num, s_num));

			}
			stack_2.push_back(result);

		}
		else { 
			for (string c : stack_1) std::cout << c;
			std::cout << "\nERROR :" << i; 
		}

	}

	return stack_2[0];
}

// Solves string expressions
double exp_solver(string exp) {
	return post_fix_exp_solver(post_fix_conv(exp));
}

double f_x(string exp, char var, double x) {
	//std::cout << x << "\n";
	//printf("%.15f\n", x);
	string point_str = std::to_string(x);
	string pattern = std::string(1, var);
	string temp_exp = std::regex_replace(exp, std::regex(pattern), point_str);
	return post_fix_exp_solver(post_fix_conv(temp_exp));

}



// Takes in a postfix expression instead. Should be much quicker.
double f_pf(vector<string> equ, string var, double x) {

	string point_str = std::to_string(x);
	for (size_t i = 0; i < equ.size(); i++) {
		if (equ[i] == var) {
			equ[i] = point_str;
		}
	}
	return post_fix_exp_solver(equ);
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
		std::cout << res;
	}
	
	return res * (upper_bound - lower_bound)/2;
}

// Calculates derivative at a point (NO CAS)
double f_prime(string exp, double point, char var) {
	
	double h = 1e-6;
	double f_prime = -1 * f_x(exp, var, point + 2 * h) + 8 * f_x(exp, var, point + h)
		- 8 * f_x(exp, var, point - h) + f_x(exp, var, point - 2 * h);
	f_prime /= (12 * h);

	return f_prime;
}

//Derivative at a point, exactly like previous function but takes in post fix.
double f_prime_pf(vector<string> equ_pf, string var, double point) {
	double h = 1e-6;

	double f_prime = -1 * f_pf(equ_pf, var, point + 2 * h) + 8 * f_pf(equ_pf, var, point + h)
		- 8 * f_pf(equ_pf, var, point - h) + f_pf(equ_pf, var, point - 2 * h);
	f_prime /= (12 * h);

	return f_prime;
}

double f_double_prime_pf(vector<string> equ_pf, string var, double point) {
	double h = 1e-6;

	double f_double_prime = -1 * f_pf(equ_pf, var, point + 2 * h) + 16 * f_pf(equ_pf, var, point + h)
		- 30 * f_pf(equ_pf, var, point) + 16 * f_pf(equ_pf, var, point - h) - f_pf(equ_pf, var, point - 2 * h);
	f_double_prime /= (12 * h * h);

	return f_double_prime;
}




// Expression solving (NO CAS)
double bisection_method(string equ, char var, double a, double b) {
	size_t small_iter = 0;
	size_t max_iter = 500;
	double precision = 1e-8;
	double c = 0;
	
	double f_a = f_x(equ, var, a);
	double f_b = f_x(equ, var, b);
	double f_c;

	// Ensure initial interval has a root TODO: and possible tangents
	while (f_a * f_b > 0) {
		a *= 2;
		b *= 2;
		f_a = f_x(equ, var, a);
		f_b = f_x(equ, var, b);
		small_iter++;
		if (small_iter >= 200) return 404.404; // No root found in a reasonable interval
	}

	for (size_t i = 0; i <= max_iter; i++) {
		c = (a + b) / 2;
		f_c = f_x(equ, var, c);

		

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



double newton_method(string equ, string var, double guess) {
	vector<string> equ_postfix = post_fix_conv(equ);

	double current = guess;
	double f_x;
	double f_x_prime;
	double new_val;
	double f_new_val;
	double precision = 1e-6;
	size_t MAX_ITER = 1000;

	for (size_t i = 0; i < MAX_ITER; i++) {
		f_x = f_pf(equ_postfix, var, current);
		f_x_prime = f_prime_pf(equ_postfix, var, current);

		if (f_x_prime != 0) {
			new_val = current - f_x / f_x_prime;
		}
		else {
			new_val = 0;
		}

		f_new_val = f_pf(equ_postfix, var, current);
		if (std::abs(new_val - current) < precision) {
			return new_val;
		}


		if (std::isnan(new_val)) {
			current += 1e-5;
		}
		else {
			current = new_val;
		}
	}

	return 404;
}
// Brent method. 
vector<double> brent_method(string equ, string var, double lim_a, double lim_b) {
	vector<string> equ_postfix = post_fix_conv(equ);
	double a = lim_a, b = lim_b;
	double f_a = f_pf(equ_postfix, var, a);
	double f_b = f_pf(equ_postfix, var, b);
	double f_s;
	double precision = 1e-3;
	size_t MAX_ITER = 100;

	// Check for initial NaN values
	//while (std::isnan(f_a)) { a += 0.05; f_a = f_pf(equ_postfix, var, a); }
	//while (std::isnan(f_b)) { b -= 0.05; f_b = f_pf(equ_postfix, var, b); }

	// Improved NaN and asymptote handling in Brent's Method
	if (std::isnan(f_a)) {
		a += (b - a) * 0.1; // Move slightly away from a
		f_a = f_pf(equ_postfix, var, a);
	}
	if (std::isnan(f_b)) {
		b -= (b - a) * 0.1; // Move slightly away from b
		f_b = f_pf(equ_postfix, var, b);
	}


	if (std::abs(f_a) < std::abs(f_b)) {
		std::swap(a, b);
		std::swap(f_a, f_b);
	}

	double c = a, f_c = f_a;
	double d = b - a, s = (a + b) / 2;
	bool mflag = true;
	f_s = f_pf(equ_postfix, var, s);
	if (std::isnan(f_s)) return { s, 1 };

	// Main Brent method loop
	for (size_t i = 0; i <= MAX_ITER; i++) {
		f_a = f_pf(equ_postfix, var, a);
		f_b = f_pf(equ_postfix, var, b);

		// **Derivative-based anomaly detection**
		double f_prime_s = f_prime_pf(equ_postfix, var, s);
		double f_prime_a = f_prime_pf(equ_postfix, var, a);
		double f_prime_b = f_prime_pf(equ_postfix, var, b);

		// Relative difference detection
		double delta_ab = std::abs((f_a - f_b) / std::max(std::abs(f_a), std::abs(f_b)));
		double delta_prime = std::max({ std::abs(f_prime_a), std::abs(f_prime_b), std::abs(f_prime_s) });

		// Detect rapid changes using derivative magnitude or a relative difference threshold


		// Detect anomalies at a, b, or s
		if (delta_ab > 0.5 && delta_ab < 1.7 && delta_prime >= 1e6) {
			printf("%lf %lf", delta_ab, delta_prime);
			if (std::abs(f_prime_a) >= 100000) {
				printf("Anomaly detected near a: %lf\n", a);
				return { a, 1 };
			}
			if (std::abs(f_prime_b) >= 100000) {
				printf("Anomaly detected near b: %lf\n", b);
				return { b, 1 };
			}
			if (std::abs(f_prime_s) >= 100000) {
				printf("Anomaly detected near s: %lf\n", s);
				return { s, 1 };
			}
		}

		// Handle NaNs
		// Improved NaN and asymptote handling in Brent's Method
		if (std::isnan(f_a)) {
			a += (b - a) * 0.1; // Move slightly away from a
			f_a = f_pf(equ_postfix, var, a);
		}
		if (std::isnan(f_b)) {
			b -= (b - a) * 0.1; // Move slightly away from b
			f_b = f_pf(equ_postfix, var, b);
		}

		f_c = f_pf(equ_postfix, var, c);

		// Secant or inverse quadratic interpolation
		if ((f_a != f_c) && (f_b != f_c)) {
			s = (a * f_b * f_c) / ((f_a - f_b) * (f_a - f_c)) +
				(b * f_a * f_c) / ((f_b - f_a) * (f_b - f_c)) +
				(c * f_a * f_b) / ((f_c - f_a) * (f_c - f_b));
		}
		else if (f_b - f_a != 0) {
			s = b - f_b * ((b - a) / (f_b - f_a));
		}

		std::vector<bool> conditions = {
			(s < (3 * a + b) / 4 || s > b),
			(mflag && std::abs(s - b) >= std::abs(b - c) / 2),
			(!mflag && std::abs(s - b) >= std::abs(c - d) / 2),
			(mflag && std::abs(b - c) < precision),
			(!mflag && std::abs(c - d) < precision)
		};

		if (condition_array(conditions)) {
			s = (a + b) / 2;
			mflag = true;
		}
		else {
			mflag = false;
		}

		f_s = f_pf(equ_postfix, var, s);
		d = c; c = b; f_c = f_b;

		if (f_a * f_s < 0) {
			b = s; f_b = f_s;
		}
		else {
			a = s; f_a = f_s;
		}

		//printf("current: %lf, %lf\n", f_s, s);
		if (std::abs(f_s) < precision) {
			if (s == lim_a || s == lim_b) return { 404.0, 0 };
			return { newton_method(equ, var, s), 0 };
		}
	}
	return { 404.0, 0 };
}




vector<double> equ_solver(string equ, string var, double a, double b) {
	vector<double> interest_points = { a, b };
	vector<double> current = { a, 0 };
	vector<double> zeros;
	size_t MAX_ITER = 50;
	size_t i = 0;

	for (size_t j = 0; j < interest_points.size()-1; j) {
		printf("iter: %zd: Checking between %lf and %lf\n", j, interest_points[j], interest_points[j + 1]);
		current = brent_method(equ, var, interest_points[j] + 0.1, interest_points[j + 1] - 0.1);
		if (current[0] != 404 && in_arr_flt(interest_points, current[0]) != true) {

			if (current[1] == 0) { // This detects an actual zero (0 for zero 1 for interest point)
				printf("Zero point %lf was found.\n", current[0]);
				zeros.push_back(current[0]);
			}
			else {
				printf("Interest point %lf was found.\n", current[0]);
			}
			interest_points.push_back(current[0]);
			//j = 0;
			interest_points = bubble_sort(interest_points);
			printf("New interest list: ");
			for (double t : interest_points) printf("%lf ", t);
			std::cout << '\n';
		}
		else {
			j++;
		}
	}

	return zeros;

}