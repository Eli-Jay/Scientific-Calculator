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
			// Combine '-' with the next number or variable
			if (i + 1 < size_exp && (std::isdigit(exp[i + 1]) || std::isalpha(exp[i + 1]))) {
				OUT.push_back("-" + std::string(1, exp[i + 1]));
				i++;  // Skip the next character since it's already added
			}
			else {
				OUT.push_back("-");
			}
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
	//std::cout << '|' << OUT[0] << '|';
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
				else {
					result = (std::pow(f_num, s_num));
					//printf("%lf %lf", f_num, s_num);
				}

			}
			stack_2.push_back(result);

		}
		else { 
			for (string c : stack_1) std::cout << c << ", \n";
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
double f_pf(vector<string> equ, const string var, double x) {
	string point_str = std::to_string(x);
	string neg_point_str = std::to_string(x);

	for (size_t i = 0; i < equ.size(); i++) {
		if (equ[i] == var) {
			equ[i] = point_str;
		}
		else if (equ[i] == "-" + var) {
			equ[i] = neg_point_str;
			//std::cout << neg_point_str << '\n';
		}
	}
	return post_fix_exp_solver(equ); // Solve the postfix expression
}


// Integration given limits (NO CAS)

double integration(double upper_bound, double lower_bound, string equ, string var) {
	vector<string> equ_pf = post_fix_conv(equ);
	string temp_exp;
	double x;
	string r_val;
	double res = 0;
	size_t search_i = equ.find(var);

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
		res += weights[i - 1] * f_pf(equ_pf, var, x);
		//std::cout << res;
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




double newton_method(string equ, string var, double guess) {
	vector<string> equ_postfix = post_fix_conv(equ);

	double current = guess;
	double f_x;
	double f_x_prime;
	double new_val;
	double f_new_val;
	double precision = 1e-6;
	size_t MAX_ITER = 20;

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

// Bisection method. Garuantees convergence.
double bisection_method(string equ, string var, double a, double b) {
	vector<string> equ_pf = post_fix_conv(equ);
	size_t small_iter = 0;
	size_t max_iter = 500;
	double precision = 0.1;
	double c = 0;

	double f_a = f_pf(equ_pf, var, a);
	double f_b = f_pf(equ_pf, var, b);
	double f_c;

	// Ensure initial interval has a root TODO: and possible tangents
	while (f_a * f_b > 0) {
		a += 0.1;
		b -= 0.1;
		f_a = f_pf(equ_pf, var, a);
		f_b = f_pf(equ_pf, var, b);
		small_iter++;
		if (small_iter >= 200) return 404.404; // No root found in a reasonable interval
	}

	for (size_t i = 0; i <= max_iter; i++) {
		c = (a + b) / 2;
		f_c = f_pf(equ_pf, var, c);



		// If f(c) is close enough to zero, return c
		if (std::abs(f_c) < precision) {
			return newton_method(equ, var, c);
			//return c;
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
		//if (std::abs(b - a) < precision) {
		//	return (a + b) / 2;
		//	//return newton_method(equ, var, a);
			
		//}

		//std::cout << c << "\n";
	}

	return 34434.404;  // Return if max iterations reached without convergence
}

vector<double> anomalies(string equ, string var, double lim_a, double lim_b) {
	vector<string> equ_postfix = post_fix_conv(equ);
	vector<double> x_anom = {lim_a, lim_b};

	double step = 0.01; // Start with a reasonable step size
	double tolerance = 0.01; // Minimum step size for refinement
	double change_threshold = 0.1; // Normalized change threshold
	double slope_threshold = 10; // Normalized slope threshold

	double prev_val = f_pf(equ_postfix, var, lim_a);
	double curr_val = prev_val;

	for (double t = lim_a + step; t <= lim_b; t += step) {
		curr_val = f_pf(equ_postfix, var, t);

		// Compute relative changes
		double delta_val = std::abs(curr_val - prev_val);
		double norm_change = delta_val / std::max(1.0, std::abs(prev_val));
		double slope = norm_change / step;

		// Detect anomalies based on normalized values
		if (norm_change > change_threshold && slope > slope_threshold) {
			x_anom.push_back(t);
			//printf("Anomaly near: %lf (norm_change: %lf, slope: %lf)\n", t, norm_change, slope);

			// Refine detection by reducing step size temporarily
			double local_step = step / 10;
			while (local_step > tolerance) {
				double refined_t = t - step;
				for (double local_t = refined_t; local_t < t; local_t += local_step) {
					double local_val = f_pf(equ_postfix, var, local_t);
					double local_delta = std::abs(local_val - prev_val);
					if (local_delta / std::max(1.0, std::abs(prev_val)) > change_threshold) {
						x_anom.push_back(local_t);
						break;
					}
					prev_val = local_val;
				}
				local_step /= 2; // Reduce step size for finer granularity
			}
			t += step; // Skip ahead to avoid redundant detections
		}

		prev_val = curr_val;
	}
	return x_anom;
}


// Brent method. 
vector<double> brent_method(string equ, string var, double lim_a, double lim_b) {
	vector<string> equ_postfix = post_fix_conv(equ);
	double a = lim_a;
	double b = lim_b;
	double f_a = f_pf(equ_postfix, var, a);
	double f_b = f_pf(equ_postfix, var, b);
	double f_s;
	double precision = 1e-3;
	size_t i = 0;
	size_t MAX_ITER = 20;

	// Checks for nan values to prevent errors near values that go to infinity.
	while (std::isnan(f_a)) {
		a += 0.05;
		f_a = f_pf(equ_postfix, var, a);
	}
	while (std::isnan(f_b)) {
		b -= 0.05;
		f_b = f_pf(equ_postfix, var, b);
	}


	if (std::abs(f_a) < std::abs(f_b)) {
		std::swap(a, b);
		std::swap(f_a, f_b);
	}

	double c = a;
	double f_c = f_a;
	double d = b - a;
	double s = (a + b) / 2;
	bool mflag = true;

	f_s = f_pf(equ_postfix, var, s);

	

	for (size_t i = 0; i <= MAX_ITER; i++) {
		//precision = std::min(precision, std::abs(b - a) * 1e-10);

		f_a = f_pf(equ_postfix, var, a);
		f_b = f_pf(equ_postfix, var, b);

		//printf("f_a, a: %lf %lf\nf_b, b: %lf %lf\nf_c, c: %lf %lf\nf_s, s: %lf, %lf\n", f_a, a, f_b, b, f_c, c, f_s, s);

		// Checks for nan values to prevent errors near values that go to infinity.
		while (std::isnan(f_a)) {
			a += precision;
			f_a = f_pf(equ_postfix, var, a);
		}

		while (std::isnan(f_b)) {
			b -= precision;
			f_b = f_pf(equ_postfix, var, b);
		}
		f_c = f_pf(equ_postfix, var, c);

		//printf("AFTER:\nf_a, a: %lf %lf\nf_b, b: %lf %lf\nf_c, c: %lf %lf\nf_s, s: %lf, %lf\n", f_a, a, f_b, b, f_c, c, f_s, s);


		// Secant or inverse quadratic interpolation
		if ((f_a != f_c) && (f_b != f_c)) {
			s = (a * f_b * f_c) / ((f_a - f_b) * (f_a - f_c)) +
				(b * f_a * f_c) / ((f_b - f_a) * (f_b - f_c)) +
				(c * f_a * f_b) / ((f_c - f_a) * (f_c - f_b));
			//std::cout << s << "S-QUAD\n";
		}
		else if (f_b - f_a != 0) {
			s = b - f_b * ((b - a) / (f_b - f_a));
			//std::cout << s << "S-NONQUAD\n";
		}
		//else s = (a + b) / 2;

		// Define conditions
		std::vector<bool> conditions = {
			(s < (3 * a + b) / 4 || s > b),               // C1
			(mflag && std::abs(s - b) >= std::abs(b - c) / 2), // C2
			(!mflag && std::abs(s - b) >= std::abs(c - d) / 2), // C3
			(mflag && std::abs(b - c) < precision),          // C4
			(!mflag && std::abs(c - d) < precision)          // C5
		};

		// Check if any condition is true; if so, use bisection
		if (condition_array(conditions)) {
			s = (a + b) / 2;
			//std::cout << s << "S-CONDITIONS\n";
			mflag = true;
		}
		else {
			mflag = false;
		}
		//std::cout << s << "S\n";
		f_s = f_pf(equ_postfix, var, s);


		d = c;
		c = b;
		f_c = f_b;

		if (f_a * f_s < 0) {
			b = s;
			f_b = f_s;
		}
		else {
			a = s;
			f_a = f_s;
		}

		if (std::abs(f_a) < std::abs(f_b)) {
			std::swap(a, b);
			std::swap(f_a, f_b);
		}

		//printf("%.19f\n", s);
		if (std::abs(f_s) < precision) {
			if (s == lim_a || s == lim_b) return { 404.0, 0 }; // Checks if it tries to converge on a value outside the limits.

			return { newton_method(equ, var, s), 0 };
		}
	}

	return { 404.0, 0 };
}





vector<double> equ_solver(string equ, string var, double a, double b) {

	vector<double> interest_points = anomalies(equ, var, a, b);
	vector<double> current = { a, 0 };
	vector<double> zeros;
	double oppos_zero;
	size_t MAX_ITER = 50;
	size_t i = 0;
	interest_points = bubble_sort(interest_points);

	//return interest_points;
	for (size_t j = 0; j < interest_points.size()-1; j) {
		printf("iter: %zd: Checking between %lf and %lf\n", j, interest_points[j], interest_points[j + 1]);
		current = brent_method(equ, var, interest_points[j], interest_points[j + 1]);
		if (current[0] != 404 && in_arr_flt(interest_points, current[0]) != true) {

			if (current[1] == 0) { // This detects an actual zero (0 for zero 1 for interest point)
				printf("Zero point %lf was found.\n Now checking opposite sign.\n", current[0]);
				zeros.push_back(current[0]);
				interest_points.push_back(current[0]);
				oppos_zero = newton_method(equ, var, current[0] * -1);
				printf("OPPO: %lf\n", oppos_zero);
				if (oppos_zero != 404 && in_arr_flt(interest_points, oppos_zero) != true ) {
					printf("Opposite sign found.\n");
					zeros.push_back(oppos_zero);
					interest_points.push_back(oppos_zero);
				}
			}
			//else {
			//	printf("Interest point %lf was found.\n", current[0]);
			//}
			
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