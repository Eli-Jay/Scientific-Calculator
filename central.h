#pragma once

#ifndef CENTRAL_H
#define CENTRAL_H

#include <vector>
#include <string>

std::vector<std::string> post_fix_conv(std::string exp);
int prec(char oper);
bool is_num(std::string num);
std::vector<std::string> post_fix_conv(std::string exp);
bool in_arr_str(std::string lst, char c);
std::vector<std::string> reverse(std::vector<std::string> list);
double post_fix_exp_solver(std::vector<std::string> stack_1);
double integration(double upper_bound, double lower_bound, std::string exp, char var = 'x');
double derivative_point(std::string exp, double point, char var);
double bisection_method(std::string equ, char var, double a, double b);
#endif // CENTRAL_H
