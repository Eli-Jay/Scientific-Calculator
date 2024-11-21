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
double integration(double upper_bound, double lower_bound, std::string equ, std::string var = "x");
double bisection_method(std::string equ, std::string var, double a, double b);
std::vector<double> brent_method(std::string equ, std::string var, double a, double b);
std::vector<double> equ_solver(std::string equ, std::string var, double a, double b);
double f_pf(std::vector<std::string> equ, std::string var, double x);
double f_double_prime_pf(std::vector<std::string> equ_pf, std::string var, double point);
double exp_solver(std::string exp);

#endif // CENTRAL_H
