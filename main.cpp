// Eli Enyart
// Main file used for user interface

#include <iostream>
#include <string>
#include "central.h"
#include <vector>
#include <iterator>
using std::string;
using std::vector;
using std::cout;
using std::cin;


int main()
{   

    //cout << '\n' << integration(1.1,-100,"(1/(x-1))",'x');
    //cout << '\n' << derivative_point("(1/(x-1))-2", 100, 'x');
    
    
    //cout << "\n" << brent_method("(x^3-6*x^2+11*x-6)/(x^2-1)+2.71^(x/5)-10/(x-0.5)", "x", 1, 10);    
 
    //return 0;
    vector<double> zeros = equ_solver("(x^3-6*x^2+11*x-6)/(x^2-1)+2.71^(x/5)-10/(x-0.5)", "x", -10, 10);
    for (double i : zeros) printf("%lf ", i);


    return 0;
    string exp;
    double solution;
    //for (char i : exp) cout << i;

    cout << "Type an expression\n";
    while (true) {
        cin >> exp;
        vector<string> p_fix = post_fix_conv(exp);
        solution = post_fix_exp_solver(p_fix);
        cout << "POST FIX: ";
        for (string i : p_fix) cout << i + " ";
        printf("\n>>\t%.8f\n\n", solution);
    
    }
    


    return 0;
}