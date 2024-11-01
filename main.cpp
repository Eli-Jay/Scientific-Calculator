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

    //cout << '\n' << integration(5,0,"2*x*x",'x');
    //cout << '\n' << derivative_point("2*x", 1, 'x');

    cout << '\n' << bisection_method("x^2", 'x', -0.02, 0.01);

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
