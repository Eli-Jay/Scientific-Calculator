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
    
    
    //cout << "\n" << brent_method("0.1*(x-1)*(x+2)*(x-3)*(x+4)*(x-5)", 'x', -4, 0);
    
    vector<long double> zeros = equ_solver("0.1*(x-1)*(x+2)*(x-3)*(x+4)*(x-17.1)", 'x', -1000, 1000);
    for (long double i : zeros) printf("%lf ", i);


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