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
    
    string exp;
    cout << "Type an expression\n";
    while (true) {
        cin >> exp;
        vector<string> p_fix = post_fix_conv(exp);
        double solution = post_fix_exp_solver(p_fix);
        cout << "POST FIX: ";
        for (string i : p_fix) cout << i + " ";
        printf("\n>>\t%.8f\n\n", solution);
    
    }
    


    return 0;
}
