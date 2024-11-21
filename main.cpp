// Eli Enyart
// Main file used for user interface

#include <iostream>
#include <string>
#include "central.h"
#include <vector>
#include <iterator>
#include <chrono>
using std::string;
using std::vector;
using std::cout;
using std::cin;


int main()
{   

    //cout << '\n' << integration(0.1,-0.1,"1/(x)","x") << '\n';
    
    
    //cout << "\n" << brent_method("(x^3-6*x^2+11*x-6)/(x^2-1)+2.71^(x/5)-10/(x-0.5)", "x", 1, 10);    
    //return 0;
    auto start = std::chrono::high_resolution_clock::now();
    vector<double> zeros = equ_solver("1/(x-3)-5", "x", -100, 100); 
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    //cout << bisection_method("1/(x)-10", "x", -10.1, 10) << '\n';
    cout << elapsed.count() << '\n';
    
    for (double i : zeros) printf("%lf ", i);

    //cout << exp_solver("-2^2");
    //for (string i : post_fix_conv("(-x)^2")) cout << i;
    return 0;
    string equ;
    double solution;
    //vector<double> zeros;
    //for (char i : exp) cout << i;

    cout << "Type an equation to find the roots of. Use variable 'x'. Type \"exit\" to close.\n";
    while (true) {
        cin >> equ;
        auto start = std::chrono::high_resolution_clock::now();
        if (equ == "exit") break;
        zeros = equ_solver(equ, "x", -5, 5);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        cout << "\n\nSolutions:  >> {";
        for (double i : zeros) printf("%lf ", i);
        cout << "}\n\n";
        cout << elapsed.count() << '\n';
    }
    


    return 0;
}