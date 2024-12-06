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

// converts user input for algorithm to understand ex 3x -> 3*x, pi->3.14.... |NOTE: -x^2 -> -1*(x^2)
string input_conversion(string input) {
    string out;

    //pi -> 3.14

    //multiplication ex: 3(1-3) -> 3*(1-3)

    //sin(x) -> equation est

    //cos(x)

    //tan(x)

    //n!

    return out;
}

// used to read in commands of the user. 
string command_reader(string input) {
    string out;
    
    if (input.substr(0, 3) == "sol") { // Command when user wants to find the zeros of an equations.
        size_t equal_location = 0;
        size_t comma_location = 0;
        string first_equ = "";
        string second_equ = "";
        vector<double> zeros;

        for (equal_location = 4; input[equal_location] != '='; equal_location++); // Finds where the '=' is located
        //for (comma_location = 4; input[comma_location] != ','; comma_location++); // Finds where the ',' is located

        first_equ = input.substr(4, equal_location - 4);
        second_equ = input.substr(equal_location + 1, input.size() - equal_location - 2);
        string comb_equ = '(' + first_equ + ")-(" + second_equ + ")"; // fx = gx -> fx - gx = 0

        zeros = equ_solver(comb_equ, "x", -300, 300);
        cout << "Solutions: {";
        for (double z : zeros) printf("%lf ", z);
        cout << "}\n\n";
     }
    

    //sol(equ_1=eq_2, var) function - equ_solver()

    //int(bottom, top, equation, var) - integrate()

    //der(equation, variable, point) - f_prime()

    //sum(equation, variable, to) - sum() (create)

    //etc...

    else {
        cout << exp_solver(input) << "\n\n"; //default (none) - exp_solver()
    }

    



    return out;
}



int main()
{   

    //cout << '\n' << integration(0.1,-0.1,"1/(x)","x") << '\n';
    
    
    //cout << "\n" << brent_method("(x^3-6*x^2+11*x-6)/(x^2-1)+2.71^(x/5)-10/(x-0.5)", "x", 1, 10); 
    //cout << find_range("x^2-3", "x")[1];
    //return 0;
    //auto start = std::chrono::high_resolution_clock::now();
    //vector<double> zeros = equ_solver("1/(x-10)-3*x", "x", -100, 100); 
    //auto end = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double> elapsed = end - start;
    //cout << bisection_method("1/(x)-10", "x", -10.1, 10) << '\n';
    //cout << elapsed.count() << '\n';
    
    //for (double i : zeros) printf("%lf ", i);

    //cout << exp_solver("-2^2");
    //for (string i : post_fix_conv("(-x)^2")) cout << i;
    //return 0;
    string command;
    double solution;
    vector<double> zeros;
    //for (char i : exp) cout << i;

    cout << "Type an equation to find the roots of. Use variable 'x'. Type \"exit\" to close.\n";
    while (true) {
        cout << ">>";
        cin >> command;
        if (command == "exit") break;
        auto start = std::chrono::high_resolution_clock::now();
        command_reader(command);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        cout << elapsed.count() << '\n';
    }
    


    return 0;
}