//This is a driver for regression.cpp, please see
//regression.h for information

#include <iostream>
#include <string>
#include"scheme.h"

using namespace std;

int main(){
    param p_input;
    //input parameters
    
    /*cout << "Number of time steps: ";
     cin >> p_input.N_t; cout << p_input.N_t;
     cout << "\nNumber of intervals for y: ";
     cin >> p_input.N_y; cout << p_input.N_y;
     cout << "\nNumber of intervals for e: ";
     cin >> p_input.N_e; cout << p_input.N_e;
     cout << "\nEnd point of domain for y:";
     cin >> p_input.L_y; cout << " [- " << p_input.L_y << ",+" << p_input.L_y << "]";
     cout << "\nEnd point of domain for e: ";
     cin >> p_input.L_e; cout << "[0," << p_input.L_e << "]";
     cout << "\nNumber of Monte Carlo samples: ";
     cin >> p_input.N_q; cout << "[0," << p_input.N_q << "]";
     cout << "\nTime horizon: ";
     cin >> p_input.T; cout << "[0," << p_input.T << "]";
     cout << "\n0(explicit) <= theta <=1(implicit) : ";
     cin >> p_input.theta; cout << p_input.theta;
     cout << "\ngamma: ";
     cin >> p_input.gamma; cout << p_input.gamma;
     cout << "\nmu: ";
     cin >> p_input.mu; cout << p_input.mu;
     cout << "\nrho = (1-1/(2*a)): ";
     cin >> p_input.rho; cout << p_input.rho;
     cout << "\nalpha: ";
     cin >> p_input.alpha; cout << p_input.alpha << "\n";
     cout << "Terminal Condition (Write \"carbon\" to choose V(T,y,e) = -alpha*e*1_{y>0}, otherwise press any key follwed by enter key to choose a constant terminal condition.): ";
     std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
     string term;
     getline(cin,term);
     
     p_input.const_term = -2;
     if (term.compare("carbon")==0){
     p_input.term_fnc = carbon;
     cout << "\nV(T,y,e) = -alpha*e*1_{y>0}";
     }
     else{
     p_input.term_fnc = constant;
     cout << "\nV(T,y,e) = ";
     cin >> p_input.const_term;
     cout << p_input.const_term << "\n";
     }
     cout << "Write \"Neumann\" choose  Neumann Artificial Boundary Condition (ABC); otherwise press any key follwed by enter key and the ABC will automatically be \"Transparent Dirichlet\":";
     std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
     string bd;
     getline(cin,bd);
     
     if (bd.compare("Neumann")==0){
     p_input.bd = Nmnn;
     cout << "\nABC: " << "Neumann" ;
     }
     else{
     p_input.bd = Drchlt;
     cout << "\nABC: " << "Transparent Dirichlet \n V(t," << -p_input.L_y << ",e) = (T-t)/(4rho) \nV(t," << p_input.L_y << ",e) = -alpha*e+(1-alpha)^2(T-t)/(4rho) \nAt (t, e," << p_input.L_e << ") ABC is obtained by linear extrapolation";
     }*/
    //Parameters preset
    p_input.N_t = 1000;
    cout << "\nNumber of time step: " << p_input.N_t;
    p_input.N_y = 100;
    cout << "\nNumber of intervals for y: " << 2*p_input.N_y;
    p_input.N_e = 100;
    cout << "\nNumber of intervals for e: " << p_input.N_e;
    p_input.L_y = 20;
    cout << "\nDomain of y: [" << -p_input.L_y << ", +" << p_input.L_y << "]";
    p_input.L_e = 20;
    cout << "\nDomain of e: [0 , +" << p_input.L_e << "]";
    p_input.N_q = 1000;
    cout << "\nNumber of samples for MC: " << p_input.N_q;
    p_input.T = 10;
    cout << "\nTime horizon: [0, " << p_input.T << "]";
    p_input.theta = .9;
    cout << "\n0(explicit) <= theta <=1(implicit) : " << p_input.theta;
    p_input.gamma = .5;
    cout << "\ngamma: " << p_input.gamma;
    p_input.mu = .1;
    cout << "\nmu: " << p_input.mu;
    p_input.rho = .9;
    cout << "\nrho=(1-1/(2*a)): " << p_input.rho;
    p_input.alpha = .1;
    cout << "\nalpha: " << p_input.alpha;
    p_input.constant = 0;
    cout << "\nTerminal Condition: " << "V(T,y,e) =" << p_input.constant << "-" << p_input.alpha << "*e*1_{y>0}";
    p_input.bd = Drchlt;
    cout << "\nArtificial Boundary Condition (ABC):";
    if (p_input.bd==Nmnn){
        cout << "Neumann" ;
    }
    else{
        cout << "Transparent Dirichlet \nV(t," << -p_input.L_y << ",e) = (T-t)/(4rho) \nV(t," << p_input.L_y << ",e) = -alpha*e+(1-alpha)^2(T-t)/(4rho) \nAt (t,e,0) ABC is obtained by some approximation of optimal control";
    }
    
    
    //Other  parameters
    p_input.log="scheme.log";
    p_input.latex_file="heatmap.tex";
    p_input.dy = p_input.L_y/static_cast<double>(p_input.N_y);
    p_input.de = p_input.L_e/static_cast<double>(p_input.N_e);
    p_input.dt = p_input.T/p_input.N_t;
    p_input.c = .5*p_input.gamma*p_input.gamma*p_input.dt/(p_input.dy*p_input.dy);
    cout << "\nc from CFL condition (<.5 necessary for explicit, i.e. whenever theta=0): " <<  p_input.c << "\n";
    try{
        implementation scheme(p_input);
        make_latex heatmap(p_input);
    }
    catch(...){
        cout << "\nException Caught\n";
    }
    cout << "\nWritten in scheme.log\n";
    return 0;
}