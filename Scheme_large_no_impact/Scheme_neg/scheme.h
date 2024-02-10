//AMERICAN PUT OPTION APPROXIMATION

//Created by Michael Angileri for the University of Michigan
//Department of Mathematics Research Experience for Undergraduates
//during Summer 2013 under the advisement of Arash Fahim

//Please see file mtrand.h for information regarding the random
//number generator. Ownership of said file is not being claimed.

//This program generates sample paths of a stock price given inputs
//from driver.cpp and returns the approximated value of an American
//put option using a least-squares algorithm



#ifndef BESSE_H
#define BESSE_H
#define _USE_MATH_DEFINES
#include <cmath>
#include<vector>
#include<random>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<string>
#include<iomanip>
#include<ctime>
#include"mtrand.h"

enum bndr_cndtn {Nmnn, Drchlt};

struct param{//Holds the parameters of the problem
	unsigned int N_t, N_y, N_e, N_q;//number of intervals for y, number of intervals for e, number of samples for Monte Carlo for optimal production
	double T, dt, dy , de, L_y, L_e , theta, gamma, mu, rho, alpha, c, constant, lambda;
    //maturity, T/N_t,  L_y/N_y, L_e/N_e, emission diffusion, emission drift, rho=(1-1/(2*a)) and a is the parameter for exponential utility, penalty, c = .5*gamma^2*dt/dy^2
    //Use double term only if terminal condition is constant
    enum bndr_cndtn bd;
    std::string log, latex_file;
    
};

struct S_i{//Holds the coordinates of ONE point and the value fucntion and its derivatives at this point
	double y, e, v, v_y, v_e ;// time, y, e and value function and derivatives wrt e and y
    S_i(){y=0; e=0; v=0; v_e=0; v_y=0;};
    S_i(double y0, double e0, double v0, double v0_y, double v0_e){y=y0; e=e0; v=v0; v_y=v0_y; v_e=v0_e;};
    S_i operator=(S_i S){
        return S_i(S.y, S.e, S.v, S.v_e, S.v_y);
    }
    S_i operator+(S_i S){
        return S_i(S.y, S.e, v+S.v, v_e+S.v_e, v_y+S.v_y);
    }
    S_i operator-(S_i S){
        return S_i(S.y, S.e, v-S.v, v_e-S.v_e, v_y-S.v_y);
    }
    S_i(const param &p, unsigned int j, unsigned int i){//assigns values to the member of S_i based on therminal condition
        y = p.dy * (static_cast<double>(j)-static_cast<double>(p.N_y));
        e = p.de * (static_cast<double>(i)-static_cast<double>(p.N_e));
        v = 0; v_e = 0; v_y = 0;
    }
};



class implementation{
	param p;
//    enum bndr_cndtn bndr_type = p.bd;
	unsigned int cur_time;
    bool const * terminal_CN_relaxed;//handles the exceptions of calculations for terminal condition
    bool const * first_block;//handles the exceptions of calculations for first block of linear system in Crank_Nicolson scheme
    // saving the mesh points and the value of the function in diffrenet methods (S1, S2)
	std::vector<S_i> S1;
//    std::vector<S_i> S2;
    std::vector<S_i> tau;//The correction term
    std::vector<S_i> q;//optimal control (production rate)
    // build a mesh with the comparison of different methods
    std::vector<S_i> D;
    std::vector<S_i> subtract(std::vector<S_i> &S_1, std::vector<S_i> &S_2);
    //vectors for solving tridiag system for implicit heat equation
    std::vector<double> y0;
    std::vector<double> f0;
    std::vector<double> l0;
    std::vector<double> u0;

    //coefficients in relaxing the nonlinear part
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> dd;
    std::vector<double> cc;
    std::vector<double> ee;
    //Ax=f: x is the value function at time n and y depends on the value function at time n+1. A is not tridiagonal, but can be handled by a tridiag system. See documentation for details.
    std::vector<double> f1;//Ax=f
    std::vector<double> y1;// Ly=f, Ux=y
    //vectors for solving non-tridiag system for relaxed Crank-Nicolson

    std::vector<double> l1;
    std::vector<double> u1;
    std::vector<double> psi;//Holds the value of optimal control (production)
    std::vector<double> approx_q;//Approximate q to construct a transparent boundary condition at e=0
    
	time_t timer;
	
	//construct the mesh
	void make_S(std::vector<S_i> &S);
    void make_B(std::vector<S_i> &S, std::vector<double> &y,unsigned int n, unsigned int i, MTRand_int32 &grand);// impose Dirichlet boundary condition
    void set_to_terminal(std::vector<S_i> &S);
    //construct the vectors c and d for tridiagonal system in implicit scheme for heat equation
    void update_lu0();
    void update_y0(unsigned int i);
    void make_f0(std::vector<S_i> &S, unsigned int i);
    void update_f0(unsigned int i);//updates f0 when the boundary condition is Dirichlet
    void reset_f0();void reset_y0();
    //making the relaxation term
    void make_coeff(std::vector<S_i> &S, unsigned int n, unsigned int i);
    // solve linear system Ax=y for relaxed Crank-Nicolson. The system is split into blocks.  Block eqns are solved recursively from the last block, whic is a tridiag system
    //then knowing the solution of a specific block, the next block becomes tridiag too. See documenty for details.
    //Coefficients needed at each block, i determines the block number
    void update_lu1();
    void update_y1(unsigned int k);
    
    void make_f1(std::vector<S_i> &S, unsigned int n, unsigned int i);
    void update_f1(std::vector<S_i> &S, unsigned int n, unsigned int i);
    void reset_y1(); void reset_f1();
    void update_c1(std::vector<S_i> &S, unsigned int i);
    void update_d1(std::vector<S_i> &S, unsigned int i);
    
    //make optimal control
    void make_psi(std::vector<S_i> &S, unsigned int n);
    
    //switching between two operators in the relaxed splitting scheme
    void implicit_scheme(std::vector<S_i> &S, unsigned int n);// one-step scheme for implicit heat eqn (S is a mesh point, n is the time step (for the purpose of time-dependent artificial boundary conditions), bd is the  type of artificial boundary condition)
    void CN_relaxed_scheme(std::vector<S_i> &S, unsigned int n);// one-step relaxed Crank-Nicolson(Bassee)(S is a mesh point, n is the time step (for the purpose of time-dependent artificial boundary conditions, bd is the  type of artificial boundary condition)
    // Make ABC for e=0
    //Buliding the derivatives of v wrt e and y
    void make_v_der(std::vector<S_i> &S);
    //correction term
    void make_tau(std::vector<S_i> &S, std::vector<S_i> &tau);
    //Markovian optimal control at all points of the mesh
    void make_q(std::vector<S_i> &S, std::vector<S_i> &q);
    //writing to file
    std::ofstream value_fnc_dat, value_fnc_csv, psi_csv, log_file;//  file to write solutions and other stuff into
    void write2file(std::vector<S_i> &S, const bool &file, std::string file_name);
    //updating the derivatives
 
	
	//returns L*
//	void transpose(std::vector<double> &L);
    
    //Is int n odd?
    bool is_odd(int n);

    //indexing of vectors
    unsigned int index(unsigned int j, unsigned int i);

	//Gaussian random
    MTRand_int32 grand;
	double gauss(MTRand_int32 &grand);

	//int power function
//	unsigned int pow(unsigned int one, unsigned int two);
    
    //double power function
//    double pow(double one, unsigned int two);

    
    //Gaussian cdf
    double cdf(double x);
    
    //max
    double max(double a , double b);

public:
	implementation(const param &p_input);
	implementation(){};
    ~implementation(){};
};


class make_latex{
    param p;
    std::ofstream latex_file;
    void write_file();
    
public:
    make_latex(const param &p_input);
    make_latex(){};
};

#endif //BESSE_H