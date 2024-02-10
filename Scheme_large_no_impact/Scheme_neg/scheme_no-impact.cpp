//Please see regression.cpp for further information

#include"scheme.h"

using namespace std;


//construct the mesh and its operators
void implementation::make_S(std::vector<S_i> &S){
    for (unsigned int j = 0; j <= 2*p.N_y; j++){
        for (unsigned int i = 0; i <= 2*p.N_e; i++){
            S.push_back(S_i(p,j,i));
        }
    }
}



void implementation::set_to_terminal(std::vector<S_i> &S){
    for(unsigned int i=0; i<S.size(); i++){
        S.at(i).v = p.constant + ((S.at(i).y>0) ? -p.alpha*S.at(i).e : 0);
    }
}

//subtrating two solutions on a mesh
vector<S_i> implementation::subtract(std::vector<S_i> &S_1, std::vector<S_i> &S_2){
    vector<S_i> D;
    S_i D0;
    unsigned long int min = std::min(S_1.size(), S_2.size());
    for(unsigned int i=0; i<min; i++){
        D0.y = S_1.at(i).y;
        D0.e = S_1.at(i).e;
        D0.v = S_1.at(i).v-S_2.at(i).v;
        D0.v_e = S_1.at(i).v_e-S_2.at(i).v_e;
        D0.v_y = S_1.at(i).v_y-S_2.at(i).v_y;
        D.push_back(D0);
    }
    if (log_file.is_open()){
        log_file << "\nDifference of two solution methods calculated!";
    }
    cout << "\nDifference of two solution methods calculated!";
    return D;
}




// start to build implicit scheme components


void implementation::make_f0(std::vector<S_i> &S, unsigned int i){
    for (unsigned int j = 1; j<=2*p.N_y-1; j++){
        //Set y to the solution of previous block *b +Value function
        f0.at(j) = (1-2*p.c*(1-p.theta))*S.at(index(j,i)).v+p.c*(1-p.theta)*(S.at(index(j+1,i)).v+S.at(index(j-1,i)).v);
    }
}




void implementation::update_f0(unsigned int i){
    f0.at(1) += p.c*p.theta*y0.at(index(0,i));
    f0.at(2*p.N_y-1) += p.c*p.theta*y0.at(index(2*p.N_y,i));
}


void implementation::update_lu0(){
    u0.at(1)= 1+p.c*p.theta*((p.bd==Nmnn)?1:2);
    for(unsigned int j=2; j<=2*p.N_y-2; j++){
        l0.at(j)=-p.c*p.theta/u0.at(j-1);
        u0.at(j)=1+2*p.c*p.theta+l0.at(j)*p.c*p.theta;
    }
    l0.at(2*p.N_y-1)=-p.c*p.theta/u0.at(2*p.N_y-2);
    u0.at(2*p.N_y-1)=1+p.c*p.theta*((p.bd==Nmnn)?1:2)+l0.at(2*p.N_y-1)*p.c*p.theta;
}

void implementation::update_y0(unsigned int i){
    y0.at(index(1,i))=f0.at(1);
    for(unsigned int j=2; j<=2*p.N_y-1; j++){
        y0.at(index(j,i))=f0.at(j)-l0.at(j)*y0.at(index(j-1,i));
    }
    y0.at(index(2*p.N_y-1,i)) /= u0.at(2*p.N_y-1);
    for (unsigned int j=2*p.N_y-2; j>=1; j--){
        y0.at(index(j,i))=(y0.at(index(j,i))+p.c*p.theta*y0.at(index(j+1,i)))/u0.at(j);
    }
}


void implementation::reset_f0(){
    for (unsigned int j=0; j<=2*p.N_y-1; j++){
        //initializing vectors for solving Besse relaxation scheme in even steps
        f0.at(j)=0;
    }
}


void implementation::reset_y0(){
    for (unsigned int j=0; j<=2*p.N_y-1; j++){
        for (unsigned int i=0; i<=2*p.N_e; i++){
            y0.at(index(j,i))=0;
        }
    }
}

// implicit scheme for diffusion part
void implementation::implicit_scheme(std::vector<S_i> &S, unsigned int n){// implicit method for heat equation in even steps
    for (unsigned int i = 0; i<=2*p.N_e; i++){
        make_f0(S,i);
        if (p.bd==Drchlt){ // Dirichlet  boundary conditions must be set before the step calculations
            make_B(S,y0,n,i,grand);
            if ((n % 200 == 1) && (i==2*p.N_e)){
                if (log_file.is_open()){
                    log_file << "\nDirichlet boundary condition in step " << n/2;
                }
                cout << "\nDirichlet boundary condition in step " << n/2;
            }
            update_f0(i);//adjusting f0 when Dirichlet boundary conditions
        }
        update_lu0();
        if ((n % 200 == 1) && (i==2*p.N_e)){
            if (log_file.is_open()){
                log_file << "\nLU decomposition for implicit(explicit) scheme for diffusion term at step " << n/2 << "!";
            }
            cout << "\nLU decomposition for implicit(explicit) scheme for diffusion term at step " << n/2 << "!";
        }
        update_y0(i);
        //Setting boundary condition for heat equation
        if (p.bd==Nmnn){//Neumann boundary condition must be set after the step calculations
            make_B(S,y0,n,i,grand);
            if ((n % 200 == 1) && (i==2*p.N_e)){
                if (log_file.is_open()){
                    log_file << "\nNeumann boundary condition in step " << n/2;
                }
                cout << "\nNeumann boundary condition in step " << n/2;
            }
        }
        reset_f0();
    }
    for (unsigned int j = 0; j <= 2*p.N_y; j++){
        for (unsigned int i=0; i<=2*p.N_e; i++){
            S.at(index(j,i)).v = y0.at(index(j,i));
        }
    }
    reset_y0();
    if (n % 200 == 1){
        if (log_file.is_open()){
            log_file << "\nSolution of diffusion term at step " << n/2 << "!";
        }
        cout << "\nSolution of diffusion term at step " << n/2 << "!";
    }
}

//Building boundary conditions in each step, shared with the implicit scheme for heat equation and CN_relaxed scheme
void implementation::make_B(std::vector<S_i> &S, std::vector<double> &y, unsigned int n, unsigned int i, MTRand_int32 &grand){
    if (p.bd==Drchlt){// boundary condition for carbon terminal
        if(is_odd(n)){
            y.at(index(0,i)) = p.constant+(p.T-p.dt*static_cast<double>((n+1)/2))/(4*p.rho);
            y.at(index(2*p.N_y,i)) = p.constant-p.alpha*p.de*(static_cast<double>(i)-static_cast<double>(p.N_e))+(1-p.alpha)*(1-p.alpha)*(p.T-p.dt*static_cast<double>((n+1)/2))/(4*p.rho);
            // n+1 is to avoid discotinuity in the boundary by setting the boundary condition to the one for last step
        }
        else{
            y.at(index(0,i)) = p.constant+(p.T-p.dt*static_cast<double>((n)/2))/(4*p.rho);
            y.at(index(2*p.N_y,i)) = p.constant-p.alpha*p.de*(static_cast<double>(i)-static_cast<double>(p.N_e))+(1-p.alpha)*(1-p.alpha)*(p.T-p.dt*static_cast<double>((n)/2))/(4*p.rho);
            if (i==1){//add Dirichlet boundary condition for e=0 in CN scheme
                y.at(index(0,0)) = p.constant+(p.T-p.dt*static_cast<double>((n)/2))/(4*p.rho);
                y.at(index(2*p.N_y,0)) = p.constant+(1-p.alpha)*(1-p.alpha)*(p.T-p.dt*static_cast<double>((n)/2))/(4*p.rho);
                double e , x, w, u;
                int r ,m;
//                cout << "\nHere1: ";
                for (unsigned int j=1; j<=2*p.N_y-1; j++){
                    e = psi.at(index(j,0))*p.dt;
                    w = floor(e/p.de);
                    m = static_cast<int>(w);
                    w = e/p.de - w;
                    x= (psi.at(index(j,0))+p.mu)*p.dt;
                    u = floor(x/p.dy);
                    r = static_cast<int>(u);
                    u = x/p.dy - u;
                    if (m>2*p.N_e)
                    cout << "\nj=" << j << ", j+r=" << j+r << ", j+r+1=" << j+r+1 << ", m=" << m << ", m+1=" << m+1 << ", psi=" << psi.at(index(j,0));
                    y.at(index(j,0)) = (j+r<2*p.N_y) ? ( (1-p.rho*psi.at(index(j,0)))*psi.at(index(j,0))*p.dt +  (1-w)*((1-u)*S.at(index(j+r,m)).v+u*S.at(index(j+r+1,m)).v) + w*((1-u)*S.at(index(j+r,m+1)).v+u*S.at(index(j+r+1,m+1)).v) ) : ( p.constant-p.alpha*p.de*static_cast<double>(m)+(1-p.alpha)*(1-p.alpha)*(p.T-p.dt*static_cast<double>((n)/2))/(4*p.rho) );
                }
//                cout << "\nHere2" ;
            }
        }
    }
    else{// homogenous Neumann boundary condition
        y.at(index(0,i)) = y.at(index(1,i));
        y.at(index(2*p.N_y,i)) = y.at(index(2*p.N_y-1,i));
        if ((!is_odd(n)) && (i == 2*p.N_e-1)){
            for (unsigned int j=0; j<=2*p.N_y; j++){
                y.at(index(j,2*p.N_e)) = y.at(index(j,2*p.N_e-1));
            }
        }
    }
}

// start to  build Crank_Nicolson with relaxation components



//making the space derivatives of v, v_y and v_e
void implementation::make_v_der(std::vector<S_i> &S){
    //V_e and V_y
    for (unsigned int j=1; j<=2*p.N_y-1; j++){
        for (unsigned int i=0; i<=2*p.N_e-1; i++){
            S.at(index(j,i)).v_e=(S.at(index(j,i+1)).v-S.at(index(j,i)).v)/p.de;
            S.at(index(j,i)).v_y=(S.at(index(j+1,i)).v-S.at(index(j-1,i)).v)/(2*p.dy);
        }
    }
    for (unsigned int i=0; i<=2*p.N_e-1; i++){
        S.at(index(0,i)).v_e=(S.at(index(0,i+1)).v-S.at(index(0,i)).v)/p.de;
        S.at(index(2*p.N_y,i)).v_e=(S.at(index(2*p.N_y,i+1)).v-S.at(index(2*p.N_y,i)).v)/p.de;
    }
    for (int j=1; j<=2*p.N_y-1; j++){
        S.at(index(j,2*p.N_e)).v_y=(S.at(index(j+1,2*p.N_e)).v-S.at(index(j-1,2*p.N_e)).v)/(2*p.dy);
    }
    //artificially setting derivatives at the boundary for plotting purposes
    for (unsigned int i=0; i<=2*p.N_e; i++){
        S.at(index(0,i)).v_y=S.at(index(1,i)).v_y;
        S.at(index(2*p.N_y,i)).v_y=S.at(index(2*p.N_y-1,i)).v_y;
    }
    for (int j=0; j<=2*p.N_y; j++){
        S.at(index(j,2*p.N_e)).v_e=S.at(index(j,2*p.N_e-1)).v_e;
    }
    //    cout << "\n" << S.at(index(2*p.N_y-1,1)).v_e;
}



//make optimal control
void implementation::make_psi(std::vector<S_i> &S, unsigned int n){
    for (unsigned int j = 0; j<=2*p.N_y; j++){
        for (unsigned int i=0; i<=2*p.N_e; i++){
            psi.at(index(j,i)) = max(0,(1 + S.at(index(j,i)).v_e + S.at(index(j,i)).v_y)/(2*p.rho));
        }
    }
}


//making the relaxation term psi for even steps;
//it turns out that relaxation term shows up as coefficient of a linear system, i.e. a, b, and c. See documentation for more details
void implementation::make_coeff(std::vector<S_i> &S, unsigned int n, unsigned int k){
    for (unsigned int j = 1; j<=2*p.N_y-1; j++){
        a.at(j) = p.mu*p.dt/(2*p.dy) + p.dt*psi.at(index(j,k))/(4*p.dy);
        b.at(j) = p.dt*psi.at(index(j,k))/(2*p.de);
        c.at(j) = p.dt*psi.at(index(j,k))/2;
    }
    if (p.bd==Nmnn){//Neumann
        if (*first_block){
            dd.at(1)=1+p.theta*a.at(1);
            ee.at(1)=-p.theta*a.at(1);
            for (unsigned j=2; j<=2*p.N_y-2; j++){
                dd.at(j)=1;
                ee.at(j)=-p.theta*a.at(j);
                cc.at(j)=p.theta*a.at(j);
            }
            dd.at(2*p.N_y-1)=1-p.theta*a.at(2*p.N_y-1);
            cc.at(2*p.N_y-1)=p.theta*a.at(2*p.N_y-1);
        }
        else{
            dd.at(1)=1+p.theta*(a.at(1)+b.at(1));
            ee.at(1)=-p.theta*a.at(1);
            for (unsigned j=2; j<=2*p.N_y-2; j++){
                dd.at(j)=1+p.theta*b.at(j);
                ee.at(j)=-p.theta*a.at(j);
                cc.at(j)=p.theta*a.at(j);
            }
            dd.at(2*p.N_y-1)=1+p.theta*(b.at(2*p.N_y-1)-a.at(2*p.N_y-1));
            cc.at(2*p.N_y-1)=p.theta*a.at(2*p.N_y-1);
        }
    }
    else{//Dirichlet
        if (*first_block){
            dd.at(1)=1-p.theta*b.at(1);
            ee.at(1)=-p.theta*a.at(1);
            for (unsigned j=2; j<=2*p.N_y-2; j++){
                dd.at(j)=1-p.theta*b.at(j);
                ee.at(j)=-p.theta*a.at(j);
                cc.at(j)=p.theta*a.at(j);
            }
            dd.at(2*p.N_y-1)=1-p.theta*b.at(2*p.N_y-1);
            cc.at(2*p.N_y-1)=p.theta*a.at(2*p.N_y-1);
        }
        else{
            dd.at(1)=1-p.theta*b.at(1);
            ee.at(1)=-p.theta*a.at(1);
            for (unsigned j=2; j<=2*p.N_y-2; j++){
                dd.at(j)=1-p.theta*b.at(j);
                ee.at(j)=-p.theta*a.at(j);
                cc.at(j)=p.theta*a.at(j);
            }
            dd.at(2*p.N_y-1)=1-p.theta*b.at(2*p.N_y-1);
            cc.at(2*p.N_y-1)=p.theta*a.at(2*p.N_y-1);
        }
    }
    if ((n % 200 ==0) && (k==0)){
        if (log_file.is_open()){
            log_file << "\nCoefficients of (implicit) Crank-Nicolson built at step " << n/2 << "!";
        }
        cout << "\nCoefficients of (implicit) Crank-Nicolson built at step " << n/2 << "!";
    }
}

void implementation::make_f1(std::vector<S_i> &S, unsigned int n, unsigned int k){
    if (p.bd==Nmnn){
        for (unsigned int j = 1; j<=2*p.N_y-1; j++){
            //Set y to the solution of previous block *b +Value function
            f1.at(j) = p.theta*b.at(j)*y1.at(index(j,k+1));//b*solution of previous block
            f1.at(j) += (1-(1-p.theta)*b.at(j))*S.at(index(j,k)).v + (1-p.theta)*a.at(j)*S.at(index(j+1,k)).v - (1-p.theta)*a.at(j)*S.at(index(j-1,k)).v + (1-p.theta)*b.at(j)*S.at(index(j,k+1)).v + c.at(j);//+Value function at current block
        }
    }
    else{//make f for Dirichlet
        for (unsigned int j = 1; j<=2*p.N_y-1; j++){
            f1.at(j) = -p.theta*b.at(j)*y1.at(index(j,k-1));
            f1.at(j) += (1+(1-p.theta)*b.at(j))*S.at(index(j,k)).v + (1-p.theta)*a.at(j)*S.at(index(j+1,k)).v - (1-p.theta)*a.at(j)*S.at(index(j-1,k)).v - (1-p.theta)*b.at(j)*S.at(index(j,k-1)).v + c.at(j);
        }
    }
}

void implementation::update_f1(std::vector<S_i> &S, unsigned int n, unsigned int k){//used only in the case of Dirichlet boundary condition
    f1.at(1) -= a.at(1)*p.theta*y1.at(index(0,k));
    f1.at(2*p.N_y-1) += a.at(2*p.N_y-1)*p.theta*y1.at(index(2*p.N_y,k));
}


void implementation::update_lu1(){
    u1.at(1)=dd.at(1);
    for(unsigned int j=2; j<=2*p.N_y-1; j++){
        l1.at(j)=cc.at(j)/u1.at(j-1);
        u1.at(j)=dd.at(j)-l1.at(j)*ee.at(j-1);
    }
}

void implementation::update_y1(unsigned int k){
    y1.at(index(1,k))=f1.at(1);
    for(unsigned int j=2; j<=2*p.N_y-1; j++){
        y1.at(index(j,k))=f1.at(j)-l1.at(j)*y1.at(index(j-1,k));
    }
    y1.at(index(2*p.N_y-1,k)) = y1.at(index(2*p.N_y-1,k)) /u1.at(2*p.N_y-1);
    for (unsigned int j=2*p.N_y-2; j>=1; j--){
        y1.at(index(j,k))=(y1.at(index(j,k))-ee.at(j)*y1.at(index(j+1,k)))/u1.at(j);
    }
}

void implementation::reset_f1(){
    for (unsigned int j=0; j<=2*p.N_y; j++){
        //initializing vector f1 for solving CN relaxation scheme in even steps
        f1.at(j)=0;
    }
}

void implementation::reset_y1(){
    for (unsigned int j=0; j<=2*p.N_y; j++){
        //initializing vector y1 for solving CN relaxation scheme in even steps
        for (unsigned int i=0;i<=2*p.N_e;i++){
            y1.at(index(j,i))=0;
        }
    }
}

//Crank_Nicolson algorithm with relaxation proposed by Besse

void implementation::CN_relaxed_scheme(std::vector<S_i> &S, unsigned int n){
    make_v_der(S);//Required to build psi
    make_psi(S,n);
    bool block1 = true;
    first_block = &block1;
    for (unsigned int i=1; i<=2*p.N_e; i++){
        unsigned int k = (p.bd==Nmnn) ? 2*p.N_e-i : i;//Neumann goes backward to solve the block matrix while Dirichlet goes forward. k switches the direction based on the choice of ABC
        make_coeff(S,n,k);
        if (p.bd==Drchlt){// Dirichlet  boundary conditions must be set before the step calculations
            make_B(S,y1,n,i,grand);
            if ((i==2*p.N_e) && (n % 200 == 0)){
                if (log_file.is_open()){
                    log_file << "\nDirichlet boundary condition in step " << n/2;
                }
                cout << "\nDirichlet boundary condition in step " << n/2;
            }
        }
        make_f1(S,n,k);
        if(p.bd==Drchlt){
            update_f1(S,n,k);
        }
        update_lu1();
        update_y1(k);
        //setting the value function to solution of linear system
        if (p.bd==Nmnn){// Dirichlet  boundary conditions must be set after the step calculations
            make_B(S,y1,n,k,grand);
            if ((k==0) && (n % 200 == 0)){
                if (log_file.is_open()){
                    log_file << "\nNeumann boundary condition in step " << n/2;
                }
                cout << "\nNeumann boundary condition in step " << n/2;
            }
        }
        if(*first_block){
            block1 = false;
            first_block = &block1;
        }
        reset_f1();
    }
    for (unsigned int j = 0; j <= 2*p.N_y; j++){
        for (unsigned int i=0; i<=2*p.N_e; i++){
            S.at(index(j,i)).v = y1.at(index(j,i));
        }
    }
    reset_y1();
    if (*terminal_CN_relaxed){
        bool terminal = false;//accomodating terminal condition
        terminal_CN_relaxed = &terminal;
    }
    if (n % 200 == 0){
        if (log_file.is_open()){
            log_file << "\nSolution of the advection-nonlinear term calculated at step " << n/2 << "!";
        }
        cout << "\nSolution of the advection-nonlinear term calculated at step " << n/2 << "!";
    }
}



// making the correction term

void implementation::make_tau(std::vector<S_i> &S, std::vector<S_i> &tau){
    S_i D0;
    for(unsigned int i=0; i<S.size(); i++){
        D0.y = S.at(i).y;
        D0.e = S.at(i).e;
        D0.v = S.at(i).v_y;
        D0.v_e = 0;
        D0.v_y = 0;
        tau.at(i).v=D0.v;
    }
}

// making optimal control at all points of the mesh

void implementation::make_q(std::vector<S_i> &S, std::vector<S_i> &q){
    for (unsigned int j=0; j<=2*p.N_y; j++){
        for (unsigned int i=0; i<=2*p.N_e; i++){
            q.at(index(j,i)).v = psi.at(index(j,i));
        }
    }
    //    S_i D0;
    //    for(unsigned int i=0; i<S.size(); i++){
    //        D0.y = S.at(i).y;
    //        D0.e = S.at(i).e;
    //        D0.v = psi.at(index(i/(p.N_e+1),i%(p.N_e+1)));
    //        D0.v_e = 0;
    //        D0.v_y = 0;
    //        q.at(i).v=D0.v;//updating marginal optimal production
    //    }
}

// test of int being odd

bool implementation::is_odd(int n){
    if ( n % 2== 0 )
        return false;
    else
        return true;
}



//max

double implementation::max(double a, double b){
    return (a<b) ?b :a;
}



//indexing of vectors


unsigned int implementation::index(unsigned int j, unsigned int i){
    unsigned int index = 0;
    index += j*(2*p.N_e+1) + i;
    return index;
}



void implementation::write2file(vector<S_i> &S, const bool &file, std::string file_name){
    if (file){
        value_fnc_csv.open(file_name+".csv");
        value_fnc_csv  << "y" << ", " << "e" << ", " <<  "y(t;y;e)" << "\n";
        for (unsigned int j=0; j<=2*p.N_y; j++){
            for (unsigned int i=0; i<=2*p.N_e; i++){
                value_fnc_csv <<  S.at(index(j,i)).y << ",  " << S.at(index(j,i)).e << ",  " <<  S.at(index(j,i)).v << "\n";
            }
        }
        value_fnc_csv.close();
    }
    else{
        value_fnc_dat.open(file_name+".dat");
        for (unsigned int j=p.N_y/2; j<=3*p.N_y/2; j++){
            for (unsigned int i=p.N_e/2; i<=3*p.N_e/2; i++){
                value_fnc_dat << S.at(index(j,i)).y << "  " << S.at(index(j,i)).e << "  " <<  S.at(index(j,i)).v << "\n";
                //                if (S.at(index(j,i)).v<0)
                //                    cout <<"\nj=" << j << ", i=" << i << ", q^*=" << S.at(index(j,i)).v;
            }
            value_fnc_dat<< "\n";
        }
        value_fnc_dat.close();
    }
}

// running implementation with a given set of parameters

implementation::implementation(const param &p_input){
    time(&timer);
    p = p_input;
    
    //declarations for wrting into files of tyoe csv and dat + file names
    const bool csv=true;
    const bool dat=false;
    const std::string term_cond="terminal";
    const std::string solution="solution";
    const std::string correction="correction";
    const std::string control="optimal_control";
    
    log_file.open(p.log.c_str());
    //Writing the inputs into the log file
    if (log_file.is_open()){
        log_file << "Number of time intervals: " << p.N_t;
        log_file  << "\nNumber of intervals for y: " << 2*p.N_y;
        log_file  << "\nNumber of intervals for e: " << 2*p.N_e;
        log_file  << "\nDomain of y: [" << -p.L_y << ", +" << p.L_y << "]";
        log_file  << "\nDomain of e: [" << -p.L_e << ", +" << p.L_e << "]";
        log_file << "\nNumber of samples for MC: " << p.N_q;
        log_file  << "\nTime horizon: [0, " << p.T << "]";
        log_file  << "\n0(explicit) <= theta <=1(implicit) : " << p.theta;
        log_file  << "\ngamma: " << p.gamma;
        log_file  << "\nmu: " << p.mu;
        log_file  << "\nrho=(1-1/(2*a)): " << p.rho;
        log_file  << "\nalpha: " << p.alpha;
        log_file << "\nTerminal Condition: " << "V(T,y,e) = -alpha*e*1_{y>0}";
        log_file << "\nArtificial Boundary Condition (ABC):";
        if (p.bd==Nmnn){
            log_file  << "Homogeneous Neumann" ;
        }
        else{
            log_file << "Transparent Dirichlet \nV(t," << -p_input.L_y << ",e) = (T-t)/(4rho) \nV(t," << p_input.L_y << ",e) = -alpha*e+(1-alpha)^2(T-t)/(4rho) \nAt (t,e,0) ABC is obtained by some approximation of optimal control";
        }
        log_file << "\nc from CFL condition (<.5 necessary for explicit, i.e. whenever theta=0): " <<  p.c << "\n";
    }
    //boolean variable to manage to isolate the calculations in the terminal step
    bool terminal = true;
    terminal_CN_relaxed = &terminal;
    //boolean variable to manage the first block of linear system in CN scheme
    //initailize structures to 0
    //initializing vectors for solving implicit heat equation in odd
    for (unsigned int j=0; j<=2*p.N_y; j++){
        //initializing vectors for solving implicit heat equation in odd steps
        l0.push_back(0);u0.push_back(0);f0.push_back(0);
        //initializing vectors for solving CN relaxation scheme in even steps
        a.push_back(0);b.push_back(0);c.push_back(0);
        dd.push_back(0);ee.push_back(0);cc.push_back(0);
        l1.push_back(0);u1.push_back(0);f1.push_back(0);
        for (unsigned int i=0; i<=2*p.N_e; i++){//psi should be saved for the next step, y1 saves the value function at current until we are done with the value function oin previous step
            y0.push_back(0);psi.push_back(0);y1.push_back(0);
        }
    }
    make_S(S1);// builds the mesh  for the scheme
    make_S(tau);
    make_psi(S1,0);
    make_S(q);
    //    make_S(S2);
    set_to_terminal(S1);//Set the solution equal to terminal condition
    //    set_to_terminal(S2);
    if (log_file.is_open()){
        log_file << "\nMesh built!";
    }
    cout << "\nMesh built!";
    //Writing terminal condition to a file
    write2file(S1,csv,term_cond.c_str());
    write2file(S1,dat,term_cond.c_str());

    //implementation of the loop-frog scheme, odd steps heat equation, even steps Crank-Nicolson
    for (int n = 2*p.N_t-1; n>=0; n--){
        if (is_odd(n)) {// solving the implicit or explicit heat equation
            //            p.bd=Drchlt;
            implicit_scheme(S1,n);
            //            p.bd=Nmnn;
            //            implicit_scheme(S2,n);
        }
        else{// solving the nonlinear (Crank_Nicolson with relaxation or puse explicit)
            //            p.bd=Drchlt;
            CN_relaxed_scheme(S1,n);
            //            p.bd=Nmnn;
            //            CN_relaxed_scheme(S2,n);
        }
    }
    //    D = subtract(S1,S2);
    make_v_der(S1);//To calculated optimal control and tau
    make_tau(S1,tau);
    make_psi(S1,0);
    make_q(S1,q);
    
    if (log_file.is_open())
        log_file << "\nWriting to files.....";
    cout << "\nWriting to files.....";
    
    write2file(S1,csv,solution.c_str());
    write2file(tau,csv,correction.c_str());
    write2file(q,csv,control.c_str());
    if (log_file.is_open()){
        log_file << "\nCSV files created!";
    }
    cout << "\nCSV files created!";
    
    write2file(S1,dat,solution.c_str());
    write2file(tau,dat,correction.c_str());
    write2file(q,dat,control.c_str());
    if (log_file.is_open()){
        log_file << "\nDAT files created!";
    }
    cout << "\nDAT files created!";
    if (log_file.is_open()){
        log_file << "\n" << difftime(time(NULL), timer) << " seconds for generation\n";
    }
    cout << "\n" << difftime(time(NULL), timer) << " seconds for generation";
    log_file.close();
}

make_latex::make_latex(const param &p_input){
    p = p_input;
    latex_file.open(p.latex_file.c_str());
    write_file();
    latex_file.close();
}

void make_latex::write_file(){
    if (latex_file.is_open()){
        latex_file << "\\documentclass[border=10pt]{article}\n\\usepackage{pgfplots}\n \\usepackage{caption}\n\\usepackage{tikz}\n\\usepackage{subcaption}                  \n\\addtolength{\\textwidth}{2cm}\n\\addtolength{\\hoffset}{-1cm}\n%\\addtolength{\\textheight}{2cm}\n%\\addtolength{\\voffset}{-1cm}\n\\pgfplotsset{compat=newest}\n\\begin{document}\n\\[  -\\partial_t V - \\frac{\\gamma^2}{2}\\partial_{yy} V -\\mu \\partial_{y} V -\\frac{1}{4\\rho}(1+\\partial_e V + \\partial_y V)_+^2 = 0  \\]   \nThe parameters of the problem:\\\\ \n";
        latex_file << "$\\gamma=" << p.gamma << "$, $\\alpha=" << p.alpha << "$, $\\beta=1$, $\\pi=q(1-q)$, $\\eta(q)=\\lambda(q)=q$, $T=" << p.T << "$, $\\mu=" << p.mu << "$ and $\\varrho=" << p.rho << "$ ($a=" << 1/(2*(1-p.rho))<< "$),  bounded domain for $(t,y,e)$ is $[0," << p.T << "]\\times[" << -p.L_y << ", " << p.L_y << "]\\times[" << -p.L_e << ", +" << p.L_e << "]$ with ";
        if (p.bd==Nmnn){
            latex_file << "homogeneous Neumann boundary condition ";
            latex_file << "\n\\[\nV_e(t," << p.L_e << ",y)=V_y(t,e," << -p.L_y << ")= V_y(t,e," << p.L_y << ")=0\n\\]";
        }
        else{
            latex_file << "transparent Dirichlet boundary condition ";
            latex_file << "\n\\[V(t,e," << -p.L_y << ")=\\frac{(T-t)}{4\\varrho},\\hspace*{.3cm} V(t,e," << p.L_y << ")=-\\alpha e+(1-\\alpha)^2\\frac{(T-t)}{4\\varrho}\n\\]";
        }
        latex_file << "\nand terminal condition $V(T,e,y)=";
        latex_file << "-\\alpha e1_{\\{y>0\\}}$.";
        latex_file << "\n\\begin{figure}";
        latex_file << "\n\\begin{subfigure}[b]{0.5\\textwidth}";
        latex_file << "\n\\centering";
        latex_file << "\n\\resizebox{\\linewidth}{!}{";
        latex_file << "\n\\begin{tikzpicture}";
        latex_file << "\n\\begin{axis}[xlabel={$y$}, ylabel={$e$}, grid=none, colormap/jet, colorbar, \n%point meta min=3, point meta max=-10, \ncolorbar style={\n /pgf/number format/precision=3,\n %ytick={3.777,3.778,3.779,3.780,3.781,3.782,3.783,3.784}\n},\nview={0}{90}]";
        latex_file << "\n\\addplot3[surf,shader=flat] file {terminal.dat};";
        latex_file << "\n\\end{axis}";
        latex_file << "\n\\end{tikzpicture}";
        latex_file << "\n}";
        latex_file << "\n\\caption{Terminal condition $V(T,\\cdot)$}";
        latex_file << "\n\\end{subfigure}";
        latex_file << "\n\\begin{subfigure}[b]{0.5\\textwidth}";
        latex_file << "\n\\centering";
        latex_file << "\n\\resizebox{\\linewidth}{!}{";
        latex_file << "\n\\begin{tikzpicture}";
        latex_file << "\n\\begin{axis}[xlabel={$y$}, ylabel={$e$}, grid=none, colormap/jet, colorbar, \n%point meta min=0, point meta max=-10, \ncolorbar style={\n /pgf/number format/precision=3, %ytick={3.777,3.778,3.779,3.780,3.781,3.782,3.783,3.784}\n},\nview={0}{90}]";
        latex_file << "\n\\addplot3[surf,shader=flat] file {solution.dat};";
        latex_file << "\n\\end{axis}";
        latex_file << "\n\\end{tikzpicture}";
        latex_file << "\n}";
        latex_file << "\n\\caption{Solution $V(0,\\cdot)$}";
        latex_file << "\n\\end{subfigure}";
        latex_file << "\n\\end{figure}";
        latex_file << "\n\\begin{figure}";
        latex_file << "\n\\begin{subfigure}[b]{0.5\\textwidth}";
        latex_file << "\n\\centering";
        latex_file << "\n\\resizebox{\\linewidth}{!}{";
        latex_file << "\n\\begin{tikzpicture}";
        latex_file << "\\begin{axis}[xlabel={$y$}, ylabel={$e$}, grid=none, colormap/jet,               colorbar,\n%point meta min=-.7, point meta max=.1, \ncolorbar style={\n/pgf/number format/precision=3,\n%ytick={.1110,.1111,.1112,.1113,.1114}\n},";
        latex_file << "\nview={0}{90}]";
        latex_file << "\n\\addplot3[surf,shader=flat] file {correction.dat};";
        latex_file << "\n\\end{axis}";
        latex_file << "\n\\end{tikzpicture}";
        latex_file << "\n}";
        latex_file << "\n\\caption{Correction term $\\tau=V_y$}";
        latex_file << "\n\\end{subfigure}";
        latex_file << "\n\\begin{subfigure}[b]{0.5\\textwidth}";
        latex_file << "\n\\centering";
        latex_file << "\n\\resizebox{\\linewidth}{!}{";
        latex_file << "\n\\begin{tikzpicture}";
        latex_file << "\\begin{axis}[xlabel={$y$}, ylabel={$e$}, grid=none, colormap/jet,               colorbar,\n%point meta min=0, point meta max=.6, \ncolorbar style={\n/pgf/number format/precision=3,\n%ytick={0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1}\n},";
        latex_file << "\nview={0}{90}]";
        latex_file << "\n\\addplot3[surf,shader=flat]  file {optimal_control.dat};";
        latex_file << "\n\\end{axis}";
        latex_file << "\n\\end{tikzpicture}";
        latex_file << "\n}";
        latex_file << "\n\\caption{Optimal Control $\\frac{1}{4}(1+V_e+V_y)_+$}";
        latex_file << "\n\\end{subfigure}";
        latex_file << "\n\\end{figure}";
        latex_file << "\n\\end{document}\n";
    }
};



//standard normal random number generation
double implementation::gauss(MTRand_int32 &grand){
    double u = (static_cast<double>(grand()) + .5)/(4294967296.);
    double v = (static_cast<double>(grand()) + .5)/(4294967296.);
    return sqrt(-2*log(u))*cos(2*3.14159265358979323846*v); //Box-Muller
}

