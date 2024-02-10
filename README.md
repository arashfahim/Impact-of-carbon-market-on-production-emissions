# Impact of carbon market on production emissions
This repository contains the solution to the optimal control problem for the research paper https://arxiv.org/abs/2312.03665
The code solves the control problem in three different cases. Each case is implemented in a different folder. For more information, see the reseach paper.
The numerical scheme is solving the HJB equation using Clark-Nicholson scheme. The boundary conditions on the computational window is chosen to accomodate the degeneracy of one of the variables.

$$0=-\partial_t V^{(2)}-(\mu_t-\lambda_t\gamma_t)V^{(2)}_y$$ 



$$-\frac12 \gamma^2 V^{(2)}_{yy}-\sup_{q\ge 0}\theta(t,q,V^{(2)}_y,V^{(2)}_e)$$

$$V^{(2)}(T,y,e)= -\alpha{1}_{\{y>0\}}e$$
