#include "utils.hpp"
#include "h_funcs.hpp"
int main() {
  namespace ublas = boost::numeric::ublas;
  utils::Params p;
  //param_init_args
  //---params,eps,lamda,T,wc,osc_c,mode_c  
  utils::params_init(p, 0, 2.5, .2, .25, 28, 2);
  std::cout << "wk" << p.osc_freq << std::endl;
  std::cout << "ck" << p.osc_c <<std::endl<<std::endl;
  utils::vector_type  psi(2* (unsigned)p.num_bath);
  utils::init_psi(psi,p);
  utils::vector_type psi_prime(2*(unsigned)p.num_bath);
  

  h_funcs::h_sb(psi, psi_prime, p);
  
  h_funcs::h_s(psi, psi_prime, p);
  
  h_funcs::h_b(psi,psi_prime, p);
  
  //std::cout << "original:" << std::endl << psi << std::endl;
  //std::cout << "prime: " << std::endl <<psi_prime << std::endl <<std::endl;
}
