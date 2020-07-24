#include "utils.hpp"
#include "h_funcs.hpp"
#include <boost/numeric/odeint.hpp>

namespace ublas = boost::numeric::ublas;
const int count_osc = 12;
const int count_mode = 4;
//read https://aip.scitation.org/doi/full/10.1063/1.5126945

void write_population( const utils::vector_type &x , const double t );
void rhs(const utils::vector_type &x, utils::vector_type &dxdt, const double /* t */);

void test_ode(const utils::vector_type &x, const utils::Params &p);

int main() {

  utils::Params p;
  //param_init_args
  //---params,eps,lamda,T,wc,osc_c,mode_c  
  utils::params_init(p, 0, 2.5, .2, .25,  count_osc,count_mode);
  // std::cout << "wk" << p.osc_freq << std::endl;
  //  std::cout << "ck" << p.osc_c <<std::endl<<std::endl;
 
  
  utils::vector_type  psi(4* (unsigned)p.num_bath);
  utils::init_psi(psi,p);
  utils::vector_type psi_prime(4*(unsigned)p.num_bath);

  //test_ode(psi,p);
  boost::numeric::odeint::integrate( rhs, psi , 0.0 , 12.0 , 0.0001 , write_population );
  
  //std::cout << "original:" << std::endl << psi[0] << std::endl;
  //std::cout << "prime: " << std::endl <<psi_prime[2*p.num_bath] << std::endl <<std::endl;
}

void test_ode(const utils::vector_type &x, const utils::Params &p){
  
  std::vector<double> dxdt(4* (unsigned)p.num_bath);

  h_funcs::h_sb(x, dxdt, p);  
  h_funcs::h_s(x, dxdt, p);
  h_funcs::h_b(x,dxdt, p);
  
  
  for (std::vector<double>::const_iterator i = x.begin(); i != x.end(); ++i)
    std::cout << *i << ' ';
  std::cout << std::endl;

  for (std::vector<double>::const_iterator i = dxdt.begin(); i != dxdt.end(); ++i)
    std::cout << *i << ' ';
  std::cout << std::endl;
  
}
  

void write_population( const utils::vector_type &x , const double t ){
  utils::Params p;
  //param_init_args
  //---params,eps,lamda,T,wc,osc_c,mode_c  
  utils::params_init(p, 0, 2.5, .2, .25,  count_osc,count_mode);

  double total;
  int n_half = 2*p.num_bath;

  for(int i = 0; i < p.num_bath; ++i){
    total += pow(x[i],2)+pow(x[i+n_half],2);
    total -= pow(x[i+p.num_bath],2)+pow(x[i + p.num_bath+n_half],2);
  }
  
  std::cout << t << '\t' << total << '\t'  << std::endl;
  
}


void rhs(const utils::vector_type &x, utils::vector_type &dxdt, const double /* t */){
  utils::Params p;
  //param_init_args
  //---params,eps,lamda,T,wc,osc_c,mode_c  
  utils::params_init(p, 0, 2.5, .2, .25, count_osc, count_mode);

#pragma omp parallel for
  for(int i = 0; i < (4*p.num_bath); ++i){
    dxdt[i] = 0;
  }
  
  h_funcs::h_sb(x, dxdt, p);  
  h_funcs::h_s(x, dxdt, p);
  h_funcs::h_b(x,dxdt, p);
    
}


