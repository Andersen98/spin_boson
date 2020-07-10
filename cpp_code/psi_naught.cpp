#include  "psi_naught.hpp"
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
namespace ublas = boost::numeric::ublas;

void init_psi(ublas::vector<double> &psi, const Params &p){
  ublas::vector<int> number_state((unsigned)p.osc_count);
  double omega_sum = ublas::norm_1(p.osc_freq)*0.5;
  
  for(int i = 0; i < p.num_bath; ++i){

    //update the number state

    
    for( int j = 0; j < p.osc_count; ++j){
      if( number_state (j) < p.mode_count){
	break;
      }
      number_state(j) = 0;
      number_state(j+1) += 1;
    }
    
    //std::cout << number_state << std::endl;

    double energy_i = ublas::inner_prod(p.osc_freq,number_state) + omega_sum;

    psi(i) = exp(-energy_i);
    
    //increment the number state
    number_state(0) += 1;
						   
  }
  psi /= ublas::norm_2(psi);
  
}
