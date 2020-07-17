#include "utils.hpp"
#include "h_funcs.hpp"
namespace h_funcs{
  void h_s(const utils::vector_type &psi_in, utils::vector_type &psi_out,const utils::Params &p){

    int num_bath = p.num_bath;
#pragma omp parallel for
    for(int i = 0; i < (p.num_bath); ++i){
      psi_out(i) += psi_in(i+num_bath);
      psi_out(i+num_bath) += psi_in(i);
      psi_out(i) += psi_in(i)*p.epsilon;
      psi_out(i+num_bath) -= psi_in(i+num_bath)*p.epsilon;
     
    }


  }

  void h_b(const utils::vector_type &psi_in, utils::vector_type &psi_out,const utils::Params &p){
    int num_bath = p.num_bath;
    double omega_sum = ublas::norm_1(p.osc_freq)*0.5;
#pragma omp parallel for
    for(int i = 0; i < num_bath; ++i){
      ublas::vector<int> num_state(p.osc_count);
      utils::update_number_state(num_state,i,p);
      double energy_i = ublas::inner_prod(num_state,p.osc_freq) + omega_sum;
      psi_out(i) += psi_in(i)*energy_i;
      psi_out(i+num_bath) += psi_in(i + num_bath)*energy_i;
    }
  }
  
  
}
