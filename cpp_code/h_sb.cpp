#include "h_funcs.hpp"
#include "utils.hpp"

namespace h_funcs{
  
  void h_sb(const utils::vector_type &psi_in, utils::vector_type &psi_out,const utils::Params &p){
    
#pragma omp parallel for
    for( int i = 0; i < p.num_bath; ++i){
      ublas::vector<int> num_state((unsigned)p.osc_count);
      //solve psi out for each ind i
      //sum over ther different values of Q_k 
      utils::update_number_state(num_state,i,p);

      for(int k = 0; k < p.osc_count; ++ k){
      
	double factor_k = p.osc_c(k)*sqrt(1.0/(2.0* p.osc_freq(k)));
	int n_k = num_state(k);
      
	//adagger_k raises the kth term in num states for psi_in
	//----this should just be a reindexing, which takes
	// [n1,n2-1, n3...] -> sqrt(n2) [ n1,n2,n3...]
	if( n_k > 0){
	  int a_dag_idx = i - pow(p.mode_count,k);
	  psi_out(i) += sqrt(n_k)*factor_k*psi_in(a_dag_idx);
	  psi_out(i+p.num_bath) -= sqrt(n_k)*factor_k*psi_in(a_dag_idx+p.num_bath);
	}

	//a_k lowers the kth term in psi_in
	//[n1,n2+1,n3..]-> sqrt(n2+1) [n1,n2,n3]
	if( (n_k+1) < p.mode_count){
	  int a_idx = i + pow(p.mode_count,k);
	  psi_out(i) += sqrt(n_k+1)*factor_k*psi_in(a_idx);
	  psi_out(i+p.num_bath) -= sqrt(n_k+1)*factor_k*psi_in(a_idx+p.num_bath);

	}
	
      }//end for k
      
    }//end for i
      

  }
}
