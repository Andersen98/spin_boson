#include "utils.hpp"

namespace utils{
  
  void update_number_state(boost::numeric::ublas::vector<int> &num_state, int idx, const Params &p){

    num_state.clear();
    int mode_count = p.mode_count;
    int digit = 0;
    while(idx !=0){
      int mod = idx%mode_count;
      num_state(digit) = mod;
      idx = idx/mode_count;
      ++digit;
    }
  }


  
  void params_init(Params &p,double epsilon,double lambda, double T, double wc, int osc_c, int mode_c){
    p.epsilon = epsilon;
    p.lambda = lambda;
    p.T = T;
    p.beta = 1.0/T;
    p.wc = wc;
    p.wf = 2.0*wc;
    p.w0 = p.wf/(2.0*osc_c);
    p.osc_count = osc_c;
    p.mode_count = mode_c;
    p.osc_c.resize(osc_c);
    p.osc_freq.resize(osc_c);
    for (int i = 0; i < osc_c; ++i){
      p.osc_freq(i) = i*(p.wf-p.w0)/((double)(osc_c-1))+p.w0;
      p.osc_c(i) = sqrt(2.0*lambda/osc_c)*p.osc_freq(i);
    }
    p.num_bath = pow(mode_c,osc_c);
  }


  
  void init_psi(vector_type &psi, const Params &p){


    double omega_sum = ublas::norm_1(p.osc_freq)*0.5;
    double psi_norm = 0;

    #pragma omp parallel for 
    for(int i = 0; i < p.num_bath; ++i){
      //update the number state
      ublas::vector<int> number_state((unsigned)p.osc_count);
      utils::update_number_state(number_state,i,p);
      //std::cout << number_state<<std::endl;
      //compute sum of energies for oscillators
      double energy_i = ublas::inner_prod(p.osc_freq,number_state) + omega_sum;
      psi(i) = exp(-p.beta*energy_i/2.0);
      
    }
    
    #pragma omp parallel for reduction(+:psi_norm)
    for(int i = 0; i < p.num_bath; ++i){
      psi_norm += pow(psi(i),2);
    }

    psi_norm = sqrt(psi_norm);


    #pragma omp parallel for
    for(int i = 0; i < p.num_bath; ++i){
	psi(i) = psi(i)/psi_norm;
    }
    
  }

}
