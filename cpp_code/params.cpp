#include "params.hpp"



void params_init(Params &p,double epsilon,double lambda, double T, double wc, unsigned osc_c, unsigned mode_c){
  p.epsilon = epsilon;
  p.lambda = lambda;
  p.T = T;
  p.wc = wc;
  p.wf = 2*wc;
  p.w0 = p.wf/(2*osc_c);
  p.osc_count = osc_c;
  p.mode_count = mode_c;

  p.osc_freq = new double[osc_c];
  p.osc_c = new double[osc_c];

  for (unsigned i = 0; i < osc_c; ++i){
    p.osc_freq[i] = (1+i)*(p.wf/osc_c);
    p.osc_c[i] = lambda/osc_c;
  }
  p.num_bath = mode_c^osc_c;
}
