#include "params.hpp"
#include <cmath>


void params_init(Params &p,double epsilon,double lambda, double T, double wc, int osc_c, int mode_c){
  p.epsilon = epsilon;
  p.lambda = lambda;
  p.T = T;
  p.wc = wc;
  p.wf = 2*wc;
  p.w0 = p.wf/(2*osc_c);
  p.osc_count = osc_c;
  p.mode_count = mode_c;

  //  boost::numeric::ublas::vector<double> p.osc_freq(osc_c);
  // boost::numeric::ublas::vector<double> p.osc_c(osc_c);
  p.osc_c.resize(osc_c);
  p.osc_freq.resize(osc_c);
  for (int i = 0; i < osc_c; ++i){
    p.osc_freq(i) = (1+i)*(p.wf/osc_c);
    p.osc_c(i) = lambda/osc_c;
  }
  p.num_bath = pow(mode_c,osc_c);
}
