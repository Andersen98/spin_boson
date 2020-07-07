#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>

//the physical and simulation parameters
struct Params{
  //inergy is in units of Delta, so we are left with the
  //physical params epsilon lambda and T
  double epsilon,lambda,T;

  //parameters for the bath spectrum
  //wc - characteristic bath frequency
  //wf - the cutoff frequency
  //w0 - the initial freqency
  double wc,wf,w0;

  //simulation parameters - the oscillator count is the number of oscillators to use for the truncated bath
  //mode count gives the count for the oscillator modes (ex. if count = 3 then we have states |0> |1> |2>
  unsigned osc_count,mode_count;

  //oscillator freqency omega
  double *osc_freq;
  //oscillator coupling array
  double *osc_c;
  //number of terms
  unsigned num_bath;

};


void params_init(Params &p,double epsilon,double lambda, double T, double wc, unsigned osc_c, unsigned mode_c);


  
int main() {
  namespace ublas = boost::numeric::ublas;
  Params p;
  params_init(p, 0, 2.5, .2, .25, 10, 2);
  ublas::vector<double> v (3);
  for (unsigned i = 0; i < v.size(); ++i){
    v (i) = i;
  }
  std::cout << v << std::endl;
}

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
