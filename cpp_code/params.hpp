#ifndef __PARAMSINIT__
#define __PARAMSINIT__ 1
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
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
  int osc_count,mode_count;

  //oscillator freqency vector
  boost::numeric::ublas::vector<double> osc_freq;
  //oscillator coupling vector
  boost::numeric::ublas::vector<double> osc_c;
  //number of terms
  int num_bath;

};


void params_init(Params &p,double epsilon,double lambda, double T, double wc, int osc_c, int mode_c);


#endif
