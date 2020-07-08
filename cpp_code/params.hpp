#ifndef __PARAMSINIT__
#define __PARAMSINIT__ 1
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


#endif
