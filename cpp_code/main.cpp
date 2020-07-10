#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "params.hpp"
#include "psi_naught.hpp"

void init_psi(boost::numeric::ublas::vector<double> *psi, const Params &p);

int main() {
  namespace ublas = boost::numeric::ublas;
  Params p;
  //param_init_args
  //---params,eps,lamda,T,wc,osc_c,mode_c
  params_init(p, 0, 2.5, .2, .25, 30, 2);
  ublas::vector<double> psi(2* (unsigned)p.num_bath);
  init_psi(psi,p);
  
}
