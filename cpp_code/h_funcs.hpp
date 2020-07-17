#ifndef __H_FUNCS__
#define __H_FUNCS__ 1


#include "utils.hpp"

namespace ublas = boost::numeric::ublas;

namespace h_funcs{
  
  void h_sb(const utils::vector_type &psi_in,  utils::vector_type &psi_out,const utils::Params &p);

  void h_s(const utils::vector_type &psi_in,  utils::vector_type &psi_out,const utils::Params &p);

  void h_b(const utils::vector_type &psi_in,  utils::vector_type &psi_out,const utils::Params &p);
}

#endif
