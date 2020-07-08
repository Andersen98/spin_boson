#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "params.hpp"


  
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
