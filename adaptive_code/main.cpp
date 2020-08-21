#include <iostream>
#include <vector>
#include <boost/program_options.hpp>


using namespace std;
namespace po = boost::program_options;

int main(int ac, const char ** av){

  //Declare supported options
  po::options_description desc("Allowed options: ");
  po::desc.add_options()
    ("help", "produce help message")
    ("epsilon", po::value<double>(), "Sets the energy bias")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

  if(vm.count("help")){
    std::cout << desc << endl;
  }

  if(vm.count("compression")) {
    cout << "essilon was set to "
	 << vm["epsilon"].as<double>() << endl;
  }else{
    cout << "epsilon not set, assuming value of zero"
	 << endl;
  }
  
  cout << "Hello world"<< endl;


}
