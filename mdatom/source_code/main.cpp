#include <MDSimulation.h>
#include <ChainMDSimulation.h>
#include <iostream>

/*!
 * main function.
 * The program expects one argument for the input file (parameters), an one optional argument for
 * the  file with the initial coordinates.
 */

int main(int argc, char* argv[]) {
  std::cout << "hi" << std::endl;

    switch (argc) {
        case 2: break;
        case 3: break;
        case 4: break;
        default:
        std::cerr << "Usage: mdatom input_file [coordinate_file] [bonds_file] > output \n";
            return 1;
    }

    std::string parameterFile = argv[1];
    std::string coordinatesFile = ""; // NB: might be empty
    std::string bondsFile = "";
    if (argc > 2)
        coordinatesFile = argv[2];
    
    if (argc > 3)
        bondsFile = argv[3];
    
    // Check whether we are running a normal simulation or a chain simulation
    if (argc <= 3) {
      MDSimulation md(std::cout);
      try {
          md.performSimulation(parameterFile, coordinatesFile);
      }
      catch (std::exception& e) {
          std::cerr << e.what();
          return 1;
      }
    } else {
      ChainMDSimulation md(std::cout);
      try {
          md.performSimulation(parameterFile, coordinatesFile, bondsFile);
      }
      catch (std::exception& e) {
          std::cerr << e.what();
          return 1;
      }
    }
    

    return 0;
}
