#ifndef COORDINATESVELOCITIESANDBONDSINITIALIZER_H
#define COORDINATESVELOCITIESANDBONDSINITIALIZER_H

#include "CoordinatesAndVelocitiesInitializer.h"

// Extends CoordinatesAndVelocitiesInitializer by an input file specifying
// the bonds between the given atoms. The file should have the following format:
//   NO_OF_BONDS
//   PAIR1_ATOM1 PAIR1_ATOM2
//   PAIR2_ATOM1 PAIR2_ATOM2
//   ...
// where we index the atoms starting at 0 as specified in the input file
// for the coordinates.
// Automatic generation of bonds not supported at the moment.
class CoordinatesVelocitiesAndBondsInitializer : public CoordinatesAndVelocitiesInitializer {
  public:
    CoordinatesVelocitiesAndBondsInitializer(MDRunOutput& mdoutput, const MDParameters& parameters, 
      std::string coordinatesFileName, std::string bondsFileName);
    void initialize(std::vector<double>& positions, std::vector<double>& velocities, 
                    std::vector<std::vector<bool>>& bonds);

  private:
    std::string fileName_bonds;
    std::ifstream fin_bonds;
};

#endif // COORDINATESVELOCITIESANDBONDSINITIALIZER_H
