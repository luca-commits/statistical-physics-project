#include "CoordinatesVelocitiesAndBondsInitializer.h"
#include "BinaryIO.h"
#include "RandomNumberGenerator.h"
#include <cmath>
#include <stdexcept>
#include <utility>
#include <string>

CoordinatesVelocitiesAndBondsInitializer::CoordinatesVelocitiesAndBondsInitializer(MDRunOutput& mdoutput, const MDParameters& parameters, std::string coordinatesFileName, std::string bondsFileName)
    : CoordinatesAndVelocitiesInitializer(mdoutput, parameters, coordinatesFileName),
      fileName_bonds(bondsFileName) {

}

void CoordinatesVelocitiesAndBondsInitializer::initialize(std::vector<double>& positions, 
    std::vector<double>& velocities, std::vector<std::vector<bool>>& bonds) {
    
    // Initialize positions and velocities as usual
    CoordinatesAndVelocitiesInitializer::initialize(positions, velocities);
    
    // Check whether we are simulating a chain and using an automatically generated mesh
    if (par.chainMdType != ChainSimType::noChains &&
        par.xvInitialization != InitialXVGenerator::generateInitialX) {
        fin_bonds.open(fileName_bonds, std::ios::in);
        if (fin_bonds.bad())
            throw std::runtime_error("can't open" + fileName_bonds);
        
        unsigned int nat = par.numberAtoms;
        
        // Read in number of bonds
        unsigned int nbonds;
        fin_bonds >> nbonds;
        std::cout << "no of bonds: " << nbonds << std::endl << std::flush;
        // Throw exception if there are too many bonds
        if (nbonds > nat * nat)
            throw std::runtime_error("NBONDS (" + fileName_bonds + ") = " + 
                                     std::to_string(nbonds) + " > NATOM^2");
        
        bonds.resize(nat);
        
        for (int i = 0; i < nat; i++) {
          bonds[i].resize(nat);
          for (int j = 0; j < nat; j++) {
            bonds[i][j] = false;
          }
        }
        
        // Fill bonds vector from bonds file
        for (int counter = 0; counter < nbonds; counter++) {
            unsigned int i, j;
            fin_bonds >> i;
            fin_bonds >> j;
            
            bonds[i][j] = true;
            bonds[j][i] = true;
        }
    } else if(par.chainMdType != ChainSimType::noChains &&
              par.xvInitialization == InitialXVGenerator::generateInitialX) {
        throw std::runtime_error("Automatic generation of chains not supported!");
    } else {
        bonds.clear(); // make sure no bond terms are calculated if not needed
    }
    
}
