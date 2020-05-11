#include "CoordinatesVelocitiesAndBondsInitializer.h"

CoordinatesVelocitiesAndBondsInitializer::CoordinatesAndVelocitiesInitializer(MDRunOutput& mdoutput,
    const MDParameters& parameters, std::string coordinatesFileName, std::string bondsFileName)
    : output(mdOutput), 
      par(parameters), 
      fileName(std::move(coordinatesFileName)),
      fileName_bonds(std::move(bondsFileName)) {

}

void CoordinatesVelocitiesAndBondsInitializer::initialize(std::vector<double>& positions, 
    std::vector<double>& velocities, std::vector<std::pair<int, int>> bonds) {
    
    // Initialize positions and velocities as usual
    this->initialize(positions, velocities);
    
    // Load bonds file if needed
    if (par.chainMdType != ChainSimType::noChains &&
        par.xvInitialization != InitialXVGenerator::generateInitialX) {
        fin_bonds.open(fileName_bonds, std::ios::in);
        if (fin_bonds.bad())
            throw std::runtime_error("can't open" + fileName_bonds);
    }
    
    if (par.chainMdType != ChainSimType::noChains &&
        par.xvInitialization != InitialXVGenerator::generateInitialX) {
        unsigned int nat = par.numberAtoms;
        
        // Read in number of bonds
        unsigned int nbonds;
        fin_bonds >> nbonds;
        
        if (nbonds > nat * nat)
            throw std::runtime_error("NBONDS (" + fileName_bonds + ") = " + 
                                     std::tostring(nbonds) + " > NATOM^2");
        
        // Fill bonds vector from bonds file
        for (int counter = 0; counter < nbonds; counter++) {
            unsigned int i, j;
            fin_bonds >> i;
            fin_bonds >> j;
            
            bonds.push_back(std::make_pair<unsigned int, unsigned int>(i, j));
        }
    } else {
        bonds.clear(); // make sure no bond terms are calculated if not needed
    }
    
}
