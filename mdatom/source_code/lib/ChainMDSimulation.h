#ifndef CHAINMDSIMULATION_H
#define CHAINMDSIMULATION_H

#include "MDSimulation.h"
#include "CoordinatesVelocitiesAndBondsInitializer.h"

/*!
 * This class launches a chain MD simulation starting from parameters and, optionally, coordinates and bonds.
 */
class ChainMDSimulation : public MDSimulation {
  public:
    explicit ChainMDSimulation(std::ostream& outputStream);
    
    // Simulation from parameters, coords and bonds in files
    void performSimulation(const std::string& parFile, const std::string& coordFile,
                           const std::string& bondsFile);
    // Simulation from already loaded parameters; coords and bonds are in files
    void performSimulation(const MDParameters& par, const std::string& coordinateFile,
                           const std::string& bondsFile);
  
  private:
    // build internal data structures of coords, velocities and bonds from files
    void initializeCoordinatesVelocitiesAndBonds(const std::string& coordinateFile, 
                                                 const std::string& bondsFile);
    
    std::vector<std::vector<bool>> bonds; //bonds[i][j] = 1 if atom i is bonded to atom j
};

#endif // CHAINMDSIMULATION_H
