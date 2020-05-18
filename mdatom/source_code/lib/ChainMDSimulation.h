#ifndef CHAINMDSIMULATION_H
#define CHAINMDSIMULATION_H

#include "MDSimulation.h"
#include "CoordinatesVelocitiesAndBondsInitializer.h"

/*!
 * This class launches a MD simulation starting from parameters and, optionally, coordinates.
 */
class ChainMDSimulation : public MDSimulation {
  public:
    explicit ChainMDSimulation(std::ostream& outputStream);
    
    void performSimulation(const std::string& parFile, const std::string& coordFile,
                           const std::string& bondsFile = "");
    void performSimulation(const MDParameters& par, const std::string& coordinateFile,
                           const std::string& bondsFile = "");
  
  private:
    void initializeCoordinatesVelocitiesAndBonds(const std::string& coordinateFile, 
                                                 const std::string& bondsFile);
    
    std::vector<std::vector<bool>> bonds;//bonds[i][j] = 1 if atom i is bonded to atom j
};

#endif // CHAINMDSIMULATION_H
