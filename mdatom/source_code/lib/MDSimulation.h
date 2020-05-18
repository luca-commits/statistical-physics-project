#ifndef MDSIMULATION_H
#define MDSIMULATION_H

#include "MDParameters.h"
#include "MDRunOutput.h"
#include "Timer.h"
#include "AveragedRadialDistribution.h"
#include <string>
#include <vector>
#include <utility>
#include <iostream>

/*!
 * This class launches a MD simulation starting from parameters and, optionally, coordinates.
 */
class MDSimulation {
  public:
    /*! Constructor; the output of the MD simulation will be redirected to outputStream. */
    explicit MDSimulation(std::ostream& outputStream);

    /*! Perform a simulation based on a parameter file and an (optional) coordinate file. */
    void performSimulation(const std::string& parFile, const std::string& coordFile = "");
    /*! Perform a simulation based parameters and an (optional) coordinate file. */
    void performSimulation(const MDParameters& par, const std::string& coordinateFile = "");
  
  protected:
    void prepareRun();
    void checkParameterValidity();
    void initializeCoordinatesAndVelocities(const std::string& coordinateFile);
    virtual void executeMDIterations();
    virtual void printRadialDistribution(const AveragedRadialDistribution& radialDistribution);
    virtual void finalizeRun();
  
    MDParameters parameters;
    MDRunOutput output;
    Timer timer;
    std::vector<double> positions, velocities;
    std::vector<std::pair<int,int>> bonds; //contains all pairs of bonded atoms in the chain
};

#endif // MDSIMULATION_H
