#include "ChainMDSimulation.h"
#include "TrajectoryFileWriter.h"
#include "ChainMDRun.h"
#include "CoordinatesAndVelocitiesInitializer.h"
#include "ParameterIO.h"
#include "ParameterValidityChecker.h"
#include "RandomNumberGenerator.h"

ChainMDSimulation::ChainMDSimulation(std::ostream& outputStream)
  : MDSimulation(outputStream) {

}

void ChainMDSimulation::performSimulation(const std::string& parameterFile, 
  const std::string& coordinateFile, const std::string& bondsFile) {
  
  MDParameters par = ParameterIO::readParameters(parameterFile);
  performSimulation(par, coordinateFile, bondsFile);
}

// copied from MDSimulation, only changed init class to accomodate for bonds
void ChainMDSimulation::performSimulation(const MDParameters& par, 
  const std::string& coordinateFile, const std::string& bondsFile) {
  
  parameters = par;
  prepareRun();
  checkParameterValidity();
  initializeCoordinatesVelocitiesAndBonds(coordinateFile, bondsFile);
  executeMDIterations();
  finalizeRun();
}

void ChainMDSimulation::initializeCoordinatesVelocitiesAndBonds(const std::string& coordinateFile,
  const std::string& bondsFile) {
  
  // Use initializer class to build coords, velocities and bonds
  CoordinatesVelocitiesAndBondsInitializer xvbInitializer(output, parameters, coordinateFile, bondsFile);
  xvbInitializer.initialize(positions, velocities, bonds);
}

void ChainMDSimulation::executeMDIterations() {
    TrajectoryFileWriter trajectoryWriter(parameters, "coords.final", "coords.traj");
    trajectoryWriter.writeBeforeRun();

    timer.mdStart();
    ChainMDRun mdRun(parameters, output, trajectoryWriter);
    mdRun.run(positions, velocities, bonds);
    timer.mdEnd();

    printRadialDistribution(mdRun.getRadialDistribution());
    trajectoryWriter.writeFinalCoordinates(positions, velocities);
}
