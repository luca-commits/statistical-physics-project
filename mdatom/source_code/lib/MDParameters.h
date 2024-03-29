#ifndef MDPARAMETERS_H
#define MDPARAMETERS_H

#include <string>

enum class FinalCoordinateFileFormat {
    binary,
    ascii
};

enum class TrajectoryFileFormat {
    binary,
    ascii
};

enum class InitialXVGenerator {
    generateInitialX,
    xFromAscii,
    xFromBinary,
    xvFromAscii,
    xvFromBinary
};

enum class SimulationType {
    constantEnergy,
    constantTemperature
};

enum class ChainSimType {
    noChains,			//normal md simulation
    noAngles,		    //only bond and angle terms
    complete			//bond, angle and torsion terms
};
/*!
 * This struct contains the parameters for a MD simulation for the Lennard-Jones model.
 */
struct MDParameters {
    static MDParameters defaultParameters();

    std::string title;
    int numberAtoms;
    double atomicMass;
    SimulationType mdType;
    ChainSimType chainMdType;
    double boxSize[3];
    int numberMDSteps;
    double initialTime;
    double timeStep;
    double initialTemperature;
    double targetTemperature;
    double temperatureCouplingTime;
    int randomSeed;
    InitialXVGenerator xvInitialization;
    double coordInitializationSpread;
    FinalCoordinateFileFormat finalXVOutput;
    int numberAtomsOnBoxEdge[3];
    double epsilonLJ;
    double sigmaLJ;
    double interactionCutoffRadius;
    int propertyPrintingInterval;
    int numberRadialDistrPoints;
    double radialDistrCutoffRadius;
    bool trajectoryOutput;
    TrajectoryFileFormat trajectoryOutputFormat;
    int trajectoryOutputInterval;
    double gamma;
    double theta;
};

FinalCoordinateFileFormat finalCoordinateFileFormatFromInt(int ntxo);
int finalCoordinateFileFormatToInt(FinalCoordinateFileFormat ntxo);

InitialXVGenerator initialXVGeneratorFromInt(int ntxi);
int initialXVGeneratorToInt(InitialXVGenerator ntxi);

TrajectoryFileFormat trajectoryFileFormatFromInt(int ntpw);
int trajectoryFileFormatToInt(TrajectoryFileFormat ntpw);

SimulationType simulationTypeFromInt(int ntt);
int simulationTypeToInt(SimulationType ntt);

ChainSimType chainSimTypeFromInt(int ntt);
int chainSimTypeToInt(ChainSimType ntt);

#endif // MDPARAMETERS_H
