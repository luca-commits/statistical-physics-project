#ifndef CHAINMDRUN_H
#define CHAINMDRUN_H

#include "MDRun.h"
#include "ChainInteractionCalculator.h"

// Inherited class from MDRun which accomodates for bonds between atoms
class ChainMDRun : public MDRun {
  public:
    ChainMDRun(const MDParameters& parameters, MDRunOutput& out, TrajectoryFileWriter& trajectoryFileWriter);
  private:
    virtual void performStep(std::vector<double>& positions, std::vector<double>& velocities, int nstep, double time);
    
    ChainInteractionCalculator chainForceCalculator;
};

#endif
