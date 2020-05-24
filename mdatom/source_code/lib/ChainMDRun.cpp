#include "ChainMDRun.h"
#include "PeriodicBoundaryConditions.h"
#include "CenterOfMassCalculator.h"
#include "TrajectoryFileWriter.h"
#include <cmath>

ChainMDRun::ChainMDRun(const MDParameters& parameters,
  MDRunOutput& out, TrajectoryFileWriter& trajectoryFileWriter)
  : MDRun(parameters, out, trajectoryFileWriter),
    chainForceCalculator(parameters) {

};

// Modified original function from MDRun to use new ChainInteractionCalculator for chains
void ChainMDRun::performStep(std::vector<double> &positions, std::vector<double> &velocities, std::vector<std::vector<bool>> &bonds, int nstep, double time) {
  /* put atoms in central periodic box */
  PeriodicBoundaryConditions::recenterAtoms(par.numberAtoms, positions, par.boxSize);

  /* calculate forces, potential energy, virial
   * and contribution to the radial distribution function
   */
  chainForceCalculator.calculate(positions, bonds, forces);
  radialDistribution.addInstantaneousDistribution(chainForceCalculator.getInstantaneousRadialDistribution());
  double vir = chainForceCalculator.getVirial();
  properties[2] = chainForceCalculator.getPotentialEnergy();
  properties[3] = vir;

  /* determine velocity scaling factor, when coupling to a bath */
  double scal = 1;
  if (par.mdType == SimulationType::constantTemperature) {
      double dtt = par.timeStep / par.temperatureCouplingTime;
      scal = std::sqrt(1 + dtt * (ekin0 / ekg - 1));
  }

  /* perform leap-frog integration step,
   * calculate kinetic energy at time t-dt/2 and at time t,
   * and calculate pressure
   */
  double oldKineticEnergy = 0.;
  double newKineticEnergy = 0.;
  for (int j3 = 0; j3 < nat3; j3++) {
      double oldVelocity = velocities[j3];
      double newVelocity = (oldVelocity + forces[j3] * dtm) * scal;
      oldKineticEnergy += newVelocity * newVelocity;
      newKineticEnergy += (oldVelocity + newVelocity) * (oldVelocity + newVelocity);
      velocities[j3] = newVelocity;
      positions[j3] += newVelocity * par.timeStep;
  }
  oldKineticEnergy *= (par.atomicMass / 2.);
  newKineticEnergy *= (par.atomicMass / 8.);
  properties[1] = newKineticEnergy;
  properties[0] = properties[1] + properties[2];
  double pres = 2. * (newKineticEnergy - vir) / (vol * 3.);
  properties[4] = pres;
  properties[5] = scal;
  if (par.mdType == SimulationType::constantTemperature) {
      ekg = oldKineticEnergy;
  }

  /* update arrays for averages and fluctuations */
  for (int m = 0; m < numberProperties; m++) {
      averages[m] += properties[m];
      fluctuations[m] += properties[m] * properties[m];
  }

  printOutputForStep(positions, velocities, nstep, time);
}

void ChainMDRun::run(std::vector<double> &x, std::vector<double> &v, std::vector<std::vector<bool>> &bonds) {
    forces.resize(x.size());
    synchronizedPositions.resize(x.size());
    radialDistribution.setZero();

    initializeVariables();
    initializeTemperature(v);

    output.printInitialTemperature(properties[1] / fac);
    output.printIterationStart();

    /* dynamics step */
    double time = par.initialTime;
    for (int nstep = 0; nstep < par.numberMDSteps; nstep++) {
        time += par.timeStep;
        performStep(x, v, bonds, nstep, time);
    }

    printAverages(time);
}
