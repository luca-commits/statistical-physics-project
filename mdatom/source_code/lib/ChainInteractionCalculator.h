#ifndef CHAININTERACTIONCALCULATOR_H
#define CHAININTERACTIONCALCULATOR_H

#include "InteractionCalculator.h"
#include <utility>


/*This class calculates the interactions between atoms considering covalent bonds. This 
 * includes bond terms, angle terms and torsion terms
 */

class ChainInteractionCalculator : public InteractionCalculator{

    public:
      explicit ChainInteractionCalculator(const MDParameters& parameters);  //constructor
      
      // called from the outside, calculates all interactions
      void ChainInteractionCalculator::calculate(const std::vector<double>& positions, const std::vector<std::vector<bool>>& bonds, std::vector<double>& forces)
      
      //loops over all atoms and adds angle contributions to E_pot by calling calculateInteractionA
      //the name doesn't make too much sense on its own, it's just to show that it is a parallel to 
      //InteractionCalculator.h 
      void calculateA (const std::vector<double>& positions, const std::vector<std::vector<bool>>& bonds);
    private:

      void calculateAngle(int i, int j, int k, const std::vector<double>& positions,
                          const std::vector<std::vector<bool>>& bonds);
      void calculateDihedral (int i, int j, int k, int l, const std::vector<double>& positions,
                              const std::vector<std::vector<bool>>& bonds); 
      //only works for chains
      //loops over all atoms and adds angle contributions to E_pot by calling calculateInteractionA
      //only need to loop over one variable since every atom has one angle (correct me if I'm wrong)
      void calculateA (const std::vector<double>& positions, const std::vector<std::vector<bool>> bonds);
      //sets the angle by calling calculateAngle, then calls calculatePotentialA
      void calculateInteractionA (int i, const std::vector<double>& positions,
                                  const std::vector<std::vector<bool>>&  bonds);
      //calculates only the potential contribution of the angle 
      void calculatePotentialA();
      //set the dihedral angle + does everything that calculateInteraction does (in the base class)
      void calculateInteraction (int i, int j, const std::vector<double>& positions,
                                 const std::vector<std::vector<bool>>& bonds, const std::vector<double> forces);
      //calculates the potential contribution from Coulomb interaction, bond terms and dihedral terms
      void calculatePotentialAndForceMagnitude () override;

      //This function is called by InteractionCalculatorA. Sice calculateInteractionA only adds to
      // the potential, only the latter has to be reset to zero
      void resetPotentialToZero();
      void initializeValues () override;


      double theta0;
      double gamma;

      double dihedral_ijkl;
      double angle_ijk;
      
      bool bond_ij;
      
      double kb, ka; // constants for bond contribution to potential
      unsigned int n;

      ChainSimType type;
};  

#endif //CHAININTERACTIONCALCULATOR_H


