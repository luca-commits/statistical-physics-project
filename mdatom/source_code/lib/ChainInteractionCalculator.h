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
      
      //loops over all atoms and adds angle contributions to E_pot by calling calculateInteractionA
      void calculateA (const std::vector<double>& positions, const std::vector<std::pair<int>> bonds,
                                const std::vector<double>& forces);
    private:

      void calculateAngle(int i, int j, int l, const std::vector<double>& positions,
                          const std::vector<std::pair<int, int>>& bonds);
      void calculateDihedral (int i, int j, int k, int l, const std::vector<double>& positions,
                              const std::vector<std::pair<int, int>>& bonds); 
      //loops over all atoms and adds angle contributions to E_pot by calling calculateInteractionA
      void calculateA (const std::vector<double>& positions, const std::vector<std::pair<int, int>> bonds,
                                const std::vector<double>& forces);
      //sets the angle by calling calculateAngle, then calls calculatePotential
      void calculateInteractionA (int i, int k, int j, const std::vector<double>& positions
                                 std::vector<double>& forces, const std::vector<std::pair<int, int>>&  bonds);
      //calculates only the potential contribution of the angle 
      void calculatePotentialA();
      //set the dihedral angle + does all calculateInteraction does in the base class
      void calculateInteraction (int i, int j, const std::vector<double>& positions,
                                 const std::vector<std::pair<int, int>>& bonds, const std::vector<double> forces);
      //calculates the potential contribution from Coulomb interaction, bond terms and dihedral terms
      void calculatePotentialAndForceMagnitude override();

      void resetVariablesToZero override();
      void initializeVariables override();


      double theta0;
      double gamma;

      double dihedral_ijkl;
      double angle_ijk;

      ChainSimulationType type;
};  

#endif //CHAININTERACTIONCALCULATOR_H


