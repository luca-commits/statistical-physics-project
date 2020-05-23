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
      void calculate(const std::vector<double>& positions, const std::vector<std::vector<bool>>& bonds, std::vector<double>& forces);
    private:

      void calculateAngle(int i, int j, int k, const std::vector<double>& positions,
                          const std::vector<std::vector<bool>>& bonds);
      void calculateDihedral (int i, int j, int k, int l, const std::vector<double>& positions,
                              const std::vector<std::vector<bool>>& bonds); 
      //only works for chains
      //loops over all atoms and adds angle contributions to E_pot by calling calculateInteractionA
      //only need to loop over one variable since every atom has one angle (correct me if I'm wrong)
      void calculateA (const std::vector<double>& positions, const std::vector<std::vector<bool>> bonds);
      //loops over all bonds and adds bond contributions to E_pot
      void calculateB (const std::vector<double>& positions, const std::vector<std::vector<bool>> bonds);
      //sets the angle by calling calculateAngle, then calls calculatePotentialA
      void calculateInteractionA (int i, const std::vector<double>& positions,
                                  const std::vector<std::vector<bool>>&  bonds);
      //calculates only the potential contribution of the angle 
      void calculatePotentialAndForceMagnitudeA();
      //set the dihedral angle + does everything that calculateInteraction does (in the base class)
      void calculateInteraction (int i, int j, const std::vector<double>& positions,
                                 const std::vector<std::vector<bool>>& bonds, std::vector<double>& forces);
      //calculates the potential contribution from Coulomb interaction, bond terms and dihedral terms
      void calculatePotentialAndForceMagnitude () override;
      
      void calculateForceAndVirialContributions(int i, int j, std::vector<double>& forces) () override;

      void calculateForceAndVirialContributionsA(int i, std::vector<double>& forces);

      void initializeValues () override;

      double Vn;
      double r0;
      double theta0;
      double gamma;

      double dihedral_ijkl;
      double angle_ijk;
      
      bool bond_ij;
      
      double kb, ka; // constants for bond contribution to potential
      unsigned int n;

      double ei; //potential due to atom i
      double dijb;//force due to bonds on the atoms i and j of a bond
      double did;//force due to dihedrals on the first atom of the dihedral quartett
      double djd;//force due to dihedrals on the second atom of the dihedral quartett
      double dkd;//force due to dihedrals on the third atom of the dihedral quartett
      double dld;//force due to dihedrals on the fourth atom of the dihedral quartett
      double dia;//force on the first element of the angle ijk 
      double dja;//force on the third element of hte angle ijk (all these forces are devided by the inter-particle vector)

      vector<int> pa, pb; //these are the orthonormal vectors I need for the calculation of force contributions

      double length_one_dihedral; //the length of the first bond in a dihedral quartett
      double length_two_dihedral; 
      double length_three_dihedral;
      
      double length_one_angle;
      double length_two_angle;
      ChainSimType type;
};  

#endif //CHAININTERACTIONCALCULATOR_H


