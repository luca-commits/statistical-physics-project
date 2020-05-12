#ifndef CHAININTERACTIONCALCULATOR_H
#define CHAININTERACTIONCALCULATOR_H

#include "IteractionCalculator.h"
#include <utility>


/*This class calculates the interactions between atoms considering covalent bonds. This 
 * includes bond terms, angle terms and torsion terms
 */

class ChainInteractionCalculator : public InteractionCalculator{
    public:
      explicit ChainInteractionCalculator(const MDParameters& parameters);  //constructor
      
      void calculateInteraction (int i, int j, const std::vector<double>& positions
                                 std::vector<double>& forces, const std::vector<std::pair<unsigned int, unsigned int>&  bonds);
      
      void calculateAngle(int i, int j, int l, const std::vector<double>& positions,
                          const std::vector<std::pair<unsigned int, unsigned int>& bonds);

      void calculateDihedral (int i, int j, int k, int l, const std::vector<double>& positions,
                          const std::vector<std::pair<unsigned int, unsigned int>& bonds); 
    
    private:
      double angle_ijk;
      double dihedral_ijkl;







#endif //CHAININTERACTIONCALCULATOR_H


