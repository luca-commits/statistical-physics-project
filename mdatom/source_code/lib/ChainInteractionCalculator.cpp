#include "ChainInteractionCalculator.h"
#include <algorithm>
#include <exception>
#include <cmath>
#include <iostream>

double dist(int i, int j, const std::vector<double>& pos) {
  return std::sqrt((pos[3*i  ]-pos[3*j  ]) * (pos[3*i  ]-pos[3*j  ]) +
                   (pos[3*i+1]-pos[3*j+1]) * (pos[3*i+1]-pos[3*j+1]) + 
                   (pos[3*i+2]-pos[3*j+2]) * (pos[3*i+2]-pos[3*j+2]));
}

double dot(int i, int j, const std::vector<double>& pos) {
  return pos[3*i  ] * pos[3*j  ] +
         pos[3*i+1] * pos[3*j+1] +
         pos[3*i+2] * pos[3*j+2];
}

double dot(const std::vector<double>& vec1, const std::vector<double>& vec2) {
  return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

std::vector<double> cross(int i, int j, const std::vector<double>& pos) {
  std::vector<double> ret;
  ret.push_back(pos[3*i+1] * pos[3*j+2] - pos[3*i+2] * pos[3*j+1]);
  ret.push_back(pos[3*i+2] * pos[3*j  ] - pos[3*i  ] * pos[3*j+2]);
  ret.push_back(pos[3*i  ] * pos[3*j+1] - pos[3*i+1] * pos[3*j  ]);
  
  return ret;
}

std::vector<double> cross(const std::vector<double>& v1, const std::vector<double>& v2) {
  std::vector<double> ret;
  ret.push_back(v1[1] * v2[2] - v1[2] * v2[1]);
  ret.push_back(v1[2] * v2[0] - v1[0] * v2[2]);
  ret.push_back(v1[0] * v2[1] - v1[1] * v2[0]);
  
  return ret;
}

double norm(int i, const std::vector<double>& pos) {
  return std::sqrt(pos[3*i  ] * pos[3*i  ] +
                   pos[3*i+1] * pos[3*i+1] +
                   pos[3*i+2] * pos[3*i+2]);
}

void ChainInteractionCalculator::calculateAngle(int i, int j, int k, const std::vector<double>& positions,
                          const std::vector<std::vector<bool>>& bonds) {
    
  if (i == j || i == k || j == k) {
    angle_ijk = 0;
    return;
  }
  
  if (!bonds[i][j] || !bonds[j][k]) {
    angle_ijk = 0;
    return;
  }
  
  double r_ij = dist(i, j, positions);
  double r_ik = dist(i, k, positions);
  double r_jk = dist(j, k, positions);
  
  angle_ijk = std::acos((r_ij * r_ij + r_jk * r_jk - r_ik * r_ik) / (2 * r_ij * r_jk));
}

ChainInteractionCalculator::ChainInteractionCalculator(const MDParameters& parameters) 
  : InteractionCalculator::InteractionCalculator(parameters) {
  type = parameters.chainMdType;
}

void ChainInteractionCalculator::calculateDihedral (int i, int j, int k, int l, const std::vector<double>& pos,
  const std::vector<std::vector<bool>>& bonds) {
  
  // TODO: Implement check whether bonds between four atoms exist.
  // Luca: if you call calculateDihedral between four atoms then there is a bond between j and k,
  //       this check has to be done in the function InteractionCalculator
  
  // Calculations taken from 
  // https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
  std::vector<double> b1;
  b1.push_back(pos[3*j  ] - pos[3*i  ]);
  b1.push_back(pos[3*j+1] - pos[3*i+1]);
  b1.push_back(pos[3*j+2] - pos[3*i+2]);
  
  std::vector<double> b2;
  b2.push_back(pos[3*k  ] - pos[3*j  ]);  
  b2.push_back(pos[3*k+1] - pos[3*j+1]);
  b2.push_back(pos[3*k+2] - pos[3*j+2]);
  
  std::vector<double> b3;
  b3.push_back(pos[3*l  ] - pos[3*k  ]);
  b3.push_back(pos[3*l+1] - pos[3*k+1]);
  b3.push_back(pos[3*l+2] - pos[3*k+2]);
  
  std::vector<double> n1 = cross(b1, b2);
  std::vector<double> n2 = cross(b2, b3);
  
  std::vector<double> m1 = cross(n1, b2);
  
  double x = dot(n1, n2) / (norm(0, n1) * norm(0, n2));
  double y = dot(m1, n2) / (norm(0, n1) * norm(0, b2) * norm(0, n2));
  dihedral_ijkl = std::atan2(y, x);
}

void ChainInteractionCalculator::calculate(const std::vector<double>& positions, const std::vector<std::vector<bool>>& bonds, std::vector<double>& forces) {
    resetVariablesToZero(forces);
    initializeValues();
    
    if(par.chainMdType == ChainSimType::complete){
        // calculateA(positions, bonds);
    }        

    for (int i = 0; i < par.numberAtoms - 1; i++) {
        for (int j = i + 1; j < par.numberAtoms; j++) {
            calculateInteraction(i, j, positions, bonds, forces);
        }
    }
    virial /= 2.;
}

void ChainInteractionCalculator::calculateA (const std::vector<double>& positions, 
                                             const std::vector<std::vector<bool>> bonds){

    for(int i = 1; i < par.numberAtoms - 1; ++i){
        calculateInteractionA(i, positions, bonds);
    }
}

void ChainInteractionCalculator::calculateInteractionA(int i, const std::vector<double>& positions, 
                                                              const std::vector<std::vector<bool>>& bonds){
    calculateAngle(i-1, i, i + 1, positions, bonds);
    calculatePotentialA();
    // std::cout << "angle contribution to energy: " << ei << " " << angle_ijk << std::endl;
    potentialEnergy += ei;
}

void ChainInteractionCalculator::calculatePotentialA(){
    ei = ka * std::pow((angle_ijk - theta0), 2);
}


void ChainInteractionCalculator::calculateInteraction(int i, int j, const std::vector<double>& positions,
      const std::vector<std::vector<bool>>& bonds, std::vector<double>& forces) {
    applyPeriodicBoundaryConditions(i, j, positions);
    calculateSquaredDistance();
    
    dij = 0;
    
    if (type != ChainSimType::noChains) {
      bond_ij = bonds[i][j];

      
      // bond and dihedral contribution to potential energy
      if (bond_ij) {
          double val1 = kb * std::pow(std::sqrt(rij2) - r0, 2);
          potentialEnergy += val1;
          dij -= 2 * kb * (std::sqrt(rij2) - r0);

          std::cout << "dij from bond  " << 2 * kb * (std::sqrt(rij2) - r0) << std::endl;
          
          std::cout << "bond contribution to energy:" << val1 << " dist to eq: " << (std::sqrt(rij2) - r0) << std::endl;
      }
      
      if (j == i+1 && j < par.numberAtoms && i > 0)
          calculateDihedral(i-1, i, j, j+1, positions, bonds);
          double val2 = Vn * 3. * (1. + std::cos(3 * dihedral_ijkl - gamma));
          // potentialEnergy += val2;
          /// std::cout << "dihedral contribution to energy: " << val2 << " dihedral angle: " << dihedral_ijkl << std::endl;
    }
    
    if (rij2 < rcutf2) {
        calculatePotentialAndForceMagnitude();
        potentialEnergy += eij;
        
    }
    
    calculateForceAndVirialContributions(i, j, forces);
    
    radialDistribution.addPairAtSquaredDistance(rij2);
}

void ChainInteractionCalculator::calculatePotentialAndForceMagnitude() {
    // E_vdW: Lennard-Jones Potential
    double riji2 = 1.0 / rij2; // inverse inter-particle distance squared
    double riji6 = riji2 * riji2 * riji2; // inverse inter-particle distance (6th power)
    double crh = c12 * riji6;
    double crhh = crh - 2 * c6; //  L-J potential work variable
    eij= crhh * riji6;
    // derived with WolframAlpha
    // dij += 12 * par.epsilonLJ * sig6 * (std::pow(rij2, 3) - sig6) / std::pow(rij2, 7);
    dij += 6. * (crh + crhh) * riji6 * riji2;

    std::cout << "dij from LJ  " << 6. * (crh + crhh) * riji6 * riji2 << std::endl;


      // std::cout << "LJ-contribution to energy: " << eij << std::endl;
}

void ChainInteractionCalculator::initializeValues() {
    sig6 = par.sigmaLJ * par.sigmaLJ;
    sig6 = sig6 * sig6 * sig6;
    c6 = par.epsilonLJ * sig6;
    c12 = c6 * sig6;
    rcutf2 = par.interactionCutoffRadius * par.interactionCutoffRadius;
    for (int m = 0; m < 3; m++)
        inverseBoxLength[m] = 1.0 / par.boxSize[m];

    Vn = 5.86; // kJ / mol
    gamma = 0; // rad
    ka = 167.36; // kJ / (mol * rad^2)
    kb = 1294.04e2; // kJ / (mol * nm^2)
    theta0 = 1.911; // rad
    r0 = 0.1526; // nm
}


