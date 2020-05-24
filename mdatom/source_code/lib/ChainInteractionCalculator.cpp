
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

double cross(const double v1[3], const double v2[3]) {
  double ret[3];
  ret[0] = (v1[1] * v2[2] - v1[2] * v2[1]);
  ret[1] = (v1[2] * v2[0] - v1[0] * v2[2]);
  ret[2] = (v1[0] * v2[1] - v1[1] * v2[0]);
  
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


void ChainInteractionCalculator::calculateSecondAngle(int i, int j, int k, const std::vector<double>& positions,
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
  
  angle_jkl = std::acos((r_ij * r_ij + r_jk * r_jk - r_ik * r_ik) / (2 * r_ij * r_jk));
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
         calculateA(positions, bonds);
    }   


    calculateD(const std::vector<double>& positions, const std::vector<bool>& bonds);    

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

void ChainInteractionCalculator::calculateD (const std::vector<double>& positions,
                                            const std::vector<std::vector<bool>> bonds){
    for(int i = 1; i < par.numberAtoms -2; ++i){
        calculateInteractionD(i-1, i, i+1, i+2, positions, bonds);
    }
}

void ChainInteractionCalculator::calculateInteractionD(int i, int j, int k, int l, 
                                                       const std::vector<double>& positions,
                                                       const std::vector<std::vector<bool>>& bonds,
                                                       const std::vector<double>& forces){
    applyPeriodicBoundaryConditions(i, j, k, l, positions);
    calculateSquaredDistance();
    calculateDihedral();
    calculateAngle(i, j, k, positions, bonds);
    calculateSecondAngle(j, k, l, positions, bonds);
    calculateForceAndVirialContributionsD(i, j, k, l, forces, positions);
}

void ChainInteractionCalculator::calculateForceAndVirialContributionsD(int i, int j, int k, int l,
                                                                       std::vector<double>& forces
                                                                       std::vector<double>& positions){
    double forcea[3];
    double forceb[3];                                                                               
    double forcec[3];
    double forced[3];
    double neg_xij[3];
    double neg_xjk[3];
    for(int i = 0; i < 3; ++i){
        neg_xij[i] = -xij[i];
        nex_xjk[i] = -xjk[i];
    }
    double p1_strich[3] = cross(neg_xij, xjk);
    double norm_p1_strich = std::sqrt(p1_strich[0] * p1_strich[0] + 
                                      p1_strich[1] * p1_strich[1] + p1_strich[2] * p1_strich[2]);
    double p2_strich[3] = cross(neg_xjk, xkl);
    double norm_p2strich = std::sqrt(p2_strich[0] * p2_strich[0] + 
                                      p2_strich[1] * p2_strich[1] + p2_strich[2] * p2_strich[2]);
    double p1[3];
    double p2[3];
    for(int i = 0; i < 3; ++i){
        p1[i] = p1_strich[i] / norm_p1_strich[i];
        p2[i] = p2_strich[i] / norm_p2strich[i];
    }
    double coeff_a = 0.5 * 1/(std::sqrt(xij2) * std::sin(angle_ijk)*Vn)*(std::sin(dihedral_ijkl) +
                                                                         2*std::sin(2*dihedral_ijkl) +
                                                                         3*std::sin(3*dihedral_ijkl));
    double coeff_b = 0.5 * 1/(std::sqrt(xkl2) * std::sin(angle_jkl)*Vn)*(std::sin(dihedral_ijkl) +
                                                                         2*std::sin(2*dihedral_ijkl) +
                                                                         3*std::sin(3*dihedral_ijkl));
    for(int i = 0; i < 3; ++i){
        forcea[i] = coeff_a * p1[i];
        forceb[i] = coeff_b * p2[i];
    }
    
    //compute force fc
    double oc[3];
    for(int i = 0; i < 3; ++i){
        oc[i] = xjk[i]/2;
    }
    
    double norm_oc2;
    for(int i; i < 3; ++i)
        norm_oc2 += oc[i] * oc[i];

    double coeff_c = 1 / norm_oc2;
    
    //now all the computations needed for tc
    double first_term[3];
    double second_term[3];
    double third_term[3];
                                                                                                           //////////////////////// 
    first_term = cross(oc, forced);
    second_term = cross(xkl, forced);
    third_term = cross(neg_xjk, forcea);
    for(int i = 0; i < 3; ++i){
        second_term[i] /= 2.;
        third_term[i] /= 2.;
    }
    //now compute tc
    for(int i = 0; i < 3; ++i){
        tc[i] = -(first_term[i] + second_term[i] + third_term[i]);
    }
    
    double tc_x_oc[3] = cross(tc, oc)
    //now compute fc
    for(int i = 0; i < 3; ++i){
        tc[i] = coeff_c * tc_x_oc;
    }
    
    for(int i = 0; i < 3; ++i){
        fc[i] = -fa[i] -fb[i] - fc[i];
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


void ChainInteractionCalculator::applyPeriodicBoundaryConditions(int i, int j, int k, int l,
                                                                 std::vector<double>& positions){
    int i3 = 3 * i;
    int j3 = 3 * j;
    int k3 = 3 * k;
    int l3 = 3 * l;
    for (int m = 0; m < 3; ++m){
          xij[m] = positions[i3 + m] - positions[j3 +m];
          xjk[m] = positions[j3 + m] - positions[k3 + m];
          xkl[m] = positions[k3 + m] - positions[l3 + m];
    }
}

void ChainInteractionCalculator::calculateSquaredDistance(){
    rij = 0;
    rjk = 0;
    rkl = 0;
    for(int m = 0; m < 3; ++m){
        rij2 += xij[m] * xij[m];
        rjk2 += xjk[m] * xjk[m];
        rkl2 += xkl[m] * xkl[m];
    }
}


