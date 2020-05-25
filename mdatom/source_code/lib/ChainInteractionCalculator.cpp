#include "ChainInteractionCalculator.h"

#include <algorithm>
#include <exception>
#include <cmath>
#include <iostream>
#include <vector>

inline int nearestInteger(double x) {
    return x > 0 ? static_cast<int>(x + 0.5) : static_cast<int>(x - 0.5);
}

double ChainInteractionCalculator::dist(int i, int j, const std::vector<double>& pos) {
  std::vector<double> vec(3);

  int i3 = 3 * i;
  int j3 = 3 * j;
  for (int m = 0; m < 3; m++) {
      vec[m] = pos[i3 + m] - pos[j3 + m];
      vec[m] = vec[m] - nearestInteger(vec[m]*inverseBoxLength[m]) * par.boxSize[m];
  }

  return std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
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

double* cross(const double v1[3], const double v2[3], double result[3]) {
  result[0] = (v1[1] * v2[2] - v1[2] * v2[1]);
  result[1] = (v1[2] * v2[0] - v1[0] * v2[2]);
  result[2] = (v1[0] * v2[1] - v1[1] * v2[0]);

  return result;
}


double norm(int i, const std::vector<double>& pos) {
  return std::sqrt(pos[3*i  ] * pos[3*i  ] +
                   pos[3*i+1] * pos[3*i+1] +
                   pos[3*i+2] * pos[3*i+2]);
}

std::vector<double> ChainInteractionCalculator::connect(int i, int j, const std::vector<double> positions) {
  std::vector<double> ret(3);

  int i3 = 3 * i;
  int j3 = 3 * j;
  for (int m = 0; m < 3; m++) {
      ret[m] = positions[i3 + m] - positions[j3 + m];
      ret[m] = ret[m] - nearestInteger(ret[m]*inverseBoxLength[m]) * par.boxSize[m];
  }

  return ret;
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

  // Kosinussatz
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
  : InteractionCalculator::InteractionCalculator(parameters),
    xjk(3), xkl(3),
    dfi(3), dfj(3), dfk(3) {
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

void ChainInteractionCalculator::calculate (std::vector<double>& positions, const std::vector<std::vector<bool>>& bonds, std::vector<double>& forces) {

    resetVariablesToZero(forces);
    initializeValues();

    if (par.chainMdType == ChainSimType::complete) {
        calculateA(positions, bonds, forces);
        calculateD(positions, bonds, forces);
        // calculateB(positions, bonds, forces);
    } else if (par.chainMdType == ChainSimType::noAngles) {
        calculateD(positions, bonds, forces);
        // calculateB(positions, bonds, forces);
    }

    for (int i = 0; i < par.numberAtoms - 1; i++) {
        for (int j = i + 1; j < par.numberAtoms; j++) {
            calculateInteraction(i, j, positions, bonds, forces);
        }
    }
    virial /= 2.;
}

void ChainInteractionCalculator::calculateA (const std::vector<double>& positions,
                                             const std::vector<std::vector<bool>> bonds,
                                             std::vector<double>& forces) {

    for(int i = 1; i < par.numberAtoms - 1; ++i){
        calculateInteractionA(i, positions, bonds, forces);
    }
}

void ChainInteractionCalculator::calculateD (std::vector<double>& positions,
                                             const std::vector<std::vector<bool>> bonds,
                                             std::vector<double>& forces) {
    for(int i = 1; i < par.numberAtoms - 2; ++i){
        calculateInteractionD(i-1, i, i+1, i+2, positions, bonds, forces);
    }
}

void ChainInteractionCalculator::calculateB (const std::vector<double>& positions,
                                             const std::vector<std::vector<bool>> bonds,
                                             std::vector<double>& forces) {
    for (int i = 1; i < par.numberAtoms - 1; i++) {
        calculateInteractionB(i, i+1, positions, bonds, forces);
    }
}

void ChainInteractionCalculator::calculateInteractionD(int i, int j, int k, int l,
                                                       std::vector<double>& positions,
                                                       const std::vector<std::vector<bool>>& bonds,
                                                       std::vector<double>& forces){
    applyPeriodicBoundaryConditions(i, j, k, l, positions);
    calculateSquaredDistance();
    calculateDihedral(i,j,k,l, positions, bonds);
    calculateAngle(i, j, k, positions, bonds);
    calculateSecondAngle(j, k, l, positions, bonds);

    if (k == j+1 && k < par.numberAtoms && j > 0 && bonds[j][k]) {
          double E_pot_dihedral = Vn * (4. +
                                        std::cos(3 * dihedral_ijkl - gamma) -
                                        std::cos(2 * dihedral_ijkl - gamma) +
                                        std::cos(1 * dihedral_ijkl - gamma));

          potentialEnergy += E_pot_dihedral;
          // std::cout << "dihedral contribution to energy: " << E_pot_dihedral << " dihedral angle: " << dihedral_ijkl << std::endl;
    }

    calculateForceAndVirialContributionsD(i, j, k, l, forces, positions);
}

void ChainInteractionCalculator::calculateForceAndVirialContributionsD(int i, int j, int k, int l,
                                                                       std::vector<double>& forces,
                                                                       std::vector<double>& positions){
    std::vector<double> forcea(3);
    std::vector<double> forceb(3);
    std::vector<double> forcec(3);
    std::vector<double> forced(3);
    std::vector<double> neg_xij(3);
    std::vector<double> neg_xjk(3);

    for(int m = 0; m < 3; ++m){
        neg_xij[m] = -xij[m];
        neg_xjk[m] = -xjk[m];
    }

    std::vector<double> p1_strich(3);
    std::vector<double> p2_strich(3);

    p1_strich = cross(neg_xij, xjk);
    double norm_p1_strich = std::sqrt(p1_strich[0] * p1_strich[0] +
                                      p1_strich[1] * p1_strich[1] +
                                      p1_strich[2] * p1_strich[2]);
    p2_strich = cross(xkl, neg_xjk);
    double norm_p2_strich = std::sqrt(p2_strich[0] * p2_strich[0] +
                                      p2_strich[1] * p2_strich[1] +
                                      p2_strich[2] * p2_strich[2]);

    std::vector<double> p1(3);
    std::vector<double> p2(3);
    for (int m = 0; m < 3; ++m){
        p1[m] = p1_strich[m] / norm_p1_strich;
        p2[m] = p2_strich[m] / norm_p2_strich;
    }

    double coeff_a = 1/(std::sqrt(rij2) * std::sin(angle_ijk))*Vn*(std::sin(dihedral_ijkl) +
                                                                   2*std::sin(2*dihedral_ijkl) +
                                                                   3*std::sin(3*dihedral_ijkl));

    double coeff_d = 1/(std::sqrt(rkl2) * std::sin(angle_jkl))*Vn*(std::sin(dihedral_ijkl) +
                                                                   2*std::sin(2*dihedral_ijkl) +
                                                                   3*std::sin(3*dihedral_ijkl));
    for (int m = 0; m < 3; ++m){
        forcea[m] = coeff_a * p1[m];
        forced[m] = coeff_d * p2[m];
    }

    //compute force fc
    std::vector<double> oc(3);
    for (int m = 0; m < 3; ++m){
        oc[m] = xjk[m]/2;
    }

    double norm_oc2;
    for (int m = 0; m < 3; ++m)
        norm_oc2 += oc[m] * oc[m];

    double coeff_c = 1 / norm_oc2;

    //now all the computations needed for tc
    std::vector<double> first_term(3);
    std::vector<double> second_term(3);
    std::vector<double> third_term(3);

    first_term = cross(oc, forced);
    second_term = cross(xkl, forced);
    third_term = cross(neg_xij, forcea);
    for(int m = 0; m < 3; ++m){
        second_term[m] /= 2.;
        third_term[m] /= 2.;
    }
    //now compute tc
    std::vector<double> tc(3);
    for(int m = 0; m < 3; ++m){
        tc[m] = -(first_term[m] + second_term[m] + third_term[m]);
    }

    std::vector<double> tc_x_oc(3);
    tc_x_oc = cross(tc, oc);
    //now compute fc
    for(int m = 0; m < 3; ++m){
        forcec[m] = coeff_c * tc_x_oc[m];
    }

    for(int m = 0; m < 3; ++m){
        forceb[m] = -forcea[m] - forcec[m] - forced[m];
    }

    //now add the forces to the right places in the forces vector
    int i3 = 3 * i;
    int j3 = 3 * j;
    int k3 = 3 * k;
    int l3 = 3 * l;

    for (int m = 0; m < 3; ++m){
        forces[i3 + m] -= forcea[m];
        forces[j3 + m] -= forceb[m];
        forces[k3 + m] -= forcec[m];
        forces[l3 + m] -= forced[m];
    }///fehlt nocht virial
}

void ChainInteractionCalculator::calculateInteractionA(int i, const std::vector<double>& positions,
                                                              const std::vector<std::vector<bool>>& bonds,
                                                              std::vector<double>& forces){
    if (bonds[i-1][i] && bonds[i][i+1]) {
      calculateAngle(i-1, i, i+1, positions, bonds);

      // Calculate forces for every atom part of angle in every direction (x, y, z)
      std::vector<double> v_ba = connect(i, i-1, positions);
      std::vector<double> v_bc = connect(i, i+1, positions);
      std::vector<double> v_cb = connect(i+1, i, positions);

      std::vector<double> v_baxbc = cross(v_ba, v_bc);

      std::vector<double> pa = cross(v_ba, v_baxbc);
      std::vector<double> pc = cross(v_cb, v_baxbc);

      for (int m = 0; m < 3; m++) {
        dfi[m] = -2 * ka * (angle_ijk - theta0) / (norm(0, pa) * dist(i-1, i, positions)) * pa[m];
        dfk[m] = -2 * ka * (angle_ijk - theta0) / (norm(0, pc) * dist(i, i+1, positions)) * pc[m];
        dfj[m] = -dfi[m] - dfk[m];

        std::cout << angle_ijk << " " << pa[m] << " " << pc[m] << std::endl;
      }

      calculatePotentialAndForceMagnitudeA();
      // std::cout << "angle contribution to pot energy: " << ei << " angle: " << angle_ijk << std::endl;
      potentialEnergy += ei;
      calculateForceAndVirialContributionsA(i-1, i, i+1, positions, forces);
    }
}

void ChainInteractionCalculator::calculatePotentialAndForceMagnitudeA(){
    ei = ka * std::pow((angle_ijk - theta0), 2);
}

void ChainInteractionCalculator::calculateForceAndVirialContributionsA(int i, int j, int k, const std::vector<double>& pos,
                                                                       std::vector<double>& forces) {
    int i3 = 3 * i;
    int j3 = 3 * j;
    int k3 = 3 * k;

    for (int m = 0; m < 3; m++) {
      forces[i3 + m] -= dfi[m];
      forces[j3 + m] -= dfj[m];
      forces[k3 + m] -= dfk[m];

      // dont know how to calculate virial in this case
      // virial -= 0.5 * xij[m] * df;
    }
}

void ChainInteractionCalculator::calculateInteractionB(int i, int j, const std::vector<double>& pos,
                                                       const std::vector<std::vector<bool>>& bonds, std::vector<double>& forces) {
  InteractionCalculator::applyPeriodicBoundaryConditions(i, j, pos);
  calculateSquaredDistance();

  dij = 0;

    if (type != ChainSimType::noChains) {
        bond_ij = bonds[i][j];

        // bond and contribution to potential energy
        if (bond_ij) {
            double val1 = kb * std::pow(std::sqrt(rij2) - r0, 2);
            potentialEnergy += val1;
            dij -= 2 * kb * (std::sqrt(rij2) - r0);

            std::cout << "dij from bond  " << 2 * kb * (std::sqrt(rij2) - r0) << std::endl;

            std::cout << "bond contribution to energy:" << val1 << " dist to eq: " << (std::sqrt(rij2) - r0) << std::endl;
        }

        int i3 = 3 * i;
        int j3 = 3 * j;

        for (int m = 0; m < 3; m++) {
            // Force increment in direction of inter-particle vector
            //(note: xij[m]/rij is unit vector in inter-particle direction.)
            double df = xij[m] * dij;
            forces[i3 + m] += df;
            forces[j3 + m] -= df;
            virial -= xij[m] * df;
        }
    }
}

void ChainInteractionCalculator::calculateInteraction(int i, int j, const std::vector<double>& positions,
      const std::vector<std::vector<bool>>& bonds, std::vector<double>& forces) {
    InteractionCalculator::applyPeriodicBoundaryConditions(i, j, positions);
    calculateSquaredDistance();

    dij = 0;

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
    rij2 = 0;
    rjk2 = 0;
    rkl2 = 0;
    for(int m = 0; m < 3; ++m){
        rij2 += xij[m] * xij[m];
        rjk2 += xjk[m] * xjk[m];
        rkl2 += xkl[m] * xkl[m];
    }
}
