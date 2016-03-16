//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <string> 
#include <utility> 
#include <vector> 
#include <map> 

#include "PROC_SA_sm_0.h"
#include "HelAmps_sm.h"
#include "process_base_classes.h"

using namespace MG5_sm; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: u u~ > e+ ve mu- vm~ WEIGHTED=8 @1
// Process: c c~ > e+ ve mu- vm~ WEIGHTED=8 @1
// Process: d d~ > e+ ve mu- vm~ WEIGHTED=8 @1
// Process: s s~ > e+ ve mu- vm~ WEIGHTED=8 @1

//--------------------------------------------------------------------------
// Initialize process.

PROC_SA_sm_0::PROC_SA_sm_0(Parameters_sm &params):
params(params)
{
  // Set external particle masses for this matrix element
  mME.push_back(params.ZERO); 
  mME.push_back(params.ZERO); 
  mME.push_back(params.ZERO); 
  mME.push_back(params.ZERO); 
  mME.push_back(params.ZERO); 
  mME.push_back(params.ZERO); 

  mapFinalStates[{-11, 12, 13, -14}] = 
  {
    {
      &PROC_SA_sm_0::matrix_1_uux_epvemumvmx, 
      true, 
      {
        std::make_pair(2, -2), std::make_pair(4, -4)
      }
      , 
      64, 
      36
    }
    , 
    {
      &PROC_SA_sm_0::matrix_1_ddx_epvemumvmx, 
      true, 
      {
        std::make_pair(1, -1), std::make_pair(3, -3)
      }
      , 
      64, 
      36
    }
  }; 

}

//--------------------------------------------------------------------------
// Evaluate |M|^2, return a map of final states

std::map < std::pair < int, int > , double > PROC_SA_sm_0::sigmaKin(const
    std::vector < std::vector<double> > &initialMomenta, const std::vector <
    std::pair < int, std::vector<double> > > &finalState)
{

  // Set initial particle momenta
  momenta[0] = (double * ) (&initialMomenta[0][0]); 
  momenta[1] = (double * ) (&initialMomenta[1][0]); 

  // Suppose final particles are passed in the "correct" order
  std::vector<int> selectedFinalState(6 - 2); 
  size_t index = 2; 
  for (auto const &finalPart: finalState)
  {
    selectedFinalState[index - 2] = finalPart.first; 
    momenta[index] = (double * ) (&finalPart.second[0]); 
    index++; 
  }

  // Set the event specific parameters
  params.setDependentParameters(); 
  params.setDependentCouplings(); 

  // Initialise result object
  std::map < std::pair < int, int > , double > result; 

  // Define permutation
  int perm[6]; 
  for(int i = 0; i < 6; i++ )
  {
    perm[i] = i; 
  }

  for(auto &me: mapFinalStates[selectedFinalState])
  {
    double me_sum = 0; 
    double me_mirror_sum = 0; 
    for(int ihel = 0; ihel < 64; ihel++ )
    {
      if(me.goodHel[ihel])
      {
        double sum = 0.; 
        calculate_wavefunctions(perm, helicities[ihel]); 
        double meTemp = (this->* (me.meCall))(); 
        sum += meTemp; 
        me_sum += meTemp/me.denominator; 

        if(me.hasMirrorProcess)
        {
          perm[0] = 1; 
          perm[1] = 0; 
          // Calculate wavefunctions
          calculate_wavefunctions(perm, helicities[ihel]); 
          // Mirror back
          perm[0] = 0; 
          perm[1] = 1; 
          meTemp = (this->* (me.meCall))(); 
          sum += meTemp; 
          me_mirror_sum += meTemp/me.denominator; 
        }

        if( !sum)
          me.goodHel[ihel] = false; 
      }
    }

    for (auto const &initialState: me.initialStates)
    {
      result[initialState] = me_sum; 
      if (me.hasMirrorProcess)
        result[std::make_pair(initialState.second, initialState.first)] =
            me_mirror_sum;
    }
  }


  return result; 
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void PROC_SA_sm_0::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  // Calculate all wavefunctions
  static std::complex<double> w[26][18]; 

  ixxxxx(&momenta[perm[0]][0], mME[0], hel[0], +1, w[0]); 
  oxxxxx(&momenta[perm[1]][0], mME[1], hel[1], -1, w[1]); 
  ixxxxx(&momenta[perm[2]][0], mME[2], hel[2], -1, w[2]); 
  oxxxxx(&momenta[perm[3]][0], mME[3], hel[3], +1, w[3]); 
  oxxxxx(&momenta[perm[4]][0], mME[4], hel[4], +1, w[4]); 
  ixxxxx(&momenta[perm[5]][0], mME[5], hel[5], -1, w[5]); 
  FFV1P0_3(w[0], w[1], params.GC_2, params.ZERO, params.ZERO, w[6]); 
  FFV2_3(w[2], w[3], params.GC_100, params.mdl_MW, params.mdl_WW, w[7]); 
  FFV1_1(w[4], w[6], params.GC_3, params.ZERO, params.ZERO, w[8]); 
  FFV2_5_3(w[0], w[1], params.GC_51, params.GC_58, params.mdl_MZ,
      params.mdl_WZ, w[9]);
  FFV2_4_1(w[4], w[9], params.GC_50, params.GC_59, params.ZERO, params.ZERO,
      w[10]);
  FFV2_2(w[5], w[9], params.GC_62, params.ZERO, params.ZERO, w[11]); 
  FFV2_3(w[5], w[4], params.GC_100, params.mdl_MW, params.mdl_WW, w[12]); 
  FFV1_2(w[2], w[6], params.GC_3, params.ZERO, params.ZERO, w[13]); 
  FFV2_4_2(w[2], w[9], params.GC_50, params.GC_59, params.ZERO, params.ZERO,
      w[14]);
  FFV2_1(w[3], w[9], params.GC_62, params.ZERO, params.ZERO, w[15]); 
  FFV2_2(w[0], w[7], params.GC_100, params.ZERO, params.ZERO, w[16]); 
  FFV1P0_3(w[0], w[1], params.GC_1, params.ZERO, params.ZERO, w[17]); 
  FFV1_1(w[4], w[17], params.GC_3, params.ZERO, params.ZERO, w[18]); 
  FFV2_3_3(w[0], w[1], params.GC_50, params.GC_58, params.mdl_MZ,
      params.mdl_WZ, w[19]);
  FFV2_4_1(w[4], w[19], params.GC_50, params.GC_59, params.ZERO, params.ZERO,
      w[20]);
  FFV2_2(w[5], w[19], params.GC_62, params.ZERO, params.ZERO, w[21]); 
  FFV1_2(w[2], w[17], params.GC_3, params.ZERO, params.ZERO, w[22]); 
  FFV2_4_2(w[2], w[19], params.GC_50, params.GC_59, params.ZERO, params.ZERO,
      w[23]);
  FFV2_1(w[3], w[19], params.GC_62, params.ZERO, params.ZERO, w[24]); 
  FFV2_2(w[0], w[12], params.GC_100, params.ZERO, params.ZERO, w[25]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV2_0(w[5], w[8], w[7], params.GC_100, amp[0]); 
  FFV2_0(w[5], w[10], w[7], params.GC_100, amp[1]); 
  FFV2_0(w[11], w[4], w[7], params.GC_100, amp[2]); 
  VVV1_0(w[6], w[12], w[7], params.GC_4, amp[3]); 
  VVV1_0(w[12], w[7], w[9], params.GC_53, amp[4]); 
  FFV2_0(w[13], w[3], w[12], params.GC_100, amp[5]); 
  FFV2_0(w[14], w[3], w[12], params.GC_100, amp[6]); 
  FFV2_0(w[2], w[15], w[12], params.GC_100, amp[7]); 
  FFV2_0(w[16], w[1], w[12], params.GC_100, amp[8]); 
  FFV2_0(w[5], w[18], w[7], params.GC_100, amp[9]); 
  FFV2_0(w[5], w[20], w[7], params.GC_100, amp[10]); 
  FFV2_0(w[21], w[4], w[7], params.GC_100, amp[11]); 
  VVV1_0(w[17], w[12], w[7], params.GC_4, amp[12]); 
  VVV1_0(w[12], w[7], w[19], params.GC_53, amp[13]); 
  FFV2_0(w[22], w[3], w[12], params.GC_100, amp[14]); 
  FFV2_0(w[23], w[3], w[12], params.GC_100, amp[15]); 
  FFV2_0(w[2], w[24], w[12], params.GC_100, amp[16]); 
  FFV2_0(w[25], w[1], w[7], params.GC_100, amp[17]); 

}
double PROC_SA_sm_0::matrix_1_uux_epvemumvmx() 
{

  static std::complex<double> ztemp; 
  static std::complex<double> jamp[1]; 
  // The color matrix
  static const double denom[1] = {1}; 
  static const double cf[1][1] = {{3}}; 

  // Calculate color flows
  jamp[0] = -amp[0] - amp[1] - amp[2] - amp[3] - amp[4] - amp[5] - amp[6] -
      amp[7] - amp[8];

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(int i = 0; i < 1; i++ )
  {
    ztemp = 0.; 
    for(int j = 0; j < 1; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  return matrix; 
}

double PROC_SA_sm_0::matrix_1_ddx_epvemumvmx() 
{

  static std::complex<double> ztemp; 
  static std::complex<double> jamp[1]; 
  // The color matrix
  static const double denom[1] = {1}; 
  static const double cf[1][1] = {{3}}; 

  // Calculate color flows
  jamp[0] = -amp[9] - amp[10] - amp[11] - amp[12] - amp[13] - amp[14] - amp[15]
      - amp[16] - amp[17];

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(int i = 0; i < 1; i++ )
  {
    ztemp = 0.; 
    for(int j = 0; j < 1; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  return matrix; 
}



