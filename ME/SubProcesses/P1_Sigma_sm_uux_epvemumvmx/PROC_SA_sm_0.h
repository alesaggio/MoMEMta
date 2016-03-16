//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_uux_epvemumvmx_H
#define MG5_Sigma_sm_uux_epvemumvmx_H

#include <complex> 
#include <vector> 
#include <utility> 
#include <map> 

#include "../../src/Parameters_sm.h"
#include "../../src/process_base_classes.h"

//==========================================================================
// A class for calculating the matrix elements for
// Process: u u~ > e+ ve mu- vm~ WEIGHTED=8 @1
// Process: c c~ > e+ ve mu- vm~ WEIGHTED=8 @1
// Process: d d~ > e+ ve mu- vm~ WEIGHTED=8 @1
// Process: s s~ > e+ ve mu- vm~ WEIGHTED=8 @1
//--------------------------------------------------------------------------

// Forward declaration needed to template correctly __MatrixElement in the class
class PROC_SA_sm_0; 

  class PROC_SA_sm_0: public CPPProcess
  {
    public:

      // Constructor & destructor
      PROC_SA_sm_0(Parameters_sm &params); 
      virtual ~PROC_SA_sm_0() {}; 

      // Calculate flavour-independent parts of cross section.
      virtual std::map < std::pair < int, int > , double > sigmaKin(const
          std::vector < std::vector<double> > &initialMomenta, const
          std::vector < std::pair < int, std::vector<double> > > &finalState);

      // Info on the subprocess.
      virtual std::string name() const {return "u u~ > e+ ve mu- vm~ (sm)";}

      const std::vector<double> & getMasses() const {return mME;} 

    private:

      // default constructor should be hidden
      PROC_SA_sm_0(); 

      // list of helicities combinations
      const int helicities[64][6] = {{-1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1,
          -1, 1}, {-1, -1, -1, -1, 1, -1}, {-1, -1, -1, -1, 1, 1}, {-1, -1, -1,
          1, -1, -1}, {-1, -1, -1, 1, -1, 1}, {-1, -1, -1, 1, 1, -1}, {-1, -1,
          -1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, -1, 1}, {-1,
          -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, 1}, {-1, -1, 1, 1, -1, -1},
          {-1, -1, 1, 1, -1, 1}, {-1, -1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1},
          {-1, 1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, 1}, {-1, 1, -1, -1, 1,
          -1}, {-1, 1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, -1}, {-1, 1, -1, 1,
          -1, 1}, {-1, 1, -1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1}, {-1, 1, 1, -1,
          -1, -1}, {-1, 1, 1, -1, -1, 1}, {-1, 1, 1, -1, 1, -1}, {-1, 1, 1, -1,
          1, 1}, {-1, 1, 1, 1, -1, -1}, {-1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1,
          -1}, {-1, 1, 1, 1, 1, 1}, {1, -1, -1, -1, -1, -1}, {1, -1, -1, -1,
          -1, 1}, {1, -1, -1, -1, 1, -1}, {1, -1, -1, -1, 1, 1}, {1, -1, -1, 1,
          -1, -1}, {1, -1, -1, 1, -1, 1}, {1, -1, -1, 1, 1, -1}, {1, -1, -1, 1,
          1, 1}, {1, -1, 1, -1, -1, -1}, {1, -1, 1, -1, -1, 1}, {1, -1, 1, -1,
          1, -1}, {1, -1, 1, -1, 1, 1}, {1, -1, 1, 1, -1, -1}, {1, -1, 1, 1,
          -1, 1}, {1, -1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1}, {1, 1, -1, -1, -1,
          -1}, {1, 1, -1, -1, -1, 1}, {1, 1, -1, -1, 1, -1}, {1, 1, -1, -1, 1,
          1}, {1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1}, {1, 1, -1, 1, 1,
          -1}, {1, 1, -1, 1, 1, 1}, {1, 1, 1, -1, -1, -1}, {1, 1, 1, -1, -1,
          1}, {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1},
          {1, 1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1}};

      // Private functions to calculate the matrix element for all subprocesses
      // Calculate wavefunctions
      void calculate_wavefunctions(const int perm[], const int hel[]); 
      std::complex<double> amp[18]; 
      double matrix_1_uux_epvemumvmx(); 
      double matrix_1_ddx_epvemumvmx(); 

      // map of final states
      std::map < std::vector<int> , std::vector < __MatrixElement <
          PROC_SA_sm_0 >> > mapFinalStates;

      // Reference to the model parameters instance passed in the constructor
      Parameters_sm& params; 

      // vector with external particle masses
      std::vector<double> mME; 

      // vector with momenta (to be changed each event)
      double * momenta[6]; 
  }; 



  #endif  // MG5_Sigma_sm_uux_epvemumvmx_H

