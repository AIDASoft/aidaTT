# v00-07

* 2017-06-23 Frank Gaede ([PR#15](https://github.com/AIDASoft/aidaTT/pull/15))
  -  adapt examples to recent namespace change in DD4hep
  - replace DDSurfaces with dd4hep::rec
  - fix -Winfinite-recursion in Vector5 (operator+)

* 2017-06-20 Andre Sailer ([PR#14](https://github.com/AIDASoft/aidaTT/pull/14))
  - Adapt to changes in namespaces and LCDD -->  Detector

* 2017-06-26 Frank Gaede ([PR#16](https://github.com/AIDASoft/aidaTT/pull/16))
  - fix some Coverity issues

* 2017-05-08 Andre Sailer ([PR#13](https://github.com/AIDASoft/aidaTT/pull/13))
  - CMake: Instead of using GENERATE_PACKAGE_CONFIG_FILES call configure_file and install directly. Macro was used from DD4hep before (and is no longer available there).

# v00-06

* 2017-04-03 Emilia Leogrande ([PR#11](https://github.com/AIDASoft/aidaTT/pull/11))
  -Replace ILDCellID0 with LCTrackerCellID

* 2017-04-12 Shaojun Lu ([PR#12](https://github.com/AIDASoft/aidaTT/pull/12))
  - Added "-DGBL_EIGEN_SUPPORT_ROOT" to allow  GBLInterface.cc still use ROOT from GBL.

# v00-05 
  
F. Gaede
- replaced GSL with Eigen for Vector5 and fiveByfiveMatrix
  - needs Eigen Headers to be installed
  - removed simpleFitter for now (depends on GSL and is used nowehre)
  
#  v00-04

F. Gaede
   - added set/getMass() methods to trajectory
   - removed default pion mass in
   - computeQMS()
   - computeEnergyLoss() mass has to be explicetly provided as argument
   - made compatible with c++11
   - removed -ansi -pedantic -Wno-long-long
   - fixed narrowing in initializer lists
  
#  v00-03
  
 
 - many changes, fixes and updates (for details, see svn -v log )
 

F.Gaede

 - added material_ntuples.cpp to create ntuple with material
   properties from model and surfaces
 - drawx0.C makes comparison plots for x0 
 
 - fixed scattering matrices 
   (setting off diagonal elements to 0.)

 - implmented cone intersection 
   using a newtonian method

  - added optional use of streamlog
  -  adapted lcio\_tracks.cpp to (re)fit
    tracks including material effects
   - fixed some compiler warnings
  - fixed geometry initialization for examples
   ( added  IGeometry& instance(const std::string& initString ) )
  - added example material\_effects/check_materials.cpp

  - cleaned up trajectory class ( documentation, code,...)
    - fixed issue w/ energy loss for jacobian
    - removed unneeded members and methods
      - surface intersections, bfield,...

  - fixed issue w/ calculation of curvilinear coordinate system (use s==0.)
  - pulls for central muons are now as good as w/ DDKalTest

 - finished first implementation of taking energy loss into account for
   initial track parameters:
     - trajectorElements now hold a track state that is valid until the next element
     - track states are computed in addMeasurment/addScatterer
   - implemented pointAt/tangentAt for trajectory with many track parameters
 - moved energy loss calculation to materialutils.cc
 - updated examples to using IGeometry instead of DD4hepGeometry
 - some code cleanup ...

  - cleaned up geometry interface:
    - made IGeometry a singleton
    - implemented access to B field map
    - added doxygen documentation
  - added materialUtils.hh
    - computeQMS ( needs debugging )
  - added a track state to trajectoryElement

  - restrucured code base:
    - renamed helixGymnastics.hh to helixUtil.hh
    - removed helixHelpers.hh
      -> combine all helix relates utilitiy functions in helixUtil.hh
      -> combine fitting releated helper functions in utilities.hh   
    - documented helix functions
    - created second version taking only Vector5 and Vector3D as
      arguments to have less overhead when used outside aidaTT
    - cleaned up definition of track parameterization 
       -> access functions defined in trackParameterizationLCIO.hh
    - moved Vector3D, Vector5 and fiveByfiveMatrix to core
    - renamed helpers to util
   - added implementations of Vector5 and fiveByFiveMatrix
      using the Eigen library for later use  ...

   - bug fix for intersection calculation w/ cylinders:
      - return solution w/ positive s if both solutions
        are ~equal
   - bug fix for calculation of phi from x,y:
     - take abs(dphi) to be less than M_PI

   - added code to treat 1D measurements (commented out for now)
   - write out Millipede-2 file for debugging in GBLInterface

   - restrict intersections to first half arc of the track

   - GBLInterface
   - calculate only intersections with positive s and sort list
   - fix calculation of qms

  - changed IFittingAlgorithm, GBLInterface and
    trajectory to be able to return more than one result
    (at different sites) for the fit

Y.Voutsinas

  - calculating energy loss only in surfaces that are considered scatterers
  - correction in  the calculation of particle's energy
  - first implementation of energy loss - still needs some debugging & testing
  - modify trajectoryElement as to obtain an element that isn't a scatterer or a measurement, 
   use IP as such an element and add it to the trajectory

  - re-implement seeding and iterative fitting in the example, removing some not necessary debugging output

  - modifying function intersectsWithinZCylinderBounds in order to return both solutions and not only the closest one

  - correcting the sign in the calculation of qop

  -  fixing a bug in the calculation of the first's element arc-length & 
   modify the lcio_tracks example in order to cope with TPC tracks

  - crosschecking the system's transformations

  - unitTests/finalTrackTest.cc
    adding test to check the creation of an lcio track from an aida track
  - adding test of covariance matrix in initial track parameters test
  - test for comparing the trajectory cov. mat with the one that is written out to lcio track
  
#  v00-02
  
F.Gaede:
 - fixed major memory leak in Vector5 and fiveByFiveMatrix 
   ( unneeded alloc in assignment operator)

 - made c'tors and accessors in Vector5 and fiveByFiveMatrix 
   more efficient by removing unneeded range checks and
   initializations

Y.Voutsinas:
 -  adding test for functions calculateCurvature & calculateLambda
 -  correcting transformation from perigee to L3

  
#  v00-01 

C. Rosemann:
- first release:
   - implementation of GBL fitting using a DD4hep model that provides DDRec::Surfaces
   - implementation of multiple scattering still under validation.

