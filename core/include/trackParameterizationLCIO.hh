#ifndef trackParameterizationLCIO_hh
#define trackParameterizationLCIO_hh

/**
 * Define access functions for LCIO track parameters (trackParameterizationLCIO.hh), as described in 
 * Kraemer, Track parameters in LCIO [http://www-flc.desy.de/lcnotes/notes/LC-DET-2006-004.pdf]
 */

namespace aidaTT {
  

  /** Enum for accessing the track parameters in L3/LCIO perigee convention,
   *  defining the order of the parameters in vectors and martrices 
   */
  enum {
    OMEGA = 0,   
    TANL,
    PHI0,
    D0,
    Z0
  } ;

  /// omega, signed curvature
  inline double calculateOmega(const trackParameters& tp) {
    return tp(OMEGA);
  }

  /// tan(Lambda)
  inline double calculateTanLambda(const trackParameters& tp){
    return tp(TANL);
  }

  /// phi0, azimuth angle of track at PCA in xy-plane
  inline double calculatePhi0(const trackParameters& tp){
    return tp(PHI0);
  }

  /// d0, distance of closest approach in xy-plane
  inline double calculateD0(const trackParameters& tp){
    return tp(D0);
  }

  /// z0, z coordinate of track at PCA in xy-plane
  inline double calculateZ0(const trackParameters& tp){
    return tp(Z0);
  }

  /// signed curvature
  inline double calculateCurvature(const trackParameters& tp){
    return tp(OMEGA);
  }

  /// dip angle in sz-plane
  inline double calculateLambda(const trackParameters& tp){
    return atan(tp(TANL));
  }

  //========= define same functions for accessing Vector5 ==============

  /// omega, signed curvature
  inline double calculateOmega(const Vector5& tp) {
    return tp(OMEGA);
  }

  /// tan(Lambda)
  inline double calculateTanLambda(const Vector5& tp){
    return tp(TANL);
  }

  /// phi0, azimuth angle of track at PCA in xy-plane
  inline double calculatePhi0(const Vector5& tp){
    return tp(PHI0);
  }

  /// d0, distance of closest approach in xy-plane
  inline double calculateD0(const Vector5& tp){
    return tp(D0);
  }

  /// z0, z coordinate of track at PCA in xy-plane
  inline double calculateZ0(const Vector5& tp){
    return tp(Z0);
  }

  /// signed curvature
  inline double calculateCurvature(const Vector5& tp){
    return tp(OMEGA);
  }

  /// dip angle in sz-plane
  inline double calculateLambda(const Vector5& tp){
    return atan(tp(TANL));
  }


}


#endif

