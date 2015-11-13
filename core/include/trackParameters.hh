#ifndef TRACKPARAMETERS_H
#define TRACKPARAMETERS_H

#include <vector>
#include <bitset>
#include <stdexcept>
#include <ostream>

#include "fiveByFiveMatrix.hh"
#include "Vector3D.hh"
#include "Vector5.hh"

namespace aidaTT
{
  /** Main class that holds the track parameters, where a track is fully characterized by 
   *  23 floating point numbers.
   *  Five values parametrize a helix with, the reference point adds three values and the 
   *  covariance matrix has 15 independent entries.
   *  The actual definition of the parameters is done in trackParametrization.hh.
   *
   *  For now we use the LCIO track parameters (trackParameterizationLCIO.hh), as described in 
   *  Kraemer, Track parameters in LCIO [http://www-flc.desy.de/lcnotes/notes/LC-DET-2006-004.pdf]
   *
   *  @version $Id
   *  @author Ch. Rosemann, F. Gaede, DESY
   **/

  class trackParameters
  {
  public:
    /// default constructor, sets everything to zero
    trackParameters();
    /// copy ctor 
    trackParameters(const trackParameters&);
    /// assignment
    trackParameters& operator=(const trackParameters&);
    
    /// access helix parameters
    inline const Vector5& parameters() const {
      return _helixparams;
    }
    
    /// write access to helix parameters
    inline Vector5& parameters() {
      return _helixparams;
    }
    
    /// direct read access to the individual helix parameter elements by index - 
    /// use enum from trackParametrization.hh
    inline double operator()(unsigned index) const{
      return _helixparams(index);
    }
    /// direct write access to the individual helix parameter elements by index - 
    /// use enum from trackParametrization.hh
    inline double& operator()(unsigned index){
      return _helixparams(index);
    }

    /// the reference point in global, cartesian coordinates (x,y,z)
    inline const Vector3D& referencePoint() const {
      return _refpoint;
    }
    /// write access to the reference point in global, cartesian coordinates (x,y,z)
    inline Vector3D& referencePoint() {
      return _refpoint;
    }

    /// read access to covariance matrix
    inline const fullCovariance& covarianceMatrix() const {
      return _covmatrix;
    };

    /// write access to covariance matrix
    inline fullCovariance& covarianceMatrix() {
      return _covmatrix;
    };


    /// set everything at once
    void setTrackParameters(const Vector5&, const fullCovariance&, const Vector3D&);

    ///~ only set the helix parameters
    void setTrackParameters(const Vector5&);


    /// set the reference point, defaults to the nominal center
    void setReferencePoint(const Vector3D&);

    ///~ set the covariance matrix
    void setCovarianceMatrix(const fiveByFiveMatrix&);

  private:
    fiveByFiveMatrix _covmatrix;
    Vector5        _helixparams;
    Vector3D       _refpoint;
  };



  /// dump the track parameters to a stream
  inline std::ostream & operator << (std::ostream & os, const trackParameters& tp) {
    os << " [parameters: " << tp.parameters() << "]" << std::endl 
       << " { covariance: " << tp.covarianceMatrix() << " }, "
       << " ( reference point: " << tp.referencePoint() << " )";
    return os ;
  }
}


// this defines the actual param,eterization that is used:
#include "trackParameterization.hh"


#endif // TRACKPARAMETERS_H
