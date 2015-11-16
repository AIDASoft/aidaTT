#include "trackParameters.hh"

#include <iostream>

namespace aidaTT
{
  trackParameters::trackParameters() : _covmatrix() , _helixparams(), _refpoint()
  {
  }



  trackParameters::trackParameters(const trackParameters& tp) :  _covmatrix(tp._covmatrix) , _helixparams(tp._helixparams), _refpoint(tp._refpoint)
  {
  }



  trackParameters& trackParameters::operator=(const trackParameters& tp)
  {
    if(this == &tp)
      return *this;

    this->_covmatrix = tp._covmatrix;
    this->_helixparams = tp._helixparams;
    this->_refpoint = tp._refpoint;

    return *this;
  }



  void trackParameters::setTrackParameters(const Vector5& parameters, const fullCovariance& covMatrix, const Vector3D& refPoint)
  {
    setTrackParameters(parameters);
    setReferencePoint(refPoint);
    setCovarianceMatrix(covMatrix);
  }



  void trackParameters::setTrackParameters(const Vector5& params)
  {
    _helixparams = params;
  }



  void trackParameters::setReferencePoint(const Vector3D& rp)
  {
    _refpoint = rp;
  }



  void trackParameters::setCovarianceMatrix(const fullCovariance& cov)
  {
    _covmatrix = cov;
  }
}
