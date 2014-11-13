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



    Vector5 trackParameters::parameters() const
    {
        return _helixparams;
    }



    double trackParameters::operator()(unsigned int index) const
    {
        if(index > 4)
            throw std::invalid_argument("[trackParameters::operator()] Wrong index when accessing track parameters for reading.");
        else
            return _helixparams(index);
    }



    double& trackParameters::operator()(unsigned int index)
    {
        if(index > 4)
            throw std::invalid_argument("[trackParameters::operator()] Wrong index when accessing track parameters for reading.");
        else
            return _helixparams(index);
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



    #ifdef USE_LCIO
    trackParameters::trackParameters(const EVENT::TrackState* const ts ) : 
		_helixparams(ts->getOmega(), ts->getTanLambda(), ts->getPhi(), ts->getD0(), ts->getZ0()) , _refpoint( ts->getReferencePoint() )
    {
		std::vector<double> covm(ts->getCovMatrix().begin(),ts->getCovMatrix().end() );
		//std::vector<double> refp(ts->getReferencePoint().begin(),ts->getReferencePoint().end() );
		//_covmatrix( covm ); 
		//;

	}
	
	

    IMPL::TrackStateImpl* trackParameters::createPersistentTrackState()
    {
		return new IMPL::TrackStateImpl();
	}
    #endif // USE_LCIO
}
