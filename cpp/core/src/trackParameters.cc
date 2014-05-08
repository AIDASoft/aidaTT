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



    Vector5 trackParameters::getTrackParameters() const
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



    void trackParameters::print() const
    {
        //~ cout << "[trackParameters]::print() <<<<<<< START >>>>>>>>" << endl;
        //~ cout << "[trackParameters] helix parameters: Omega = " << _omega << " , TanLambda = " << _tanlambda << " , Phi0 = " << _phizero << " , D0 = " << _dzero << " , Z0 = " << _zzero << endl;
        //~ cout << "[trackParameters] reference point: [x = " << _xref << " , y = " << _yref << " , z = " << _zref << "]." << endl;
        //~ cout << "[trackParameters] covariance matrix: "
        //~ << "  [ " << _covmatrix.at(0)  <<  " , " << _covmatrix.at(1) <<  " , " << _covmatrix.at(2)
        //~ <<  " , " << _covmatrix.at(3)  <<  " , " << _covmatrix.at(4) <<  " , " << _covmatrix.at(5) <<  " , " << _covmatrix.at(6)
        //~ <<  " , " << _covmatrix.at(7)  <<  " , " << _covmatrix.at(8) <<  " , " << _covmatrix.at(9) <<  " , " << _covmatrix.at(10)
        //~ <<  " , " << _covmatrix.at(11) << " , " << _covmatrix.at(12) <<  " , " << _covmatrix.at(13) <<  " , " << _covmatrix.at(14) << "]"  << endl;

        ///~ [  " << _covmatrix.at(0) <<  " "   " << _covmatrix.at(1) <<  " "   " << _covmatrix.at(3) <<  " "   " << _covmatrix.at(6) <<  " "  " << _covmatrix.at(10) <<  " " ]
        ///~ [  " << _covmatrix.at(1) <<  " "   " << _covmatrix.at(2) <<  " "   " << _covmatrix.at(4) <<  " "   " << _covmatrix.at(7) <<  " "  " << _covmatrix.at(11) <<  " " ]
        ///~ [  " << _covmatrix.at(3) <<  " "   " << _covmatrix.at(4) <<  " "   " << _covmatrix.at(5) <<  " "   " << _covmatrix.at(8) <<  " "  " << _covmatrix.at(12) <<  " " ]
        ///~ [  " << _covmatrix.at(6) <<  " "   " << _covmatrix.at(7) <<  " "   " << _covmatrix.at(8) <<  " "   " << _covmatrix.at(9) <<  " "  " << _covmatrix.at(13) <<  " " ]
        ///~ [ " << _covmatrix.at(10) <<  " "  " << _covmatrix.at(11) <<  " "  " << _covmatrix.at(12) <<  " "  " << _covmatrix.at(13) <<  " "  " << _covmatrix.at(14) <<  " " ]


        std::cout << "[trackParameters]::print() >>>>>>>> END <<<<<<<" << std::endl;
    }
}
