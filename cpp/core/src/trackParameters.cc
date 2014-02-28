#include "trackParameters.hh"

#include <iostream>

using namespace aidaTT;
using namespace std;



trackParameters::trackParameters() : _omega(0.), _tanlambda(0.), _phizero(0.), _dzero(0.), _zzero(0.), _xref(0.), _yref(0.), _zref(0.)
{
    _covmatrix.resize(15);
    _helixparams.resize(5);
    _refpoint.resize(3);
    for(unsigned int i = 0; i < 15; ++i)
        _covmatrix[i] = 0.;
}



trackParameters::trackParameters(const trackParameters& tp)
{
    _covmatrix.resize(15);
    _helixparams.resize(5);
    _refpoint.resize(3);
    for(unsigned int i = 0; i < 15; ++i)
        _covmatrix[i] = tp.getCovarianceMatrix().at(i);

    ////! INCOMPLETE
    //~ tp.getTrackParameters();
    //~ tp.getReferencePoint();
}



trackParameters trackParameters::operator=(const trackParameters& tp)
{
    return trackParameters();
}



vector<double> trackParameters::getTrackParameters() const
{
    return _helixparams;
}



void trackParameters::setTrackParameters(const vector<double>& parameters, const vector<double>& refPoint, const vector<double>& covMatrix)
{
    setTrackParameters(parameters);
    setReferencePoint(refPoint);
    setCovarianceMatrix(covMatrix);
}



void trackParameters::setTrackParameters(const vector<double>& p)
{
    _omega     = p.at(0);
    _tanlambda = p.at(1);
    _phizero   = p.at(2);
    _dzero     = p.at(3);
    _zzero     = p.at(4);
    _helixparams = p;
}




void trackParameters::setReferencePoint(const vector<double>& rp)
{
    _xref = rp.at(0);
    _yref = rp.at(1);
    _zref = rp.at(2);
    _refpoint = rp;
}



void trackParameters::print() const
{
    cout << "[trackParameters]::print() <<<<<<< START >>>>>>>>" << endl;
    cout << "[trackParameters] helix parameters: Omega = " << _omega << " , TanLambda = " << _tanlambda << " , Phi0 = " << _phizero << " , D0 = " << _dzero << " , Z0 = " << _zzero << endl;
    cout << "[trackParameters] reference point: [x = " << _xref << " , y = " << _yref << " , z = " << _zref << "]." << endl;
    cout << "[trackParameters] covariance matrix: "
         << "  [ " << _covmatrix.at(0)  <<  " , " << _covmatrix.at(1) <<  " , " << _covmatrix.at(2)
         <<  " , " << _covmatrix.at(3)  <<  " , " << _covmatrix.at(4) <<  " , " << _covmatrix.at(5) <<  " , " << _covmatrix.at(6)
         <<  " , " << _covmatrix.at(7)  <<  " , " << _covmatrix.at(8) <<  " , " << _covmatrix.at(9) <<  " , " << _covmatrix.at(10)
         <<  " , " << _covmatrix.at(11) << " , " << _covmatrix.at(12) <<  " , " << _covmatrix.at(13) <<  " , " << _covmatrix.at(14) << "]"  << endl;

    ///~ [  " << _covmatrix.at(0) <<  " "   " << _covmatrix.at(1) <<  " "   " << _covmatrix.at(3) <<  " "   " << _covmatrix.at(6) <<  " "  " << _covmatrix.at(10) <<  " " ]
    ///~ [  " << _covmatrix.at(1) <<  " "   " << _covmatrix.at(2) <<  " "   " << _covmatrix.at(4) <<  " "   " << _covmatrix.at(7) <<  " "  " << _covmatrix.at(11) <<  " " ]
    ///~ [  " << _covmatrix.at(3) <<  " "   " << _covmatrix.at(4) <<  " "   " << _covmatrix.at(5) <<  " "   " << _covmatrix.at(8) <<  " "  " << _covmatrix.at(12) <<  " " ]
    ///~ [  " << _covmatrix.at(6) <<  " "   " << _covmatrix.at(7) <<  " "   " << _covmatrix.at(8) <<  " "   " << _covmatrix.at(9) <<  " "  " << _covmatrix.at(13) <<  " " ]
    ///~ [ " << _covmatrix.at(10) <<  " "  " << _covmatrix.at(11) <<  " "  " << _covmatrix.at(12) <<  " "  " << _covmatrix.at(13) <<  " "  " << _covmatrix.at(14) <<  " " ]


    cout << "[trackParameters]::print() >>>>>>>> END <<<<<<<" << endl;


}
