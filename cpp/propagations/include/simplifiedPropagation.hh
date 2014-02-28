#ifndef SIMPLIFIEDPROPAGATION_HH
#define SIMPLIFIEDPROPAGATION_HH

#include "IPropagation.hh"

namespace aidaTT
{

class simplifiedPropagation : public IPropagation
{
    fiveByFiveMatrix  getJacobian(double dw, double qop, const Vector3D& tstart, const Vector3D& tend, const Vector3D& bfield);
};

}

#endif //SIMPLIFIEDPROPAGATION_HH
