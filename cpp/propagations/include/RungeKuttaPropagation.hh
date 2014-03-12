#ifndef RUNGEKUTTAPROPAGATION_H
#define RUNGEKUTTAPROPAGATION_H

#include "IPropagation.hh"

namespace aidaTT
{

class RungeKuttaPropagation : public Ipropagation
{
    fiveByFiveMatrix  getJacobian(double dw, double qop, const Vector3D& tstart, const Vector3D& tend, const Vector3D& bfield);
};

}

#endif // RUNGEKUTTAPROPAGATION_H
