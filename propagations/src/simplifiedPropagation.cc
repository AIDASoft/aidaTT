#include "simplifiedPropagation.hh"
#include <cmath>

namespace aidaTT
{
  bool simplifiedPropagation::getJacobian(fiveByFiveMatrix& jac, double dw, double qop, const Vector3D& tStart, const Vector3D& tEnd, const Vector3D& bfield, double)
    {

      std::cout << " what about this one????? " << std::endl;

        jac.Unit();
        const double dz = tStart.z() - tEnd.z();
        const double cosLambda = 1. / sqrt(1. + dz * dz);
        jac(2, 0) = -bfield.z() * dw;
        jac(3, 0) = -0.5 * bfield.z() * dw * dw * cosLambda;
        jac(3, 2) = dw * cosLambda;
        jac(4, 1) = dw;

        return true;
    }
}
