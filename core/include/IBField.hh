#ifndef IBFIELD_HH
#define IBFIELD_HH

#include "Vector3D.hh"

namespace aidaTT
{


    /* simple interface to a bfield implementation
     *
     * provides the magnetic field in three components,
     * taken as global cartesian coordinates.
     *
     */


    class IBField
    {
        public:
            virtual Vector3D BField(const Vector3D&) const = 0;

            virtual Vector3D BField(const double x, const double y, const double z) const = 0;

            virtual double Bx(const Vector3D&) const = 0;
            virtual double By(const Vector3D&) const = 0;
            virtual double Bz(const Vector3D&) const = 0;
    };
}
#endif // IBFIELD_HH
