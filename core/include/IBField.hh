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
            virtual ~IBField()
            {
                ;
            };

            Vector3D getBField(const Vector3D&) const;

            Vector3D getBField(const double x, const double y, const double z) const;

            double Bx(const Vector3D&) const;
            double By(const Vector3D&) const;
            double Bz(const Vector3D&) const;
    };
}
#endif // IBFIELD_HH
