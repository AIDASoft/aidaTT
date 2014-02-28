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

            Vector3D const getBField(const Vector3D&);

            Vector3D const getBField(const double x, const double y, const double z);

            double const getBx(const Vector3D&);
            double const getBy(const Vector3D&);
            double const getBz(const Vector3D&);

    };
}
#endif // IBFIELD_HH
