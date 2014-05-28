#ifndef CONSTANTSOLENOIDBFIELD_HH
#define CONSTANTSOLENOIDBFIELD_HH

/* the most simple magnetic field
 *
 * constant in value, only Bz != 0
 */

#include "IBField.hh"

namespace aidaTT
{

    class ConstantSolenoidBField : public IBField
    {
        public:
            ConstantSolenoidBField(const double field) : _bz(field) {};
            virtual ~ConstantSolenoidBField()
            {
                ;
            };

            Vector3D getBField(const Vector3D&) const
            {
                return Vector3D(0., 0., _bz);
            };
            Vector3D getBField(const double x, const double y, const double z) const
            {
                return Vector3D(0., 0., _bz);
            };

            double Bx(const Vector3D&) const
            {
                return 0.;
            };
            double By(const Vector3D&) const
            {
                return 0.;
            };
            double Bz(const Vector3D&) const
            {
                return _bz;
            };

        private:
            ConstantSolenoidBField();

            const double _bz;
    };
}
#endif // CONSTANTSOLENOIDBFIELD_HH
