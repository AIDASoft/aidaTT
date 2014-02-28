#ifndef IMATERIAL_HH
#define IMATERIAL_HH

#include <utility>
#include "fiveByFiveMatrix.hh"

namespace aidaTT
{

    class IMaterial
    {
        public:
            ///~ get the integrated radiation length between two points
            ///~ characterised by their arc length and the method to get from one to another
            double getIntegralX0(double, double, const fiveByFiveMatrix&);

            ///~ calculate the average nuclear charge Z and the effective length between the two points
            std::pair<double, double> getAverageZ(double, double, const fiveByFiveMatrix&);
    };

}

#endif // IMATERIAL_HH
