#include "utilities.hh"

namespace aidaTT
{

    /// SUPER DUMMY !!! only L3 type: [ Omega, tan(lambda), phi_0, d_0m z_0, ]
    double calculateX0(const trackParameters& tp)
    {
        /// L3 type only
        return cos(tp.getTrackParameters()(2)) * tp.getTrackParameters()(3);
    }



    double calculateY0(const trackParameters& tp)
    {
        /// L3 type only
        return sin(tp.getTrackParameters()(2)) * tp.getTrackParameters()(3);
    }



    double calculatePhi0(const trackParameters& tp)
    {
        /// L3 type only
        return tp.getTrackParameters()(2);
    }



    double  calculateTanLambda(const trackParameters& tp)
    {
        /// L3 type only
        return tp.getTrackParameters()(1);
    }



    double calculateZ0(const trackParameters& tp)
    {
        /// L3 type only
        return tp.getTrackParameters()(4);
    }
}
