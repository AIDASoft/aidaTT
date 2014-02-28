#ifndef IMEASUREMENTSURFACE_H
#define IMEASUREMENTSURFACE_H

#include "trackParameters.hh"
#include "fiveByFiveMatrix.hh"

namespace aidaTT
{

    /*** Abstract measurement surfaces:: the interface to a measurement layer through a geometry API.
     *
     *  @version $Id
     *  @author Ch. Rosemann, DESY
     *
     * A measurement is characterized by being associated to a certain detector.
     * A detector/a measurement surface has to be uniquely identifiable.
     *
     * To construct the \chi^2  for parameter estimation two parts are needed:
     *      - the uncorrelated measurements m (in n dimensions)
     *      - the measurement directions M (in the according n dimensions)
     *
     * There are only a limited number of cases to consider:
     *      - 2D measurements are most common
     *      - 1D measurements are also used
     *
     * A trajectoryElement that is actually a measurement has three ingrdients:
     *      - a detector id (object)
     *      - the measurement values (in global coordinates)
     *      - the measurement covariance (in the same global coordinates)
     *
     * Required from the geometry is:
     *      - the orthogonal measurement directions
     *      - the orthogonality is defined such as the covariance becomes diagonal
     *
     * Internally there is a local, orthogonal system describing the surface.
     * For any point on the surface, this system can be accessed and then rotated.
     *
     *
     * TODO: what about navigation to a surface -
     *      extrapolation or propagation
     *
     *  proposal: inherit from ISurface
     *      inherit first layer to BoundedPlane, BoundedCylinder
     *      inherit from first to second layer special/analytical propagation solution: Zplane, ZYPlane, ZCylinder
     *
     ***/

    class IMeasurementSurface
    {
        public:
            IMeasurementSurface(const unsigned int, void* = 0);  // dimension plus whatever id is wanted

            ///~ a measurement is characterized by its dimensions
            virtual void* getCharacterization() const = 0;
            virtual unsigned int getDimensionality() const = 0;

            ///~ or rather:
            //~ nDimensionalDirectionVector  getMeasurementDirections(const trajectoryElement& aHit) const;
            //~ nDimensionalMeasurementVector getMeasurementValues( const trajectoryElement& aHit) const;

            ///~ the main functionality is the calculation of the intersection
            ///~ the dimension of the result depends on the dimension of the measurement
            std::vector<double> calculateIntersection(const trackParameters&, const fiveByFiveMatrix&);

    };
}
#endif // IMEASUREMENTSURFACE_H
