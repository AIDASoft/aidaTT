#ifndef TRAJECTORYELEMENT_H
#define TRAJECTORYELEMENT_H

#include "IGeometry.hh"
#include "trackParameters.hh"
#include "fiveByFiveMatrix.hh"

#include <utility>
#include <vector>

namespace aidaTT
{


    /*
     * this is the main class that provides info on the elements that make up a track
     * - the measurements
     * - the material
     * - the information how to get to the next/previous element
     *
     * the definitions of:
     *  what a measurementsurface is
     *  what material it is
     * are external!
     *
     * an element is uniquely identified and identifiable by its arc length
     * for a fixed point of reference (!)
     *  ?? how can this be made explicit
     *
     * An element is anything that is of interest in the track fitting or extrapolation:
     *  - the key elements are the measurements, without them there is no track
     *  - also the material between the measurements or extrapolation points is crucial
     *  - also points of interest, e.g. the calorimeter surface, the beam pipe, the nominal IP, etc.
     *
     * If the trajectoryElement is nothing but a point on the trajectory, then only the jacobian to the next element will be known.
     *
     * If a trajectoryElement is a measurement, then several things must be set:
     *  - the connection to the measurement surface
     *  - the residual(s) associated with the measurement
     *  - the associated uncertainty (or -ties) matching the residual(s)
     *  - the dimensionality of the measurement (although this is clear from the surface)
     *
     * If a trajectoryElement contains material, then the following must be initialized:
     *  - the connection to the material information that allows updating
     *  - the current estimate of the radiation length and atomic charge Z, possibly averaged.
     * [ note: this follows the idea, that any thick scatterer can be approximated by two thin scatters. ]
     *
     * All information that is put into the trajectoryElement can also be retrieved.
     *
     */


    class trajectoryElement
    {
        public:
            ///~ standard constructor A for measurements: arc length is given, the surface it belongs to and some identification
            trajectoryElement(double, const Surface&, const std::vector<double>&, void* = NULL);

            ///~ constructor B: only the arc length is given and some identification
            trajectoryElement(double, void* = NULL);

            ///~ constructor C: only arc length and the jacobian to the next element plus some identification
            trajectoryElement(double, const fiveByFiveMatrix&, void* = NULL);

            ///~ constructor D: everything is already known: arc length, surface, the jacobian to the next element and some identification
            trajectoryElement(double,  const Surface&, const std::vector<double>&, const fiveByFiveMatrix&, void* = NULL);

            /// the getting routines
            const Surface& getSurface() const
            {
                return *_surface;
            };

            const fiveByFiveMatrix& getJacobianToNextElement() const
            {
                return _jacobianToNext;
            };

            std::pair<trackParameters, fullCovariance> getFullState() const;
            trackParameters  getStateVector() const;

            bool hasMeasurement() const
            {
                return (_measDim > 0);
            };

            // the following depend on the type of element:
            unsigned int getMeasurementDimension() const
            {
                return _measDim;
            };
            const std::vector<double>& getMeasurementResiduals() const
            {
                return _residuals;
            } ;
            const std::vector<double>& getMeasurementErrors() const
            {
                return _resolutions;
            } ;

            /// the setting routines

            ///~ set the jacobian to the next element
            void setJacobianToNextElement(const fiveByFiveMatrix&);

        private:
            ///~ no construction without the arc length!
            trajectoryElement();

            ///~ no assigment, no copying
            trajectoryElement(const trajectoryElement&);
            trajectoryElement operator=(const trajectoryElement&);

            void _calcResiduals();
            void _calcMaterial();

            double _arclength;
            const Surface* const _surface;
            unsigned int _measDim;
            std::vector<double> _residuals;
            std::vector<double> _resolutions;
            fiveByFiveMatrix _jacobianToNext;
            const void* const _id; // just store
    };

}

#endif // TRAJECTORYELEMENT_H
