#ifndef TRAJECTORYELEMENT_H
#define TRAJECTORYELEMENT_H

#include "IMaterial.hh"
#include "IMeasurementSurface.hh"

#include "trackParameters.hh"
#include "fiveByFiveMatrix.hh"
#include <utility>

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
            ///~ constructor A: only arc length plus some identification
            trajectoryElement(double, void*);

            ///~ constructor B: only arc length and the jacobian to the next element plus some identification
            trajectoryElement(double, const fiveByFiveMatrix&, void*);

            /// the getting routines

            const IMeasurementSurface& getMeasurementSurface();
            const IMaterial& getMaterial();

            const fiveByFiveMatrix& getJacobianToNextElement();

            std::pair<trackParameters, fullCovariance> const getFullState();
            trackParameters const getStateVector();

            // CR: unsure, if this is needed here....
            //~ const trajectoryElement& getNextElement();
            //~ const trajectoryElement& getPreviousElement();

            bool hasMaterial() const
            {
                return _material;
            };
            bool hasMeasurement() const
            {
                return _measurement;
            };

            // the following depend on the type of element:
            unsigned int getMeasurementDimension() const;
            const double* getMeasurementResiduals() const;
            const double* getMeasurementErrors() const;

            /// the setting routines

            ///~ set the pointer to the material information
            ///~ this can be used or discarded, if provided
            void setMaterial(const IMaterial&);

            ///~ this sets the most important material properties (e.g. averaged):
            ///~  radiation length and Z
            void addMaterial(double, double);

            ///~ associate the measurement with a geometrical entity
            ///~ this provides access to the local measurement directions
            void setMeasurementSurface(const IMeasurementSurface&);

            ///~ set the actual values of a measurement
            ///~ they are characterized by the values of the residuals and their standard deviation
            void addMeasurement(double* residuals, double* error, unsigned int dimension);

            ///~ set the jacobian to the next element
            void setJacobianToNextElement(const fiveByFiveMatrix&) ;

        private:
            ///~ no construction without the arc length!
            trajectoryElement();
            trajectoryElement(const trajectoryElement&);
            trajectoryElement operator=(const trajectoryElement&);

            bool _material;
            bool _measurement;


            double _arcLength;
            fiveByFiveMatrix _jacobianToNext;
    };

}

#endif // TRAJECTORYELEMENT_H
