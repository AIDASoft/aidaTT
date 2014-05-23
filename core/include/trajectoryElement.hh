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
            trajectoryElement(double, const ISurface&, const std::vector<double>&, void* = NULL);

            ///~ constructor B: only the arc length is given and some identification
            trajectoryElement(double, void* = NULL);

            //~ ///~ constructor C: only arc length and the jacobian from the previous element plus some identification
            //~ trajectoryElement(double, const fiveByFiveMatrix&, void* = NULL);
//~ 
            //~ ///~ constructor D: everything is already known: arc length, surface, the jacobian to the next element and some identification
            //~ trajectoryElement(double,  const ISurface&, const std::vector<double>&, const fiveByFiveMatrix&, void* = NULL);

            /// the getting routines
            const ISurface& surface() const
            {
                return *_surface;
            };

            const fiveByFiveMatrix& jacobian() const
            {
                return _jacobianFromPrevious;
            };

            const fiveByFiveMatrix& jacobianFromPrevious() const
            {
                return _jacobianFromPrevious;
            };

            trackParameters  fullState() const;

            bool hasMeasurement() const
            {
               return  _measurement;
            };

            // the following depend on the type of element:
            unsigned int measurementDimension() const;
            
            const std::vector<double>& measurementResiduals() const;
            
            const std::vector<double>& measurementErrors() const;
            
            ///~ access the local curvilinear system at the given point
            const std::pair<Vector3D,Vector3D>& localCurvilinearSystem() const;
            
            ///~ access the measurement direction (if available) for the corresponding surface
            const std::vector<Vector3D>& measurementDirections() const;
            
            ///~ and finally: the projection matrix from the local track frame to the measurement system
            const std::vector<Vector3D>& localToMeasurementProjection() const;

            /// the setting routines

            ///~ set the jacobian from the previous element
            void setJacobian(const fiveByFiveMatrix&);

        private:
            ///~ no construction without the arc length!
            trajectoryElement();

            ///~ no assigment, no copying
            trajectoryElement(const trajectoryElement&);
            trajectoryElement operator=(const trajectoryElement&);

            void _calcResiduals();
            void _calcMaterial();

            ///~ 
            double _arclength;
            const ISurface* const _surface;
            fiveByFiveMatrix _jacobianFromPrevious;
            
            ///~ measurement variables:
            bool _measurement;
            
            std::vector<double> _residuals;
            std::vector<double> _resolutions;
            
            std::vector<Vector3D>* _measDirections;
      
            const void* const _id; // just store
    };

}

#endif // TRAJECTORYELEMENT_H
