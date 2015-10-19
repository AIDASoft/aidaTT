#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <vector>
#include <list>
#include <TMath.h>

#include "trackParameters.hh"
#include "trajectoryElement.hh"
#include "helixGymnastics.hh"
#include "utilities.hh"
#include "helixHelpers.hh"

#include "IPropagation.hh"
#include "IGeometry.hh"
#include "IBField.hh"
#include "IFittingAlgorithm.hh"
#include "Vector3D.hh"

#include "fitResults.hh"

namespace aidaTT
{
    /// The central class that provides functionality.
    /** The main class to provide functionality for track parameter estimation.
     *   This includes propagation and extrapolation to any specific point/surface.
     *
     *   Created with the class aidaTT; uses a specific set of
     *        - fitting algorithm
     *        - propagation method
     *        - geometry interface
     *
     *  @version $Rev$
     *  @author Ch. Rosemann, DESY
     **/

    class IFittingAlgorithm;

    class trajectory
    {
        public:
            /// create an empty trajectory
            trajectory();

            /// copy construct a trajectory -- NOT the internals
            trajectory(const trajectory&);

            /// constructor (A) trajectory that needs to be fitted from a reference track
            trajectory(const trackParameters&, IFittingAlgorithm*, const IBField*, IPropagation*, const IGeometry*);

            //~ the minimal useful constructor
            trajectory(const trackParameters&, const IGeometry*);

            ~trajectory();

            trackParameters getInitialTrackParameters() const
            {
                return _referenceParameters;
            };

            /// methods available before fitting
            const std::vector<trajectoryElement*>& trajectoryElements() const;

            ///~ add a measurement/hit to trajectory; identified by:
            ///~     a position, the resolution, the surface and some id
            void addMeasurement(const Vector3D&, const std::vector<double>&, const ISurface&, void*);
            //void addMeasurement(const Vector3D&, const TMatrixDSym&, const ISurface&, void*);

            ///~ TODO:: needs more thought!
            ///~ manually add an element to the trajectory; e.g. a point of interest
            ///~ optional: add a surface; e.g. which only contains material
            void addElement(const Vector3D&, void* id);
            void addElement(const Vector3D&, const ISurface&, void* id);
            // add a scatterer
      void addScatterer(const Vector3D& position, std::vector<double>& precision, const ISurface& surface, const trackParameters& seed_tp, void* id);
           //void addScatterer(const Vector3D& position, TMatrixDSym& precision, const ISurface& surface, const trackParameters& seed_tp, void* id); // alternative constructor for addScatterer that utilises a symmetric matrix for precision

            ///~ test whether/where a surface is intersected
            bool intersectWithSurface(const ISurface* surface, Vector3D& intersect);

            const std::vector<std::pair<double, const ISurface*> >& getIntersectionsWithSurfaces(const std::list<const ISurface*>&);

            IFittingAlgorithm* getFittingAlgorithm() const;
            IPropagation* getPropagation() const;
            IBField* getBField() const;

            /// helper method: compute the point on the trajectory for a given value of s
            Vector3D pointAt( double s) ;


            ///~ prepare: add scattering material, sort elements, calculate and add the jacobians to all elements
            void prepareForFitting();

            ///~ the actual fit call
            bool fit();

            const fitResults& getFitResults();

            /// TODO:: placeholder; only z component in constant bfield
            double Bz() const
            {
                return _bfieldZ;
            };

            // methods after fitting
            //~ std::vector<trajectoryElement*> getFittedTrajectoryElements() const;
            //~ std::vector<trajectoryElement*> getOutliers() const;

            //~ /// quantify the results
            //~ double getChiSquare() const;
            //~ unsigned int getNDF() const;

        private:

            /// disable assignment
            trajectory operator=(const trajectory&);

            // the  internal parts
            trackParameters          _referenceParameters;
            std::vector<trajectoryElement*>  _initialTrajectoryElements;

            std::vector<std::pair<double, const ISurface*> > _intersectionsList;

            IFittingAlgorithm* _fittingAlgorithm;
            const IBField* const _bfield;
            IPropagation* _propagation;
            const IGeometry* const _geometry;

            double _bfieldZ;

    public:

      bool _calculateIntersectionWithSurface(const ISurface*, double&, Vector2D* = NULL, Vector3D* = NULL);
      void _calculateLocalCoordinates(const ISurface*, const Vector3D&, Vector2D*, Vector3D* = NULL);
      
      bool _intersectsWithinZCylinderBounds(const ISurface*, double&, Vector2D* = NULL, Vector3D* = NULL);
      bool _intersectWithinZPlaneBounds(const ISurface*, double&, Vector2D* = NULL, Vector3D* = NULL);
      bool _intersectWithinZDiskBounds(const ISurface*, double&, Vector2D* = NULL, Vector3D* = NULL);

    };
}

#endif // TRAJECTORY_H
