#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <vector>
#include <list>

#include "trackParameters.hh"
#include "trajectoryElement.hh"
#include "helixGymnastics.hh"

#include "IPropagation.hh"
#include "IGeometry.hh"
#include "IBField.hh"
#include "IFittingAlgorithm.hh"
#include "Vector3D.hh"

#include "fitResults.hh"

#ifdef USE_LCIO
#include "lcio.h"
#include "IMPL/TrackImpl.h"
#endif // USE_LCIO

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

            //~ MarlinTrk:: virtual int addHit(trajectoryElement* hit) = 0 ;
            //~ MarlinTrk:: virtual int addAndFit( trajectoryElement* hit, double& chi2increment, double maxChi2Increment=DBL_MAX ) = 0 ;
            //~ MarlinTrk:: virtual int testChi2Increment( trajectoryElement* hit, double& chi2increment ) = 0 ;

            //~ MarlinTrk:: virtual int getHitsInFit( std::vector<std::pair<trajectoryElement*, double> >& hits ) = 0 ;
            //~ MarlinTrk:: virtual int getOutliers( std::vector<std::pair<trajectoryElement*, double> >& hits ) = 0 ;
            //~ MarlinTrk:: virtual int getNDF( int& ndf ) = 0 ;

            //~ MarlinTrk:: virtual int extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
            //~ MarlinTrk:: virtual int extrapolate( const gear::Vector3D& point, trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;

            //~ MarlinTrk:: virtual int propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
            //~ MarlinTrk:: virtual int propagate( const gear::Vector3D& point, trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;

            //~ MarlinTrk:: virtual int getTrackState( IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
            //~ MarlinTrk:: virtual int getTrackState( trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;

            /// methods available before fitting
            const std::vector<trajectoryElement*>& trajectoryElements() const;

            ///~ add a measurement/hit to trajectory; identified by:
            ///~     a position, the resolution, the surface and some id
            void addMeasurement(const Vector3D&, const std::vector<double>&, const ISurface&, void*);

            ///~ TODO:: needs more thought!
            ///~ manually add an element to the trajectory; e.g. a point of interest
            ///~ optional: add a surface; e.g. which only contains material
            void addElement(const Vector3D&, void* id);
            void addElement(const Vector3D&, const ISurface&, void* id);

            ///~ test whether/where a surface is intersected
            bool intersectWithSurface(const ISurface* surface, Vector3D& intersect);

            const std::vector<std::pair<double, const ISurface*> >& getIntersectionsWithSurfaces(const std::list<const ISurface*>&);

            IFittingAlgorithm* getFittingAlgorithm() const;
            IPropagation* getPropagation() const;
            IBField* getBField() const;

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


            #ifdef USE_LCIO
			/// persistency actions 
            /// create a persistent trajectory
            IMPL::TrackImpl* createPersistentTrajectory();

            ///~ create from persistent object
            trajectory(const EVENT::Track* const);
            #endif // USE_LCIO


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

            bool _calculateIntersectionWithSurface(const ISurface*, double&, Vector2D* = NULL);
            void _calculateLocalCoordinates(const ISurface*, const Vector3D&, Vector2D*);

            bool _intersectsWithinZCylinderBounds(const ISurface*, double&, Vector2D* = NULL);
            bool _intersectWithinZPlaneBounds(const ISurface*, double&, Vector2D* = NULL);
            bool _intersectWithinZDiskBounds(const ISurface*, double&, Vector2D* = NULL);



            /*
            // ?
            //~ MarlinTrk::   virtual int initialise( bool fitDirection ) = 0 ;
            //~ MarlinTrk:: virtual int initialise(  const EVENT::TrackState& ts, double bfield_z, bool fitDirection ) = 0 ;

            BREAKS THE CONCEPT

            // kalmanning
            //~ MarlinTrk:: virtual int smooth( trajectoryElement* hit ) = 0 ;

            // no layers available -- can't assume
            //~ MarlinTrk:: virtual int propagateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;
            //~ MarlinTrk:: virtual int propagateToLayer( int layerID, trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;
            //~ MarlinTrk:: virtual int intersectionWithLayer( int layerID, gear::Vector3D& point, int& detElementID, int mode=modeClosest ) = 0  ;
            //~ MarlinTrk:: virtual int intersectionWithLayer( int layerID, trajectoryElement* hit, gear::Vector3D& point, int& detElementID, int mode=modeClosest ) = 0  ;
            //~ MarlinTrk:: virtual int extrapolateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;
            //~ MarlinTrk:: virtual int extrapolateToLayer( int layerID, trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;



            /// DON'T UNDERSTAND
            * Steve: "not needed"
            //~ MarlinTrk:: virtual int getTrackerHitAtPositiveNDF( trajectoryElement*& trkhit ) = 0 ;

            * what's a detElement?
            //~ MarlinTrk:: virtual int propagateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
            //~ MarlinTrk:: virtual int propagateToDetElement( int detEementID, trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
            //~ MarlinTrk:: virtual int extrapolateToDetElement( int detEementID, trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
            //~ MarlinTrk:: virtual int extrapolateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
            //~ MarlinTrk:: virtual int intersectionWithDetElement( int detElementID, gear::Vector3D& point, int mode=modeClosest ) = 0  ;
            //~ MarlinTrk:: virtual int intersectionWithDetElement( int detEementID, trajectoryElement* hit, gear::Vector3D& point, int mode=modeClosest ) = 0  ;
            */
    };
}

#endif // TRAJECTORY_H
