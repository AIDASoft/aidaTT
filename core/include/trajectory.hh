#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <vector>
#include <list>

#include "trackParameters.hh"
#include "trajectoryElement.hh"
#include "IPropagation.hh"
#include "IGeometry.hh"
#include "IFittingAlgorithm.hh"

namespace aidaTT
{

    /*** The main class to provide functionality for track parameter estimation.
     *   This includes propagation and extrapolation to any specific point/surface.
     *
     *   Created with the class aidaTT; uses a specific set of
     *        - fitting algorithm
     *        - propagation method
     *        - geometry interface
     *
     *  @version $Rev$
     *  @author Ch. Rosemann, DESY
     ***/

    class IFittingAlgorithm;

    class trajectory
    {
        public:
            /// create an empty trajectory
            trajectory();

            /// copy construct a trajectory
            trajectory(const trajectory&);

            /// trajectory that are to be fitted need full information at initialization
            trajectory(const trackParameters&, const std::vector<trajectoryElement>&,
                       const IFittingAlgorithm*, const IPropagation*, const IGeometry*);

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
            const std::vector<trajectoryElement>& trajectoryElements() const;

            ///~ add a measurement/hit to trajectory
            void addMeasurement(const Vector3D& position, const ISurface& surface, void* id);

            ///~ test whether/where a surface is intersected
            bool intersectWithSurface(const ISurface* surface, Vector3D& intersect );

            const std::vector<std::pair<double, const ISurface*> >& getIntersectionsWithSurfaces(const std::list<const ISurface*>&);
  

            IFittingAlgorithm* getFittingAlgorithm() const;
            IPropagation* getPropagationMethod() const;
            
            
            ///~ prepare / fit
            void prepareForFitting();
            
            unsigned int fit();

            /// methods after fitting
            std::vector<trajectoryElement> getFittedTrajectoryElements() const;
            std::vector<trajectoryElement> getOutliers() const;

            /// quantify the results
            double getChiSquare() const;
            unsigned int getNDF() const;

        private:
            /* disable assignment */
            trajectory operator=(const trajectory&);

            // the  internal parts
            trackParameters          _referenceParameters;
            std::vector<trajectoryElement>         _initialTrajectoryElements;

            std::vector<std::pair<double, const ISurface*> > _intersectionsList;

            const IFittingAlgorithm* const _fittingAlgorithm;
            const IPropagation* const _propagation;
            const IGeometry* const _geometry;

            bool _intersectsWithinZCylinderBounds(const ISurface*, double&);
            bool _intersectWithinZPlaneBounds(const ISurface*, double&);
            bool _intersectWithinZDiskBounds(const ISurface*, double&);

            double _calculateRadius() const;
            double _calculateXCenter() const;
            double _calculateYCenter() const;

            double _calculatePhifromXY(double, double) const;
            double _calculateSfromXY(double, double) const;
            double _calculateSfromXY(std::pair<double, double> k) const
            {
                return _calculateSfromXY(k.first, k.second);
            }
            double _calculateZfromS(double) const;

            std::pair<Vector3D, Vector3D>* _calculateLocalCurvilinearSystem(double);

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
