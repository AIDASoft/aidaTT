#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <vector>

#include "trackParameters.hh"
#include "trajectoryElement.hh"
#include "IPropagation.hh"
#include "IGeometryAPI.hh"
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
     *  @version $Id
     *  @author Ch. Rosemann, DESY
     ***/


    class trajectory
    {
        public:
            /// A trajectory can only be constructed with full information
            trajectory(const trackParameters&, const std::vector<trajectoryElement>&,
                       IFittingAlgorithm*, IPropagation*, IGeometryAPI*);

            ~trajectory();

            /// methods available before fitting
            std::vector<trajectoryElement> getTrajectoryElements();
            std::vector<trajectoryElement> getMeasurements();

            IFittingAlgorithm* getFittingAlgorithm();
            IPropagation* getPropagationMethod();

            bool fit();

            /// methods after fitting
            std::vector<trajectoryElement> getFittedTrajectoryElements();
            std::vector<trajectoryElement> getFittedMeasurements();
            std::vector<trajectoryElement> getOutliers();

/// quantify the results
            double getChiSquare();
            unsigned int getNDF();

        private:
            /* disable default and copy constructor; no assignment */
            trajectory();
            trajectory(const trajectory&);
            trajectory operator=(const trajectory&);

            // the  internal parts
            const trackParameters          _referenceParameters;
            std::vector<trajectoryElement> _initialTrajectoryElements;

            const IFittingAlgorithm* _fittingAlgorithm;
            const IPropagation* _propagation;
            const IGeometryAPI* _geometry;


            // and go on with the legacy code:

            /*

            // adding hits is not so trivial => need to know where they are located:
            *   front & end is easy; other positions are potentially expensive
              virtual int addHit(trajectoryElement* hit) = 0 ;

            // unsure about the initialisation functions. what for?
              virtual int initialise( bool fitDirection ) = 0 ;
              virtual int initialise(  const EVENT::TrackState& ts, double bfield_z, bool fitDirection ) = 0 ;


            // straightforward call
                  virtual int fit( double maxChi2Increment=DBL_MAX ) = 0 ;


            // combination of adding (same problems) and fitting (easy)
              virtual int addAndFit( trajectoryElement* hit, double& chi2increment, double maxChi2Increment=DBL_MAX ) = 0 ;

            // different name, same function
              virtual int testChi2Increment( trajectoryElement* hit, double& chi2increment ) = 0 ;

            // kalmanning
            *   virtual int smooth() = 0 ;
              virtual int smooth( trajectoryElement* hit ) = 0 ;

            // two questions: a) strange concept of a single return b) what's the defition of track state?
              virtual int getTrackState( IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
              virtual int getTrackState( trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;

            // straightforward
              virtual int getHitsInFit( std::vector<std::pair<trajectoryElement*, double> >& hits ) = 0 ;
              virtual int getOutliers( std::vector<std::pair<trajectoryElement*, double> >& hits ) = 0 ;
              virtual int getNDF( int& ndf ) = 0 ;

            /// had this before -- wtf??
              virtual int getTrackerHitAtPositiveNDF( trajectoryElement*& trkhit ) = 0 ;

            // propagation: extrapolation + material fx // method still unclear
            * expensive!? material caching? -> integration along trajectory
            *   simpler if intersection is known
              virtual int propagate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
              virtual int propagate( const gear::Vector3D& point, trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;

             /// oops, this needs a layered geometry. assume?
              virtual int propagateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;
              virtual int propagateToLayer( int layerID, trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;

              // this i don't understand...?
              virtual int propagateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
              virtual int propagateToDetElement( int detEementID, trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;

            ///// same thing, different colour -- extrapolation
              virtual int extrapolate( const gear::Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
              virtual int extrapolate( const gear::Vector3D& point, trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) = 0 ;
              virtual int extrapolateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;
              virtual int extrapolateToLayer( int layerID, trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest ) = 0  ;
              virtual int extrapolateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;
              virtual int extrapolateToDetElement( int detEementID, trajectoryElement* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) = 0  ;

              // again layer construction
              virtual int intersectionWithLayer( int layerID, gear::Vector3D& point, int& detElementID, int mode=modeClosest ) = 0  ;
              virtual int intersectionWithLayer( int layerID, trajectoryElement* hit, gear::Vector3D& point, int& detElementID, int mode=modeClosest ) = 0  ;


              // detelement == surface?
              virtual int intersectionWithDetElement( int detElementID, gear::Vector3D& point, int mode=modeClosest ) = 0  ;
              virtual int intersectionWithDetElement( int detEementID, trajectoryElement* hit, gear::Vector3D& point, int mode=modeClosest ) = 0  ;
            */


    };
}

#endif // TRAJECTORY_H
