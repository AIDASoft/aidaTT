#include "api2gbl.h"

using namespace std;
using namespace aidaTT;
using namespace gbl;


api2gbl::api2gbl()
{
}

/* the routine
 * 
 * loop over all elements in the trajectory
 * 1) create GblPoint with p2pJacobian
 * 2) - (a) if measurement add with [projection2Measurement, residuals, precision]
 * 		(b) if scatterer 
 * 		(c) if neither, then:
 * 3) create trajectory from list of GBL points & fit
 * 4) [if success]: get & prepare results

 * 
 * 			/// Create vector of GBL points.
			vector < GblPoint > theListOfPoints;

			for vector<helperHit*>::const_iterator aHit:
			* 
			* 
				// create a point from the jacobian
				GblPoint point(point2pointJacobian);
				// now add the measurement point
				if ((*aHit)->isMeasurement()) {
					point.addMeasurement((*aHit)->getLocalToMeasurementProjection(), (*aHit)->getResiduals(), (*aHit)->getPrecision());
					if ((*aHit)->getNumGlobals() > 0) point.addGlobals((*aHit)->getGlobalLabels(), (*aHit)->getGlobalDerivatives());
				} else {
					// this is not so nice...? the index is counted 1-based (for technical reasons, not because of fortran)
					_refPointIndex = (aHit - theListOfHits.begin()) + 1;
				}

				// store the point in the list that will be handed to the trajectory
				theListOfPoints.push_back(point);

				// calculate next jacobian, if it's not the last point!
				if (aHit < (endHit - 1))// guard against running over the end
				{
					// which then finally is the 3d path length difference ! don't forget the cos lambda
					double deltaW = ((*(aHit + 1))->getS() - (*aHit)->getS()) / cos(seedLambda);
					double phi1 = (*aHit)->getPhi();
					// now assign the next p2p jacobian in local curvilinear parameters
					point2pointJacobian = hlx.analyticalHelixJacobian(phi1, deltaW);
				}
				delete *aHit;
			}
			theListOfHits.clear();

			/// Create trajectory from list of GBL points.
			_trajectory = new GblTrajectory(theListOfPoints, _curvature);

			/// Fit trajectory.
			if (!_trajectory->fit(_chisquare, _ndf, _lostweight, _fitOptions)) {
				if (_writeMillepedeOut) {
					trackerAlcaSelector sel(*seedTrack);
					if (sel.isSelected()) _trajectory->milleOut(*_milleBinary);
				}
				TVectorD trackparameters(5);
				TMatrixDSym trackcovariance(5);
				if (_ndf > 0) // check if the fit was successful
				{
					// get the results at a given label in local cl track parameters
					// the track parameters are corrections to the curvilinear track parameters
					_trajectory->getResults(_refPointIndex, trackparameters, trackcovariance);
					// the results need to be transformed to the LCIO perigee representation ( d0, phi0, Omega, z0, tanLambda)
					TMatrixD jacCL2LCIO = hlx.perigeeToLCIOJacobian() * hlx.curvilinearToPerigeeJacobian();
					TVectorD correctionVec = jacCL2LCIO * trackparameters;
					TMatrixDSym covarianceMatrix = trackcovariance.Similarity(jacCL2LCIO);
					/* prepare LCIO output; unknown parts first:
					 *
					 * theTrack->setTypeBit (int index, bool val=true) // ?
					 * theTrack->setdEdx (float dEdx)
					 * theTrack->setdEdxError (float dEdxError)
					 * theTrack->setReferencePoint (float *rPnt) */

					TrackImpl* theTrack = new TrackImpl();
					theTrack->setOmega(seedOmega + correctionVec[2]);
					theTrack->setTanLambda(seedTanLambda + correctionVec[4]);
					theTrack->setPhi(seedPhi + correctionVec[1]);
					theTrack->setD0(seedD0 + correctionVec[0]);
					theTrack->setZ0(seedZ0 + correctionVec[3]);
					theTrack->setReferencePoint(refPoint);
					// compressed covariance matrix
					float compressedCovariance[15];
					int ind = 0;
					for (int i = 0; i < 5; ++i)
					for (int j = 0; j <= i; ++j)
					compressedCovariance[ind++] = covarianceMatrix[i][j];
					theTrack->setCovMatrix(compressedCovariance);

					theTrack->setChi2(_chisquare);
					theTrack->setNdf(_ndf);

					// and finally add the hits to the track
					for (vector<TrackerHit*>::const_iterator theHit = seedTrack->getTrackerHits().begin(), hitEnd =
							seedTrack->getTrackerHits().end(); theHit < hitEnd; ++theHit)
					theTrack->addHit(*theHit);

					m_out(DEBUG) << " +++ the fitted track has a chi^2 value of " << _chisquare << " for ndf = " << _ndf
					<< " and the parameters [omega, tanLambda, phi0, d0, z0] = [" << theTrack->getOmega() << "," << theTrack->getTanLambda()
					<< "," << theTrack->getPhi() << "," << theTrack->getD0() << "," << theTrack->getZ0() << "]." << endl;
					/// Fill output collection.
					outputTrackCollection->addElement(theTrack);
				}
			}
			// cleanup
			if (_trajectory) delete _trajectory;

