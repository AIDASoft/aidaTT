#ifdef USE_GBL
#include "GBLInterface.hh"

namespace aidaTT
{
    GBLInterface::GBLInterface() : _trajectory(NULL), _correctionVector(NULL), _covarianceMatrix(NULL)
    {}



    GBLInterface::~GBLInterface()
    {
        if(_trajectory != NULL)
            delete _trajectory;
    }



    bool GBLInterface::fit(const trajectory& TRAJ)
    {
        /* several bits of information are needed to initialize the gbl:
         *  - a vector of GblPoints, which in turn need a p2p jacobian to be instantiated
         *  -- a measurement GblPoint needs the projection matrix, the residual vector and the precision
         *  -- a scattering GblPoint needs the precision (expected inverse standard deviation)
         */

        /// create vector of GBL points
        std::vector < gbl::GblPoint > theListOfPoints;

        const std::vector<trajectoryElement*>& elements = TRAJ.trajectoryElements();

        for(std::vector<trajectoryElement*>::const_iterator element = elements.begin(), last = elements.end(); element < last; ++element)
            {

                const fiveByFiveMatrix& jac = (*element)->jacobian();

                ///~ initialise point with jacobian from last to the current element
                gbl::GblPoint point(TMatrixD(5, 5, jac.array()));

                if((*element)->hasMeasurement())
                    {
                        /// three elements are needed to add a measurement to a gblpoint:
                        /// 1) the projection matrix from the local track frame to the measurement system
                        /// 2) the residuals in the measurement directions
                        /// 3) the inverse resolution = precision of the measurements

                        const unsigned int mDim = (*element)->measurementDimension();

                        //~ 1) projection matrix -- the basis change matrix
                        const std::vector<Vector3D>& projL2M = (*element)->localToMeasurementProjection();

                        ///~ convert the data from the vector into an array in two steps:
                        std::vector<double> projElements;
                        for(std::vector<Vector3D>::const_iterator it = projL2M.begin(), last = projL2M.end(); it < last; ++it)
                            {
                                projElements.push_back(it->x());
                                projElements.push_back(it->y());
                                projElements.push_back(it->z());
                            }
                        const double* pL2M = projElements.data();

                        //~ 2) the residuals in the measurement direction
                        const std::vector<double>& residuals = (*element)->measurementResiduals();

                        //~ convert the vector to an array:
                        const double* res = residuals.data();

                        //~ 3) the precision of the measurements -- the inverse of the resolution
                        const std::vector<double>& errors = (*element)->measurementErrors();

                        std::vector<double> precision;
                        for(unsigned int errorid = 0; errorid < errors.size(); ++errorid)
                            {
                                if(errors.at(errorid) <= 1e-18)
                                    precision.push_back(1e+99);
                                else
                                    precision.push_back(1. / errors.at(errorid));
                            }
                        //~ convert the vector to an array:
                        const double* prec = precision.data();

                        point.addMeasurement(TMatrixD(3, mDim, pL2M), TVectorD(mDim, res), TVectorD(mDim, prec));
                    }
                if((*element)->isScatterer())
                    {
                        //~ now evaluate the scattering material
                        //~ the addScatterer routine will add a thin scatterer at the given position
                        //~ a thin scatter only changes the local direction, no offset. Multiple step approach:
                        //~ -- 1.) determine whether scatterer or not
                        //~ -- 2.) if scatterer: thin or thick
                        //~ -- 3.) if thick scatterer: calculate two positions and material properties at the points!
                        //~ the scattering info enters through the inverse covariance matrix, doesn't need to be diagonalized before invocation
                        //~ point.gbl::addScatterer ( TVectorD notNeededHere, TMatrixDSym aPrecision   );


                    }

                // store the point in the list that will be handed to the trajectory
                theListOfPoints.push_back(point);
            }

//~ if( !theListOfPoints.isValid() )
        //~ throw something;

        _trajectory = new gbl::GblTrajectory(theListOfPoints, true); /// TODO: pass info about magnetic field


        _trajectory->fit(_chisquare, _ndf, _lostweight);

        return true;
    }
}
/*  ORIGINAL


 (1) the field information
        Vector3D bfield = getBField();
        double Bz = bfield.z();// rho component is accessible with field.rho()


  (2)    // loop over the seed tracks
        for inputIter in seedTrackCollectionVector {

          (a) get values from seed track --

            Track* seedTrack = dynamic_cast<Track *>(*inputIter);
            double seedD0 = seedTrack->getD0();
            double seedPhi = seedTrack->getPhi();
            double seedOmega = seedTrack->getOmega();
            double seedZ0 = seedTrack->getZ0();
            double seedTanLambda = seedTrack->getTanLambda();
            const double seedLambda = atan(seedTanLambda);

            //CHK track parameters are defined at (distance of closest approach to) reference point !?
            const float* refPoint = seedTrack->getReferencePoint();
            double refPos[] = {refPoint[0], refPoint[1], refPoint[2]};

            //CHK arc-length is measured from reference point
            float rInner = seedTrack->getRadiusOfInnermostHit();

            vector<TrackerHit*> hitList = seedTrack->getTrackerHits();


            if (Bz == NULL) {
                _curvature = false;
            } else {
                _curvature = true;
            }


         (b) build reference helix
            simpleHelix hlx(-seedOmega, seedPhi, -seedD0, seedTanLambda, seedZ0, refPos, Bz * 0.0002998);

         (c) order hits: first create list with arc length
            /// Order hits by arc-length
            double offset = 0.;
            vector<indexArcPair> hitWithArcList;
            for (unsigned int i=0; i< hitList.size(); ++i) {
                double position[] = {hitList[i]->getPosition()[0], hitList[i]->getPosition()[1]};
                double sArc = hlx.getArcLengthXY(position) - offset;
                hitWithArcList.push_back( make_pair(i, sArc));
            }
            typedef std::pair<int, double> indexArcPair;

            // find smallest arc-length (and first/last hit)
            indexArcPair closest(hitWithArcList.front()), first(closest), last(closest);

            for (vector<indexArcPair>::const_iterator thePair = hitWithArcList.begin()+1; thePair < hitWithArcList.end(); ++thePair) {
                if (abs(thePair->second) < abs(closest.second)) closest = *thePair;
                if (thePair->second < first.second) first = *thePair;
                if (thePair->second > last.second) last = *thePair;
            }

            /// Add reference point if not at a hit (or close by)
            const double averageStep = (last.second-first.second)/(hitWithArcList.size()-1);
            if (abs(closest.second/averageStep) > 1.0E-5) {
                closest = make_pair(-1, 0.);
                hitWithArcList.push_back( closest );
            }

            sort(hitWithArcList.begin(), hitWithArcList.end(), compare);


            /// Create vector of GBL points. from helper hits.
            vector < GblPoint > theListOfPoints;

            // number of valid hits (with measurement)
            unsigned int numValid = 0;

            // point 2 point jacobian in local curvilinear parameters
            TMatrixD point2pointJacobian(5, 5);
            double sOld = 0.;

            // Energy loss estimators.
            vector <double> eLossEstimators;

            // number of hits lost (due to threshold)
            unsigned int numLow = 0;

            for (vector<indexArcPair>::const_iterator thePair = hitWithArcList.begin(); thePair < hitWithArcList.end(); ++thePair) {
                int index = thePair->first;
                gblHelperHit aHit = index >= 0 ? gblHelperHit(*hitList[index], hlx, _writeMillepedeOut, _encodedModuleID) : gblHelperHit(seedPhi);
                if (aHit.isValid()) {
                    // jacobian from previous point
                    double deltaW = (aHit.getS() - sOld) / cos(seedLambda);
                    double phi1 = aHit.getPhi();
                    point2pointJacobian = hlx.analyticalHelixJacobian(phi1, deltaW);
                    sOld = aHit.getS();

                    // create a point from the jacobian
                    GblPoint point(point2pointJacobian);
                    // now add the measurement point
                    if (aHit.isMeasurement()) {
                        numValid++;
                        point.addMeasurement(aHit.getLocalToMeasurementProjection(), aHit.getResiduals(), aHit.getPrecision());
                        if (aHit.getNumGlobals() > 0) point.addGlobals(aHit.getGlobalLabels(), aHit.getGlobalDerivatives());
                        if (aHit.getEScaled() > 0.) eLossEstimators.push_back(1./sqrt(aHit.getEScaled()));
                        else numLow++;
                    }

                    // store the point in the list that will be handed to the trajectory
                    theListOfPoints.push_back(point);

                    // is reference point?
                    if (index == closest.first) _refPointIndex = theListOfPoints.size();
                } else m_out(DEBUG) << " invalid hit: " << thePair->first << ", " << thePair->second << endl;
            }
            if (numValid < 3) {
                m_out(DEBUG) << " ++ track number " << endIter - inputIter << " skipped - too few valid hits " << endl;
                continue;
            }



            /// Create trajectory from list of GBL points.
            _trajectory = new GblTrajectory(theListOfPoints, _curvature);

            /// Fit trajectory.
            if (!_trajectory->fit(_chisquare, _ndf, _lostweight, _fitOptions)) {

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

                    TrackImpl* theTrack = new TrackImpl();
                    theTrack->setOmega(seedOmega + correctionVec[2]);
                    theTrack->setTanLambda(seedTanLambda + correctionVec[4]);
                    theTrack->setPhi(seedPhi + correctionVec[1]);
                    theTrack->setD0(seedD0 + correctionVec[0]);
                    theTrack->setZ0(seedZ0 + correctionVec[3]);
                    theTrack->setReferencePoint(refPoint);
                    theTrack->setRadiusOfInnermostHit(rInner);

                    // compressed covariance matrix
                    float compressedCovariance[15];
                    int ind = 0;
                    for (int i = 0; i < 5; ++i)
                    for (int j = 0; j <= i; ++j)
                    compressedCovariance[ind++] = covarianceMatrix[i][j];
                    theTrack->setCovMatrix(compressedCovariance);

                    theTrack->setChi2(_chisquare);
                    theTrack->setNdf(_ndf);

                    // simple dEdx from median(1/sqrt(EScaled))
                    const unsigned int numEstimators = eLossEstimators.size();

                    if (numEstimators > numLow) {
                        sort(eLossEstimators.begin(), eLossEstimators.end());
                        const int nTot = numEstimators + numLow;
                        const int indMed = nTot/2;
                        const double eMed = (nTot % 2 == 0) ? 0.5*(eLossEstimators[indMed-1]+eLossEstimators[indMed]) :eLossEstimators[indMed]; // median
                        vector <double> absDev;
                        for (unsigned int i=0; i<numEstimators; ++i) absDev.push_back(abs(eLossEstimators[i]-eMed));
                        sort(absDev.begin(), absDev.end());
                        const double eMAD = (nTot % 2 == 0) ? 0.5*(absDev[indMed-1]+absDev[indMed]) :absDev[indMed];// median of absolute deviations, sigma(gauss) = 1.4826 * MAD
                        const float ndEdx = nTot;
                        const float dEdx = 1.0/(eMed*eMed);
                        const float dEdxError = 2. * 1.4826 * eMAD * dEdx / eMed / sqrt(ndEdx);
                        m_out(DEBUG) << " dEdx " << ndEdx << " " << dEdx << " " << dEdxError << endl;
                        theTrack->setdEdx(dEdx);
                        theTrack->setdEdxError(dEdxError);
                    }
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
        }
        // final storage to the event
        if (_outputIsPersistent) evt->addCollection(outputTrackCollection, _outputTrackCollectionName);

    */

#endif
