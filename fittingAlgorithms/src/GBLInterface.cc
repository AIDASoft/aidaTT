#ifdef USE_GBL
#include "GBLInterface.hh"
#include "utilities.hh"

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

                std::cout << " ---  GBLInterface::fit - element : " <<  **element << std::endl ;


                if((*element)->hasMeasurement())
                    {
                        /// three elements are needed to add a measurement to a gblpoint:
                        /// 1) the projection matrix from the local track frame to the measurement system
                        /// 2) the residuals in the measurement directions
                        /// 3) the inverse resolution = precision of the measurements

                        const unsigned int mDim = (*element)->measurementDimension();

                        if(mDim > 2)
                            throw std::invalid_argument("Error: Currently only 1D or 2D measurements are implemented.");

                        //~ 1) projection matrix -- the basis change matrix
                        const std::vector<double>& projLocal2Meas = (*element)->localToMeasurementProjection();

                        /// convert the data from the vector into the matrix:
                        /// convention is that the first row comes first in the data
                        TMatrixD pL2M(2, 2);
                        pL2M(0, 0) = projLocal2Meas.at(0);
                        pL2M(0, 1) = projLocal2Meas.at(1);
                        pL2M(1, 0) = projLocal2Meas.at(2);
                        pL2M(1, 1) = projLocal2Meas.at(3);


                        //~ 2) the residuals in the measurement direction
                        const std::vector<double>& residuals = (*element)->measurementResiduals();

                        //~ convert the vector to an array:
                        const double* resid = residuals.data();

                        //~ 3) the precision of the measurements -- the inverse of the resolution
                        const std::vector<double>& precision = (*element)->precisions();

                        //~ convert the vector to an array:
                        const double* prec = precision.data();

                        // fixed size of arguments: 2D in measurements!
                        point.addMeasurement(pL2M, TVectorD(2, resid), TVectorD(2, prec));
                    }
                if((*element)->isScatterer())
                    {
                        // TODO :: implement!
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
        /// TODO :: check validity before continuing!

        //~ this is not elegant -- delete previous data before starting again
        if(_trajectory != NULL)
            delete _trajectory;

        _trajectory = new gbl::GblTrajectory(theListOfPoints, true); /// TODO: pass info about magnetic field

        unsigned int returnValue = _trajectory->fit(_chisquare, _ndf, _lostweight);

        if(returnValue == 0)
            {
                _fillResults(TRAJ);
                return true;
            }
        else
            return false;
    }



    void GBLInterface::_fillResults(const trajectory& TRAJ)
    {
        bool v = true;
        double chs = _chisquare;
        unsigned int n = _ndf;
        double wl = _lostweight;

        TVectorD tpCorr(5);
        TMatrixDSym trackcovariance(5);

        //~ get the results at a given label in local cl track parameters
        //~ the track parameters are corrections to the curvilinear track parameters
        _trajectory->getResults(0, tpCorr, trackcovariance);

        Vector5 clCorrections(tpCorr[0], tpCorr[1], tpCorr[2], tpCorr[3], tpCorr[4]);

        fiveByFiveMatrix cl2L3Jacobian  = curvilinearToL3Jacobian(TRAJ.getInitialTrackParameters(), Vector3D(0., 0., TRAJ.Bz())) ;
        Vector5 L3corrections           =  cl2L3Jacobian * clCorrections;

        Vector5 fittedParameters = TRAJ.getInitialTrackParameters().parameters() + L3corrections;


        trackParameters tp;
        tp.setTrackParameters(fittedParameters);

        fiveByFiveMatrix testCovMat;
        for(int i = 0 ; i < 5 ; i++)
            {
                for(int j = 0 ; j < 5 ; j++)
                    {
                        testCovMat(i, j) = trackcovariance(i, j);
                    }
            }

        //fg --- need to transform the covariance matrix from CL to perigee ----
        fiveByFiveMatrix cl2L3JacobianT(cl2L3Jacobian) ;
        cl2L3JacobianT.Transpose() ;

        // std::cout << " initial covariance: " << testCovMat << std::endl ;
        // std::cout << " Jacobian :          " <<  cl2L3Jacobian << std::endl ;
        // std::cout << " JacobianT:          " <<  cl2L3JacobianT << std::endl ;

        fiveByFiveMatrix finalCov  = testCovMat * cl2L3JacobianT ;

        // std::cout << " half transformed:  " << finalCov   << std::endl ;

        finalCov = cl2L3Jacobian * finalCov ;

        // std::cout << " final   covariance: " << finalCov   << std::endl ;

        tp.setCovarianceMatrix(finalCov);

        _theResults.setResults(v, chs, n, wl, tp);

    }

}


#endif
