/*
 * this is an example of the usage
 * it's based on the example1.cpp in GBL by Claus Kleinwort
 */

#include <iostream>
#include "AidaTT.hh"
#include "ConstantSolenoidBField.hh"
#include "analyticalPropagation.hh"
//#include "SurfaceGeometry.hh"
#include "GBLInterface.hh"

using namespace std;
//using namespace aidaTT;


int main()
{

    // create an AidaTT master object
//  aidaTT::AidaTT *a = new aidaTT::AidaTT();

// create the different objects needed for fitting
// first a constant field parallel to z, 1T
    aidaTT::ConstantSolenoidBField bfield(1.);

    // create the propagation object
    aidaTT::analyticalPropagation* propagation = new aidaTT::analyticalPropagation();

    // create the geometry object
    // aidaTT::SurfaceGeometry* = new SurfaceGeometry();

    // create the fitter object
    aidaTT::GBLInterface* fitter = new aidaTT::GBLInterface();

    // now create a list of trajectory objects for the fitter



}



/// Simple example.
/**
 * Create points on initial trajectory, create trajectory from points,
 * fit and write trajectory to MP-II binary file,
 * get track parameter corrections and covariance matrix at points.
 *
 * Equidistant measurement layers and thin scatterers, propagation
 * with simple jacobian (quadratic in arc length differences).
 * Curvilinear system (U,V,T) as local coordinate system.
 *
 * This example simulates and refits tracks in a system of planar detectors
 * with 2D measurements in a constant magnet field in Z direction using
 * the curvilinear system as local system and (Q/P, slopes, offsets) as
 * local track parameters. The true track parameters are
 * randomly smeared with respect to a (constant and straight) reference
 * trajectory with direction (lambda, phi) and are used (only) for the
 * on-the-fly simulation of the measurements and scatterers. The predictions
 * from the reference trajectory are therefore always zero and the residuals
 * needed (by addMeasurement) are equal to the measurements.
 */
/*


    unsigned int nTry = 10000; //: number of tries
    unsigned int nLayer = 10; //: number of detector layers
    std::cout << " Gbltst $Rev: 93 $ " << nTry << ", " << nLayer << std::endl;

    TRandom *r = new TRandom3();

    clock_t startTime = clock();


// track direction - no magnetic field
    double sinLambda = 0.3;
    double cosLambda = sqrt(1.0 - sinLambda * sinLambda);
    double sinPhi = 0.;
    double cosPhi = sqrt(1.0 - sinPhi * sinPhi);

// define curvilinear coordinate system (UVT)
// tDir = (cosLambda * cosPhi, cosLambda * sinPhi, sinLambda)
// U = Z x T / |Z x T|, V = T x U
    TMatrixD uvDir(2, 3);
    uvDir[0][0] = -sinPhi;
    uvDir[0][1] = cosPhi;
    uvDir[0][2] = 0.;
    uvDir[1][0] = -sinLambda * cosPhi;
    uvDir[1][1] = -sinLambda * sinPhi;
    uvDir[1][2] = cosLambda;


// measurement resolution
    TVectorD measErr(2);
    measErr[0] = 0.001;
    measErr[1] = 0.001;

    TVectorD measPrec(2); // (independent) precisions
    measPrec[0] = 1.0 / (measErr[0] * measErr[0]);
    measPrec[1] = 1.0 / (measErr[1] * measErr[1]);

    TMatrixDSym measInvCov(2); // inverse covariance matrix
    measInvCov.Zero();
    measInvCov[0][0] = measPrec[0];
    measInvCov[1][1] = measPrec[1];

// scattering error
    TVectorD scatErr(2);
    scatErr[0] = 0.001;
    scatErr[1] = 0.001;
    TVectorD scatPrec(2);
    scatPrec[0] = 1.0 / (scatErr[0] * scatErr[0]);
    scatPrec[1] = 1.0 / (scatErr[1] * scatErr[1]);

// (RMS of) CurviLinear track parameters (Q/P, slopes, offsets)
    TVectorD clPar(5);
    TVectorD clErr(5);
    clErr[0] = 0.001;
    clErr[1] = -0.1;
    clErr[2] = 0.2;
    clErr[3] = -0.15;
    clErr[4] = 0.25;

    TMatrixDSym clCov(5), clSeed(5);
    unsigned int seedLabel = 0;

    double bfac = 0.2998; // Bz*c for Bz=1
    double step = 1.5 / cosLambda; // constant steps in RPhi

    double Chi2Sum = 0.;
    int NdfSum = 0;
    double LostSum = 0.;
    int numFit = 0;

    for (unsigned int iTry = 1; iTry <= nTry; ++iTry) {
        // curvilinear track parameters
        for (unsigned int i = 0; i < 5; ++i) {
            clPar[i] = clErr[i] * r->Gaus();
        }
        clCov.Zero();
        for (unsigned int i = 0; i < 5; ++i) {
            clCov[i][i] = 1.0 * (clErr[i] * clErr[i]);
        }
//      std::cout << " Try " << iTry << ":" << clPar << std::endl;
        TMatrixD addDer(2, 2);
        addDer.Zero();
        addDer[0][0] = 1.;
        addDer[1][1] = 1.;
// arclength
        double s = 0.;
        TMatrixD jacPointToPoint(5, 5);
        jacPointToPoint.UnitMatrix();
// create list of points
        std::vector<GblPoint> listOfPoints;
        listOfPoints.reserve(2 * nLayer);

        for (unsigned int iLayer = 0; iLayer < nLayer; ++iLayer) {
//          std::cout << " Layer " << iLayer << ", " << s << std::endl;
//     measurement directions
            double sinStereo = (iLayer % 2 == 0) ? 0. : 0.1;
            double cosStereo = sqrt(1.0 - sinStereo * sinStereo);
            TMatrixD mDirT(3, 2);
            mDirT.Zero();
            mDirT[1][0] = cosStereo;
            mDirT[2][0] = sinStereo;
            mDirT[1][1] = -sinStereo;
            mDirT[2][1] = cosStereo;
// projection measurement to local (curvilinear uv) directions (duv/dm)
            TMatrixD proM2l = uvDir * mDirT;
// projection local (uv) to measurement directions (dm/duv)
            TMatrixD proL2m = proM2l;
            proL2m.Invert();
// point with (independent) measurements (in measurement system)
            GblPoint point(jacPointToPoint);
// measurement - prediction in measurement system with error
            TVectorD meas = proL2m * clPar.GetSub(3, 4);
//MP            meas += addDer * addPar; // additional parameters
            for (unsigned int i = 0; i < 2; ++i) {
                meas[i] += measErr[i] * r->Gaus();
            }
            point.addMeasurement(proL2m, meas, measPrec);
            /* point with (correlated) measurements (in local system)
             GblPoint point(jacPointToPoint);
             // measurement - prediction in local system with error
             TVectorD meas(2);
             for (unsigned int i = 0; i < 2; ++i) {
             meas[i] = measErr[i] * r->Gaus() + addDer(i,0) * 0.0075;
             }
             meas = proM2l * meas + clPar.GetSub(3, 4);
             TMatrixDSym localInvCov = measInvCov;
             localInvCov.SimilarityT(proL2m);
             point.addMeasurement(meas, localInvCov);

             // additional local parameters? * /
//          point.addLocals(addDer);
//MP            point.addGlobals(globalLabels, addDer);
            addDer *= -1.; // Der flips sign every measurement
// add point to trajectory
            listOfPoints.push_back(point);
            unsigned int iLabel = listOfPoints.size();
            if (iLabel == seedLabel) {
                clSeed = clCov.Invert();
            }
// propagate to scatterer
            jacPointToPoint = gblSimpleJacobian(step, cosLambda, bfac);
            clPar = jacPointToPoint * clPar;
            clCov = clCov.Similarity(jacPointToPoint);
            s += step;
            if (iLayer < nLayer - 1) {
                TVectorD scat(2);
                scat.Zero();
                // point with scatterer
                GblPoint point(jacPointToPoint);
                point.addScatterer(scat, scatPrec);
                listOfPoints.push_back(point);
                iLabel = listOfPoints.size();
                if (iLabel == seedLabel) {
                    clSeed = clCov.Invert();
                }
                // scatter a little
                for (unsigned int i = 0; i < 2; ++i) {
                    clPar[i + 1] += scatErr[i] * r->Gaus();
                    clCov[i + 1] += scatErr[i] * scatErr[i];
                }
                // propagate to next measurement layer
                clPar = jacPointToPoint * clPar;
                clCov = clCov.Similarity(jacPointToPoint);
                s += step;
            }
        }
//
        // create trajectory
        //GblTrajectory traj(listOfPoints);
        GblTrajectory traj(listOfPoints, seedLabel, clSeed); // with external seed
        //traj.printPoints();
        if (not traj.isValid()) {
            std::cout << " Invalid GblTrajectory -> skip" << std::endl;
            continue;
        }
// fit trajectory
        double Chi2;
        int Ndf;
        double lostWeight;
        traj.fit(Chi2, Ndf, lostWeight);
//      std::cout << " Fit: " << Chi2 << ", " << Ndf << ", " << lostWeight << std::endl;
        /*       TVectorD aCorrection(5);
         TMatrixDSym aCovariance(5);
         traj.getResults(1, aCorrection, aCovariance);
         std::cout << " cor " << std::endl;
         aCorrection.Print();
         std::cout << " cov " << std::endl;
         aCovariance.Print(); * /
    /* look at residuals
     for (unsigned int label=1; label<=listOfPoints.size(); ++label) {
       unsigned int numData=0;
       std::cout << " measResults, label " << label << std::endl;
       TVectorD residuals(2), measErr(2), resErr(2), downWeights(2);
       traj.getMeasResults(label, numData, residuals, measErr, resErr, downWeights);
       std::cout << " measResults, numData " << numData << std::endl;
       // residuals.Print(); measErr.Print(); resErr.Print();
     }  /
// debug printout
        //traj.printTrajectory();
        //traj.printPoints();
        //traj.printData();
// write to MP binary file
//MP        traj.milleOut(mille);
        Chi2Sum += Chi2;
        NdfSum += Ndf;
        LostSum += lostWeight;
        numFit++;
    }
    clock_t endTime = clock();
    double diff = endTime - startTime;
    double cps = CLOCKS_PER_SEC;
    std::cout << " Time elapsed " << diff / cps << " s" << std::endl;
    std::cout << " Chi2/Ndf = " << Chi2Sum / NdfSum << std::endl;
    std::cout << " Tracks fitted " << numFit << std::endl;
*/

