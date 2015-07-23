#include "initialTrackTest.hh"

using namespace std;
using namespace aidaTT;
using namespace lcio;

initialTrackTest::initialTrackTest() : UnitTest("InitialTrackTest", __FILE__)
{

}



void initialTrackTest::_test()
{

  IMPL::TrackImpl* initialTrack = new IMPL::TrackImpl();

  IMPL::TrackStateImpl* trkState = new IMPL::TrackStateImpl() ;

  double chi2 = 225.256 ;
  int ndf = 64 ;

  const float rp[3] = {0., 0., 0.};


  const float cm[25] = { 1., 0. , 0., 0., 0.,
			 0., 1. , 0., 0., 0.,
			 0., 0. , 1., 0., 0.,
			 0., 0. , 0., 1., 0.,
			 0., 0. , 0., 0., 1.
  };
  
  std::vector<double> cmvec (cm, cm + sizeof(cm) / sizeof(cm[0]) );
  
  trkState->setOmega( -0.00021121 );
  trkState->setTanLambda( 00.874685 );
  trkState->setPhi( -0.411434 );
  trkState->setD0( 9.3253e-06 );
  trkState->setZ0( -0.000176128 );

  trkState->setReferencePoint( rp );

  trkState->setCovMatrix( cm ) ;

  trkState->setLocation(lcio::TrackState::AtIP);

  initialTrack->addTrackState(trkState);
  initialTrack->setChi2(chi2) ;
  initialTrack->setNdf(ndf) ;

  const double hp[5] = { trkState->getOmega() / dd4hep::mm, trkState->getTanLambda(), trkState->getPhi(), trkState->getD0()  * dd4hep::mm, trkState->getZ0() * dd4hep::mm };

  aidaTT::trackParameters iTP(  aidaTT::readLCIO( initialTrack->getTrackState( lcio::TrackState::AtIP) )    );
  iTP.setCovarianceMatrix(cmvec);

  fiveByFiveMatrix  CovMat = iTP.covarianceMatrix();

  test_(floatCompare(iTP.parameters()(0), hp[0]));
  test_(floatCompare(iTP.parameters()(1), hp[1]));
  test_(floatCompare(iTP.parameters()(2), hp[2]));
  test_(floatCompare(iTP.parameters()(3), hp[3]));
  test_(floatCompare(iTP.parameters()(4), hp[4]));

  test_(floatCompare(CovMat(0,0), cm[0]));
  test_(floatCompare(CovMat(0,1), cm[1]));
  test_(floatCompare(CovMat(0,2), cm[2]));
  test_(floatCompare(CovMat(0,3), cm[3]));
  test_(floatCompare(CovMat(0,4), cm[4]));
  test_(floatCompare(CovMat(1,0), cm[5]));
  test_(floatCompare(CovMat(1,1), cm[6]));
  test_(floatCompare(CovMat(1,2), cm[7]));
  test_(floatCompare(CovMat(1,3), cm[8]));
  test_(floatCompare(CovMat(1,4), cm[9]));
  test_(floatCompare(CovMat(2,0), cm[10]));
  test_(floatCompare(CovMat(2,1), cm[11]));
  test_(floatCompare(CovMat(2,2), cm[12]));
  test_(floatCompare(CovMat(2,3), cm[13]));
  test_(floatCompare(CovMat(2,4), cm[14]));
  test_(floatCompare(CovMat(3,0), cm[15]));
  test_(floatCompare(CovMat(3,1), cm[16]));
  test_(floatCompare(CovMat(3,2), cm[17]));
  test_(floatCompare(CovMat(3,3), cm[18]));
  test_(floatCompare(CovMat(3,4), cm[19]));
  test_(floatCompare(CovMat(4,0), cm[20]));
  test_(floatCompare(CovMat(4,1), cm[21]));
  test_(floatCompare(CovMat(4,2), cm[22]));
  test_(floatCompare(CovMat(4,3), cm[23]));
  test_(floatCompare(CovMat(4,4), cm[24]));


  double calcLambda = calculateLambda(iTP);
  double calcCurv = calculateCurvature(iTP);


  test_(floatCompare(iTP.parameters()(0), calcCurv));
  test_(floatCompare(atan(iTP.parameters()(1)), calcLambda));

}



void initialTrackTest::run()
{
    _test();
}
