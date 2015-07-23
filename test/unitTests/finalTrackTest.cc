#include "finalTrackTest.hh"

using namespace std;
using namespace aidaTT;
using namespace lcio;

finalTrackTest::finalTrackTest() : UnitTest("FinalTrackTest", __FILE__)
{

}



void finalTrackTest::_test()
{

  IMPL::TrackImpl* finalTrack = new IMPL::TrackImpl() ;

  IMPL::TrackStateImpl* trkState = new IMPL::TrackStateImpl() ;

  aidaTT::trackParameters fTP ;  


  Vector5 paramVec(-0.00021121, 00.874685, -0.411434, 9.3253e-06, -0.000176128 ) ;

  fTP.setTrackParameters(paramVec) ;

  const float covariancematrix[25] = { 1., 0. , 0., 0., 0.,
			 0., 1. , 0., 0., 0.,
			 0., 0. , 1., 0., 0.,
			 0., 0. , 0., 1., 0.,
			 0., 0. , 0., 0., 1.
  };
  
  std::vector<double> cmvec (covariancematrix, covariancematrix + sizeof(covariancematrix) / sizeof(covariancematrix[0]) );

  fTP.setCovarianceMatrix(cmvec);


  trkState = aidaTT::createLCIO( fTP );

  const double hp[5] = { trkState->getOmega() / dd4hep::mm, trkState->getTanLambda(), trkState->getPhi(), trkState->getD0()  * dd4hep::mm, trkState->getZ0() * dd4hep::mm };

  test_(floatCompare(fTP.parameters()(0), hp[0]));
  test_(floatCompare(fTP.parameters()(1), hp[1]));
  test_(floatCompare(fTP.parameters()(2), hp[2]));
  test_(floatCompare(fTP.parameters()(3), hp[3]));
  test_(floatCompare(fTP.parameters()(4), hp[4]));


  fiveByFiveMatrix  CovMat = fTP.covarianceMatrix();
  std::vector<float> cm  =  trkState->getCovMatrix();

  test_(floatCompare(CovMat(0,0), cm[5] / (dd4hep::mm * dd4hep::mm)));
  test_(floatCompare(CovMat(0,1), cm[12] / dd4hep::mm));
  test_(floatCompare(CovMat(0,2), cm[4] / dd4hep::mm));
  test_(floatCompare(CovMat(0,3), cm[3]));
  test_(floatCompare(CovMat(0,4), cm[8]));
  test_(floatCompare(CovMat(1,0), cm[12] / dd4hep::mm));
  test_(floatCompare(CovMat(1,1), cm[14]));
  test_(floatCompare(CovMat(1,2), cm[11]));
  test_(floatCompare(CovMat(1,3), cm[10] * dd4hep::mm));
  test_(floatCompare(CovMat(1,4), cm[13] * dd4hep::mm));
  test_(floatCompare(CovMat(2,0), cm[4] / dd4hep::mm));
  test_(floatCompare(CovMat(2,1), cm[11]));
  test_(floatCompare(CovMat(2,2), cm[2]));
  test_(floatCompare(CovMat(2,3), cm[1] * dd4hep::mm));
  test_(floatCompare(CovMat(2,4), cm[7] * dd4hep::mm));
  test_(floatCompare(CovMat(3,0), cm[3]));
  test_(floatCompare(CovMat(3,1), cm[10] * dd4hep::mm));
  test_(floatCompare(CovMat(3,2), cm[1] * dd4hep::mm));
  test_(floatCompare(CovMat(3,3), cm[0] * dd4hep::mm * dd4hep::mm));
  test_(floatCompare(CovMat(3,4), cm[6] * dd4hep::mm * dd4hep::mm));
  test_(floatCompare(CovMat(4,0), cm[8]));
  test_(floatCompare(CovMat(4,1), cm[13] * dd4hep::mm));
  test_(floatCompare(CovMat(4,2), cm[7] * dd4hep::mm ));
  test_(floatCompare(CovMat(4,3), cm[6] * dd4hep::mm * dd4hep::mm));
  test_(floatCompare(CovMat(4,4), cm[9] * dd4hep::mm * dd4hep::mm));



}



void finalTrackTest::run()
{
    _test();
}
