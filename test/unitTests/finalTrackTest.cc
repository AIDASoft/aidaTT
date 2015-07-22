#include "finalTrackTest.hh"

using namespace std;
using namespace aidaTT;
using namespace lcio;

finalTrackTest::finalTrackTest() : UnitTest("FinalTrackTest", __FILE__)
{

}



void finalTrackTest::_test()
{

  IMPL::TrackImpl* finalTrack = new IMPL::TrackImpl();

  IMPL::TrackStateImpl* trkState = new IMPL::TrackStateImpl() ;

  aidaTT::trackParameters fTP ;  


  const Vector5& paramVec(-0.00021121, 00.874685, -0.411434, 9.3253e-06, -0.000176128 ) ;

  fTP.setTrackParameters(paramVec) ;


  trkState = aidaTT::createLCIO( fTP );

  const double hp[5] = { trkState->getOmega() / dd4hep::mm, trkState->getTanLambda(), trkState->getPhi(), trkState->getD0()  * dd4hep::mm, trkState->getZ0() * dd4hep::mm };

  test_(floatCompare(fTP.parameters()(0), hp[0]));
  test_(floatCompare(fTP.parameters()(1), hp[1]));
  test_(floatCompare(fTP.parameters()(2), hp[2]));
  test_(floatCompare(fTP.parameters()(3), hp[3]));
  test_(floatCompare(fTP.parameters()(4), hp[4]));


}



void finalTrackTest::run()
{
    _test();
}
