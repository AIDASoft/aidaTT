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

  test_(floatCompare(iTP.parameters()(0), hp[0]));
  test_(floatCompare(iTP.parameters()(1), hp[1]));
  test_(floatCompare(iTP.parameters()(2), hp[2]));
  test_(floatCompare(iTP.parameters()(3), hp[3]));
  test_(floatCompare(iTP.parameters()(4), hp[4]));



}



void initialTrackTest::run()
{
    _test();
}
