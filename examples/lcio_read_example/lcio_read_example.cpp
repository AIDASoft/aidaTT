#ifdef USE_LCIO

#include "lcio.h"
#include "IO/LCReader.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimTrackerHit.h"
#include "UTIL/LCTrackerConf.h"

#include "IGeometry.hh"

#include <map>

using namespace std ;
using namespace lcio;


int main(int argc, char** argv)
{
  if(argc < 2)
    {
      std::cout << " usage: ./lcio_read_example ILDEx.xml ILDExSimu.slcio" << std::endl ;
      return 1;
    }

  std::string inFile =  argv[1] ;

  const aidaTT::IGeometry& geom = aidaTT::IGeometry::instance(inFile) ;

  const std::vector<const aidaTT::ISurface*>& surfaces = geom.getSurfaces() ;


  // create map of surfaces
  std::map< long64, const aidaTT::ISurface* > surfMap ;


  for(std::vector<const aidaTT::ISurface*>::const_iterator surf = surfaces.begin() ; surf != surfaces.end() ; ++surf)
    {
      surfMap[(*surf)->id() ] = (*surf) ;
    }



  /// lcio stuff
  std::string lcioFileName = argv[2] ;

  LCReader* rdr = LCFactory::getInstance()->createLCReader() ;
  rdr->open(lcioFileName) ;

  LCEvent* evt = 0 ;

  std::vector< std::string > colNames ;
  colNames.push_back("VXDCollection") ;

  UTIL::BitField64 idDecoder(LCTrackerCellID::encoding_string()) ;


  while((evt = rdr->readNextEvent()) != 0)
    {
      for(unsigned icol = 0, ncol = colNames.size() ; icol < ncol ; ++icol)
	{

	  LCCollection* col = evt->getCollection(colNames[ icol ]) ;

	  int nHit = col->getNumberOfElements() ;

	  for(int i = 0 ; i < nHit ; ++i)
	    {

	      SimTrackerHit* sHit = (SimTrackerHit*) col->getElementAt(i) ;

	      long64 id = sHit->getCellID0() ;

	      idDecoder.setValue(id) ;
	      //      std::cout << " simhit with cellid : " << idDecoder << std::endl ;

	      const aidaTT::ISurface* surf = surfMap[ id ] ;

	      std::cout << " surface " << (*surf) << " found for id : " << std::hex << id << std::dec  ;
	    }




	}
    }


  return 0;
}


#endif // USE_LCIO
