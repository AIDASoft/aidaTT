#ifdef USE_LCIO

#ifndef LCIOPERSISTENCY_HH
#define LCIOPERSISTENCY_HH

#include "lcio.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/TrackStateImpl.h"

#include "trackParameters.hh"
#include "trajectory.hh"

namespace aidaTT
{
    trackParameters readLCIO(const EVENT::TrackState* const ts) ;

    IMPL::TrackStateImpl* createLCIO(const trackParameters& tp);


//~ trajectory readLCIO(const EVENT::Track* const t)
//~ {
//~ }
//~
//~
//~
//~ IMPL::TrackImpl* createLCIO(const trajectory& traj)
//~ {
//~ }
}
#endif // LCIOPERSISTENCY_HH
#endif // USE_LCIO
