#ifndef __OXSX_DIST_FILLER__
#define __OXSX_DIST_FILLER__
#include <CutCollection.h>
#include <EventSystematicManager.h>

class DataSet;
class BinnedED;
class CutLog;
class DistFiller
{
public:
  // default is to take all the events
  static void FillDist(BinnedED &, const DataSet &, const CutCollection &cuts_ = CutCollection(),
                       EventSystematicManager sysMan_ = EventSystematicManager(), int nEv_ = -1);

  static void FillDist(BinnedED &, const DataSet &, const CutCollection &cuts_, CutLog &log_,
                       EventSystematicManager sysMan_ = EventSystematicManager(), int nEv_ = -1);
};

#endif
