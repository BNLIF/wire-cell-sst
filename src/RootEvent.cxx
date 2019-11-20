#include "WCPSst/RootEvent.h"

#include "TClonesArray.h"
#include <vector>

WCPSst::RootEvent::RootEvent()
    : number(-1)
    , run(-1)
    , subrun(-1)
    , nchannels(0)
    , channelid(new std::vector<int>)
    , signal(new TClonesArray)
{
}

WCPSst::RootEvent::~RootEvent()
{
    delete channelid;
    signal->Delete();
    delete signal;
}
