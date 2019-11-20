#ifndef WIRECELLSST_ROOTEVENT_H
#define WIRECELLSST_ROOTEVENT_H

#include "TClonesArray.h"
#include <vector>

namespace WCPSst{

    struct RootEvent {
	RootEvent();
	~RootEvent();
	

	int number;
	int run;
	int subrun;
	int nchannels;
	std::vector<int> *channelid;
	TClonesArray* signal;
    };
}

#endif
