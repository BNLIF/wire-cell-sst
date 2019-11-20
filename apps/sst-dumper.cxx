#include "WCPNav/FrameDataSource.h"
#include "WCPNav/SimDataSource.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/Util.h"
#include "WCPData/Trace.h"
#include "WCPData/SimTruth.h"

#include <iostream>
using namespace std;

int main(int argc, char** argv)
{
    if (argc < 1) {
	cerr << "usage: test_dump celltree_file.root [TDirectory/path]" << endl;
	return 1;
    }
    const char* root_file_name = argv[1];

    const char* tpath = "/Event/Sim";
    if (argc > 2) {
	tpath = argv[2];
    }

    // Note, we could cheat and just make an SST version of the FDS
    // and then get direct access to its SimDS mixing facet, but let's
    // pretend to do things proper.
    //WCPSst::FrameDataSource* fds = 0;

    //
    WCP::FrameDataSource* fds = 0;

    fds = WCPSst::make_fds(root_file_name);
    if (!fds) {
	cerr << "ERROR: failed to get FDS from " << root_file_name << endl;
	return 1;
    }

    WCP::SimDataSource* sim = dynamic_cast<WCP::SimDataSource*>(fds);
    if (!sim) {
	cerr << "ERROR: the FDS is not also an SimDS " << endl;
	return 2;
    }

    cerr << "Got " << fds->size() 
	 << " frames from " << tpath 
	 << " in " << root_file_name << endl;

    for (size_t iframe=0; iframe<fds->size(); ++iframe) {

	if (fds->jump(iframe) < 0) {
	    cerr << "Failed to go to frame " << iframe << endl;
	    return 3;
	}

	const WCP::Frame& frame = fds->get();
	cerr << "Frame " << iframe << " has " << frame.traces.size() << " traces" << endl;

	for (size_t itrace=0; itrace < frame.traces.size(); ++itrace) {
	    const WCP::Trace& trace = frame.traces[itrace];
	    
	    cerr << "Frame: " << iframe
		 << " Trace:" << itrace
		 << " chid:" << trace.chid
		 << " tbin:" << trace.tbin
		 << " #bins:" << trace.charge.size()
		 << endl;

	}
	
	WCP::SliceDataSource sds(*fds);
	for (int islice=0; islice < sds.size(); ++islice) {
	    if (sds.jump(islice) < 0) {
		cerr << "Slice " << islice << " failed" << endl;
		return 2;
	    }
	    const WCP::Slice& slice = sds.get();
	    // vector<pair<int,float>>
	    const WCP::Channel::Group& cg = slice.group();
	    float qtot=0;
	    for (int ich=0; ich<cg.size(); ++ich) {
		qtot += cg[ich].second;
	    }
	    cerr << "slice " << islice << " @" << slice.tbin()
		 << " qtot=" << qtot << " in " << cg.size() << " channels" << endl;

	}



	WCP::SimTruthSelection sts = sim->truth();
	cerr << "Got " << sts.size() << " true hits" << endl;
	for (size_t itruth = 0; itruth < sts.size(); ++itruth) {
	    const WCP::SimTruth* st = sts[itruth];
	    cerr << "Hit: "
		 << " @ (" << st->x() << " " << st->y() << " " << st->z() << ")"
		 << " q=" << st->charge()
		 << " tdc=" << st->tdc()
		 << endl;
	}

    } // frames


    return 0;

}
