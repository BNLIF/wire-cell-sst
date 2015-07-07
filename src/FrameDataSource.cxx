#include "WireCellSst/FrameDataSource.h"

#include "TClonesArray.h"
#include "TH1F.h"

#include <iostream>		// debug


WireCellSst::FrameDataSource::FrameDataSource(TTree& ttree)
    : WireCell::FrameDataSource()
    , event_tree(&ttree)
    , sim_tree(0)
    , event()
    , simtruth()

{
    // sigh, we can't do things this simply because the ttree does not
    // have a single branch.  
    // tree->SetBranchAddress(name, &event);

    event_tree->SetBranchAddress("eventNo" , &event.number);
    event_tree->SetBranchAddress("runNo"   , &event.run);
    event_tree->SetBranchAddress("subRunNo", &event.subrun);

    event_tree->SetBranchAddress("calib_nChannel", &event.nchannels);
    event_tree->SetBranchAddress("calib_channelId", &event.channelid);
    event_tree->SetBranchAddress("calib_wf", &event.signal);

    
}

WireCellSst::FrameDataSource::~FrameDataSource()
{
  frame.clear();
}

int WireCellSst::FrameDataSource::size() const
{
    return event_tree->GetEntries();
}

int WireCellSst::FrameDataSource::jump(int frame_number)
{
    if (frame.index == frame_number) {
	return frame_number;
    }

    frame.clear();		// win or lose, we start anew

    if (frame_number < 0) {	// underflow
	return frame_number;
    }

    int siz = event_tree->GetEntry(frame_number);
    if (siz <= 0 ) {
	return -1;
    }

    if (frame_number >= siz) {
	return -1;
    }

    // load into frame
    int nchannels = event.channelid->size();
    for (size_t ind=0; ind < nchannels; ++ind) {
	TH1F* signal = dynamic_cast<TH1F*>(event.signal->At(ind));
	if (!signal) {
	    return -1;
	}

	WireCell::Trace trace;
	trace.chid = event.channelid->at(ind);

	trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
	for (int ibin=1; ibin <= signal->GetNbinsX(); ++ibin) {
	    trace.charge.push_back(signal->GetBinContent(ibin));
	    //std::cerr << ibin << " " << signal->GetBinContent(ibin) << std::endl;
	}
	frame.traces.push_back(trace);
    }
    frame.index = frame_number;
    return frame.index;
}

void WireCellSst::FrameDataSource::set_sim_tree(TTree& tree)
{
    sim_tree = &tree;

    tree.SetBranchAddress("simide_size", &rootsimtruth.size);
    tree.SetBranchAddress("simide_channelIdY", &rootsimtruth.channelIdY);
    tree.SetBranchAddress("simide_trackId", &rootsimtruth.trackId);
    tree.SetBranchAddress("simide_tdc", &rootsimtruth.tdc);
    tree.SetBranchAddress("simide_x", &rootsimtruth.x);
    tree.SetBranchAddress("simide_y", &rootsimtruth.y);
    tree.SetBranchAddress("simide_z", &rootsimtruth.z);
    tree.SetBranchAddress("simide_numElectrons", &rootsimtruth.q);
}

WireCell::SimTruthSelection WireCellSst::FrameDataSource::truth() const
{
    WireCell::SimTruthSelection ret;
    if (! sim_tree) {
	return ret;
    }

    if (frame.index < 0) {	// underflow
	return ret;
    }

    simtruth.clear();
    rootsimtruth.clear();
    int siz = sim_tree->GetEntry(frame.index);
    if (siz <= 0 || frame.index >= siz) {
	return ret;
    }


    for (int ind=0; ind < rootsimtruth.size; ++ind) {
	WireCell::SimTruth st((*rootsimtruth.x)[ind],
			      (*rootsimtruth.y)[ind],
			      (*rootsimtruth.z)[ind],
			      (*rootsimtruth.q)[ind],
			      (*rootsimtruth.tdc)[ind],
			      (*rootsimtruth.trackId)[ind]);
	ret.push_back(&(*simtruth.insert(st).first));
    }
    return ret;
}


