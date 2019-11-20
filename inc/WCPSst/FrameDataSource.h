#ifndef WIRECELLSST_FRAMEDATASOURCE_H
#define  WIRECELLSST_FRAMEDATASOURCE_H

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/SimDataSource.h"
#include "WCPData/SimTruth.h"
#include "WCPSst/RootEvent.h"
#include "WCPSst/RootSimTruth.h"

#include "TTree.h"

namespace WCPSst {

    /**
       
     */
    class FrameDataSource : public WCP::FrameDataSource, virtual public WCP::SimDataSource {
	mutable TTree *event_tree, *sim_tree;	// or TChain
	WCPSst::RootEvent event;
	mutable WCPSst::RootSimTruth rootsimtruth;
	mutable WCP::SimTruthSet simtruth;

      public:
	FrameDataSource(TTree& event_tree, const char* br="calib");
	virtual ~FrameDataSource();

	/// Return the number of frames this data source knows about.  Return -1 if unlimited.
	virtual int size() const;

	/// Explicitly set the "frame" (event) to process.  Frame number returned or -1 on error.
	virtual int jump(int frame_number);

	/// Must explicitly set sim tree if truth() is going to be useful
	void set_sim_tree(TTree& sim_tree);

	/// Access to the sim truth objects
	WCP::SimTruthSelection truth() const; 
	
	

    };

}

#endif
