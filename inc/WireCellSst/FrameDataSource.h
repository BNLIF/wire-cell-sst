#ifndef WIRECELLSST_FRAMEDATASOURCE_H
#define  WIRECELLSST_FRAMEDATASOURCE_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/SimDataSource.h"
#include "WireCellData/SimTruth.h"
#include "WireCellSst/RootEvent.h"
#include "WireCellSst/RootSimTruth.h"

#include "TTree.h"

namespace WireCellSst {

    /**
       
     */
    class FrameDataSource : public WireCell::FrameDataSource, virtual public WireCell::SimDataSource {
	mutable TTree *event_tree, *sim_tree;	// or TChain
	WireCellSst::RootEvent event;
	mutable WireCellSst::RootSimTruth rootsimtruth;
	mutable WireCell::SimTruthSet simtruth;

      public:
	FrameDataSource(TTree& event_tree);
	virtual ~FrameDataSource();

	/// Return the number of frames this data source knows about.  Return -1 if unlimited.
	virtual int size() const;

	/// Explicitly set the "frame" (event) to process.  Frame number returned or -1 on error.
	virtual int jump(int frame_number);

	/// Must explicitly set sim tree if truth() is going to be useful
	void set_sim_tree(TTree& sim_tree);

	/// Access to the sim truth objects
	WireCell::SimTruthSelection truth() const; 
	
	

    };

}

#endif
