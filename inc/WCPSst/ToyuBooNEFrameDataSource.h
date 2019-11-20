#ifndef WIRECELLSST_TOYUBOONEFRAMEDATASOURCE_H
#define  WIRECELLSST_TOYUBOONEFRAMEDATASOURCE_H

#include "WCPNav/FrameDataSource.h"
#include "WCPSst/RootEvent.h"
#include "WCPNav/GeomDataSource.h"

#include "TTree.h"
#include "TH1F.h"

namespace WCPSst {

    /**
       
     */
    class ToyuBooNEFrameDataSource : public WCP::FrameDataSource {
	mutable TTree* tree;	// or TChain
	WCPSst::RootEvent event;

      public:
	ToyuBooNEFrameDataSource(TTree& tree, const WCP::GeomDataSource& gds,int bins_per_frame1 = 9600);
	virtual ~ToyuBooNEFrameDataSource();

	void Save();

	/// Return the number of frames this data source knows about.  Return -1 if unlimited.
	virtual int size() const;

	/// Explicitly set the "frame" (event) to process.  Frame number returned or -1 on error.
	virtual int jump(int frame_number);

    private:
	const WCP::GeomDataSource& gds;
	int nwire_u, nwire_v, nwire_w;

	TH1F **hu;
	TH1F **hv;
	TH1F **hw;	
    };

}

#endif
