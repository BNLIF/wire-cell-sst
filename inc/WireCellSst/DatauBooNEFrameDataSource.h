#ifndef WIRECELLSST_DATAUBOONEFRAMEDATASOURCE_H
#define  WIRECELLSST_DATAUBOONEFRAMEDATASOURCE_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellSst/RootEvent.h"
#include "WireCellNav/GeomDataSource.h"
#include "WireCellData/GeomWire.h"

#include "TTree.h"
#include "TH1F.h"

namespace WireCellSst {

 

    /**
       
     */
    class DatauBooNEFrameDataSource : public WireCell::FrameDataSource {
	mutable TTree* tree;	// or TChain
	WireCellSst::RootEvent event;

      public:
	DatauBooNEFrameDataSource(TTree& tree, const WireCell::GeomDataSource& gds,int bins_per_frame1 = 9600);
	virtual ~DatauBooNEFrameDataSource();

	void Save();

	/// Return the number of frames this data source knows about.  Return -1 if unlimited.
	virtual int size() const;

	/// Explicitly set the "frame" (event) to process.  Frame number returned or -1 on error.
	virtual int jump(int frame_number);

	bool chirp_check(double rms, int plane, int channel);
	double correlation1(TH1F *h1, TH1F *h2);
	
	WireCell::WireMap& get_u_map(){return uplane_map;};
	WireCell::WireMap& get_v_map(){return vplane_map;};
	WireCell::WireMap& get_w_map(){return wplane_map;};

    private:
	const WireCell::GeomDataSource& gds;
	int nwire_u, nwire_v, nwire_w;

	 WireCell::WireMap uplane_map;
	 WireCell::WireMap vplane_map;
	 WireCell::WireMap wplane_map;
	 
	TH1F **hu;
	TH1F **hv;
	TH1F **hw;	
    };

    

}

#endif
