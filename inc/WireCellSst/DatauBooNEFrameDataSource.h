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
	
	void zigzag_removal(TH1F *h1, int plane, int channel_no);
	void chirp_id(TH1F *h1, int plane, int channel_no);
	void SignalFilter(TH1F *h1);
	double CalcRMSWithFlags(TH1F *hist);
	void RemoveFilterFlags(TH1F *hist);
	void RawAdaptiveBaselineAlg(TH1F *hist);
	void NoisyFilterAlg(TH1F *hist, int plane, int channel_no);

	WireCell::WireMap& get_u_map(){return uplane_map;};
	WireCell::WireMap& get_v_map(){return vplane_map;};
	WireCell::WireMap& get_w_map(){return wplane_map;};

	WireCell::ChirpMap& get_u_cmap(){return uchirp_map;};
	WireCell::ChirpMap& get_v_cmap(){return vchirp_map;};
	WireCell::ChirpMap& get_w_cmap(){return wchirp_map;};


    private:
	const WireCell::GeomDataSource& gds;
	int nwire_u, nwire_v, nwire_w;

	WireCell::WireMap uplane_map;
	WireCell::WireMap vplane_map;
	WireCell::WireMap wplane_map;

	WireCell::ChirpMap uchirp_map;
	WireCell::ChirpMap vchirp_map;
	WireCell::ChirpMap wchirp_map;
	
	std::map<int, float> urms_map;
	std::map<int, float> vrms_map;
	std::map<int, float> wrms_map;
	
	TH1F **hu;
	TH1F **hv;
	TH1F **hw;	
    };

    

}

#endif
