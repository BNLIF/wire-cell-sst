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
      //mutable TTree* tree;	// or TChain
      //WireCellSst::RootEvent event;

      public:
	DatauBooNEFrameDataSource(const char* root_file, const WireCell::GeomDataSource& gds,int bins_per_frame1 = 9600);
	virtual ~DatauBooNEFrameDataSource();

	void Save();
	void Clear();

	/// Return the number of frames this data source knows about.  Return -1 if unlimited.
	virtual int size() const;

	/// Explicitly set the "frame" (event) to process.  Frame number returned or -1 on error.
	virtual int jump(int frame_number);

	bool chirp_check(double rms, int plane, int channel);
	double correlation1(TH1F *h1, TH1F *h2);
	
 	void GetChannelStatus(TH1F *h1, int plane, int chan, bool& isCut, double &rmsOut);
	void zigzag_removal(TH1F *h1, int plane, int channel_no, int flag_RC = 1);
	bool ID_RC(TH1F *h1, int plane, int channel_no);
	 

	void chirp_id(TH1F *h1, int plane, int channel_no);

	void chirp_raise_baseline(TH1F *h1, int bin1, int bin2);

	void SignalFilter(TH1F *h1);
	double CalcRMSWithFlags(TH1F *hist);
	void RemoveFilterFlags(TH1F *hist);
	void RawAdaptiveBaselineAlg(TH1F *hist);
	void NoisyFilterAlg(TH1F *hist, int plane, int channel_no);
	void RemovePMTSignalCollection(TH1F *hist, float rms);

	int get_run_no(){return run_no;};
	int get_subrun_no(){return subrun_no;};
	int get_event_no(){return event_no;};

	WireCell::WireMap& get_u_map(){return uplane_map;};
	WireCell::WireMap& get_v_map(){return vplane_map;};
	WireCell::WireMap& get_w_map(){return wplane_map;};

	WireCell::ChirpMap& get_u_cmap(){return uchirp_map;};
	WireCell::ChirpMap& get_v_cmap(){return vchirp_map;};
	WireCell::ChirpMap& get_w_cmap(){return wchirp_map;};


    private:
	const WireCell::GeomDataSource& gds;
	const char* root_file;
	int nwire_u, nwire_v, nwire_w;

	int nevents;

	int run_no, subrun_no, event_no;

	WireCell::WireMap uplane_map;
	WireCell::WireMap vplane_map;
	WireCell::WireMap wplane_map;

	WireCell::ChirpMap uchirp_map;
	WireCell::ChirpMap vchirp_map;
	WireCell::ChirpMap wchirp_map;
	
	std::map<int, float> urms_map;
	std::map<int, float> vrms_map;
	std::map<int, float> wrms_map;

	// RC+RC 
	TH1F *h_rc;
	TH1 *hm_rc, *hp_rc;
	TH1F *h_1us;
	TH1 *hm_1us, *hp_1us;
	TH1F *h_2us;
	TH1 *hm_2us, *hp_2us;
    };

    

}

#endif
