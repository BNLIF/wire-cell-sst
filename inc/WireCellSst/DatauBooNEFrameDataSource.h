#ifndef WIRECELLSST_DATAUBOONEFRAMEDATASOURCE_H
#define  WIRECELLSST_DATAUBOONEFRAMEDATASOURCE_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellSst/RootEvent.h"
#include "WireCellNav/GeomDataSource.h"
#include "WireCellData/GeomWire.h"
#include "WireCellData/PMTNoiseROI.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

namespace WireCellSst {

 

    /**
       
     */
    class DatauBooNEFrameDataSource : public WireCell::FrameDataSource {
      //mutable TTree* tree;	// or TChain
      //WireCellSst::RootEvent event;

      public:
      DatauBooNEFrameDataSource(const char* root_file, const WireCell::GeomDataSource& gds,int bins_per_frame1 = 9600, int inactiveNF = 0, int flag_add_noise = 0);
	DatauBooNEFrameDataSource(const TH2F *hu_raw, const TH2F *hv_raw, const TH2F *hw_raw, TTree *T_bad, TTree *T_lf, TTree *Trun, const WireCell::GeomDataSource& gds);

	virtual ~DatauBooNEFrameDataSource();

	void Save();
	void Clear();

	/// Return the number of frames this data source knows about.  Return -1 if unlimited.
	virtual int size() const;

	/// Explicitly set the "frame" (event) to process.  Frame number returned or -1 on error.
	virtual int jump(int frame_number);

    // Dummy noise filter for processing simulation w/o noise
	virtual int jump_no_noise(int frame_number);

	void fix_ADC_shift(int chid, TH1F *h1);
	int shift_right(int value, int n, int filling, int totalBit);
	int lowest_bits(int value, int n);


	bool chirp_check(double rms, int plane, int channel);
	double correlation1(TH1F *h1, TH1F *h2);
	
 	void GetChannelStatus(TH1F *h1, int plane, int chan, bool& isCut, double &rmsOut);
	void zigzag_removal(TH1F *h1, int plane, int channel_no, int flag_RC = 1, int flag_restore=0);
	bool ID_RC(TH1F *h1, int plane, int channel_no);

	bool ID_lf_noisy(TH1F *h1);

	void chirp_id(TH1F *h1, int plane, int channel_no);

	void chirp_raise_baseline(TH1F *h1, int bin1, int bin2);

	void SignalFilter(TH1F *h1);
	double CalcRMSWithFlags(TH1F *hist);
	void RemoveFilterFlags(TH1F *hist);
	void RawAdaptiveBaselineAlg(TH1F *hist);
	void NoisyFilterAlg(TH1F *hist, int plane, int channel_no);
	
	void IDPMTSignalCollection(TH1F *hist, float rms, int channel);
	void IDPMTSignalInduction(TH1F *hist, float rms, int plane, int channel);
	void RemovePMTSignal(TH1F *hist, int start_bin, int end_bin, int flag=0);
	
	void Simu_Noise_uBooNE_Empirical(TH1F *h1, Int_t plane, Int_t channel);

	int get_run_no(){return run_no;};
	int get_subrun_no(){return subrun_no;};
	int get_event_no(){return event_no;};

	WireCell::WireMap& get_u_map(){return uplane_map;};
	WireCell::WireMap& get_v_map(){return vplane_map;};
	WireCell::WireMap& get_w_map(){return wplane_map;};

	WireCell::ChirpMap& get_u_cmap(){return uchirp_map;};
	WireCell::ChirpMap& get_v_cmap(){return vchirp_map;};
	WireCell::ChirpMap& get_w_cmap(){return wchirp_map;};

	std::set<int>& get_lf_noisy_channels(){return lf_noisy_channels;};

    private:
	std::vector<WireCell::PMTNoiseROI*> PMT_ROIs;
	
	int flag_add_noise;
	int inactiveNF;

	const WireCell::GeomDataSource& gds;
	const char* root_file;
	int nwire_u, nwire_v, nwire_w;

	int nevents;
	
	int flag_mis_config;
	
	bool load_results_from_file;

	int run_no, subrun_no, event_no;

	WireCell::WireMap uplane_map;
	WireCell::WireMap vplane_map;
	WireCell::WireMap wplane_map;

	WireCell::ChirpMap uchirp_map;
	WireCell::ChirpMap vchirp_map;
	WireCell::ChirpMap wchirp_map;
	
	std::set<int> lf_noisy_channels;

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

	/* TH1F *hu_resp; */
	/* TH1F *hv_resp; */
    };

    

}

#endif
