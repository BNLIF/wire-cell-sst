#ifndef WIRECELLSST_NEWDATAUBOONEFRAMEDATASOURCE_H
#define WIRECELLSST_NEWDATAUBOONEFRAMEDATASOURCE_H

#include "WCPNav/FrameDataSource.h"
#include "WCPSst/RootEvent.h"
#include "WCPNav/GeomDataSource.h"
#include "WCPData/GeomWire.h"

#include "TTree.h"
#include "TH1F.h"

namespace WCPSst {

 

    /**
       
     */
    class NewDatauBooNEFrameDataSource : public WCP::FrameDataSource {
	mutable TTree* tree;	// or TChain
	WCPSst::RootEvent event;

      public:
	NewDatauBooNEFrameDataSource(TTree& tree, const WCP::GeomDataSource& gds,int bins_per_frame1 = 9600);
	virtual ~NewDatauBooNEFrameDataSource();

        Int_t fPlaneNum;

	void Save();

	/// Return the number of frames this data source knows about.  Return -1 if unlimited.
	virtual int size() const;

	/// Explicitly set the "frame" (event) to process.  Frame number returned or -1 on error.
	virtual int jump(int frame_number);

        Double_t CalcRMSWithFlags(TH1F *hist);
        void ChirpFilterAlg(TH1F *hist);
        void ZigzagFilterAlg(TH1F *hist);
        void SignalFilterAlg(TH1F *hist);
        void NoisyFilterAlg(TH1F *hist,int planeNum);
        void WaveFilterAlg(TH1F **filtHists);
        void RawAdaptiveBaselineAlg(TH1F *filtHist) ;
        void RemoveFilterFlags(TH1F *filtHist);


    private:
	const WCP::GeomDataSource& gds;
	int nwire_u, nwire_v, nwire_w;

	TH1F **hu;
        TH1F **hv;
        TH1F **hw;
    };

    

}

#endif
