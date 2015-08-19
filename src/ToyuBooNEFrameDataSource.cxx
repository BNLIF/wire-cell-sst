#include "WireCellSst/ToyuBooNEFrameDataSource.h"

#include "TClonesArray.h"
#include "TH1F.h"
#include "TFile.h"
#include "WireCellData/GeomWire.h"

using namespace WireCell;

WireCellSst::ToyuBooNEFrameDataSource::ToyuBooNEFrameDataSource(TTree& ttree, const WireCell::GeomDataSource& gds,int bins_per_frame1)
    : WireCell::FrameDataSource()
    , tree(&ttree)
    , gds(gds)
    , event()
{
  bins_per_frame = bins_per_frame1;
    // sigh, we can't do things this simply because the ttree does not
    // have a single branch.  
    // tree->SetBranchAddress(name, &event);

    tree->SetBranchAddress("eventNo" , &event.number);
    tree->SetBranchAddress("runNo"   , &event.run);
    tree->SetBranchAddress("subRunNo", &event.subrun);

    // tree->SetBranchAddress("calib_nChannel", &event.nchannels);
    // tree->SetBranchAddress("calib_channelId", &event.channelid);
    // tree->SetBranchAddress("calib_wf", &event.signal);

    tree->SetBranchAddress("raw_nChannel", &event.nchannels);
    tree->SetBranchAddress("raw_channelId", &event.channelid);
    tree->SetBranchAddress("raw_wf", &event.signal);

    GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
    GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
    GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));
    
    nwire_u = wires_u.size();
    nwire_v = wires_v.size();
    nwire_w = wires_w.size();
    
    hu = new TH1F*[nwire_u];
    hv = new TH1F*[nwire_v];
    hw = new TH1F*[nwire_w];
    
    for (int i=0;i!=nwire_u;i++){
      hu[i] = new TH1F(Form("U2_%d",i),Form("U2_%d",i),bins_per_frame,0,bins_per_frame);
    }
    for (int i=0;i!=nwire_v;i++){
      hv[i] = new TH1F(Form("V2_%d",i),Form("V2_%d",i),bins_per_frame,0,bins_per_frame);
    }
    for (int i=0;i!=nwire_w;i++){
      hw[i] = new TH1F(Form("W2_%d",i),Form("W2_%d",i),bins_per_frame,0,bins_per_frame);
    }

}

WireCellSst::ToyuBooNEFrameDataSource::~ToyuBooNEFrameDataSource()
{
  for (int i=0;i!=nwire_u;i++){
    delete hu[i] ;
  }
  delete hu;
  for (int i=0;i!=nwire_v;i++){
    delete hv[i] ;
  }
  delete hv;
  for (int i=0;i!=nwire_w;i++){
    delete hw[i] ;
  }
  delete hw;
}

int WireCellSst::ToyuBooNEFrameDataSource::size() const
{
    return tree->GetEntries();
}

void WireCellSst::ToyuBooNEFrameDataSource::Save(){
  TFile *file = new TFile("temp_data.root","RECREATE");
  for (int i=0;i!=nwire_u;i++){
    TH1F *huu = (TH1F*)hu[i]->Clone(Form("U1_%d",i));
  }
  for (int i=0;i!=nwire_v;i++){
    TH1F *hvv = (TH1F*)hv[i]->Clone(Form("V1_%d",i));
  }
  for (int i=0;i!=nwire_w;i++){
    TH1F *hww = (TH1F*)hw[i]->Clone(Form("W1_%d",i));
  }
  file->Write();
  file->Close();
}


int WireCellSst::ToyuBooNEFrameDataSource::jump(int frame_number)
{
 


    if (frame.index == frame_number) {
	return frame_number;
    }

    frame.clear();		// win or lose, we start anew

    if (frame_number < 0) {	// underflow
	return frame_number;
    }

    int siz = tree->GetEntry(frame_number);
    if (siz <= 0 ) {
	return -1;
    }

    if (frame_number >= siz) {
	return -1;
    }

    for (int i=0;i!=nwire_u;i++){
      hu[i]->Reset();
    }
    for (int i=0;i!=nwire_v;i++){
      hv[i]->Reset();
    }
    for (int i=0;i!=nwire_w;i++){
      hw[i]->Reset();
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
	trace.charge.resize(bins_per_frame, 0.0);
	// for (int ibin=1; ibin <= signal->GetNbinsX()/4.; ++ibin) {
	//   trace.charge.push_back(signal->GetBinContent(4*(ibin-1)+1)+
	// 			 signal->GetBinContent(4*(ibin-1)+2)+
	// 			 signal->GetBinContent(4*(ibin-1)+3)+
	// 			 signal->GetBinContent(4*(ibin-1)+4)
	// 			   );
	// }
	//std::cout << signal->GetNbinsX() << std::endl;

	TH1F *htemp;
	float threshold;
	if (trace.chid < nwire_u){
	  htemp = hu[trace.chid];
	  threshold = 2048;
	}else if (trace.chid < nwire_u + nwire_v){
	  htemp = hv[trace.chid - nwire_u];
	  threshold = 2048;
	}else{
	  htemp = hw[trace.chid - nwire_u - nwire_v];
	  threshold = 400;
	}
	
	//std::cout << signal->GetNbinsX() << std::endl;
	
	if (bins_per_frame == 9600){
	  for (int ibin=0; ibin != signal->GetNbinsX()-1; ibin++) {
	    trace.charge.at(ibin)=(signal->GetBinContent(ibin+1)-threshold);
	    htemp->SetBinContent(ibin+1,signal->GetBinContent(ibin+1)-threshold);
	  }
	  htemp->SetBinContent(signal->GetNbinsX()-1,trace.charge.at(signal->GetNbinsX()-2));
	  trace.charge.at(signal->GetNbinsX()-1) = trace.charge.at(signal->GetNbinsX()-2);
	}else{
	  for (int ibin=0; ibin != bins_per_frame; ibin++) {
	    trace.charge.at(ibin)=(signal->GetBinContent(ibin+1)-threshold);
	    htemp->SetBinContent(ibin+1,signal->GetBinContent(ibin+1)-threshold);
	  }	
	}


	frame.traces.push_back(trace);
    }
    frame.index = frame_number;
    return frame.index;
}


