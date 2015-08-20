#include "WireCellSst/DatauBooNEFrameDataSource.h"

#include "TClonesArray.h"
#include "TH1F.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "WireCellData/GeomWire.h"
#include "TF1.h"

using namespace WireCell;

WireCellSst::DatauBooNEFrameDataSource::DatauBooNEFrameDataSource(TTree& ttree, const WireCell::GeomDataSource& gds,int bins_per_frame1)
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

WireCellSst::DatauBooNEFrameDataSource::~DatauBooNEFrameDataSource()
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

int WireCellSst::DatauBooNEFrameDataSource::size() const
{
    return tree->GetEntries();
}

void WireCellSst::DatauBooNEFrameDataSource::Save(){
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


int WireCellSst::DatauBooNEFrameDataSource::jump(int frame_number)
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
	
	for (int ibin=0; ibin != bins_per_frame; ibin++) {
	  htemp->SetBinContent(ibin+1,signal->GetBinContent(ibin+1)-threshold);
	}	
    }

    int nu = 2400, nv = 2400, nw = 3456;
    int ntotal = nu + nv + nw;
    
    double value_re[9600],value_im[9600];
    TVirtualFFT *ifft;

    int n = bins_per_frame;
    int nbin = bins_per_frame;

    TF1 *f1 = new TF1("f1","gaus");
    Double_t par[3];

    // U-plane
    for (int i=0;i!=nu;i++){
      TH1 *hm = 0;
      TH1 *hp = 0;
      TH1 *fb = 0;
      
      TH1F *h1 = hu[i];
      
      hm = h1->FFT(0,"MAG");
      hp = h1->FFT(0,"PH");
    
      for (int j=0;j!=nbin;j++){
	double rho = hm->GetBinContent(j+1);
	double phi = hp->GetBinContent(j+1);
	
	if (j==0) rho = 0;
	
	if (j<=3500 || j> nbin-3500){
	  value_re[j] = rho*cos(phi)/nbin;
	  value_im[j] = rho*sin(phi)/nbin;
	}else{
	  value_re[j] = 0;
	  value_im[j] = 0;
	}
      }
      
      ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
      ifft->SetPointsComplex(value_re,value_im);
      ifft->Transform();
      fb = TH1::TransformHisto(ifft,0,"Re");
    
      for (int j=0;j!=nbin;j++){
	h1->SetBinContent(j+1,fb->GetBinContent(j+1));
      }
      
      
      //remove pedestal
      float mean = h1->GetSum()/h1->GetNbinsX();
      float rms = 0;
      for (int j=0;j!=nbin;j++){
	rms += pow(h1->GetBinContent(j+1)-mean,2);
      }
      rms = sqrt(rms/h1->GetNbinsX());
      TH1F h2("h2","h2",100,mean-10*rms,mean+10*rms);
      for (int j=0;j!=nbin;j++){
	h2.Fill(h1->GetBinContent(j+1));
      }
      h2.Fit(f1,"Q0","");
      f1->GetParameters(par);

      for (int j=0;j!=nbin;j++){
	h1->SetBinContent(j+1,h1->GetBinContent(j+1)-par[1]);
      }
      

      
      delete hp;
      delete hm;
      delete ifft;
      delete fb;
    }



     // V-plane
    for (int i=0;i!=nv;i++){
      TH1 *hm = 0;
      TH1 *hp = 0;
      TH1 *fb = 0;
      
      TH1F *h1 = hv[i];
      
      hm = h1->FFT(0,"MAG");
      hp = h1->FFT(0,"PH");
    
      for (int j=0;j!=nbin;j++){
	double rho = hm->GetBinContent(j+1);
	double phi = hp->GetBinContent(j+1);
	
	if (j==0) rho = 0;
	
	if (j<=3500 || j> nbin-3500){
	  value_re[j] = rho*cos(phi)/nbin;
	  value_im[j] = rho*sin(phi)/nbin;
	}else{
	  value_re[j] = 0;
	  value_im[j] = 0;
	}
      }
      
      ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
      ifft->SetPointsComplex(value_re,value_im);
      ifft->Transform();
      fb = TH1::TransformHisto(ifft,0,"Re");
    
      for (int j=0;j!=nbin;j++){
	h1->SetBinContent(j+1,fb->GetBinContent(j+1));
      }
      
      
      //remove pedestal
      float mean = h1->GetSum()/h1->GetNbinsX();
      float rms = 0;
      for (int j=0;j!=nbin;j++){
	rms += pow(h1->GetBinContent(j+1)-mean,2);
      }
      rms = sqrt(rms/h1->GetNbinsX());
      TH1F h2("h2","h2",100,mean-10*rms,mean+10*rms);
      for (int j=0;j!=nbin;j++){
	h2.Fill(h1->GetBinContent(j+1));
      }
      h2.Fit(f1,"Q0","");
      f1->GetParameters(par);

      for (int j=0;j!=nbin;j++){
	h1->SetBinContent(j+1,h1->GetBinContent(j+1)-par[1]);
      }
      

      
      delete hp;
      delete hm;
      delete ifft;
      delete fb;
    }




     // W-plane
    for (int i=0;i!=nw;i++){
      TH1 *hm = 0;
      TH1 *hp = 0;
      TH1 *fb = 0;
      
      TH1F *h1 = hw[i];
      
      hm = h1->FFT(0,"MAG");
      hp = h1->FFT(0,"PH");
    
      for (int j=0;j!=nbin;j++){
	double rho = hm->GetBinContent(j+1);
	double phi = hp->GetBinContent(j+1);
	
	if (j==0) rho = 0;
	
	if (j<=3500 || j> nbin-3500){
	  value_re[j] = rho*cos(phi)/nbin;
	  value_im[j] = rho*sin(phi)/nbin;
	}else{
	  value_re[j] = 0;
	  value_im[j] = 0;
	}
      }
      
      ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
      ifft->SetPointsComplex(value_re,value_im);
      ifft->Transform();
      fb = TH1::TransformHisto(ifft,0,"Re");
    
      for (int j=0;j!=nbin;j++){
	h1->SetBinContent(j+1,fb->GetBinContent(j+1));
      }
      
      
      //remove pedestal
      float mean = h1->GetSum()/h1->GetNbinsX();
      float rms = 0;
      for (int j=0;j!=nbin;j++){
	rms += pow(h1->GetBinContent(j+1)-mean,2);
      }
      rms = sqrt(rms/h1->GetNbinsX());
      TH1F h2("h2","h2",100,mean-10*rms,mean+10*rms);
      for (int j=0;j!=nbin;j++){
	h2.Fill(h1->GetBinContent(j+1));
      }
      h2.Fit(f1,"Q0","");
      f1->GetParameters(par);

      for (int j=0;j!=nbin;j++){
	h1->SetBinContent(j+1,h1->GetBinContent(j+1)-par[1]);
      }
      

      
      delete hp;
      delete hm;
      delete ifft;
      delete fb;
    }

    
    // start to remove low frequency noise ... 
    std::vector<int> used_num;
    std::vector<int> bad_num;
    int plane;

    //deal with u plane first
    plane = 0;
    used_num.clear();
    bad_num.clear();
    for (int i=0;i!=nu;i++){
      used_num.push_back(0);
      bad_num.push_back(0);
    }
    //check RMS
    for (int i=0;i!=nu;i++){
      float mean = hu[i]->GetSum()/hu[i]->GetNbinsX();
      float rms = 0;
      int nbin1 = 0;
      for (int j=0;j!=nbin;j++){
	if (fabs(hu[i]->GetBinContent(j+1)-mean)<50){ 
	  rms += pow(hu[i]->GetBinContent(j+1)-mean,2);
	  nbin1 ++;
	}
      }
      rms = sqrt(rms/nbin1);
      if (chirp_check(rms, plane, i)){
	bad_num.at(i) = 1;
      }
      //      std::cout << i << " " << nbin << " " << mean << " " << rms << std::endl;
    }
    
    WireSelectionV uplane_all;
   
    
    
    for (int i=0;i!=nu;i++){
      WireSelection uplane;
      if (used_num.at(i)==0 && bad_num.at(i)==0){
	//std::cout << i << std::endl;
	used_num.at(i)=1;
	uplane.push_back(i);
	for (int j=i;j!=nu;j++){
	  if (used_num.at(j)==0 && bad_num.at(j)==0){
	    double corr = correlation1(hu[i],hu[j]);
	    if (corr > 0.2) {
	      used_num.at(j)=1;
	      uplane.push_back(j);
	    }
	  }
	}
	uplane_all.push_back(uplane);
	for (int j=0;j!=uplane.size();j++){
	  uplane_map[uplane.at(j)] = uplane;
	}
      }
    }



    //deal with v plane first
    plane = 1;
    used_num.clear();
    bad_num.clear();
    for (int i=0;i!=nv;i++){
      used_num.push_back(0);
      bad_num.push_back(0);
    }
    //check RMS
    for (int i=0;i!=nv;i++){
      float mean = hv[i]->GetSum()/hv[i]->GetNbinsX();
      float rms = 0;
      int nbin1 = 0;
      for (int j=0;j!=nbin;j++){
	if (fabs(hv[i]->GetBinContent(j+1)-mean)<50){ 
	  rms += pow(hv[i]->GetBinContent(j+1)-mean,2);
	  nbin1 ++;
	}
      }
      rms = sqrt(rms/nbin1);
      if (chirp_check(rms, plane, i)){
	bad_num.at(i) = 1;
      }
      //  cout << rms << endl;
    }
    
    WireSelectionV vplane_all;
    
    
    
    for (int i=0;i!=nv;i++){
      WireSelection vplane;
      if (used_num.at(i)==0 && bad_num.at(i)==0){
	//std::cout << i << std::endl;
	used_num.at(i)=1;
	vplane.push_back(i);
	for (int j=i;j!=nv;j++){
	  if (used_num.at(j)==0 && bad_num.at(j)==0){
	    double corr = correlation1(hv[i],hv[j]);
	    if (corr > 0.2) {
	      used_num.at(j)=1;
	      vplane.push_back(j);
	    }
	  }
	}
	vplane_all.push_back(vplane);
	for (int j=0;j!=vplane.size();j++){
	  vplane_map[vplane.at(j)] = vplane;
	}
      }
      
      
    }
    

    //deal with w plane first
    plane = 2;
    used_num.clear();
    bad_num.clear();
    for (int i=0;i!=nw;i++){
      used_num.push_back(0);
      bad_num.push_back(0);
    }
    //check RMS
    for (int i=0;i!=nw;i++){
      float mean = hw[i]->GetSum()/hw[i]->GetNbinsX();
      float rms = 0;
      int nbin1 = 0;
      for (int j=0;j!=nbin;j++){
	if (fabs(hw[i]->GetBinContent(j+1)-mean)<50){ 
	  rms += pow(hw[i]->GetBinContent(j+1)-mean,2);
	  nbin1 ++;
	}
      }
      rms = sqrt(rms/nbin1);
      if (chirp_check(rms, plane, i)){
	bad_num.at(i) = 1;
      }
      //  cout << rms << endl;
    }
    
    WireSelectionV wplane_all;
    
    
    
    for (int i=0;i!=nw;i++){
      WireSelection wplane;
      if (used_num.at(i)==0 && bad_num.at(i)==0){
	//std::cout << i << std::endl;
	used_num.at(i)=1;
	wplane.push_back(i);
	for (int j=i;j!=nw;j++){
	  if (used_num.at(j)==0 && bad_num.at(j)==0){
	    double corr = correlation1(hw[i],hw[j]);
	    if (corr > 0.2) {
	      used_num.at(j)=1;
	      wplane.push_back(j);
	    }
	  }
	}
	wplane_all.push_back(wplane);
	for (int j=0;j!=wplane.size();j++){
	  wplane_map[wplane.at(j)] = wplane;
	}
      }
      
      
    }


    // int test = 0;
    // for (int i=0;i!=uplane_all.size();i++){
    //   test += uplane_all.at(i).size();
    // }
    // std::cout << test << " " << uplane_map.size() << std::endl;



    for (int i=0;i!=uplane_all.size();i++){
      std::cout << "U " << i << " " << uplane_all.size() << std::endl;
      TH1F *h3 = new TH1F("h3","h3",100,-50,50);
      if (uplane_all.at(i).size()>30){
	for (int j=0;j!=nbin;j++){
	  h3->Reset();
	  for (int k=0;k!=uplane_all.at(i).size();k++){
	    if (fabs(hu[uplane_all.at(i).at(k)]->GetBinContent(j+1))<50)
	      h3->Fill(hu[uplane_all.at(i).at(k)]->GetBinContent(j+1));
	  }
	  // if (h3->GetSum()>20){
	  //   par[1] = 0;
	  //   par[2] = 10;
	  //   par[0] = 10;
	  //   f1->SetParameters(par);
	  //   h3->Fit("f1","Q0","");
	  //   f1->GetParameters(par);
	  //   if (fabs(par[1])>50) par[1] = h3->GetMean();
	  // }else{
	  //   par[1] = h3->GetMean();
	  // }
	  
	  if (h3->GetSum()>10){
	    Double_t xq = 0.5;
	    h3->GetQuantiles(1,&par[1],&xq);
	  }else{
	    par[1] = 0;
	  }

	  for (int k=0;k!=uplane_all.at(i).size();k++){
	    hu[uplane_all.at(i).at(k)]->SetBinContent(j+1,hu[uplane_all.at(i).at(k)]->GetBinContent(j+1)-par[1]);
	  }
	}
      }
      delete h3;
    }
    
    for (int i=0;i!=vplane_all.size();i++){
      std::cout << "V " << i << " " << vplane_all.size() << std::endl;
      TH1F *h3 = new TH1F("h3","h3",100,-50,50);
      if (vplane_all.at(i).size()>30){
	for (int j=0;j!=nbin;j++){
	  h3->Reset();
	  for (int k=0;k!=vplane_all.at(i).size();k++){
	    if (fabs(hv[vplane_all.at(i).at(k)]->GetBinContent(j+1))<50)
	      h3->Fill(hv[vplane_all.at(i).at(k)]->GetBinContent(j+1));
	  }
	  // if (h3->GetSum()>20){
	  //   par[1] = 0;
	  //   par[2] = 10;
	  //   par[0] = 10;
	  //   f1->SetParameters(par);
	  //   h3->Fit("f1","Q0","");
	  //   f1->GetParameters(par);
	  //   if (fabs(par[1])>50) par[1] = h3->GetMean();
	  // }else{
	  //   par[1] = h3->GetMean();
	  // }
	  
	  if (h3->GetSum()>10){
	    Double_t xq = 0.5;
	    h3->GetQuantiles(1,&par[1],&xq);
	  }else{
	    par[1] = 0;
	  }
	  
	  for (int k=0;k!=vplane_all.at(i).size();k++){
	    hv[vplane_all.at(i).at(k)]->SetBinContent(j+1,hv[vplane_all.at(i).at(k)]->GetBinContent(j+1)-par[1]);
	  }
	}
      }
      delete h3;
    }
    
    
    for (int i=0;i!=wplane_all.size();i++){
      std::cout << "W " << i << " " << wplane_all.size() << std::endl;
      TH1F *h3 = new TH1F("h3","h3",100,-50,50);
      if (wplane_all.at(i).size()>30){
	for (int j=0;j!=nbin;j++){
	  h3->Reset();
	  for (int k=0;k!=wplane_all.at(i).size();k++){
	    if (fabs(hw[wplane_all.at(i).at(k)]->GetBinContent(j+1))<50)
	      h3->Fill(hw[wplane_all.at(i).at(k)]->GetBinContent(j+1));
	  }
	  // if (h3->GetSum()>20){
	  //   par[1] = 0;
	  //   par[2] = 10;
	  //   par[0] = 10;
	  //   f1->SetParameters(par);
	  //   h3->Fit("f1","Q0","");
	  //   f1->GetParameters(par);
	  //   if (fabs(par[1])>50) par[1] = h3->GetMean();
	  // }else{
	  //   par[1] = h3->GetMean();
	  // }
	  
	  if (h3->GetSum()>10){
	    Double_t xq = 0.5;
	    h3->GetQuantiles(1,&par[1],&xq);
	  }else{
	    par[1] = 0;
	  }
	  
	  for (int k=0;k!=wplane_all.at(i).size();k++){
	    hw[wplane_all.at(i).at(k)]->SetBinContent(j+1,hw[wplane_all.at(i).at(k)]->GetBinContent(j+1)-par[1]);
	  }
	}
      }
      delete h3;
    }



    
    for (size_t ind=0; ind < nchannels; ++ind) {
      TH1F* signal;
      WireCell::Trace trace;
      trace.chid = event.channelid->at(ind);

      trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
      trace.charge.resize(bins_per_frame, 0.0);

      int flag = 0;
      if (trace.chid < nwire_u){
	signal = hu[trace.chid];
	if (uplane_map.find(trace.chid)==uplane_map.end()) flag = 1;
      }else if (trace.chid < nwire_u + nwire_v){
	signal = hv[trace.chid - nwire_u];
	if (vplane_map.find(trace.chid-nwire_u)==vplane_map.end()) flag = 1;
      }else{
	signal = hw[trace.chid - nwire_u - nwire_v];
	if (wplane_map.find(trace.chid-nwire_u-nwire_v)==wplane_map.end()) flag = 1;
      }
      
      

      if (flag ==0){
	for (int ibin=0; ibin != bins_per_frame; ibin++) {
	  trace.charge.at(ibin) = signal->GetBinContent(ibin+1);
	}
      }else{
	//crazy stuff
	//flag = 1
	for (int ibin=0; ibin != bins_per_frame; ibin++) {
	  trace.charge.at(ibin) = 0.;//signal->GetBinContent(ibin+1);
	}
      }
      
      frame.traces.push_back(trace);
    }

    frame.index = frame_number;
    return frame.index;
}


bool WireCellSst::DatauBooNEFrameDataSource::chirp_check(double rms, int plane, int channel){
  if (plane == 0){
    if (channel < 100){
      if (rms >1 && rms < 5)
	return false;      
    }else if (channel >= 100 && channel<2000){
      if (rms > 2 && rms < 10)
	return false;
    }else if (channel >= 2000 && channel < 2400){
      if (rms >1.5&& rms < 5) 
	return false;
    }
  }else if (plane == 1){
    if (channel <290){
      if (rms > 1 && rms < 5)
	return false;
    }else if (channel>=290 && channel < 2200){
      if (rms > 2 && rms < 10)
	return false;
    }else if (channel >=2200){
      if (rms > 1 && rms < 5)
	return false;
    }
  }else if (plane == 2){
    if (rms > 2 && rms < 8)
      return false;
  }
  return true;
}


double WireCellSst::DatauBooNEFrameDataSource::correlation1(TH1F *h1, TH1F *h2){
 Double_t sxx=0,syy=0,sxy=0;
  
  Double_t mean1 = h1->GetSum()/h1->GetNbinsX();
  Double_t mean2 = h2->GetSum()/h2->GetNbinsX();
  

  for (int i=0;i!=h1->GetNbinsX();i++){
    if (fabs(h1->GetBinContent(i+1)-mean1)<30 &&
	fabs(h2->GetBinContent(i+1)-mean2)<30){
      sxx += pow(h1->GetBinContent(i+1),2) - pow(mean1,2);
      syy += pow(h2->GetBinContent(i+1),2) - pow(mean2,2);
      sxy += h1->GetBinContent(i+1) * h2->GetBinContent(i+1) - mean1*mean2;
    }
  }
  
  // file->Close();
  
  Double_t r = sxy/sqrt(sxx*syy);
  
  return r;
}
