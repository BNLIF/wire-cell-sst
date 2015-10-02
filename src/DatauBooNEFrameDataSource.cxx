#include "WireCellSst/DatauBooNEFrameDataSource.h"

#include "TClonesArray.h"
#include "TH1F.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "WireCellData/GeomWire.h"
#include "TF1.h"
#include "TMath.h"

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
  std::cout << "Saved file" << std::endl;
}


void WireCellSst::DatauBooNEFrameDataSource::zigzag_removal(TH1F *h1, int plane, int channel_no){
  TVirtualFFT *ifft;
  double value_re[9600],value_im[9600];
  
  int n = bins_per_frame;
  int nbin = bins_per_frame;
  
  TF1 *f1 = new TF1("f1","gaus");
  double par[3];

  TH1 *hm = 0;
  TH1 *hp = 0;
  TH1 *fb = 0;
  
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
  
  //save rms for later use
  if (plane == 0 ){
    urms_map[channel_no] = rms;
  }else if (plane == 1){
    vrms_map[channel_no] = rms;
  }else if (plane == 2){
    wrms_map[channel_no] = rms;
  }
  // std::cout << plane << " " << channel_no << " " << rms << std::endl;

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
  delete f1;

}

void WireCellSst::DatauBooNEFrameDataSource::chirp_id(TH1F *hist, int plane, int channel_no){
  const int windowSize = 20;
  const double chirpMinRMS = 0.9;
  const double maxNormalNeighborFrac = 0.20;
  const int maxTicks =bins_per_frame ;

  int counter = 0;
  double ADCval;
  double runningAmpMean = 0.0;
  double runningAmpRMS = 0.0;
  int numLowRMS = 0;
  int firstLowRMSBin = -1;
  int lastLowRMSBin = -1;
  bool lowRMSFlag = false;
  double RMSfirst = 0.0;
  double RMSsecond = 0.0;
  double RMSthird = 0.0;
  int numNormalNeighbors = 0;
  int numBins = hist->GetNbinsX();
   
  for(int i = 0; i < numBins; i++)
    {
      ADCval = hist->GetBinContent(i+1);
      runningAmpMean += ADCval;
      runningAmpRMS += TMath::Power(ADCval,2.0);
      
      counter++;
      if(counter == windowSize)
	{
	  runningAmpMean /= (double)windowSize;
	  runningAmpRMS /= (double)windowSize;
	  runningAmpRMS = TMath::Sqrt(runningAmpRMS-TMath::Power(runningAmpMean,2.0));
	  
	  RMSfirst = RMSsecond;
	  RMSsecond = RMSthird;
	  RMSthird = runningAmpRMS;
	  
	  if(runningAmpRMS < chirpMinRMS)
	    {
	      numLowRMS++;
	    }
	  
	  if(i >= 3*windowSize-1)
	    {
	      if((RMSsecond < chirpMinRMS) && ((RMSfirst > chirpMinRMS) || (RMSthird > chirpMinRMS)))
		{
		  numNormalNeighbors++;
		}
	      
	      if(lowRMSFlag == false)
		{
		  if((RMSsecond < chirpMinRMS) && (RMSthird < chirpMinRMS))
		    {
		      lowRMSFlag = true;
		      firstLowRMSBin = i-2*windowSize+1;
		      lastLowRMSBin = i-windowSize+1;
		    }
		  
		  if((i == 3*windowSize-1) && (RMSfirst < chirpMinRMS) && (RMSsecond < chirpMinRMS))
		    {
		      lowRMSFlag = true;
		      firstLowRMSBin = i-3*windowSize+1;
		      lastLowRMSBin = i-2*windowSize+1;
		    }
		}
	      else
		{
		  if((RMSsecond < chirpMinRMS) && (RMSthird < chirpMinRMS))
		    {
		      lastLowRMSBin = i-windowSize+1;
		    }
		}
	    }
	  
	  counter = 0;
	  runningAmpMean = 0.0;
	  runningAmpRMS = 0.0;
	}
    }
  
  double chirpFrac = ((double) numLowRMS)/(((double) maxTicks)/((double) windowSize));
  double normalNeighborFrac = ((double) numNormalNeighbors)/((double) numLowRMS);
  
  if(((normalNeighborFrac < maxNormalNeighborFrac) || ((numLowRMS < 2.0/maxNormalNeighborFrac) && (lastLowRMSBin-firstLowRMSBin == numLowRMS*windowSize))) && (numLowRMS > 4))
    {
      firstLowRMSBin = TMath::Max(1,firstLowRMSBin-windowSize);
      lastLowRMSBin = TMath::Min(numBins,lastLowRMSBin+2*windowSize);
      
      if((numBins-lastLowRMSBin) < windowSize)
	{
	  lastLowRMSBin = numBins;
	}
      
      if(chirpFrac > 0.99)
	{
	  firstLowRMSBin = 1;
	  lastLowRMSBin = numBins;
	}
      

      for(int i = 0; i < numBins; i++)
	{
	  if((i+1 >= firstLowRMSBin) && (i+1 <= lastLowRMSBin))
	    {
	      hist->SetBinContent(i+1,10000.0);
	    }
	}

      //How to save the chirp results ... 
      
      if (plane == 0){
	//u-plane
	if (uchirp_map.find(channel_no) == uchirp_map.end()){
	  std::pair<int,int> abc(firstLowRMSBin-1, lastLowRMSBin-1);
	  uchirp_map[channel_no] = abc;
	}else{
	  if (firstLowRMSBin < uchirp_map[channel_no].first )
	    uchirp_map[channel_no].first = firstLowRMSBin-1;
	  if (lastLowRMSBin > uchirp_map[channel_no].second)
	    uchirp_map[channel_no].second = lastLowRMSBin-1;
	}
      }else if (plane == 1){
	//v-plane
	if (vchirp_map.find(channel_no) == vchirp_map.end()){
	  std::pair<int,int> abc(firstLowRMSBin-1, lastLowRMSBin-1);
	  vchirp_map[channel_no] = abc;
	}else{
	  if (firstLowRMSBin < vchirp_map[channel_no].first )
	    vchirp_map[channel_no].first = firstLowRMSBin-1;
	  if (lastLowRMSBin > vchirp_map[channel_no].second)
	    vchirp_map[channel_no].second = lastLowRMSBin-1;
	}
      }else if (plane == 2){
	//w-plane
	if (wchirp_map.find(channel_no) == wchirp_map.end()){
	  std::pair<int,int> abc(firstLowRMSBin-1, lastLowRMSBin-1);
	  wchirp_map[channel_no] = abc;
	}else{
	  if (firstLowRMSBin < wchirp_map[channel_no].first )
	    wchirp_map[channel_no].first = firstLowRMSBin-1;
	  if (lastLowRMSBin > wchirp_map[channel_no].second)
	    wchirp_map[channel_no].second = lastLowRMSBin-1;
	}
      }
    }
}

void WireCellSst::DatauBooNEFrameDataSource::SignalFilter(TH1F *hist){
  const double sigFactor = 4.0;
  const int padBins = 8;
  
  double rmsVal = CalcRMSWithFlags(hist);
  double sigThreshold = sigFactor*rmsVal;
  
  double ADCval;
  std::vector<bool> signalRegions;
  int numBins = hist->GetNbinsX();

  for(int i = 0; i < numBins; i++)
    {
      ADCval = hist->GetBinContent(i+1);
      
      if(((ADCval > sigThreshold) || (ADCval < -1.0*sigThreshold)) && (ADCval < 4096.0))
	{
	  signalRegions.push_back(true);
	}
      else
	{
	  signalRegions.push_back(false);
	}
    }
  
  for(int i = 0; i < numBins; i++)
    {
      if(signalRegions[i] == true)
	{
	  for(int j = TMath::Max(0,i-padBins); j < TMath::Min(numBins,i+padBins); j++)
	    {
	      ADCval = hist->GetBinContent(j+1);
	      if(ADCval < 4096.0)
		{
		  hist->SetBinContent(j+1,ADCval+20000.0);
		}
	    }
	}
    }  
  
}

double WireCellSst::DatauBooNEFrameDataSource::CalcRMSWithFlags(TH1F *hist){
  double ADCval;
  double theMean = 0.0;
  double theRMS = 0.0;
  int waveformSize = hist->GetNbinsX();
  int counter = 0;

  for(int i = 0; i < waveformSize; i++)
    {
      ADCval = hist->GetBinContent(i+1);
      
      if(ADCval < 4096.0)
	{
	  theMean += ADCval;
	  theRMS += TMath::Power(ADCval,2.0);
	  counter++;
	}
    }
  
  if(counter == 0)
    {
      theMean = 0.0;
      theRMS = 0.0;
    }
  else
    {
      theMean /= (double)counter;
      theRMS /= (double)counter;
    theRMS = TMath::Sqrt(theRMS-TMath::Power(theMean,2.0));
    }
  
  return theRMS;

}

void WireCellSst::DatauBooNEFrameDataSource::RemoveFilterFlags(TH1F *filtHist)
{
  double ADCval;
  int numBins = filtHist->GetNbinsX();
  for(int i = 0; i < numBins; i++)
    {
      ADCval = filtHist->GetBinContent(i+1);
      
      if(ADCval > 4096.0)
	{
	  if(ADCval > 10000.0)
	    filtHist->SetBinContent(i+1,ADCval-20000.0);
	  else
	    filtHist->SetBinContent(i+1,0.0);
	}
    }
  
  return;
}

void  WireCellSst::DatauBooNEFrameDataSource::RawAdaptiveBaselineAlg(TH1F *filtHist)
{
  const int windowSize = 20;
  
  int numBins = filtHist->GetNbinsX();
  int minWindowBins = windowSize/2;
  
  double baselineVec[numBins];
  bool isFilledVec[numBins];
  
  int numFlaggedBins = 0;
  for(int j = 0; j < numBins; j++)
    {
      if(filtHist->GetBinContent(j+1) == 10000.0)
	{
	  numFlaggedBins++;
	}
    }
  if(numFlaggedBins == numBins) return; // Eventually replace this with flag check
  
  double baselineVal = 0.0;
  int windowBins = 0;
  int index;
  double ADCval;
  for(int j = 0; j <= windowSize/2; j++)
    {
      ADCval = filtHist->GetBinContent(j+1);
      if(ADCval < 4096.0)
	{
	  baselineVal += ADCval;
	  windowBins++;
	}
    }
  
  if(windowBins == 0)
    baselineVec[0] = 0.0;
  else
    baselineVec[0] = baselineVal/((double) windowBins);
  
  if(windowBins < minWindowBins)
    isFilledVec[0] = false;
  else
    isFilledVec[0] = true;
  
  int oldIndex;
  int newIndex;
  for(int j = 1; j < numBins; j++)
    {
      oldIndex = j-windowSize/2-1;
      newIndex = j+windowSize/2;
      
      if(oldIndex >= 0)
	{
	  ADCval = filtHist->GetBinContent(oldIndex+1);
	  if(ADCval < 4096.0)
	    {
	      baselineVal -= filtHist->GetBinContent(oldIndex+1);
	      windowBins--;
	    }
	}
      if(newIndex < numBins)
	{  
	  ADCval = filtHist->GetBinContent(newIndex+1);
	  if(ADCval < 4096)
	    {
	      baselineVal += filtHist->GetBinContent(newIndex+1);
	      windowBins++;
	    }
	}
      
      if(windowBins == 0)
	baselineVec[j] = 0.0;
      else
	baselineVec[j] = baselineVal/windowBins;
      
      if(windowBins < minWindowBins)
	isFilledVec[j] = false;
      else
	isFilledVec[j] = true;
    }
  
  int downIndex;
  int upIndex;
  bool downFlag;
  bool upFlag;
  for(int j = 0; j < numBins; j++)
    {
      downFlag = false;
      upFlag = false;
      
      ADCval = filtHist->GetBinContent(j+1);
      if(ADCval != 10000.0)
	{
	  if(isFilledVec[j] == false)
	    {
	      downIndex = j;
	      while((isFilledVec[downIndex] == false) && (downIndex > 0) && (filtHist->GetBinContent(downIndex+1) != 10000.0))
		{
		  downIndex--;
		}
	      
	      if(isFilledVec[downIndex] == false)
		downFlag = true;
	      
	      upIndex = j;
	      while((isFilledVec[upIndex] == false) && (upIndex < numBins-1) && (filtHist->GetBinContent(upIndex+1) != 10000.0))
		{
		  upIndex++;
		}
	      
	      if(isFilledVec[upIndex] == false)
		upFlag = true;
	      
	      if((downFlag == false) && (upFlag == false))
		baselineVec[j] = ((j-downIndex)*baselineVec[downIndex]+(upIndex-j)*baselineVec[upIndex])/((double) upIndex-downIndex);
	      else if((downFlag == true) && (upFlag == false))
		baselineVec[j] = baselineVec[upIndex];
	      else if((downFlag == false) && (upFlag == true))
		baselineVec[j] = baselineVec[downIndex];
	      else
		baselineVec[j] = 0.0;
	    }
	  
	  filtHist->SetBinContent(j+1,ADCval-baselineVec[j]);
	}
    }
}

void WireCellSst::DatauBooNEFrameDataSource::NoisyFilterAlg(TH1F *hist, int planeNum, int channel_no)
{
  double rmsVal = CalcRMSWithFlags(hist);
  // const double maxRMSCut[3] = {10.0,10.0,5.0};
  double maxRMSCut[3] = {10.0,10.0,5.0};
  double minRMSCut[3] = {2,2,2};

  if (planeNum == 0){
    if (channel_no < 100){
      maxRMSCut[0] = 5;
      minRMSCut[0] = 1;
    }else if (channel_no >= 100 && channel_no<2000){
      maxRMSCut[0] = 10;
      minRMSCut[0] = 2;
    }else if (channel_no >= 2000 && channel_no < 2400){
      maxRMSCut[0] = 5;
      minRMSCut[0] = 1;
    }
  }else if (planeNum == 1){
    if (channel_no <290){
      maxRMSCut[1] = 5;
      minRMSCut[1] = 1;
    }else if (channel_no>=290 && channel_no < 2200){
      maxRMSCut[1] = 10;
      minRMSCut[1] = 2;
    }else if (channel_no >=2200){
      maxRMSCut[1] = 5;
      minRMSCut[1] = 1;
    }
  }else if (planeNum == 2){
    maxRMSCut[2] = 8;
    minRMSCut[2] = 2;
  }
  
  // std::cout << "Xin: " << planeNum << " " << channel_no << " " << rmsVal << std::endl;
   
  if(rmsVal > maxRMSCut[planeNum] || rmsVal < minRMSCut[planeNum])
    //if(rmsVal > maxRMSCut[planeNum] )
    //if( rmsVal < minRMSCut[planeNum])
    {
      Int_t numBins = hist->GetNbinsX();
      for(Int_t i = 0; i < numBins; i++)
	{
	  hist->SetBinContent(i+1,10000.0);
	}         
      
      if (planeNum == 0){
	//u-plane
	if (uchirp_map.find(channel_no) == uchirp_map.end()){
	  std::pair<int,int> abc(0, numBins-1);
	  uchirp_map[channel_no] = abc;
	}else{
	  uchirp_map[channel_no].first = 0;
	  uchirp_map[channel_no].second = numBins-1;
	}
      }else if (planeNum == 1){
	if (vchirp_map.find(channel_no) == vchirp_map.end()){
	  std::pair<int,int> abc(0, numBins-1);
	  vchirp_map[channel_no] = abc;
	}else{
	  vchirp_map[channel_no].first = 0;
	  vchirp_map[channel_no].second = numBins-1;
	}
      }else if (planeNum == 2){
	if (wchirp_map.find(channel_no) == wchirp_map.end()){
	  std::pair<int,int> abc(0, numBins-1);
	  wchirp_map[channel_no] = abc;
	}else{
	  wchirp_map[channel_no].first = 0;
	  wchirp_map[channel_no].second = numBins-1;
	}
      }
                 
    }
  
  return;
}


int WireCellSst::DatauBooNEFrameDataSource::jump(int frame_number)
{
 
  // return frame_number;

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

    
    std::cout << "Load Data " << std::endl;
    // load into frame
    int nchannels = event.channelid->size();
    for (size_t ind=0; ind < nchannels; ++ind) {
	TH1F* signal = dynamic_cast<TH1F*>(event.signal->At(ind));
	if (!signal) {
	    return -1;
	}

	WireCell::Trace trace;
	trace.chid = event.channelid->at(ind);

	//	trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
	//trace.charge.resize(bins_per_frame, 0.0);

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


    std::cout << "Remove ZigZag " << std::endl;
    //deal with the zig zag noise
    for (int i=0;i!=nu;i++){
      zigzag_removal(hu[i],0,i);
    }
    for (int i=0;i!=nv;i++){
      zigzag_removal(hv[i],1,i);
    }
    for (int i=0;i!=nw;i++){
      zigzag_removal(hw[i],2,i);
    }

    std::cout << "Identify Chirping" << std::endl;
    // deal with the chirping ... set chirping part > 10000
    for (int i=0;i!=nu;i++){
      chirp_id(hu[i],0,i);
    }
    for (int i=0;i!=nv;i++){
      chirp_id(hv[i],1,i);
    }
    for (int i=0;i!=nw;i++){
      chirp_id(hw[i],2,i);
    }

    
    std::cout << "Adaptive Baseline " << uchirp_map.size() << " " << 
      vchirp_map.size() << " " << wchirp_map.size() << std::endl;
    // do the adaptive baseline ... 
    for (auto it = uchirp_map.begin(); it!= uchirp_map.end(); it++){
      SignalFilter(hu[it->first]);
      RawAdaptiveBaselineAlg(hu[it->first]);
      RemoveFilterFlags(hu[it->first]);
    }
    for (auto it = vchirp_map.begin(); it!= vchirp_map.end(); it++){
      SignalFilter(hv[it->first]);
      RawAdaptiveBaselineAlg(hv[it->first]);
      RemoveFilterFlags(hv[it->first]);
    }
    for (auto it = wchirp_map.begin(); it!= wchirp_map.end(); it++){
      SignalFilter(hw[it->first]);
      RawAdaptiveBaselineAlg(hw[it->first]);
      RemoveFilterFlags(hw[it->first]);
    }
    
    std::cout << "Noisy Channel " << std::endl;
    // deal with the noisy signal, and put them into chirping map 
    for (int i=0;i!=nu;i++){
      SignalFilter(hu[i]);
      NoisyFilterAlg(hu[i],0,i);
      RemoveFilterFlags(hu[i]);
    }
    for (int i=0;i!=nv;i++){
      SignalFilter(hv[i]);
      NoisyFilterAlg(hv[i],1,i);
      RemoveFilterFlags(hv[i]);
    }
    for (int i=0;i!=nw;i++){
      SignalFilter(hw[i]);
      NoisyFilterAlg(hw[i],2,i);
      RemoveFilterFlags(hw[i]);
    }

    
    // deal with coherent noise removal 
    int n = bins_per_frame;
    int nbin = bins_per_frame;
    double par[3];
    
    WireSelectionV uplane_all;
    for (int i=0;i!=53;i++){
      WireSelection uplane;
      for (int j=0;j!=48;j++){
    	int num = i*48+j;
    	if (num < nu){
    	  uplane.push_back(num);
    	}
      }
      if (uplane.size() !=0) {
    	uplane_all.push_back(uplane);
    	for (int j=0;j!=uplane.size();j++){
    	  uplane_map[uplane.at(j)] = uplane;
    	}
      }
    }
    
    
    WireSelectionV vplane_all;
    for (int i=0;i!=53;i++){
      WireSelection vplane;
      for (int j=0;j!=48;j++){
    	int num = i*48+j;
    	if (num < nv){
    	  vplane.push_back(num);
    	}
      }
      if (vplane.size() !=0) {
    	vplane_all.push_back(vplane);
    	for (int j=0;j!=vplane.size();j++){
    	  vplane_map[vplane.at(j)] = vplane;
    	}
      }
    }
    
    WireSelectionV wplane_all;
    for (int i=0;i!=76;i++){
      WireSelection wplane;
      for (int j=0;j!=48;j++){
    	int num = i*48+j;
    	if (num < nw){
    	  wplane.push_back(num);
    	}
      }
      if (wplane.size() !=0) {
    	wplane_all.push_back(wplane);
    	for (int j=0;j!=wplane.size();j++){
    	  wplane_map[wplane.at(j)] = wplane;
    	}
      }
    }
    


    for (int i=0;i!=uplane_all.size();i++){
      //std::cout << "U " << i << " " << uplane_all.size() << std::endl;
      //calculate maximum_rms;
      float max_rms = 0;
      for (int j=0;j!=uplane_all.at(i).size();j++){
	if (urms_map[uplane_all.at(i).at(j)] > max_rms) max_rms = urms_map[uplane_all.at(i).at(j)];
      }
      // if (max_rms<10) max_rms = 10;
      // std::cout << i << " " << max_rms << std::endl;

      
      TH1F *h3 = new TH1F("h3","h3",int(12*max_rms),-6*max_rms,6*max_rms);
      if (uplane_all.at(i).size()>0){
    	for (int j=0;j!=nbin;j++){
    	  h3->Reset();
    	  for (int k=0;k!=uplane_all.at(i).size();k++){
    	    if (fabs(hu[uplane_all.at(i).at(k)]->GetBinContent(j+1))<5*max_rms && 
    		fabs(hu[uplane_all.at(i).at(k)]->GetBinContent(j+1))>0.001)
    	      h3->Fill(hu[uplane_all.at(i).at(k)]->GetBinContent(j+1));
    	  }
	  
    	  if (h3->GetSum()>0){
    	    double xq = 0.5;
    	    h3->GetQuantiles(1,&par[1],&xq);
    	  }else{
    	    par[1] = 0;
    	  }

	  // if (j == 100) std::cout << par[1] << " " << h3->GetSum() << std::endl;
	  
    	  for (int k=0;k!=uplane_all.at(i).size();k++){
	    if (fabs(hu[uplane_all.at(i).at(k)]->GetBinContent(j+1))>0.001)
	      hu[uplane_all.at(i).at(k)]->SetBinContent(j+1,hu[uplane_all.at(i).at(k)]->GetBinContent(j+1)-par[1]);
    	  }
    	}
      }
      delete h3;
    }
    
    for (int i=0;i!=vplane_all.size();i++){
      //std::cout << "V " << i << " " << vplane_all.size() << std::endl;
      float max_rms = 0;
      for (int j=0;j!=vplane_all.at(i).size();j++){
	if (vrms_map[vplane_all.at(i).at(j)] > max_rms) max_rms = vrms_map[vplane_all.at(i).at(j)];
      }
      // if (max_rms<10) max_rms = 10;
       

      TH1F *h3 = new TH1F("h3","h3",int(12*max_rms),-6*max_rms,6*max_rms);
      if (vplane_all.at(i).size()>0){
    	for (int j=0;j!=nbin;j++){
    	  h3->Reset();
    	  for (int k=0;k!=vplane_all.at(i).size();k++){
    	    if (fabs(hv[vplane_all.at(i).at(k)]->GetBinContent(j+1))<5*max_rms && 
    		fabs(hv[vplane_all.at(i).at(k)]->GetBinContent(j+1))>0.001)
    	      h3->Fill(hv[vplane_all.at(i).at(k)]->GetBinContent(j+1));
    	  }
	  
    	  if (h3->GetSum()>0){
    	    double xq = 0.5;
    	    h3->GetQuantiles(1,&par[1],&xq);
    	  }else{
    	    par[1] = 0;
    	  }
	  
    	  for (int k=0;k!=vplane_all.at(i).size();k++){
	    if (fabs(hv[vplane_all.at(i).at(k)]->GetBinContent(j+1))>0.001)
	      hv[vplane_all.at(i).at(k)]->SetBinContent(j+1,hv[vplane_all.at(i).at(k)]->GetBinContent(j+1)-par[1]);
    	  }
    	}
      }
      delete h3;
    }

    for (int i=0;i!=wplane_all.size();i++){
      float max_rms = 0;
      for (int j=0;j!=wplane_all.at(i).size();j++){
	if (wrms_map[wplane_all.at(i).at(j)] > max_rms) max_rms = wrms_map[wplane_all.at(i).at(j)];
      }
      // if (max_rms<10) max_rms = 10;
     
      //std::cout << "W " << i << " " << wplane_all.size() << std::endl;
      TH1F *h3 = new TH1F("h3","h3",int(12*max_rms),-6*max_rms,6*max_rms);
      if (wplane_all.at(i).size()>0){
    	for (int j=0;j!=nbin;j++){
    	  h3->Reset();
    	  for (int k=0;k!=wplane_all.at(i).size();k++){
    	    if (fabs(hw[wplane_all.at(i).at(k)]->GetBinContent(j+1))<5*max_rms && 
    		fabs(hw[wplane_all.at(i).at(k)]->GetBinContent(j+1))>0.001)
    	      h3->Fill(hw[wplane_all.at(i).at(k)]->GetBinContent(j+1));
    	  }
	  
    	  if (h3->GetSum()>0){
    	    double xq = 0.5;
    	    h3->GetQuantiles(1,&par[1],&xq);
    	  }else{
    	    par[1] = 0;
    	  }
	  
    	  for (int k=0;k!=wplane_all.at(i).size();k++){
	    if (fabs(hw[wplane_all.at(i).at(k)]->GetBinContent(j+1))>0.001)
	      hw[wplane_all.at(i).at(k)]->SetBinContent(j+1,hw[wplane_all.at(i).at(k)]->GetBinContent(j+1)-par[1]);
    	  }
    	}
      }
      delete h3;
    }



    double sum_u = 0, sum_v = 0, sum_w = 0;
    for (auto it = uchirp_map.begin(); it!= uchirp_map.end();it++){
      sum_u += it->second.second - it->second.first  + 1;
      //std::cout << "U: " << it->first << " " << it->second.first << " " << it->second.second << std::endl;
    }
    for (auto it = vchirp_map.begin(); it!= vchirp_map.end();it++){
      sum_v += it->second.second - it->second.first  + 1;
      //std::cout << "V: " << it->first << " " << it->second.first << " " << it->second.second << std::endl;
    }
    for (auto it = wchirp_map.begin(); it!= wchirp_map.end();it++){
      sum_w += it->second.second - it->second.first  + 1;
      //std::cout << "W: " << it->first << " " << it->second.first << " " << it->second.second << std::endl;
    }
    sum_u /= bins_per_frame * nwire_u;
    sum_v /= bins_per_frame * nwire_v;
    sum_w /= bins_per_frame * nwire_w;
    
    std::cout << "Final Count " << uchirp_map.size() << " " << sum_u << " " 
	      << vchirp_map.size() << " " << sum_v << " " 
	      << wchirp_map.size() << " " << sum_w << " " << std::endl;
    double eff_3plane, eff_2plane;
    eff_3plane = (1-sum_u)*(1-sum_v)*(1-sum_w);
    eff_2plane = eff_3plane 
      + (1-sum_u)*(1-sum_v) * sum_w
      + (1-sum_v)*(1-sum_w) * sum_u
      + (1-sum_w)*(1-sum_u) * sum_v;
    std::cout << "Efficiency: " << eff_3plane << " " << eff_2plane << std::endl;



    
    // load the stuff back to the frame ... 
    
    for (size_t ind=0; ind < nchannels; ++ind) {
      TH1F* signal;
      WireCell::Trace trace;
      trace.chid = event.channelid->at(ind);

      trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
      trace.charge.resize(bins_per_frame, 0.0);

      //int flag = 0;
      if (trace.chid < nwire_u){
	signal = hu[trace.chid];
	//	if (uplane_map.find(trace.chid)==uplane_map.end()) flag = 1;
      }else if (trace.chid < nwire_u + nwire_v){
	signal = hv[trace.chid - nwire_u];
	//if (vplane_map.find(trace.chid-nwire_u)==vplane_map.end()) flag = 1;
      }else{
	signal = hw[trace.chid - nwire_u - nwire_v];
	//if (wplane_map.find(trace.chid-nwire_u-nwire_v)==wplane_map.end()) flag = 1;
      }
      
      

      //      if (flag ==0){
      for (int ibin=0; ibin != bins_per_frame; ibin++) {
	trace.charge.at(ibin) = signal->GetBinContent(ibin+1);
      }
      //}else{
	//crazy stuff
	//flag = 1
      //	for (int ibin=0; ibin != bins_per_frame; ibin++) {
      //	  trace.charge.at(ibin) = 0.;//signal->GetBinContent(ibin+1);
      //	}
      //}
      
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
      if (rms >0.65&& rms < 5) 
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
 double sxx=0,syy=0,sxy=0;
  
  double mean1 = h1->GetSum()/h1->GetNbinsX();
  double mean2 = h2->GetSum()/h2->GetNbinsX();
  

  for (int i=0;i!=h1->GetNbinsX();i++){
    if (fabs(h1->GetBinContent(i+1)-mean1)<30 &&
	fabs(h2->GetBinContent(i+1)-mean2)<30){
      sxx += pow(h1->GetBinContent(i+1),2) - pow(mean1,2);
      syy += pow(h2->GetBinContent(i+1),2) - pow(mean2,2);
      sxy += h1->GetBinContent(i+1) * h2->GetBinContent(i+1) - mean1*mean2;
    }
  }
  
  // file->Close();
  
  double r = sxy/sqrt(sxx*syy);
  
  return r;
}
