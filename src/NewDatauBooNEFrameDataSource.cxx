#include "WCPSst/NewDatauBooNEFrameDataSource.h"

#include "TClonesArray.h"
#include "TH1F.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "WCPData/GeomWire.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TMath.h"
#include <vector>

using namespace WCP;

const Bool_t runChirpFilterAlg = true;
const Bool_t runZigzagFilterAlg = true;
const Bool_t runSignalFilterAlg = true;
const Bool_t runNoisyFilterAlg = true;
const Bool_t runWaveFilterAlg = true;
const Bool_t runRawAdaptiveBaselineAlg = true;
const Bool_t recoverChirpingWaveforms = true;

const Double_t maxRMSCut[3] = {10.0,10.0,5.0};
const Int_t waveNoiseGroupNum = 48;

const Int_t maxTicks = 9594;

WCPSst::NewDatauBooNEFrameDataSource::NewDatauBooNEFrameDataSource(TTree& ttree, const WCP::GeomDataSource& gds,int bins_per_frame1)
    : WCP::FrameDataSource()
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

WCPSst::NewDatauBooNEFrameDataSource::~NewDatauBooNEFrameDataSource()
{
  for (int i=0;i!=nwire_u;i++){
    delete hu[i];
  }
  delete hu;
  for (int i=0;i!=nwire_v;i++){
    delete hv[i];
  }
  delete hv;
  for (int i=0;i!=nwire_w;i++){
    delete hw[i];
  }
  delete hw;
}

int WCPSst::NewDatauBooNEFrameDataSource::size() const
{
    return tree->GetEntries();
}

void WCPSst::NewDatauBooNEFrameDataSource::Save(){
  TFile *file = new TFile("temp_data_new.root","RECREATE");
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

int WCPSst::NewDatauBooNEFrameDataSource::jump(int frame_number)
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

    WCP::Trace trace;
    trace.chid = event.channelid->at(ind);
    
    //trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
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
      htemp->SetBinContent(ibin+1,signal->GetBinContent(ibin+1));
    }	
  }
   
  //  return frame_number; 

  // Setup histogram set for correlated wave noise filtering
  TH1F** waveNoiseHistsU = new TH1F * [waveNoiseGroupNum];
  TH1F** waveNoiseHistsV = new TH1F * [waveNoiseGroupNum];
  TH1F** waveNoiseHistsW = new TH1F * [waveNoiseGroupNum];

  for(Int_t i = 0; i < waveNoiseGroupNum; i++)
    {
      waveNoiseHistsU[i] = new TH1F(Form("waveNoiseHistU_number%d",i),";Tick;ADC Value",maxTicks,-0.5,maxTicks-0.5);
      waveNoiseHistsV[i] = new TH1F(Form("waveNoiseHistV_number%d",i),";Tick;ADC Value",maxTicks,-0.5,maxTicks-0.5);
      waveNoiseHistsW[i] = new TH1F(Form("waveNoiseHistW_number%d",i),";Tick;ADC Value",maxTicks,-0.5,maxTicks-0.5);
    }
  
  // Initialize helper variables
  Int_t channelNumU = -1;
  Int_t channelNumV = -1;
  Int_t channelNumW = -1;
  Int_t waveNoiseCounterU = -1;
  Int_t waveNoiseCounterV = -1;
  Int_t waveNoiseCounterW = -1;

  Int_t nu = nwire_u;
  Int_t nv = nwire_v;
  Int_t nw = nwire_w;
  
  ////Start with U-Plane
  for (Int_t i=0; i!=nu; i++){
  
    TH1F *h1 = hu[i];

    // Increment counters and do setup
    channelNumU++;
    waveNoiseCounterU++;
    waveNoiseHistsU[waveNoiseCounterU]->Reset();

    // Filter out chirping part of waveforms (and ASIC's in bad state) in time domain
    if(runChirpFilterAlg == true)
      {
  	ChirpFilterAlg(h1);
      }

    // Filter out zig-zag noise in frequency domain
    if(runZigzagFilterAlg == true)
      {
  	ZigzagFilterAlg(h1);
      }
    
    // Filter out signal regions in time domain
    if(runSignalFilterAlg == true)
      {
  	SignalFilterAlg(h1);
      }

    // Filter out noisy wires in time domain (post-filter RMS cut on non-chirping, non-signal part of waveform)
    if(runNoisyFilterAlg == true)
      {
  	NoisyFilterAlg(h1,0);
      }

    // Check if full wire group processed
    for(Int_t i=0; i < h1->GetNbinsX(); i++)
      {
  	waveNoiseHistsU[waveNoiseCounterU]->SetBinContent(i+1, h1->GetBinContent(i+1));
      }
    
    if(waveNoiseCounterU == waveNoiseGroupNum-1)
      {
    	// Filter out correlated wave noise in time domain using median correction across wires in group
    	if(runWaveFilterAlg == true)
    	  {
    	    WaveFilterAlg(waveNoiseHistsU);
    	  }

    	// Loop over wires in wire group
    	for(Int_t k = 0; k < waveNoiseGroupNum; k++)
    	  {
    	    // Run adaptive baseline method to smooth chirping waveforms near transition region and further improve S/N
    	    if(runRawAdaptiveBaselineAlg == true)
    	      {
    		RawAdaptiveBaselineAlg(waveNoiseHistsU[k]);
    	      }
	    
    	    // Remove temporary flags for chirping bins, noisy bins, and signal bins
    	    RemoveFilterFlags(waveNoiseHistsU[k]);

  	    Int_t n_samp = waveNoiseHistsU[k]->GetNbinsX();

  	    hu[channelNumU-waveNoiseGroupNum+1+k]->Reset();
	    
    	    for(Int_t i = 0; i < n_samp; i++)
    	      {
    	    	hu[channelNumU-waveNoiseGroupNum+1+k]->SetBinContent(i+1,waveNoiseHistsU[k]->GetBinContent(i+1));
  	      }
    	  }	
    	waveNoiseCounterU = -1;
      }
  }
    

  ////Second, V-Plane
  
  for (Int_t i=0; i!=nv; i++){
    
    // Increment counters and do setup
    channelNumV++;
    waveNoiseCounterV++;
    waveNoiseHistsV[waveNoiseCounterV]->Reset();
    
    TH1F *h1 = hv[i];

    // Filter out chirping part of waveforms (and ASIC's in bad state) in time domain
    if(runChirpFilterAlg == true)
      {
  	ChirpFilterAlg(h1);
      }
    
    // Filter out zig-zag noise in frequency domain
    if(runZigzagFilterAlg == true)
      {
  	ZigzagFilterAlg(h1);
      }

    // Filter out signal regions in time domain
    if(runSignalFilterAlg == true)
      {
  	SignalFilterAlg(h1);
      }

    // Filter out noisy wires in time domain (post-filter RMS cut on non-chirping, non-signal part of waveform)
    if(runNoisyFilterAlg == true)
      {
  	NoisyFilterAlg(h1,1);
      }

    // Check if full wire group processed
    for(Int_t i=0; i < h1->GetNbinsX(); i++)
      {
  	waveNoiseHistsV[waveNoiseCounterV]->SetBinContent(i+1, h1->GetBinContent(i+1));
      }
    
    if(waveNoiseCounterV == waveNoiseGroupNum-1)
      {
  	// Filter out correlated wave noise in time domain using median correction across wires in group
  	if(runWaveFilterAlg == true)
  	  {
  	    WaveFilterAlg(waveNoiseHistsV);
  	  }

  	// Loop over wires in wire group
  	for(Int_t k = 0; k < waveNoiseGroupNum; k++)
  	  {
  	    // Run adaptive baseline method to smooth chirping waveforms near transition region and further improve S/N
  	    if(runRawAdaptiveBaselineAlg == true)
  	      {
  		RawAdaptiveBaselineAlg(waveNoiseHistsV[k]);
  	      }
	    
  	    // Remove temporary flags for chirping bins, noisy bins, and signal bins
  	    RemoveFilterFlags(waveNoiseHistsV[k]);

	    Int_t n_samp = waveNoiseHistsV[k]->GetNbinsX();
	    
  	    hv[channelNumV-waveNoiseGroupNum+1+k]->Reset();

  	    // Fill filtered waveforms to be saved to file
  	    for(Int_t i = 0; i < n_samp; i++)
  	      {
  		hv[channelNumV-waveNoiseGroupNum+1+k]->SetBinContent(i+1,waveNoiseHistsV[k]->GetBinContent(i+1));
  	      }
  	  }	
  	waveNoiseCounterV = -1;
      }
  }


  ////Finally, W-Plane
  for (int i=0;i!=nw;i++){
    
    // Increment counters and do setup
    channelNumW++;
    waveNoiseCounterW++;
    waveNoiseHistsW[waveNoiseCounterW]->Reset();

    TH1F *h1 = hw[i];
    
    // Filter out chirping part of waveforms (and ASIC's in bad state) in time domain
    if(runChirpFilterAlg == true)
      {
  	ChirpFilterAlg(h1);
      }
    
    // Filter out zig-zag noise in frequency domain
    if(runZigzagFilterAlg == true)
      {
  	ZigzagFilterAlg(h1);
      }

    // Filter out signal regions in time domain
    if(runSignalFilterAlg == true)
      {
  	SignalFilterAlg(h1);
      }

    // Filter out noisy wires in time domain (post-filter RMS cut on non-chirping, non-signal part of waveform)
    if(runNoisyFilterAlg == true)
      {
  	NoisyFilterAlg(h1,2);
      }

    // Check if full wire group processed
    for(Int_t i=0; i < h1->GetNbinsX(); i++)
      {
  	waveNoiseHistsW[waveNoiseCounterW]->SetBinContent(i+1, h1->GetBinContent(i+1));
      }
    
    if(waveNoiseCounterW == waveNoiseGroupNum-1)
      {
  	// Filter out correlated wave noise in time domain using median correction across wires in group
  	if(runWaveFilterAlg == true)
  	  {
  	    WaveFilterAlg(waveNoiseHistsW);
  	  }

  	// Loop over wires in wire group
  	for(Int_t k = 0; k < waveNoiseGroupNum; k++)
  	  {
  	    // Run adaptive baseline method to smooth chirping waveforms near transition region and further improve S/N
  	    if(runRawAdaptiveBaselineAlg == true)
  	      {
  		RawAdaptiveBaselineAlg(waveNoiseHistsW[k]);
  	      }
	    
  	    // Remove temporary flags for chirping bins, noisy bins, and signal bins
  	    RemoveFilterFlags(waveNoiseHistsW[k]);

  	    hw[channelNumW-waveNoiseGroupNum+1+k]->Reset();
	    
  	    // Fill filtered waveforms to be saved to file
  	    for(Int_t i = 0; i < h1->GetNbinsX(); i++)
  	      {
  		hw[channelNumW-waveNoiseGroupNum+1+k]->SetBinContent(i+1,waveNoiseHistsW[k]->GetBinContent(i+1));
  	      }
  	  }
	
  	waveNoiseCounterW = -1;
      }
  }
 

  for(Int_t ind=0; ind < nchannels; ++ind)
    {
      TH1F *signal;
      WCP::Trace trace;
      trace.chid = event.channelid->at(ind);
      
      trace.tbin = 0;  //full readout, if zero suppress this would be non-zero
      
      trace.charge.resize(bins_per_frame, 0.0);
      
      if(trace.chid < nwire_u)
  	{
  	  signal = hu[trace.chid];
  	}else if(trace.chid < nwire_u + nwire_v)  
  	{
  	  signal = hv[trace.chid - nwire_u];
  	}else
  	{
  	  signal = hw[trace.chid - nwire_u - nwire_v];
  	}
      
      for(Int_t ibin=0; ibin!=bins_per_frame; ibin++ )
  	{
  	  trace.charge.at(ibin) = signal->GetBinContent(ibin+1);
  	}
      frame.traces.push_back(trace);
    }
  
  // Free memory associated with temporary histograms used in filtering out correlated wave noise
  for(Int_t i = 0; i < waveNoiseGroupNum; i++)
    {
      delete waveNoiseHistsU[i];
      delete waveNoiseHistsV[i];
      delete waveNoiseHistsW[i];
    }
  delete[] waveNoiseHistsU;
  delete[] waveNoiseHistsV;
  delete[] waveNoiseHistsW;

  frame.index = frame_number;
  return frame.index;
}


Double_t WCPSst::NewDatauBooNEFrameDataSource::CalcRMSWithFlags(TH1F *hist)
{
  Double_t ADCval;
  Double_t theMean = 0.0;
  Double_t theRMS = 0.0;
  Int_t waveformSize = hist->GetNbinsX();
  Int_t counter = 0;

  for(Int_t i = 0; i < waveformSize; i++)
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
      theMean /= (Double_t)counter;
      theRMS /= (Double_t)counter;
    theRMS = TMath::Sqrt(theRMS-TMath::Power(theMean,2.0));
    }
  
  return theRMS;
}


void WCPSst::NewDatauBooNEFrameDataSource::ChirpFilterAlg(TH1F *hist)
{
  const Int_t windowSize = 20;
  const Double_t chirpMinRMS = 0.9;
  const Double_t maxNormalNeighborFrac = 0.20;
  
  Int_t counter = 0;
  Double_t ADCval;
  Double_t runningAmpMean = 0.0;
  Double_t runningAmpRMS = 0.0;
  Int_t numLowRMS = 0;
  Int_t firstLowRMSBin = -1;
  Int_t lastLowRMSBin = -1;
  Bool_t lowRMSFlag = false;
  Double_t RMSfirst = 0.0;
  Double_t RMSsecond = 0.0;
  Double_t RMSthird = 0.0;
  Int_t numNormalNeighbors = 0;
  Int_t numBins = hist->GetNbinsX();
   
  for(Int_t i = 0; i < numBins; i++)
    {
      ADCval = hist->GetBinContent(i+1);
      runningAmpMean += ADCval;
      runningAmpRMS += TMath::Power(ADCval,2.0);
      
      counter++;
      if(counter == windowSize)
	{
	  runningAmpMean /= (Double_t)windowSize;
	  runningAmpRMS /= (Double_t)windowSize;
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
  
  Double_t chirpFrac = ((Double_t) numLowRMS)/(((Double_t) maxTicks)/((Double_t) windowSize));
  Double_t normalNeighborFrac = ((Double_t) numNormalNeighbors)/((Double_t) numLowRMS);
  
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
      
      if(recoverChirpingWaveforms == true)
	{
	  for(Int_t i = 0; i < numBins; i++)
	    {
	      if((i+1 >= firstLowRMSBin) && (i+1 <= lastLowRMSBin))
		{
		  hist->SetBinContent(i+1,10000.0);
		}
	    }
	}
      else
	{
	  for(Int_t i = 0; i < numBins; i++)
	    {
	      hist->SetBinContent(i+1,10000.0);
	    }
	}
    }
  
  return;
}


void WCPSst::NewDatauBooNEFrameDataSource::ZigzagFilterAlg(TH1F *hist)
{
  const Int_t startFiltBin = -1;
  const Int_t endFiltBin = 3700;
  
  Double_t rmsVal = CalcRMSWithFlags(hist);
  if(rmsVal < 0.5) return;
  
  Double_t waveformMean = 0.0;
  Double_t ADCval;
  Int_t counter = 0;
  Int_t numBins = hist->GetNbinsX();

  for(Int_t i = 0; i < numBins; i++)
    {    
      ADCval = hist->GetBinContent(i+1);
      if(ADCval != 10000.0)
	{
	  waveformMean += (Double_t) ADCval;
	  counter++;
	}
    }
  if(counter > 0)
    {
      waveformMean /= ((Double_t) counter);
    }
  
  TH1F *currentHist = new TH1F("","",numBins,-0.5,numBins-0.5);
  for(Int_t i = 0; i < numBins; i++)
    {
      ADCval = hist->GetBinContent(i+1);
      if(ADCval != 10000.0)
	{
	  currentHist->SetBinContent(i+1,ADCval);
	}
      else
	{
	  currentHist->SetBinContent(i+1,waveformMean);
	}
    }
  
  Double_t *reFilt = new Double_t[maxTicks];
  Double_t *imFilt = new Double_t[maxTicks];
  
  TH1 *hm = 0;
  hm = currentHist->FFT(0,"MAG");
  
  TH1 *hp = 0;
  hp = currentHist->FFT(0,"PH");
  
  for(Int_t i = 0; i < numBins; i++)
    {
      Double_t rho = hm->GetBinContent(i+1);
      Double_t phi = hp->GetBinContent(i+1);
      
      if((TMath::Min(i+1,numBins-i) > startFiltBin) && (TMath::Min(i+1,numBins-i) < endFiltBin))
	{
	  reFilt[i] = (rho*cos(phi))/numBins;
	  imFilt[i] = (rho*sin(phi))/numBins;
	}
      else
	{
	  reFilt[i] = 0.0;
	  imFilt[i] = 0.0;
	}
    }
  reFilt[0] = 0.0;
  imFilt[0] = 0.0;
  
  Int_t nFreqBins = numBins;
  TVirtualFFT *invCurrentFFTObject = TVirtualFFT::FFT(1,&nFreqBins,"C2R M K");
  invCurrentFFTObject->SetPointsComplex(reFilt,imFilt);
  invCurrentFFTObject->Transform();
  TH1F *newHist = new TH1F("","",numBins,-0.5,numBins-0.5);
  newHist = (TH1F*)TH1::TransformHisto(invCurrentFFTObject,0,"Re");
  
  for(Int_t i = 0; i < numBins; i++)
    {
      ADCval = hist->GetBinContent(i+1);
      
      if(ADCval != 10000.0)
	{
	  hist->SetBinContent(i+1,newHist->GetBinContent(i+1));
	}
    }
  
  delete hm;
  delete hp;
  delete currentHist;
  delete newHist;
  delete invCurrentFFTObject;
  delete[] reFilt;
  delete[] imFilt;
  
  return;
}


void WCPSst::NewDatauBooNEFrameDataSource::SignalFilterAlg(TH1F *hist)
{
  const Double_t sigFactor = 4.0;
  const Int_t padBins = 8;
  
  Double_t rmsVal = CalcRMSWithFlags(hist);
  Double_t sigThreshold = sigFactor*rmsVal;
  
  Double_t ADCval;
  std::vector<Bool_t> signalRegions;
  Int_t numBins = hist->GetNbinsX();

  for(Int_t i = 0; i < numBins; i++)
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
  
  for(Int_t i = 0; i < numBins; i++)
    {
      if(signalRegions[i] == true)
	{
	  for(Int_t j = TMath::Max(0,i-padBins); j < TMath::Min(numBins,i+padBins); j++)
	    {
	      ADCval = hist->GetBinContent(j+1);
	      if(ADCval < 4096.0)
		{
		  hist->SetBinContent(j+1,ADCval+20000.0);
		}
	    }
	}
    }  
  
  return;
}


void WCPSst::NewDatauBooNEFrameDataSource::NoisyFilterAlg(TH1F *hist, int planeNum)
{
  Double_t rmsVal = CalcRMSWithFlags(hist);
  const Double_t maxRMSCut[3] = {10.0,10.0,5.0};
 
  fPlaneNum = planeNum; 
  
  if(rmsVal > maxRMSCut[fPlaneNum])
    {
      Int_t numBins = hist->GetNbinsX();
      for(Int_t i = 0; i < numBins; i++)
	{
	  hist->SetBinContent(i+1,10000.0);
	}                          
    }
  
  return;
}


void WCPSst::NewDatauBooNEFrameDataSource::WaveFilterAlg(TH1F **filtHists)
{
  Int_t numBins = filtHists[0]->GetNbinsX();
  Double_t ADCval;
  std::vector<Double_t> corrVals;
  Double_t correction;
  Int_t corrValSize;

  for(Int_t i = 0; i < numBins; i++)
    {
      corrVals.clear();
      
      for(Int_t j = 0; j < waveNoiseGroupNum; j++)
	{
	  ADCval = filtHists[j]->GetBinContent(i+1);
	  if(ADCval < 4096.0)
	    {
	      corrVals.push_back(ADCval);
	    }
	}
      
      corrValSize = corrVals.size();
      sort(corrVals.begin(),corrVals.end());
      
      if(corrValSize < 2)
	{
	  correction = 0.0;
	}
      else if((corrValSize % 2) == 0)
	{
	  correction = (corrVals[corrValSize/2] + corrVals[(corrValSize/2)-1])/2.0;
	}
      else
	{
	  correction = corrVals[(corrValSize-1)/2];
	}

      for(Int_t j = 0; j < waveNoiseGroupNum; j++)
	{
	  ADCval = filtHists[j]->GetBinContent(i+1);
	  if(ADCval != 10000.0)
	    {
	      filtHists[j]->SetBinContent(i+1,TMath::Nint(ADCval-correction));
	    }
	}
    }
  
  return;
}


void  WCPSst::NewDatauBooNEFrameDataSource::RawAdaptiveBaselineAlg(TH1F *filtHist)
{
  const Int_t windowSize = 20;
  
  Int_t numBins = filtHist->GetNbinsX();
  Int_t minWindowBins = windowSize/2;
  
  Double_t baselineVec[numBins];
  Bool_t isFilledVec[numBins];
  
  Int_t numFlaggedBins = 0;
  for(Int_t j = 0; j < numBins; j++)
    {
      if(filtHist->GetBinContent(j+1) == 10000.0)
	{
	  numFlaggedBins++;
	}
    }
  if(numFlaggedBins == numBins) return; // Eventually replace this with flag check
  
  Double_t baselineVal = 0.0;
  Int_t windowBins = 0;
  Int_t index;
  Double_t ADCval;
  for(Int_t j = 0; j <= windowSize/2; j++)
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
    baselineVec[0] = baselineVal/((Double_t) windowBins);
  
  if(windowBins < minWindowBins)
    isFilledVec[0] = false;
  else
    isFilledVec[0] = true;
  
  Int_t oldIndex;
  Int_t newIndex;
  for(Int_t j = 1; j < numBins; j++)
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
  
  Int_t downIndex;
  Int_t upIndex;
  Bool_t downFlag;
  Bool_t upFlag;
  for(Int_t j = 0; j < numBins; j++)
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
		baselineVec[j] = ((j-downIndex)*baselineVec[downIndex]+(upIndex-j)*baselineVec[upIndex])/((Double_t) upIndex-downIndex);
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


void WCPSst::NewDatauBooNEFrameDataSource::RemoveFilterFlags(TH1F *filtHist)
{
  Double_t ADCval;
  Int_t numBins = filtHist->GetNbinsX();
  for(Int_t i = 0; i < numBins; i++)
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
