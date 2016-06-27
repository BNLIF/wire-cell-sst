#include "WireCellSst/DatauBooNEFrameDataSource.h"

#include "TClonesArray.h"
#include "TH1F.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "WireCellData/GeomWire.h"
#include "TF1.h"
#include "TMath.h"

using namespace WireCell;


double response(double *x, double *par){
double f = 4.31054*exp(-2.94809*x[0]/par[1])*par[0]-2.6202*exp(-2.82833*x[0]/par[1])*cos(1.19361*x[0]/par[1])*par[0]
              -2.6202*exp(-2.82833*x[0]/par[1])*cos(1.19361*x[0]/par[1])*cos(2.38722*x[0]/par[1])*par[0]
                    +0.464924*exp(-2.40318*x[0]/par[1])*cos(2.5928*x[0]/par[1])*par[0]
                          +0.464924*exp(-2.40318*x[0]/par[1])*cos(2.5928*x[0]/par[1])*cos(5.18561*x[0]/par[1])*par[0]
                                +0.762456*exp(-2.82833*x[0]/par[1])*sin(1.19361*x[0]/par[1])*par[0]
                                      -0.762456*exp(-2.82833*x[0]/par[1])*cos(2.38722*x[0]/par[1])*sin(1.19361*x[0]/par[1])*par[0]
                                            +0.762456*exp(-2.82833*x[0]/par[1])*cos(1.19361*x[0]/par[1])*sin(2.38722*x[0]/par[1])*par[0]
                                                  -2.6202*exp(-2.82833*x[0]/par[1])*sin(1.19361*x[0]/par[1])*sin(2.38722*x[0]/par[1])*par[0] 
                                                        -0.327684*exp(-2.40318*x[0]/par[1])*sin(2.5928*x[0]/par[1])*par[0] + 
                                                              +0.327684*exp(-2.40318*x[0]/par[1])*cos(5.18561*x[0]/par[1])*sin(2.5928*x[0]/par[1])*par[0]
                                                                    -0.327684*exp(-2.40318*x[0]/par[1])*cos(2.5928*x[0]/par[1])*sin(5.18561*x[0]/par[1])*par[0]
                                                                          +0.464924*exp(-2.40318*x[0]/par[1])*sin(2.5928*x[0]/par[1])*sin(5.18561*x[0]/par[1])*par[0];

 if (x[0] >0&&x[0] < 10){
    return f;
 }else{
   return 0;
 }
}


WireCellSst::DatauBooNEFrameDataSource::DatauBooNEFrameDataSource(const char* root_file, const WireCell::GeomDataSource& gds,int bins_per_frame1)
    : WireCell::FrameDataSource()
    , root_file(root_file)
    , gds(gds)
{
  bins_per_frame = bins_per_frame1;
  
  nevents = 0;
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));
  
  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();
  
  h_rc = new TH1F("h_rc","h_rc",bins_per_frame,0,bins_per_frame);
  h_1us = new TH1F("h_1us","h_1us",bins_per_frame,0,bins_per_frame);
  h_2us = new TH1F("h_2us","h_2us",bins_per_frame,0,bins_per_frame);

  TF1 *f1 = new TF1("func1",response,0,10,2);
  TF1 *f2 = new TF1("func2",response,0,10,2);
  f1->SetParameters(78.*1.012,1.0);
  f2->SetParameters(140.*1.012,2.0);


  double x0 = h_rc->GetBinCenter(1)/2.;
  for (Int_t i=0;i!=bins_per_frame;i++){
    double x = h_rc->GetBinCenter(i+1)/2.; // convert to us, assume 2 MHz digitization
    Double_t content = -1.0/2./1000 * exp(-(x-x0)/1000.); // 1 ms RC time
    if (x==x0) content +=1;
    h_rc->SetBinContent(i+1,content);
    h_1us->SetBinContent(i+1,f1->Eval(x));
    h_2us->SetBinContent(i+1,f2->Eval(x));
  }
  hm_rc = h_rc->FFT(0,"MAG");
  hp_rc = h_rc->FFT(0,"PH");
  
  hm_1us = h_1us->FFT(0,"MAG");
  hp_1us = h_1us->FFT(0,"PH");

  hm_2us = h_2us->FFT(0,"MAG");
  hp_2us = h_2us->FFT(0,"PH");

  delete f1;
  delete f2;
}

void WireCellSst::DatauBooNEFrameDataSource::Clear(){
  frame.clear();
}

WireCellSst::DatauBooNEFrameDataSource::~DatauBooNEFrameDataSource()
{
  delete h_rc;
  delete hm_rc;
  delete hp_rc;

  delete h_1us;
  delete hm_1us;
  delete hp_1us;

  delete h_2us;
  delete hm_2us;
  delete hp_2us;
}

int WireCellSst::DatauBooNEFrameDataSource::size() const
{
  //return tree->GetEntries();
  return nevents;
}

void WireCellSst::DatauBooNEFrameDataSource::Save(){

}


void WireCellSst::DatauBooNEFrameDataSource::zigzag_removal(TH1F *h1, int plane, int channel_no, int flag_RC){
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
  

  // need to restore the incorrectly set ASICs gain and shaping time ... 
  int flag_restore = 0;
  if (plane == 0){ // U-plane only    could be time-dependent 
    if (channel_no >=2016 && channel_no <= 2095 
	|| channel_no >=2192 && channel_no <=2303 
	|| channel_no >= 2352 && channel_no <2400)
      {
	flag_restore = 1;
      }
  }


  for (int j=0;j!=nbin;j++){
    
    double rho = hm->GetBinContent(j+1);
    double phi = hp->GetBinContent(j+1);
    
    if (flag_RC == 1){
      // need to remove RC+RC shapings
      if (hm_rc->GetBinContent(j+1)>0){
	rho = rho/pow(hm_rc->GetBinContent(j+1),2);
      }else{
	rho = 0;
      }
      phi = phi - 2*hp_rc->GetBinContent(j+1);
    }    

    // need to restore the incorrectly set ASICs gain and shaping time ... 
    if (flag_restore == 1){
      if (hm_1us->GetBinContent(j+1)>0){
	rho = rho / hm_1us->GetBinContent(j+1) * hm_2us->GetBinContent(j+1);
      }else{
	rho = 0;
      }
      phi = phi - hp_1us->GetBinContent(j+1) + hp_2us->GetBinContent(j+1);
    }
    // done ...

    if (plane == 0 && channel_no ==2240){
      if (j>=17&& j<=19) rho = 0;  // filter out 3.6 kHz noise for this single channel
    }

    // filter the zero frequency
    if (j==0) rho = 0; 
    
    // if (j<=3500 || j> nbin-3500){ // filter out the zigzag noise, >730 kHz noise
    value_re[j] = rho*cos(phi)/nbin;
    value_im[j] = rho*sin(phi)/nbin;
    // }else{
    //   value_re[j] = 0;
    //   value_im[j] = 0;
    // }

    // test ... 
    if (plane==0 || plane==1 ){
      if (j>=169&&j<=173) { // filter out the 36 kHz noise
	value_re[j] = 0; 
	value_im[j] = 0.;
      }
      if (plane == 0){
	if (j>=513 && j<=516){ // for U-plane, filter out the 110 kHz noise
	  value_re[j] = 0; 
	  value_im[j] = 0.;
	}

	// //another test
	// if (j==168 || j==167 || j==174 || j==175 ||
	//     j==511 || j==512 || j==517 || j==518){
	//   value_re[j] = 0.;
	//   value_im[j] = 0.;
	// }

	// //test the high-frequency component
	// if (j==857 || j==858|| j==1028 || j==1029 || j==1030 ||
	//     j==1200 || j== 1371 || j==1372 || j==1543 || j==1714 || j==1715 
	//     || j==1716 || j==1886 || j==2057){
	//   value_re[j] = 0; 
	//   value_im[j] = 0.;
	// }
	
      }
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
  
  // double xq = 0.5;
  // h2.GetQuantiles(1,&par[1],&xq);
  // xq = 0.5 + 0.34;
  // h2.GetQuantiles(1,&par[0],&xq);
  // xq = 0.5 - 0.34;
  // h2.GetQuantiels(1,&par[2],&xq);
  
  // par[2] = sqrt((pow(par[0]-par[1],2)+pow(par[2]-par[1],2))/2.);
  // par[0] = h2.GetSum();
  // f1->SetParameters(par);
  // h2.Fit(f1,"Q0","");
  // f1->GetParameters(par);
  
  // a different way to calculate the mean ... 
  double xq = 0.5;
  h2.GetQuantiles(1,&par[1],&xq);
  
  for (int j=0;j!=nbin;j++){
    h1->SetBinContent(j+1,h1->GetBinContent(j+1)-par[1]);
  }
  
  delete hp;
  delete hm;
  delete ifft;
  delete fb;
  delete f1;

}

void WireCellSst::DatauBooNEFrameDataSource::chirp_raise_baseline(TH1F *hist, int bin1, int bin2){
  for (int i=bin1;i<=bin2;i++){
    hist->SetBinContent(i+1,10000.0);
  }
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
      
      // remove this part and put the function inside the a different function ... 
      // for(int i = 0; i < numBins; i++)
      // 	{
      // 	  if((i+1 >= firstLowRMSBin) && (i+1 <= lastLowRMSBin))
      // 	    {
      // 	      hist->SetBinContent(i+1,10000.0);
      // 	    }
      // 	}

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

  // if (channel_no== 623 && planeNum==0) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  
  // if (channel_no== 678 && planeNum==1) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 1234 && planeNum==1) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 1633 && planeNum==1) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 1634 && planeNum==1) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 1673 && planeNum==1) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 1684 && planeNum==1) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 1736 && planeNum==1) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 1741 && planeNum==1) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 1792 && planeNum==1) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 1797 && planeNum==1) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 1820 && planeNum==1) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;

  // if (channel_no== 2584 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2586 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2585 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2612 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2616 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2618 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2619 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2620 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2621 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2622 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2623 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2630 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2642 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2679 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 2680 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 3264 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no== 3265 && planeNum==2) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  
   
  
  // if (channel_no==673 && planeNum==0) {
  //   std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  //   std::cout << hist->GetBinContent(100) << std::endl;
  // }
  // if (channel_no==675 && planeNum==0) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no==946 && planeNum==0) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  // if (channel_no==959 && planeNum==0) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;
  //  if (channel_no>2080 && planeNum==0) std::cout << "Xin: " << channel_no << " " << rmsVal << std::endl;

  if (planeNum == 0){
    if (channel_no < 100){
      maxRMSCut[0] = 5;
      minRMSCut[0] = 1;
    }else if (channel_no >= 100 && channel_no<2000){
      maxRMSCut[0] = 11; // increase the threshold slightly ... 
      minRMSCut[0] = 1.9;
    }else if (channel_no >= 2000 && channel_no < 2400){
      maxRMSCut[0] = 5;
      minRMSCut[0] = 0.9; // take into account FT-1 channel ... 
    }
  }else if (planeNum == 1){
    if (channel_no <290){
      maxRMSCut[1] = 5;
      minRMSCut[1] = 1;
    }else if (channel_no>=290 && channel_no < 2200){
      maxRMSCut[1] = 11;
      minRMSCut[1] = 1.9;
    }else if (channel_no >=2200){
      maxRMSCut[1] = 5;
      minRMSCut[1] = 1;
    }
  }else if (planeNum == 2){
    maxRMSCut[2] = 8;
    minRMSCut[2] = 1.3; // reduce threshold to take into account the adaptive baseline ... 
  }
  
  //test ... 
 

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

bool WireCellSst::DatauBooNEFrameDataSource::ID_RC(TH1F *h1, int plane, int channel_no){
  bool flag = false;
  TH1 *htemp_m = h1->FFT(0,"MAG");

  Double_t content[5];
  for (int i=0;i!=5;i++){
    content[i] = htemp_m->GetBinContent(i+2);
  }
  if (content[0] > content[1] && 
      content[0] > content[2] &&
      content[0] > content[3] &&
      content[0] > content[4] &&
      (content[0] + content[1] + content[2] + content[3] + content[4])/5. > 6000.){
    flag = true;
  }

  // if (flag){
  //   std::cout << content[0] << " " << content[1] << " " << content[2] << " " << 
  //     content[3] << " " << content[4] << std::endl;
  //   //remove baseline 
    

  // }
  

  delete htemp_m;

  return flag;
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

    const char* tpath = "/Event/Sim";
    TFile tfile(root_file,"read");
    TTree* tree = dynamic_cast<TTree*>(tfile.Get(tpath));

    // sigh, we can't do things this simply because the ttree does not
    // have a single branch.  
    // tree->SetBranchAddress(name, &event);
    
    tree->SetBranchStatus("*",0);
    
    tree->SetBranchStatus("eventNo",1);
    tree->SetBranchAddress("eventNo" , &event_no);
    tree->SetBranchStatus("runNo",1);
    tree->SetBranchAddress("runNo"   , &run_no);
    tree->SetBranchStatus("subRunNo",1);
    tree->SetBranchAddress("subRunNo", &subrun_no);
    
    std::vector<int> *channelid = new std::vector<int>;
    TClonesArray* esignal = new TClonesArray;
          
    tree->SetBranchStatus("raw_channelId",1);
    tree->SetBranchAddress("raw_channelId", &channelid);
    tree->SetBranchStatus("raw_wf",1);
    tree->SetBranchAddress("raw_wf", &esignal);
    
    int siz = tree->GetEntry(frame_number);
    
    
    if (siz > 0 && frame_number < siz) {
      
      TH1F **hu;
      TH1F **hv;
      TH1F **hw;
      
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
      
      
      // for (int i=0;i!=nwire_u;i++){
      //   hu[i]->Reset();
      // }
      // for (int i=0;i!=nwire_v;i++){
      //   hv[i]->Reset();
      // }
      // for (int i=0;i!=nwire_w;i++){
      //   hw[i]->Reset();
      // }
      
      
      std::cout << "Load Data " << std::endl;
      // load into frame
      int nchannels = channelid->size();
      for (size_t ind=0; ind < nchannels; ++ind) {
	TH1F* signal = dynamic_cast<TH1F*>(esignal->At(ind));
	if (!signal) continue;
	// {
	//     return -1;
	// }
	
	WireCell::Trace trace;
	trace.chid = channelid->at(ind);
	
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
      
      // //channel status check
      // if(0){ // turn this off for now ... 
      // 	std::cout << "Check Channel Status " << std::endl;
      // 	double rmsOut = 0;
      // 	for (int i=0;i!=nu;i++){
      // 		bool isCut = 0;
      // 		GetChannelStatus(hu[i], 0, i, isCut, rmsOut);
      // 		if( isCut){
      // 			if (uchirp_map.find(i) == uchirp_map.end()){
      // 	 			std::pair<int,int> abc(0, 9592);
      // 	  			uchirp_map[i] = abc;
      // 			}else{
      // 	  			uchirp_map[i].first = 0;
      // 	    			uchirp_map[i].second = 9592;
      // 			}
      // 		}
      // 	}

      // 	for (int i=0;i!=nv;i++){
      // 		bool isCut = 0;
      // 		GetChannelStatus(hv[i], 1, i, isCut, rmsOut);
      // 		if( isCut){
      // 			if (vchirp_map.find(i) == vchirp_map.end()){
      // 	 			std::pair<int,int> abc(0, 9592);
      // 	  			vchirp_map[i] = abc;
      // 			}else{
      // 	  			vchirp_map[i].first = 0;
      // 	    			vchirp_map[i].second = 9592;
      // 			}
      // 		}
      // 	}
	
      // 	for (int i=0;i!=nw;i++){
      // 		bool isCut = 0;
      // 		GetChannelStatus(hw[i], 2, i, isCut, rmsOut);
      // 		if( isCut){
      // 			if (wchirp_map.find(i) == wchirp_map.end()){
      // 	 			std::pair<int,int> abc(0, 9592);
      // 	  			wchirp_map[i] = abc;
      // 			}else{
      // 	  			wchirp_map[i].first = 0;
      // 	    			wchirp_map[i].second = 9592;
      // 			}
      // 		}
      // 	}
      // }


      //special cut on W-plane ... 
      for (int i=0;i!=nw;i++){
	bool isCut = 0;
	int chan = i;
	if (chan >=7136 - 4800 && chan <=7263 - 4800){
	  if (chan != 7200- 4800 && chan!=7215 - 4800)
	    isCut = 1;
	}
	if( isCut){
	  if (wchirp_map.find(i) == wchirp_map.end()){
	    std::pair<int,int> abc(0, 9591);
	    wchirp_map[i] = abc;
	  }else{
	    wchirp_map[i].first = 0;
	    wchirp_map[i].second = 9591;
	  }
	}
      }


      // if (uchirp_map.find(1517)!=uchirp_map.end()){
      // 	std::cout << "Xin Channel Status!" << std::endl;
      // }
      //hack for now to deal with FT-1 channel ... 
      // for (int i=2192;i!=nu;i++){
      // 	hu[i]->Scale(2);
      // } // removed, as we implement the special treatment of FT-1 channels
      
      
      // std::cout << "Before Chirp ID " << uchirp_map.size() << " " << 
      // 	vchirp_map.size() << " " << wchirp_map.size() << std::endl;

      for (int i=2016; i<2400;i++){
	int channel_no = i;
	if (channel_no >=2016 && channel_no <= 2095 
	    || channel_no >=2192 && channel_no <=2303 
	    || channel_no >= 2352 && channel_no <2400){
	  hu[i]->Scale(14./7.8); // assume 7.8 mV/fC gain
	}
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

      
      for (int i=2016; i<2400;i++){
	int channel_no = i;
	if (channel_no >=2016 && channel_no <= 2095 
	    || channel_no >=2192 && channel_no <=2303 
	    || channel_no >= 2352 && channel_no <2400){
	  hu[i]->Scale(7.8/14.); // assume 7.8 mV/fC gain
	}
      }

      
      // if (uchirp_map.find(2240)!=uchirp_map.end()){
      // 	std::cout << "1: " << 2240 << std::endl; 
      // }
      


      // for (auto it = uchirp_map.begin(); it!= uchirp_map.end(); it++){
      // 	if (it->second.first!=0 || it->second.second!=9591){
      // 	  std::cout << "U: " << it->first << " " << it->second.first << " " << it->second.second << std::endl;
      // 	}
      // }
      // for (auto it = vchirp_map.begin(); it!= vchirp_map.end(); it++){
      // 	if (it->second.first!=0 || it->second.second!=9591){
      // 	  std::cout << "V: " << it->first << " " << it->second.first << " " << it->second.second << std::endl;
      // 	}
      // }
      // for (auto it = wchirp_map.begin(); it!= wchirp_map.end(); it++){
      // 	if (it->second.first!=0 || it->second.second!=9591){
      // 	  std::cout << "W: " << it->first << " " << it->second.first << " " << it->second.second << std::endl;
      // 	}
      // }

      std::cout << "ID RC channels!" << std::endl;
      
      std::vector<int> ided_rc_uplane;
      std::vector<int> ided_rc_vplane;
      std::vector<int> ided_rc_wplane;

      for (int i=0;i!=nu;i++){
      	if (ID_RC(hu[i],0,i)){
	  ided_rc_uplane.push_back(i);
	  //std::cout << "U: " << i << std::endl;
	}
      }
      for (int i=0;i!=nv;i++){
	if (ID_RC(hv[i],0,i)){
	  ided_rc_vplane.push_back(i);
	  //std::cout << "V: " << i << std::endl;
	}
      }
      for (int i=0;i!=nw;i++){
	if (ID_RC(hw[i],0,i)){
	  ided_rc_wplane.push_back(i);
	  //std::cout << "W: " << i << std::endl;
	}
      }

      // std::cout << ided_rc_uplane.size() << " " << ided_rc_vplane.size()
      // 		<< " " << ided_rc_wplane.size() << std::endl;

     


      std::cout << "Remove ZigZag " << std::endl;
      // deal with the zig zag noise
      // filter single frequency 36 and 110 kHz
      // correct RC+RC 
      // correct misconfigured channel (need a database ...)
      for (int i=0;i!=nu;i++){
	auto it = find(ided_rc_uplane.begin(),ided_rc_uplane.end(),i);
	if (it == ided_rc_uplane.end()){
	  zigzag_removal(hu[i],0,i);
	}else{
	  zigzag_removal(hu[i],0,i,0);
	}
      }
      for (int i=0;i!=nv;i++){
	auto it = find(ided_rc_vplane.begin(),ided_rc_vplane.end(),i);
	if (it == ided_rc_vplane.end()){
	  zigzag_removal(hv[i],1,i);
	}else{
	  zigzag_removal(hv[i],1,i,0);
	}
      }
      for (int i=0;i!=nw;i++){
	auto it = find(ided_rc_wplane.begin(),ided_rc_wplane.end(),i);
	if (it == ided_rc_wplane.end()){
	  zigzag_removal(hw[i],2,i);
	}else{
	  zigzag_removal(hw[i],2,i,0);
	}
      }

      
      // // put RC into adaptive baseline group  ... 
      // for (int i=0;i!=ided_rc_uplane.size();i++){
      // 	if (uchirp_map.find(ided_rc_uplane.at(i)) == uchirp_map.end()){
      // 	  std::pair<int,int> abc(-1, -1);
      // 	  uchirp_map[ided_rc_uplane.at(i)] = abc;
      // 	}
      // }
      // for (int i=0;i!=ided_rc_vplane.size();i++){
      // 	if (vchirp_map.find(ided_rc_vplane.at(i)) == vchirp_map.end()){
      // 	  std::pair<int,int> abc(-1, -1);
      // 	  vchirp_map[ided_rc_vplane.at(i)] = abc;
      // 	}
      // }
      // for (int i=0;i!=ided_rc_wplane.size();i++){
      // 	if (wchirp_map.find(ided_rc_wplane.at(i)) == wchirp_map.end()){
      // 	  std::pair<int,int> abc(-1, -1);
      // 	  wchirp_map[ided_rc_wplane.at(i)] = abc;
      // 	}
      // }
      


      std::cout << "Adaptive Baseline " << uchirp_map.size() << " " << 
      	vchirp_map.size() << " " << wchirp_map.size() << std::endl;
      // do the adaptive baseline ... 
      for (auto it = uchirp_map.begin(); it!= uchirp_map.end(); it++){
      	chirp_raise_baseline(hu[it->first],it->second.first,it->second.second);
	SignalFilter(hu[it->first]);
      	RawAdaptiveBaselineAlg(hu[it->first]);
      	//RemoveFilterFlags(hu[it->first]);
      }
      for (auto it = vchirp_map.begin(); it!= vchirp_map.end(); it++){
	chirp_raise_baseline(hv[it->first],it->second.first,it->second.second);
      	SignalFilter(hv[it->first]);
      	RawAdaptiveBaselineAlg(hv[it->first]);
      	//RemoveFilterFlags(hv[it->first]);
	
      }
      for (auto it = wchirp_map.begin(); it!= wchirp_map.end(); it++){
	chirp_raise_baseline(hw[it->first],it->second.first,it->second.second);
      	SignalFilter(hw[it->first]);
      	RawAdaptiveBaselineAlg(hw[it->first]);
      	//RemoveFilterFlags(hw[it->first]);
      }


      // do the adaptive baseline for the bad RC channels ... 
      for (int i=0; i!=ided_rc_uplane.size();i++){
      	if (uchirp_map.find(ided_rc_uplane.at(i)) == uchirp_map.end()){
      	  SignalFilter(hu[ided_rc_uplane.at(i)]);
      	  RawAdaptiveBaselineAlg(hu[ided_rc_uplane.at(i)]);
      	}
      }
      for (int i=0; i!=ided_rc_vplane.size();i++){
      	if (vchirp_map.find(ided_rc_vplane.at(i)) == vchirp_map.end()){
      	  SignalFilter(hv[ided_rc_vplane.at(i)]);
      	  RawAdaptiveBaselineAlg(hv[ided_rc_vplane.at(i)]);
      	}
      }
      for (int i=0; i!=ided_rc_wplane.size();i++){
      	if (wchirp_map.find(ided_rc_wplane.at(i)) == wchirp_map.end()){
      	  SignalFilter(hw[ided_rc_wplane.at(i)]);
      	  RawAdaptiveBaselineAlg(hw[ided_rc_wplane.at(i)]);
      	}
      }


      // std::cout << "2" << std::endl;
      //  for (auto it = uchirp_map.begin(); it!= uchirp_map.end(); it++){
      // 	if (it->second.first!=0 || it->second.second!=9591){
      // 	  std::cout << "U: " << it->first << " " << it->second.first << " " << it->second.second << std::endl;
      // 	}
      // }
      // for (auto it = vchirp_map.begin(); it!= vchirp_map.end(); it++){
      // 	if (it->second.first!=0 || it->second.second!=9591){
      // 	  std::cout << "V: " << it->first << " " << it->second.first << " " << it->second.second << std::endl;
      // 	}
      // }
      // for (auto it = wchirp_map.begin(); it!= wchirp_map.end(); it++){
      // 	if (it->second.first!=0 || it->second.second!=9591){
      // 	  std::cout << "W: " << it->first << " " << it->second.first << " " << it->second.second << std::endl;
      // 	}
      // }
      // if (uchirp_map.find(1517)!=uchirp_map.end()){
      // 	std::cout << "Xin: Chirping !" << std::endl;
      // }

      
      
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


      // test for Brian ... 
      int pad_window = 0;


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
	float max_rms = 0;
      	for (int j=0;j!=uplane_all.at(i).size();j++){
      	  if (urms_map[uplane_all.at(i).at(j)] > max_rms) max_rms = urms_map[uplane_all.at(i).at(j)];
      	}
	TH1F *h3 = new TH1F("h3","h3",int(12*max_rms),-6*max_rms,6*max_rms);
	TH1F *h44 = new TH1F("h44","h44",nbin,0,nbin);
	
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
	    h44->SetBinContent(j+1,par[1]);
      	  }
	  
	  // do the RMS
	  float min = h44->GetMinimum();
	  float max = h44->GetMaximum();
	  TH1F *h55 = new TH1F("h55","h55",int(max-min+1),min,max+1);
	  // std::cout << max << " " << min << " " << int(max-min+1) << " " << nbin << std::endl;
	  for (int j=0;j!=nbin;j++){
	    h55->Fill(h44->GetBinContent(j+1));
	  }
	  float mean = h55->GetMean();
	  float rms = h55->GetRMS();
	  for (int j=0;j!=h55->GetNbinsX();j++){
	    int bin_center = h55->GetBinCenter(j+1);
	    if (bin_center < mean - 4.5*rms || bin_center > mean + 4.5*rms){
	      h55->SetBinContent(j+1,0);
	    }
	  }
	 
	  mean = h55->GetMean();
	  rms = h55->GetRMS();
	  delete h55;
	  // remove +- 3sigma one
	  std::vector<int> signals;
	  std::vector<bool> signalsBool;
    	  for (int j=0;j!=nbin;j++)
		signalsBool.push_back(0);

	  for (int j=0;j!=nbin;j++){
	    float content = h44->GetBinContent(j+1);
	    if (fabs(content-mean)>3.0*rms){
	      h44->SetBinContent(j+1,0);
	      //signals.push_back(j);
	      signalsBool.at(j) = 1;
	    
	      // add the front and back padding
	      for (int k=0;k!=pad_window;k++){
	        int bin = j+k+1;
	        if (bin > nbin-1) bin = nbin-1;
	        signalsBool.at(bin) = 1;
	        //auto it = find(signals.begin(),signals.end(),bin);
	        //if (it == signals.end())
		  //signals.push_back(bin);
	        bin = j-k-1;
	        if (bin <0) bin = 0;
	        signalsBool.at(bin) = 1;
	        //it = find(signals.begin(),signals.end(),bin);
	        //if (it == signals.end())
		//signals.push_back(bin);
	      }
	    }
	  }
	  for (int j=0;j!=nbin;j++)
		if( signalsBool.at(j) == 1 )
			signals.push_back(j);
	  
	  // adaptive baseline 
	  for (int j=0;j!=signals.size();j++){
	    int bin = signals.at(j);
	    int prev_bin=bin;
	    int next_bin=bin;
	    
	    int flag = 1;
	    while(flag){
	      prev_bin--;
	      if (find(signals.begin(),signals.end(),prev_bin)==signals.end() || prev_bin <=0){
	  	flag = 0;
	      }
	    }
	    
	    // prev_bin = prev_bin - pad_window;
	    // if (prev_bin <0) prev_bin = 0;



	    flag =1;
	    while(flag){
	      next_bin++;
	      if (find(signals.begin(),signals.end(),next_bin)==signals.end() || next_bin >=nbin-1){
	  	flag = 0;
	      }
	    }
	    
	    // next_bin = next_bin + pad_window;
	    // if (next_bin > nbin-1) next_bin = nbin-1; 
	    
	    // if (prev_bin>0) prev_bin --;
	    // if (next_bin<nbin-1) next_bin++; 


	    float prev_content, next_content;
	    // if (prev_bin >=4){
	    //   prev_content = (h44->GetBinContent(prev_bin+1) + h44->GetBinContent(prev_bin) + h44->GetBinContent(prev_bin-1) + 
	    // 		      h44->GetBinContent(prev_bin-2) + h44->GetBinContent(prev_bin-3))/5.;
	    // }else{
	    prev_content = h44->GetBinContent(prev_bin+1);
	    // }
	    // if (next_bin <= nbin-5){
	    //   next_content = (h44->GetBinContent(next_bin+1) + h44->GetBinContent(next_bin+2) + h44->GetBinContent(next_bin+3)+
	    // 		      h44->GetBinContent(next_bin+4) + h44->GetBinContent(next_bin+5))/5.;
	    // }else{
	    next_content = h44->GetBinContent(next_bin+1);
	    // }
	    

	    //std::cout << prev_bin << " " << bin << " " << next_bin << " " << signals.size() << std::endl;
	    float content = prev_content + (bin - prev_bin)/ (next_bin - prev_bin*1.0) 
	      * (next_content - prev_content);
	    h44->SetBinContent(bin+1,content);
	  }

	  
	  // calculate scaling coefficient ... 
	  double ave_coef = 0;
	  double_t ave_coef1 = 0;
	  std::vector<double> coef_all;
	  for (int k=0;k!=uplane_all.at(i).size();k++){
	    double sum2 = 0;
	    double sum3 = 0;
	    double coef = 0;
	    
	    for (int j=0;j!=nbin;j++){
	      if (fabs(hu[uplane_all.at(i).at(k)]->GetBinContent(j+1)) < 4 * urms_map[uplane_all.at(i).at(k)] ){
		sum2 += hu[uplane_all.at(i).at(k)]->GetBinContent(j+1) * h44->GetBinContent(j+1);
		sum3 += h44->GetBinContent(j+1) * h44->GetBinContent(j+1);
	      }
	    }
	    if (sum3 >0)
	      coef = sum2/sum3;
	    coef_all.push_back(coef);
	    
	    if (coef !=0){
	      ave_coef += coef;
	      ave_coef1 ++;
	    }
	    
	  }
	  if (ave_coef1>0)
	    ave_coef = ave_coef / ave_coef1;
	  
	  // for (int k=0;k!=uplane_all.at(i).size();k++){
	  //   std::cout << "U " << uplane_all.at(i).at(k) << " " <<  coef_all.at(k) << " " << ave_coef << " " << coef_all.at(k)/ave_coef << std::endl;
	  //   //std::cout << "U Ave: " << ave_coef << std::endl;
	  // }
	  //  

	  //h44->Reset();

	  for (int j=0;j!=nbin;j++){
	    for (int k=0;k!=uplane_all.at(i).size();k++){
      	      if (fabs(hu[uplane_all.at(i).at(k)]->GetBinContent(j+1))>0.001)
      		hu[uplane_all.at(i).at(k)]->SetBinContent(j+1,hu[uplane_all.at(i).at(k)]->GetBinContent(j+1)-h44->GetBinContent(j+1) * coef_all.at(k)/ave_coef);
	      //hu[uplane_all.at(i).at(k)]->SetBinContent(j+1,hu[uplane_all.at(i).at(k)]->GetBinContent(j+1)-h44->GetBinContent(j+1));
      	    }
	  }


      	}

      	delete h3;
	delete h44;
      }
      
      for (int i=0;i!=vplane_all.size();i++){
      	//std::cout << "V " << i << " " << vplane_all.size() << std::endl;
      	float max_rms = 0;
      	for (int j=0;j!=vplane_all.at(i).size();j++){
      	  if (vrms_map[vplane_all.at(i).at(j)] > max_rms) max_rms = vrms_map[vplane_all.at(i).at(j)];
      	}
      	// if (max_rms<10) max_rms = 10;
	
	
      	TH1F *h3 = new TH1F("h3","h3",int(12*max_rms),-6*max_rms,6*max_rms);
	TH1F *h44 = new TH1F("h44","h44",nbin,0,nbin);

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
	    h44->SetBinContent(j+1,par[1]);
      	  }
	  
	   // do the RMS
	  float min = h44->GetMinimum();
	  float max = h44->GetMaximum();
	  TH1F *h55 = new TH1F("h55","h55",int(max-min+1),min,max+1);
	  //std::cout << max << " " << min << " " << int(max-min+1) << std::endl;
	  for (int j=0;j!=nbin;j++){
	    h55->Fill(h44->GetBinContent(j+1));
	  }
	  float mean = h55->GetMean();
	  float rms = h55->GetRMS();
	  for (int j=0;j!=h55->GetNbinsX();j++){
	    int bin_center = h55->GetBinCenter(j+1);
	    if (bin_center < mean - 4.5*rms || bin_center > mean + 4.5*rms){
	      h55->SetBinContent(j+1,0);
	    }
	  }
	  
	  mean = h55->GetMean();
	  rms = h55->GetRMS();
	  delete h55;
	  // remove +- 3sigma one
	  std::vector<int> signals;
	  std::vector<bool> signalsBool;
    	  for (int j=0;j!=nbin;j++)
		signalsBool.push_back(0);

	  for (int j=0;j!=nbin;j++){
	    float content = h44->GetBinContent(j+1);
	    if (fabs(content-mean)>3.0*rms){
	      h44->SetBinContent(j+1,0);
	      //signals.push_back(j);
	      signalsBool.at(j) = 1;
	    
	      // add the front and back padding
	      for (int k=0;k!=pad_window;k++){
	        int bin = j+k+1;
	        if (bin > nbin-1) bin = nbin-1;
	        signalsBool.at(bin) = 1;
	        //auto it = find(signals.begin(),signals.end(),bin);
	        //if (it == signals.end())
		  //signals.push_back(bin);
	        bin = j-k-1;
	        if (bin <0) bin = 0;
	        signalsBool.at(bin) = 1;
	        //it = find(signals.begin(),signals.end(),bin);
	        //if (it == signals.end())
		//signals.push_back(bin);
	      }
	    }
	  }
	  for (int j=0;j!=nbin;j++)
		if( signalsBool.at(j) == 1 )
			signals.push_back(j);
	  
	  // adaptive baseline 
	  for (int j=0;j!=signals.size();j++){
	    int bin = signals.at(j);
	    int prev_bin=bin;
	    int next_bin=bin;
	    
	    int flag = 1;
	    while(flag){
	      prev_bin--;
	      if (find(signals.begin(),signals.end(),prev_bin)==signals.end() || prev_bin <=0){
	  	flag = 0;
	      }
	    }

	    // prev_bin = prev_bin - pad_window;
	    // if (prev_bin <0) prev_bin = 0;

	    flag =1;
	    while(flag){
	      next_bin++;
	      if (find(signals.begin(),signals.end(),next_bin)==signals.end() || next_bin >=nbin-1){
	  	flag = 0;
	      }
	    }
	    
	    //  next_bin = next_bin + pad_window;
	    // if (next_bin > nbin-1) next_bin = nbin-1; 

	    //  if (prev_bin>0) prev_bin --;
	    // if (next_bin<nbin-1) next_bin++; 

	    float prev_content, next_content;
	    // if (prev_bin >=4){
	    //   prev_content = (h44->GetBinContent(prev_bin+1) + h44->GetBinContent(prev_bin) + h44->GetBinContent(prev_bin-1) + 
	    // 		      h44->GetBinContent(prev_bin-2) + h44->GetBinContent(prev_bin-3))/5.;
	    // }else{
	      prev_content = h44->GetBinContent(prev_bin+1);
	    // }
	    // if (next_bin <= nbin-5){
	    //   next_content = (h44->GetBinContent(next_bin+1) + h44->GetBinContent(next_bin+2) + h44->GetBinContent(next_bin+3)+
	    // 		      h44->GetBinContent(next_bin+4) + h44->GetBinContent(next_bin+5))/5.;
	    // }else{
	      next_content = h44->GetBinContent(next_bin+1);
	    // }
	    

	    //std::cout << prev_bin << " " << bin << " " << next_bin << " " << signals.size() << std::endl;
	    float content = prev_content + (bin - prev_bin)/ (next_bin - prev_bin*1.0) 
	      * (next_content - prev_content);

	    h44->SetBinContent(bin+1,content);
	  }


	  // calculate scaling coefficient ... 
	  double ave_coef = 0;
	  double ave_coef1 = 0;
	  std::vector<double> coef_all;
	  for (int k=0;k!=vplane_all.at(i).size();k++){
	    double sum2 = 0;
	    double sum3 = 0;
	    double coef = 0;
	    
	    for (int j=0;j!=nbin;j++){
	      if (fabs(hv[vplane_all.at(i).at(k)]->GetBinContent(j+1)) < 4 * vrms_map[vplane_all.at(i).at(k)] ){
		sum2 += hv[vplane_all.at(i).at(k)]->GetBinContent(j+1) * h44->GetBinContent(j+1);
		sum3 += h44->GetBinContent(j+1) * h44->GetBinContent(j+1);
	      }
	    }
	    if (sum3 >0)
	      coef = sum2/sum3;

	    coef_all.push_back(coef);
	    
	    //	    std::cout << k << " " << coef << std::endl;

	    if (coef!=0){
	      ave_coef += coef;
	      ave_coef1 ++;
	    }
	  }
	  if (ave_coef1>0)
	    ave_coef = ave_coef / ave_coef1;
	  //  
	  //	  std::cout << "VAve: " << ave_coef << std::endl;
	  
	  // for (int k=0;k!=vplane_all.at(i).size();k++){
	  //   std::cout << "V " <<  vplane_all.at(i).at(k) << " " <<  coef_all.at(k) << " " << ave_coef << " " << coef_all.at(k)/ave_coef << std::endl;
	  //   //std::cout << "U Ave: " << ave_coef << std::endl;
	  // }

	  //h44->Reset();

	  for (int j=0;j!=nbin;j++){
	    for (int k=0;k!=vplane_all.at(i).size();k++){
      	      if (fabs(hv[vplane_all.at(i).at(k)]->GetBinContent(j+1))>0.001)
      		hv[vplane_all.at(i).at(k)]->SetBinContent(j+1,hv[vplane_all.at(i).at(k)]->GetBinContent(j+1)-h44->GetBinContent(j+1)* coef_all.at(k)/ave_coef);
	      //hv[vplane_all.at(i).at(k)]->SetBinContent(j+1,hv[vplane_all.at(i).at(k)]->GetBinContent(j+1)-h44->GetBinContent(j+1));
      	    }
	  }

      	}

      	delete h3;
	delete h44;
      }
      
      for (int i=0;i!=wplane_all.size();i++){
      	float max_rms = 0;
      	for (int j=0;j!=wplane_all.at(i).size();j++){
      	  if (wrms_map[wplane_all.at(i).at(j)] > max_rms) max_rms = wrms_map[wplane_all.at(i).at(j)];
      	}
      	// if (max_rms<10) max_rms = 10;
	
      	//std::cout << "W " << i << " " << wplane_all.size() << std::endl;
      	TH1F *h3 = new TH1F("h3","h3",int(12*max_rms),-6*max_rms,6*max_rms);
	TH1F *h44 = new TH1F("h44","h44",nbin,0,nbin);
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
	    h44->SetBinContent(j+1,par[1]);
      	  }

	   // do the RMS
	  float min = h44->GetMinimum();
	  float max = h44->GetMaximum();
	  TH1F *h55 = new TH1F("h55","h55",int(max-min+1),min,max+1);
	  //std::cout << max << " " << min << " " << int(max-min+1) << std::endl;
	  for (int j=0;j!=nbin;j++){
	    h55->Fill(h44->GetBinContent(j+1));
	  }
	  float mean = h55->GetMean();
	  float rms = h55->GetRMS();
	  for (int j=0;j!=h55->GetNbinsX();j++){
	    int bin_center = h55->GetBinCenter(j+1);
	    if (bin_center < mean - 4.5*rms || bin_center > mean + 4.5*rms){
	      h55->SetBinContent(j+1,0);
	    }
	  }
	  
	  mean = h55->GetMean();
	  rms = h55->GetRMS();
	  delete h55;
	  // remove +- 3sigma one
	  std::vector<int> signals;
	  std::vector<bool> signalsBool;
    	  for (int j=0;j!=nbin;j++)
		signalsBool.push_back(0);

	  for (int j=0;j!=nbin;j++){
	    float content = h44->GetBinContent(j+1);
	    if (fabs(content-mean)>3.0*rms){
	      h44->SetBinContent(j+1,0);
	      //signals.push_back(j);
	      signalsBool.at(j) = 1;
	    
	      // add the front and back padding
	      for (int k=0;k!=pad_window;k++){
	        int bin = j+k+1;
	        if (bin > nbin-1) bin = nbin-1;
	        signalsBool.at(bin) = 1;
	        //auto it = find(signals.begin(),signals.end(),bin);
	        //if (it == signals.end())
		  //signals.push_back(bin);
	        bin = j-k-1;
	        if (bin <0) bin = 0;
	        signalsBool.at(bin) = 1;
	        //it = find(signals.begin(),signals.end(),bin);
	        //if (it == signals.end())
		//signals.push_back(bin);
	      }
	    }
	  }
	  for (int j=0;j!=nbin;j++)
		if( signalsBool.at(j) == 1 )
			signals.push_back(j);
	  
	  // adaptive baseline 
	  for (int j=0;j!=signals.size();j++){
	    int bin = signals.at(j);
	    int prev_bin=bin;
	    int next_bin=bin;
	    
	    int flag = 1;
	    while(flag){
	      prev_bin--;
	      if (find(signals.begin(),signals.end(),prev_bin)==signals.end() || prev_bin <=0){
	  	flag = 0;
	      }
	    }

	    // prev_bin = prev_bin - pad_window;
	    // if (prev_bin <0) prev_bin = 0;

	    flag =1;
	    while(flag){
	      next_bin++;
	      if (find(signals.begin(),signals.end(),next_bin)==signals.end() || next_bin >=nbin-1){
	  	flag = 0;
	      }
	    }

	    // next_bin = next_bin + pad_window;
	    // if (next_bin > nbin-1) next_bin = nbin-1; 


	    // if (prev_bin>0) prev_bin --;
	    // if (next_bin<nbin-1) next_bin++; 

	    float prev_content, next_content;
	    // if (prev_bin >=4){
	    //   prev_content = (h44->GetBinContent(prev_bin+1) + h44->GetBinContent(prev_bin) + h44->GetBinContent(prev_bin-1) + 
	    // 		      h44->GetBinContent(prev_bin-2) + h44->GetBinContent(prev_bin-3))/5.;
	    // }else{
	      prev_content = h44->GetBinContent(prev_bin+1);
	    // }
	    // if (next_bin <= nbin-5){
	    //   next_content = (h44->GetBinContent(next_bin+1) + h44->GetBinContent(next_bin+2) + h44->GetBinContent(next_bin+3)+
	    // 		      h44->GetBinContent(next_bin+4) + h44->GetBinContent(next_bin+5))/5.;
	    // }else{
	      next_content = h44->GetBinContent(next_bin+1);
	    // }
	    

	    //std::cout << prev_bin << " " << bin << " " << next_bin << " " << signals.size() << std::endl;
	    float content = prev_content + (bin - prev_bin)/ (next_bin - prev_bin*1.0) 
	      * (next_content - prev_content);

	    h44->SetBinContent(bin+1,content);
	  }


	  // calculate scaling coefficient ... 
	  double ave_coef = 0;
	  double ave_coef1 = 0;
	  std::vector<double> coef_all;
	  for (int k=0;k!=wplane_all.at(i).size();k++){
	    double sum2 = 0;
	    double sum3 = 0;
	    double coef = 0;
	    
	    for (int j=0;j!=nbin;j++){
	      if (fabs(hw[wplane_all.at(i).at(k)]->GetBinContent(j+1)) < 4 * wrms_map[wplane_all.at(i).at(k)] ){
		sum2 += hw[wplane_all.at(i).at(k)]->GetBinContent(j+1) * h44->GetBinContent(j+1);
		sum3 += h44->GetBinContent(j+1) * h44->GetBinContent(j+1);
	      }
	    }
	    if (sum3 >0)
	      coef = sum2/sum3;
	    // std::cout << k << " " << coef << std::endl;
	    coef_all.push_back(coef);
	    if (coef!=0){
	      ave_coef += coef;
	      ave_coef1 ++;
	    }
	  }
	  if (ave_coef1)
	    ave_coef = ave_coef / ave_coef1;
	  //  
	  //	  std::cout << "WAve: " << ave_coef << std::endl;
	  // for (int k=0;k!=wplane_all.at(i).size();k++){
	  //   std::cout << "W " << wplane_all.at(i).at(k) << " " <<  coef_all.at(k) << " " << ave_coef << " " << coef_all.at(k)/ave_coef << std::endl;
	  //   //std::cout << "U Ave: " << ave_coef << std::endl;
	  // }

	  //h44->Reset();
	   
	  for (int j=0;j!=nbin;j++){
	    for (int k=0;k!=wplane_all.at(i).size();k++){
      	      if (fabs(hw[wplane_all.at(i).at(k)]->GetBinContent(j+1))>0.001)
      		hw[wplane_all.at(i).at(k)]->SetBinContent(j+1,hw[wplane_all.at(i).at(k)]->GetBinContent(j+1)-h44->GetBinContent(j+1)* coef_all.at(k)/ave_coef);
		//hw[wplane_all.at(i).at(k)]->SetBinContent(j+1,hw[wplane_all.at(i).at(k)]->GetBinContent(j+1)-h44->GetBinContent(j+1));
      	    }
	  }

      	}
      	delete h3;
	delete h44;
      }
      
      
      // put down all u 
      // for (int i=0;i!=nu;i++){
      // 	uchirp_map[i] = std::pair<int,int> (0,bins_per_frame);
      // }
      // for (int i=0;i!=nw;i++){
      // 	wchirp_map[i] = std::pair<int,int> (0,bins_per_frame);
      // }
      // vchirp_map.clear();


      // deal with the w plane to remove the PMT signal (negative pulse ...)
      for (int i=0;i!=nwire_w;i++){
	//RemovePMTSignalCollection(hw[i],wrms_map[i]);
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
	trace.chid = channelid->at(ind);
	
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
      
      // TFile *file = new TFile("temp_data.root","RECREATE");
      // for (int i=0;i!=nwire_u;i++){
      //   TH1F *huu = (TH1F*)hu[i]->Clone(Form("U1_%d",i));
      // }
      // for (int i=0;i!=nwire_v;i++){
      //   TH1F *hvv = (TH1F*)hv[i]->Clone(Form("V1_%d",i));
      // }
      // for (int i=0;i!=nwire_w;i++){
      //   TH1F *hww = (TH1F*)hw[i]->Clone(Form("W1_%d",i));
      // }
      // file->Write();
      // file->Close();
      // std::cout << "Saved file" << std::endl;
      
      for (int i=0;i!=nwire_u;i++){
	delete hu[i] ;
      }
      delete [] hu;
      for (int i=0;i!=nwire_v;i++){
	delete hv[i] ;
      }
      delete [] hv;
      for (int i=0;i!=nwire_w;i++){
	delete hw[i] ;
      }
      delete [] hw;
      
      
      tfile.Close();
      esignal->Clear("D");
      delete channelid;
      delete esignal;
      
      frame.index = frame_number;
      return frame.index;
    }else{
      
      tfile.Close();
      esignal->Clear("D");
      delete channelid;
      delete esignal;

      return -1;
    }
}


void WireCellSst::DatauBooNEFrameDataSource::RemovePMTSignalCollection(TH1F* hist, float rms){
  float rms1 = 0;
  float rms2 = 0;
  
  for (int i=0;i!=hist->GetNbinsX();i++){
    float content = hist->GetBinContent(i+1);
    if (fabs(content) < 3*rms){
      rms1 += content*content;
      rms2 ++;
    }
  }
  
  if (rms2 >0){
    rms1 = sqrt(rms1/rms2);
    for (int i=0;i!=hist->GetNbinsX();i++){
      float content = hist->GetBinContent(i+1);
      if (content < -3 *rms1){
	hist->SetBinContent(i+1,0);
      }
    }
  }
  
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

void WireCellSst::DatauBooNEFrameDataSource::GetChannelStatus(TH1F *h1, int plane, int chan, bool& isCut, double &rmsOut){
	//calculate mean
    	double mean = 0;
    	int count = 0;
    	for(unsigned int it = 0 ; it < h1->GetNbinsX() ; it++ ){
		mean += h1->GetBinContent(it+1);
		count++;
    	}
    	if( count > 0)
		mean = mean / (double) count;

	//calculate RMS
	count = 0;
    	double rms = 0;
    	for(unsigned int it = 0 ; it < h1->GetNbinsX() ; it++ ){
		int adc = h1->GetBinContent(it+1);
		rms += (adc - mean)*(adc - mean);
		count++;
    	}
    	if(count - 1 > 0)
		rms = TMath::Sqrt( rms /(double)(  count - 1  ) );

	rmsOut = rms;

	//	if (fabs(chan-1517)<5 && plane == 0) std::cout << "Xin: " << chan << " " << rmsOut << std::endl;

	//apply selection
	bool isCut1 = 0;
	bool isCut2 = 0;
	bool isCut3 = 0;

	//selection cut 1 : pedestal mean
	//if( chan < 4800 && (mean < 2046 - 50 || mean > 2046 + 50) )
	//	isCut1 = 1;
	//if( chan >= 4800 && (mean < 474 - 50 || mean > 474 + 50) )
	//	isCut1 = 1;

	//selection cut 2: noisy channel
	//if( rms > 30 )
	if( rms > 30 )
		isCut2 = 1;

	//selection cut 3: low RMS
	double limit = -1;
	double uplane_limit = 8.6;
	if( plane == 0 ){
		if( chan < 680 )
			limit = (chan - 0.)*(uplane_limit - 1.8 )/(680. - 0.) + 1.8;
		if( chan >= 680 && chan < 1728 )
			limit = uplane_limit;
		if( chan >= 1728 ) 
			limit = (chan - 1728)*( 1.7 - uplane_limit )/(2400. - 1728.) + uplane_limit;
		if( chan > 2000 && chan < 2100 ) limit = 0;
		if( chan > 2180 && chan < 2400 ) limit = 0;
	}
	if( plane == 1 ){
		if( chan < 576 )
			limit = (chan - 0.)*(4.0 - 1.3 )/(576.-0.) + 1.3;
		if( chan >= 576 && chan < 1921 )
			limit = 7.0;
		if( chan >= 1921 ) 
			limit = (chan - 1921)*( 1.5 - 4.5 )/(2400. - 1921.) + 4.5;
	}
	if( plane == 2  )
		limit = 3;
	limit = limit - 1;
	if( limit < 0 )	
		limit = 0;
	if( rms < limit )
		isCut3 = 1;

	//isCut = isCut1 + isCut2;
	isCut = isCut1+isCut2+isCut3;

	// Add special projection for W-plane signal ... 
	if (plane == 2 && chan >=7136 - 4800 && chan <=7263 - 4800){
	  if (chan != 7200- 4800 && chan!=7215 - 4800)
	    isCut = isCut + 1;
	}


	return;
  }
