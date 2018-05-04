#include "WireCellSst/DatauBooNEFrameDataSource.h"

#include "TClonesArray.h"
#include "TH1F.h"
#include "TFile.h"
#include "TRandom.h"
#include "TVirtualFFT.h"
#include "TGraph.h"
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

Double_t PoissonReal(const Double_t *k,  const Double_t *lambda) {
  return TMath::Exp(k[0]*TMath::Log(lambda[0])-lambda[0]) / TMath::Gamma(k[0]+1.);
}

void WireCellSst::DatauBooNEFrameDataSource::Simu_Noise_uBooNE_Empirical(TH1F *h1, Int_t plane, Int_t channel){
  
  Double_t L = 0;
  if (plane == 2) L = 233; // cm
  if (plane ==0 || plane == 1){
    if (channel <= 672){
      L = channel / 672. * (465.9-0.59) + 0.59;
    }else if (channel <=1726){
      L = 465.9; // cm
    }else{
      L = (2399-channel) / 672. * (465.9-0.59) + 0.59;
    }
  }
  Double_t RMS = sqrt(0.9*0.9 + pow(0.79 + 0.22*L/100.,2));

  Double_t MaxPoissArg = 100.; 
  TF1 *MyPoisson = new TF1("MyPoisson", PoissonReal,0.,MaxPoissArg,1);

  // ADC to parameter conversion
  Double_t x[14]={0.5,  0.6,  0.7,  0.8,      1.0,1.3 , 1.6, 2.0,  2.5, 3.0,3.5,4.5,9,16};
  Double_t y[14]={2.43326,2.211, 2.04155,1.9036,1.72,1.55, 1.44, 1.34,1.27, 1.22,1.19,1.15,1.076,1.048};
  TGraph *g4 = new TGraph(14,y,x);

  Double_t params[1]={0};
  
  
  //to be updated ... 
  params[0] = g4->Eval(RMS);
  
  MyPoisson->SetParameters(params);
  
  //std::cout << params[0] << " " << L << std::endl;

  TF1 *f1 = new TF1("f","[0]+([1]+[2]*x*9592./2.)*exp(-[3]*pow(x*9592./2.,[4]))"); // x in MHz
  Double_t fitpar[5]={4.74264e+01,1.35048e+02,3.25547e-01,1.61809e-06,1.93925e+00};
  f1->SetParameters(fitpar);

  Double_t value_re[10000]={0.0}, value_im[10000]={0.0};
  
  Int_t n = h1->GetNbinsX();
  
  for (Int_t i=0;i!=n;i++){
    Double_t freq=0;
    if (i<n/2.){
      freq = (i)*2./n; // assume 2 MHz digitization 
    }else{
      freq = (n-i)*2./n;
    }
    
    Double_t rho = f1->Eval(freq) *sqrt(n/9592.)* MyPoisson->GetRandom()/params[0];
    //h1->SetBinContent(i+1, rho);
    if (i==0) rho = 0;
    Double_t phi = gRandom->Uniform(0,2*3.1415926);
    value_re[i] = rho*cos(phi)/n;
    value_im[i] = rho*sin(phi)/n;
  }
  
  TVirtualFFT *ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
  ifft->SetPointsComplex(value_re,value_im);
  ifft->Transform();
  TH1 *fb = 0;
  fb = TH1::TransformHisto(ifft,fb,"Re");

  for (Int_t i=0;i!=n;i++){
    h1->SetBinContent(i+1,fb->GetBinContent(i+1));
  }

  delete f1;
  delete MyPoisson;
  delete g4;
  delete ifft;
  delete fb;
}


WireCellSst::DatauBooNEFrameDataSource::DatauBooNEFrameDataSource(const TH2F *hu_raw, const TH2F *hv_raw, const TH2F *hw_raw, TTree *T_bad, TTree *T_lf, TTree *Trun, const WireCell::GeomDataSource& gds)
    : WireCell::FrameDataSource()
    , root_file(0)
    , gds(gds)
    , nwire_u(0)
    , nwire_v(0)
    , nwire_w(0)
    , flag_mis_config(0)
    , run_no(0)
    , subrun_no(0)
    , event_no(0)
    , h_rc(0)
    , hm_rc(0)
    , hp_rc(0)
    , h_1us(0)
    , hm_1us(0)
    , hp_1us(0)
    , h_2us(0)
    , hp_2us(0)
    , hm_2us(0)
{
  load_results_from_file = true;

  // fill basic information
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));
  
  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();
  
  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);

  Trun->GetEntry(0);
  
  // fill chirp_map 
  Int_t chid=0, plane=0;
  Int_t start_time=0, end_time=0;
  T_bad->SetBranchAddress("chid",&chid);
  T_bad->SetBranchAddress("plane",&plane);
  T_bad->SetBranchAddress("start_time",&start_time);
  T_bad->SetBranchAddress("end_time",&end_time);
  for (Int_t i=0;i!=T_bad->GetEntries();i++){
    T_bad->GetEntry(i);
    
    std::pair<int,int> abc(start_time,end_time);
    if (plane == 0){
      uchirp_map[chid] = abc;
    }else if (plane == 1){
      vchirp_map[chid-nwire_u] = abc;
    }else if (plane == 2){
      wchirp_map[chid-nwire_u-nwire_v] = abc;
    }
    
  }
  
  T_lf->SetBranchAddress("channel",&chid);
  for (Int_t i=0;i!=T_lf->GetEntries();i++){
    T_lf->GetEntry(i);
    lf_noisy_channels.insert(chid);
  }

  
  

  // fill frame
  if (hu_raw->GetNbinsX() != nwire_u) 
    std::cout << "U plane channels mismatched!" << std::endl;
  if (hv_raw->GetNbinsX() != nwire_v) 
    std::cout << "V plane channels mismatched!" << std::endl;
  if (hw_raw->GetNbinsX() != nwire_w) 
    std::cout << "W plane channels mismatched!" << std::endl;
  
  
  frame.index =0;
  frame.clear();		// win or lose, we start anew
  
  bins_per_frame = hu_raw->GetNbinsY();

  // U plane
  for (size_t ind=0; ind < hu_raw->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
   
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hu_raw->GetBinContent(ind+1,ibin+1);
    }
    
    frame.traces.push_back(trace);
    
  }
  
  // V plane
  for (size_t ind=0; ind < hv_raw->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind + nwire_u;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hv_raw->GetBinContent(ind+1,ibin+1);
    }
    
    frame.traces.push_back(trace);
  }

  // W plane
  for (size_t ind=0; ind < hw_raw->GetNbinsX(); ++ind) {
    WireCell::Trace trace;
    trace.chid = ind + nwire_u + nwire_v;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hw_raw->GetBinContent(ind+1,ibin+1);
    }
    
    frame.traces.push_back(trace);
  }


 
}


WireCellSst::DatauBooNEFrameDataSource::DatauBooNEFrameDataSource(const char* root_file, const WireCell::GeomDataSource& gds,int bins_per_frame1, int flag_add_noise)
    : WireCell::FrameDataSource()
    , root_file(root_file)
    , gds(gds)
    , load_results_from_file(false)
    , flag_add_noise(flag_add_noise)
    , nwire_u(0)
    , nwire_v(0)
    , nwire_w(0)
    , run_no(0)
    , subrun_no(0)
    , event_no(0)
    , h_rc(0)
    , hm_rc(0)
    , hp_rc(0)
    , h_1us(0)
    , hm_1us(0)
    , hp_1us(0)
    , h_2us(0)
    , hp_2us(0)
    , hm_2us(0)
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
  f1->SetParameters(47.*1.012,1.1); // correction 1.10 us
  f2->SetParameters(140.*1.012,2.2); // correction 2.2 us


  double x0 = h_rc->GetBinCenter(1)/2.;
  for (Int_t i=0;i!=bins_per_frame;i++){
    double x = h_rc->GetBinCenter(i+1)/2.; // convert to us, assume 2 MHz digitization
    Double_t content = -1.0/2./1000 * exp(-(x-x0)/1000.); // 1 ms RC time
    // double a = 0.5 * (3-sqrt(5.));
    // double b = 0.5 * (3+sqrt(5.));
    // Double_t content = +0.5 /sqrt(5.) / 1000. * (pow(a,2)*exp(-a*x/1000.)
    // 						 - pow(b,2)*exp(-b*x/1000.));// RC in series
    if (x==x0) content +=1;
    h_rc->SetBinContent(i+1,content);
    h_1us->SetBinContent(i+1,f1->Eval(x));
    h_2us->SetBinContent(i+1,f2->Eval(x));
  }

 
  
  hm_rc = h_rc->FFT(0,"MAG");
  hp_rc = h_rc->FFT(0,"PH");
  //std::cout << "test1 " << std::endl;
  

  hm_1us = h_1us->FFT(0,"MAG");
  hp_1us = h_1us->FFT(0,"PH");
  //std::cout << "test2 " << std::endl;

  
  hm_2us = h_2us->FFT(0,"MAG");
  hp_2us = h_2us->FFT(0,"PH");
  //std::cout << "test3 " << std::endl;
  delete f1;
  delete f2;
}

void WireCellSst::DatauBooNEFrameDataSource::Clear(){
  frame.clear();
}

WireCellSst::DatauBooNEFrameDataSource::~DatauBooNEFrameDataSource()
{
  if (!load_results_from_file){
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
}

int WireCellSst::DatauBooNEFrameDataSource::size() const
{
  //return tree->GetEntries();
  return nevents;
}

void WireCellSst::DatauBooNEFrameDataSource::Save(){

}


void WireCellSst::DatauBooNEFrameDataSource::zigzag_removal(TH1F *h1, int plane, int channel_no, int flag_RC, int flag_restore){
  TVirtualFFT *ifft=0;
  double value_re[9600]={0.0},value_im[9600]={0.0};
  
  int n = bins_per_frame;
  int nbin = bins_per_frame;
  
  TF1 *f1 = new TF1("f1","gaus");
  double par[3];

  TH1 *hm = 0;
  TH1 *hp = 0;
  TH1 *fb = 0;
  
  hm = h1->FFT(0,"MAG");
  hp = h1->FFT(0,"PH");
  

  // // need to restore the incorrectly set ASICs gain and shaping time ... 
  // int flag_restore = 0;
  // if (plane == 0){ // U-plane only    could be time-dependent 
  //   if ((channel_no >=2016 && channel_no <= 2095 
  // 	 || channel_no >=2192 && channel_no <=2303 
  // 	 || channel_no >= 2352 && channel_no <2400))
  //     {
  // 	if (flag_mis_config)
  // 	  flag_restore = 1;
  //     }
  // }

  //  std::cout << plane << " " << channel_no << " " << flag_restore << std::endl;

  for (int j=0;j!=nbin;j++){
    
    double rho = hm->GetBinContent(j+1);
    double phi = hp->GetBinContent(j+1);
    
    if (flag_RC == 1){
      // need to remove RC+RC shapings
      if (hm_rc->GetBinContent(j+1)>0){
	rho = rho/pow(hm_rc->GetBinContent(j+1),2);
	//rho = rho/pow(hm_rc->GetBinContent(j+1),1);
      }else{
	rho = 0;
      }
      phi = phi - 2*hp_rc->GetBinContent(j+1);
      //phi = phi - hp_rc->GetBinContent(j+1);
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
      if (plane == 0 || plane == 1){
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
    if (fabs(h1->GetBinContent(j+1)-mean)<6*rms)
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

  
  // std::cout << maxTicks << " " << numBins << std::endl;
  
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

      // if (plane ==0 && channel_no==1943)
      // 	std::cout << chirpFrac << " " << firstLowRMSBin << " " << lastLowRMSBin << std::endl;
      
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
  
  double ADCval=0;
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
  double ADCval=0;
  double theMean = 0.0;
  double theRMS = 0.0;
  int waveformSize = hist->GetNbinsX();
  int counter = 0;

  double min = hist->GetMinimum();
  double max = hist->GetMaximum();
  if (max > 4096) max = 4096;
  if (min > 4096) min = 4096;
  TH1F *hrmsc = new TH1F("hrmsc","hrmsc",int(max-min+1),min,max+1);
  for(int i = 0; i < waveformSize; i++)
    {
      ADCval = hist->GetBinContent(i+1);
      
      if(ADCval < 4096.0)
   	{
	  hrmsc->Fill(ADCval);
	}
    }

  double par[3];
  if (hrmsc->GetSum()>0){
    double xq = 0.5-0.34;
    hrmsc->GetQuantiles(1,&par[0],&xq);
    xq = 0.5;
    hrmsc->GetQuantiles(1,&par[1],&xq);
    xq = 0.5+0.34;
    hrmsc->GetQuantiles(1,&par[2],&xq);
    
    theRMS = sqrt((pow(par[2]-par[1],2)+pow(par[1]-par[0],2))/2.);

  }else{
    theRMS = 0;
  }
  
  delete hrmsc;

  // for(int i = 0; i < waveformSize; i++)
  //   {
  //     ADCval = hist->GetBinContent(i+1);
      
  //     if(ADCval < 4096.0)
  // 	{
  // 	  theMean += ADCval;
  // 	  theRMS += TMath::Power(ADCval,2.0);
  // 	  counter++;
  // 	}
  //   }
  
  // if(counter == 0)
  //   {
  //     theMean = 0.0;
  //     theRMS = 0.0;
  //   }
  // else
  //   {
  //     theMean /= (double)counter;
  //     theRMS /= (double)counter;
  //   theRMS = TMath::Sqrt(theRMS-TMath::Power(theMean,2.0));
  //   }
  
 

  return theRMS;

}

void WireCellSst::DatauBooNEFrameDataSource::RemoveFilterFlags(TH1F *filtHist)
{
  double ADCval=0;
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
  for (int i=0;i!=numBins;i++){
    baselineVec[i] = 0;
    isFilledVec[i] = false;
  }
  
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
  int index=0;
  double ADCval=0.0;
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
  
  int downIndex=0;
  int upIndex=0;
  bool downFlag=false;
  bool upFlag=false;
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

  //  std::cout << rmsVal << " " << planeNum << " " << channel_no << std::endl;

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
  if (flag_add_noise){
    minRMSCut[0] = 0;
    minRMSCut[1] = 0.;
    minRMSCut[2] = 0;
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

bool WireCellSst::DatauBooNEFrameDataSource::ID_lf_noisy(TH1F *h1){
  TH1F *h1_copy = (TH1F*)h1->Clone("h1_copy");
  Double_t min = h1->GetMinimum();
  Double_t max = h1->GetMaximum();
  TH1F *h2 = new TH1F("h2","h2",int(10*(max-min)+1),min,max+1);
  Double_t valid = 0;

  for (Int_t j=0;j!=h1->GetNbinsX();j++){
    if (h1->GetBinContent(j+1)!=0){
      h2->Fill(h1->GetBinContent(j+1));
      valid ++;
    }
  }
  Double_t rms, mean;
  if (h2->GetSum()>0){
    Double_t par[3]={0.0};
    double xq = 0.5;
    h2->GetQuantiles(1,&par[0],&xq);
    xq = 0.5 - 0.34;
    h2->GetQuantiles(1,&par[1],&xq);
    xq = 0.5 + 0.34;
    h2->GetQuantiles(1,&par[2],&xq);
    rms = sqrt((pow(par[2]-par[0],2) + pow(par[1]-par[0],2))/2.);
    mean = par[0];

    //std::cout << mean << " " << rms << " ";
    
    for (Int_t j=0;j!=h1->GetNbinsX();j++){
      if (fabs(h1->GetBinContent(j+1)-mean)>3.5*rms){
	for (int k=-5;k!=5;k++){
	  h1_copy->SetBinContent(k+j+1,mean);
	}
      }
    }
  }
  delete h2;
  
  Double_t content = 0;
  if (valid >0){
    TH1 *hm = h1_copy->FFT(0,"MAG");
    for (int freq=0;freq!=250;freq++){
      content += pow(hm->GetBinContent(freq+2),1);
    }
  delete hm;
  }

  // std::cout << content << " " << valid << std::endl;
  
  delete h1_copy;
  
  if (valid >0){
    // std::cout << mean << " " << rms << " " << content << " " << valid << std::endl;
    if (content/valid>14) return true;
  }
  
  return false;
}

bool WireCellSst::DatauBooNEFrameDataSource::ID_RC(TH1F *h1, int plane, int channel_no){
  bool flag = false;
  TH1 *htemp_m = h1->FFT(0,"MAG");

  Double_t content[5]={0.0};
  for (int i=0;i!=5;i++){
    content[i] = htemp_m->GetBinContent(i+2);
  }
  // if (plane ==0 && channel_no ==384){
  //   std::cout << content[0] << " " << content[1] << " " << content[2] << " " << content[3] << " " << content[4] << std::endl;
  // }

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
  
  //  frame.clear();
  //return frame_number;

  
  if (load_results_from_file) return frame_number;

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
    
    // default is misconfigued
    flag_mis_config = 1; 

    //when is it configued correctly
    if (run_no>=5112 && run_no <= 5281 ){
      flag_mis_config = 0;
    }

    if (flag_add_noise)
      flag_mis_config = 0;

    
    if (siz > 0 && frame_number < siz) {
      
      TH1F **hu=0;
      TH1F **hv=0;
      TH1F **hw=0;
      
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
      
      // initialize the response function
      Double_t hu_res_array[120]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.364382, 0.387949, 0.411053, 0.433979, 0.456863, 0.479746, 0.502641, 0.52554, 0.548441, 0.57134, 0.591765, 0.609448, 0.626848, 0.644094, 0.661364, 0.678859, 0.695231, 0.710462, 0.726147, 0.742373, 0.761332, 0.783313, 0.806325, 0.830412, 0.857676, 0.888412, 0.920705, 0.954624, 0.990242, 1.02766, 1.06121, 1.09027, 1.12037, 1.15157, 1.18392, 1.21748, 1.25229, 1.28824, 1.32509, 1.36256, 1.40051, 1.43907, 1.47857, 1.51933, 1.56134, 1.60404, 1.72665, 1.94005, 2.16994, 2.42041, 2.69475, 3.07222, 3.67375, 4.60766, 5.91864, 7.30178, 8.3715, 8.94736, 8.93705, 8.40339, 7.2212, 5.76382, 3.8931, 1.07893, -3.52481, -11.4593, -20.4011, -29.1259, -34.9544, -36.9358, -35.3303, -31.2068, -25.8614, -20.3613, -15.3794, -11.2266, -7.96091, -5.50138, -3.71143, -2.44637, -1.57662, -0.99733, -0.62554, -0.393562, -0.249715, -0.15914, -0.100771, -0.062443, -0.037283, -0.0211508, -0.0112448, -0.00552085, -0.00245133, -0.000957821, -0.000316912, -8.51679e-05, -2.21299e-05, -1.37496e-05, -1.49806e-05, -1.36935e-05, -9.66758e-06, -5.20773e-06, -7.4787e-07, 3.71199e-06, 8.17184e-06, 1.26317e-05, 1.70916e-05, 2.15514e-05, 2.60113e-05, 3.04711e-05};
      Double_t hv_res_array[120]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0865303, 0.0925559, 0.0983619, 0.104068, 0.109739, 0.115403, 0.121068, 0.126735, 0.132403, 0.138072, 0.143739, 0.149408, 0.155085, 0.160791, 0.166565, 0.172454, 0.178514, 0.184795, 0.191341, 0.198192, 0.205382, 0.212944, 0.220905, 0.229292, 0.238129, 0.247441, 0.257256, 0.267601, 0.278502, 0.28999, 0.298745, 0.304378, 0.310105, 0.315921, 0.321818, 0.327796, 0.333852, 0.339967, 0.346098, 0.352169, 0.358103, 0.363859, 0.36945, 0.374915, 0.380261, 0.385401, 0.39016, 0.394378, 0.39804, 0.401394, 0.405145, 0.410714, 0.4205, 0.437951, 0.467841, 0.516042, 0.587738, 0.694157, 0.840763, 1.01966, 1.22894, 1.5612, 2.12348, 3.31455, 5.59355, 9.10709, 14.1756, 18.4603, 19.9517, 17.4166, 10.6683, 1.40656, -10.0638, -19.034, -23.654, -24.0558, -21.4418, -17.3229, -12.9485, -9.08912, -6.05941, -3.86946, -2.38669, -1.43678, -0.853335, -0.503951, -0.296551, -0.173029, -0.0990099, -0.0547172, -0.0287882, -0.0142758, -0.00661815, -0.00284757, -0.00115702, -0.000456456, -0.000183439, -8.04214e-05, -4.20533e-05, -2.62903e-05, -1.64098e-05, -6.68039e-06, 3.04903e-06, 1.27784e-05, 2.25079e-05, 3.22373e-05, 4.19667e-05, 5.16961e-05, 6.14255e-05, 7.11549e-05};
      
      TH1F *hu_resp = new TH1F("hu_resp","hu_resp",bins_per_frame,0,bins_per_frame);
      TH1F *hv_resp = new TH1F("hv_resp","hv_resp",bins_per_frame,0,bins_per_frame);
      for (Int_t i=0;i!=120;i++){
	hu_resp->SetBinContent(i+1,hu_res_array[i]);
	hv_resp->SetBinContent(i+1,hv_res_array[i]);
      }

     
      TH1 *hmr_u = hu_resp->FFT(0,"MAG");
      TH1 *hpr_u = hu_resp->FFT(0,"PH");
      
      TH1 *hmr_v = hv_resp->FFT(0,"MAG");
      TH1 *hpr_v = hv_resp->FFT(0,"PH");
      double value_re[9600]={0.0};
      double value_im[9600]={0.0};
      
      TF1 *filter_time = new TF1("filter_time","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
      double par[2]={1.43555e+01/200.*2.,4.95096e+00};
      filter_time->SetParameters(par);
      // original one
      TF1 *filter_low = new TF1("filter_low","(1-exp(-pow(x/0.08,8)))*(x<=0.177||x>=0.18)*(x<=0.2143||x>=0.215)*(x<=0.106||x>=0.109)*(x<=0.25||x>=0.251)"); // Xiangpan's filter .... 
      TF1 *filter_low_loose = new TF1("filter_low_loose","1-exp(-pow(x/0.005,2))"); // loose filter ... 
      

      std::cout << "Load Data " << std::endl;
      // load into frame
      int nchannels = channelid->size();

      TH1F *hnoise = new TH1F("hnoise","hnoise",bins_per_frame,0,bins_per_frame);
      
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

	TH1F *htemp=0;
	float threshold=0;

	
	if (trace.chid < nwire_u){
	  htemp = hu[trace.chid];
	  threshold = 2048;
	  if (flag_add_noise){
	    Simu_Noise_uBooNE_Empirical(hnoise,0,trace.chid);
	  }
	}else if (trace.chid < nwire_u + nwire_v){
	  htemp = hv[trace.chid - nwire_u];
	  threshold = 2048;
	  if (flag_add_noise){
	    Simu_Noise_uBooNE_Empirical(hnoise,1,trace.chid-nwire_u);
	  }
	}else{
	  htemp = hw[trace.chid - nwire_u - nwire_v];
	  threshold = 400;
	  if (flag_add_noise){
	    Simu_Noise_uBooNE_Empirical(hnoise,2,trace.chid-nwire_u-nwire_v);
	  }
	}
	
	fix_ADC_shift(trace.chid,signal);
	
	if (flag_add_noise){
	  for (int ibin=0; ibin != bins_per_frame; ibin++) {
	    // pure noise to be fixed
	    htemp->SetBinContent(ibin+1,hnoise->GetBinContent(ibin+1));
	  }
	}else{
	  for (int ibin=0; ibin != bins_per_frame; ibin++) {
	    htemp->SetBinContent(ibin+1,signal->GetBinContent(ibin+1)-threshold);
	  }
	}	
	
      }
      delete hnoise;
      
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

      if (!flag_add_noise){
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
	      std::pair<int,int> abc(0, bins_per_frame-1);
	      wchirp_map[i] = abc;
	    }else{
	      wchirp_map[i].first = 0;
	      wchirp_map[i].second = bins_per_frame-1;
	    }
	  }
	}
      }


      std::vector<int> ided_rc_uplane;
      std::vector<int> ided_rc_vplane;
      std::vector<int> ided_rc_wplane;

      if (1){
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
	if (((channel_no >=2016 && channel_no <= 2095 
	      || channel_no >=2192 && channel_no <=2303 
	      || channel_no >= 2352 && channel_no <2400)&&run_no<5112)||
	    ((channel_no>=2016&&channel_no<=2111 || 
	      channel_no>=2176&&channel_no<=2303 ||
	      channel_no>=2352&&channel_no<=2383)&&run_no>=5282&&run_no<=5810) ||
	    ((channel_no>=2016&&channel_no<=2111 ||
	      channel_no>=2128&&channel_no<=2303 ||
	      channel_no>=2320&&channel_no<=2383)&&run_no>=5811&&run_no<=6699) ||
	    ((channel_no>=2240&&channel_no<=2255)&&run_no>=6700&&run_no<=6998) ||
	    ((channel_no>=2048&&channel_no<=2079||channel_no>=2240&&channel_no<=2255)&&run_no>6998)
	    ){
	  if (flag_mis_config)
	    hu[i]->Scale(14./4.7); // assume 4.7 mV/fC gain
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
      // std::cout << "Before Chirp ID " << uchirp_map.size() << " " << 
      //  	vchirp_map.size() << " " << wchirp_map.size() << std::endl;
      
      for (int i=2016; i<2400;i++){
	int channel_no = i;
	if (((channel_no >=2016 && channel_no <= 2095 
	      || channel_no >=2192 && channel_no <=2303 
	      || channel_no >= 2352 && channel_no <2400)&&run_no<5112)||
	    ((channel_no>=2016&&channel_no<=2111 || 
	      channel_no>=2176&&channel_no<=2303 ||
	      channel_no>=2352&&channel_no>=2383)&&run_no>=5282&&run_no<=5810) ||
	    ((channel_no>=2016&&channel_no<=2111 ||
	      channel_no>=2128&&channel_no<=2303 ||
	      channel_no>=2320&&channel_no<=2383)&&run_no>=5811&&run_no<=6699) ||
	    ((channel_no>=2240&&channel_no<=2255)&&run_no>=6700&&run_no<=6998) ||
	    (channel_no>=2048&&channel_no<=2079||channel_no>=2240&&channel_no<=2255)&&run_no>6998){
	  if (flag_mis_config)
	    hu[i]->Scale(4.7/14.); // assume 4.7 mV/fC gain
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

      // if (find(ided_rc_wplane.begin(),ided_rc_wplane.end(),1140)!=ided_rc_wplane.end()) std::cout <<"Xin1: " << 1140 << std::endl;

      // std::cout << ided_rc_uplane.size() << " " << ided_rc_vplane.size()
      //  		<< " " << ided_rc_wplane.size() << std::endl;

     

      if (!flag_add_noise){
	std::cout << "Remove ZigZag " << std::endl;
	// deal with the zig zag noise
	// filter single frequency 36 and 110 kHz
	// correct RC+RC 
	// correct misconfigured channel (need a database ...)
	for (int i=0;i!=nu;i++){
	  auto it = find(ided_rc_uplane.begin(),ided_rc_uplane.end(),i);

	  int flag_restore = 0;
	  int channel_no = i;
	  if (((channel_no >=2016 && channel_no <= 2095 
	      || channel_no >=2192 && channel_no <=2303 
	      || channel_no >= 2352 && channel_no <2400)&&run_no<5112)||
	    ((channel_no>=2016&&channel_no<=2111 || 
	      channel_no>=2176&&channel_no<=2303 ||
	      channel_no>=2352&&channel_no>=2383)&&run_no>=5282&&run_no<=5810) ||
	    ((channel_no>=2016&&channel_no<=2111 ||
	      channel_no>=2128&&channel_no<=2303 ||
	      channel_no>=2320&&channel_no<=2383)&&run_no>=5811&&run_no<=6699) ||
	      ((channel_no>=2240&&channel_no<=2255)&&run_no>=6700&&run_no<=6998) ||
	      (channel_no>=2048&&channel_no<=2079||channel_no>=2240&&channel_no<=2255)&&run_no>6998){
	    flag_restore = 1;
	  }

	  //std::cout << run_no << " " << flag_mis_config << " " << channel_no << " " << flag_restore << std::endl;

	  if (it == ided_rc_uplane.end()){
	    zigzag_removal(hu[i],0,i,1,flag_restore);
	  }else{
	    zigzag_removal(hu[i],0,i,0,flag_restore);
	  }
	}
	for (int i=0;i!=nv;i++){
	  auto it = find(ided_rc_vplane.begin(),ided_rc_vplane.end(),i);
	  if (it == ided_rc_vplane.end()){
	    zigzag_removal(hv[i],1,i,1);
	  }else{
	    zigzag_removal(hv[i],1,i,0);
	  }
	}
	for (int i=0;i!=nw;i++){
	  auto it = find(ided_rc_wplane.begin(),ided_rc_wplane.end(),i);
	  if (it == ided_rc_wplane.end()){
	    zigzag_removal(hw[i],2,i,1);
	  }else{
	    zigzag_removal(hw[i],2,i,0);
	  }
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

      // if (wchirp_map.find(1140)!=wchirp_map.end()) std::cout <<"Xin2: " << 1140 << std::endl;


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


      // for (int i=2016; i<2400;i++){
      // 	int channel_no = i;
      // 	if (channel_no >=2016 && channel_no <= 2095 
      // 	    || channel_no >=2192 && channel_no <=2303 
      // 	    || channel_no >= 2352 && channel_no <2400){
      // 	  if (flag_mis_config)
      // 	    hu[i]->Scale(7.8/4.7); // assume 4.7 mV/fC gain
      // 	}
      // }
      
      

      }
      //if (wchirp_map.find(1140)!=wchirp_map.end()) std::cout <<"Xin3: " << 1140 << std::endl;


      //test the code ...
      // for (int i = 0;i!=nwire_u;i++){
      // 	for (Int_t j=0;j!=120;j++){
      // 	  hu[i]->SetBinContent(2000+j,hu[i]->GetBinContent(2000+j)+hu_resp->GetBinContent(j+1));
      // 	}
      // }
      // for (int i=0;i!=nwire_v;i++){
      // 	for (Int_t j=0;j!=120;j++){
      // 	  hv[i]->SetBinContent(2000+j,hv[i]->GetBinContent(2000+j)+hv_resp->GetBinContent(j+1));
      // 	}
      // }


      if (1){
	// test for Brian ... 
	int protection_factor = 5.0;
	float min_adc_limit = 50;

	float upper_adc_limit = 15;
	float upper_decon_limit = 0.02; // same for U and V???
	float upper_decon_limit1 = 0.09; // same for U and V???
	// new cut possibility ...  
	// float upper_decon_limit = 0.02;
	
	
	int pad_window_uf = 20;
	int pad_window_ub = 10;
	
	int pad_window_vf = 10;
	int pad_window_vb = 10;
	
	int pad_window_w = 10;

	int uplane_time_shift = 79;
	int vplane_time_shift = 82;
	
	
	// deal with coherent noise removal 
	int n = bins_per_frame;
	int nbin = bins_per_frame;
	double par[3]={0.0};
	
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
	    
	    
	    // add the code to 
	    std::vector<int> signals;
	    std::vector<bool> signalsBool;
	    for (int j=0;j!=nbin;j++)
	      signalsBool.push_back(0);
	    
	    TH1 *hm = h44->FFT(0,"MAG");
	    TH1 *hp = h44->FFT(0,"PH");
	    
	    for (Int_t j=0;j!=nbin;j++){
	      double freq=0;
	      if (j < nbin/2.){
		freq = j/(1.*nbin)*2.;
	      }else{
		freq = (nbin - j)/(1.*nbin)*2.;
	      }
	      
	      double rho = hm->GetBinContent(j+1)/hmr_u->GetBinContent(j+1) *filter_time->Eval(freq)*filter_low->Eval(freq);
	      double phi = hp->GetBinContent(j+1) - hpr_u->GetBinContent(j+1);
	      
	      
	      if (j==0) rho = 0;
	      value_re[j] = rho*cos(phi)/nbin;
	      value_im[j] = rho*sin(phi)/nbin;
	    }
	    Int_t n = nbin;
	    TVirtualFFT *ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
	    ifft->SetPointsComplex(value_re,value_im);
	    ifft->Transform();
	    TH1 *fb = TH1::TransformHisto(ifft,0,"Re");
	    
	    
	    // do the RMS
	    float min = fb->GetMinimum();
	    float max = fb->GetMaximum();
	    TH1F *h55 = new TH1F("h55","h55",int((max-min)*100+1),min,max+1e-5);
	    // std::cout << max << " " << min << " " << int(max-min+1) << " " << nbin << std::endl;
	    for (int j=0;j!=nbin;j++){
	      h55->Fill(fb->GetBinContent(j+1));
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
	    if (h55->GetRMS()>0)
	      rms = h55->GetRMS();
	    //  std::cout << 0 << " " << i << " " << rms << std::endl;
	    delete h55;
	    
	    for (int j=0;j!=nbin;j++){
	      float content = fb->GetBinContent(j+1);
	      if ((content-mean)>std::max(protection_factor*rms,upper_decon_limit)){
	    	int time_bin = j + uplane_time_shift;
	    	if (time_bin >= nbin) time_bin -= nbin;
		//	h44->SetBinContent(time_bin+1,0);
	    	signalsBool.at(time_bin) = 1;
	     	// add the front and back padding
	     	for (int k=0;k!=pad_window_ub;k++){
	    	  int bin = time_bin+k+1;
	     	  if (bin > nbin-1) bin = nbin-1;
	     	  signalsBool.at(bin) = 1;
	    	}
	     	for (int k=0;k!=pad_window_uf;k++){
	     	  int bin = time_bin-k-1;
	     	  if (bin <0) bin = 0;
	     	  signalsBool.at(bin) = 1;
	    	}
	      }
	    }



	    delete ifft;
	    delete fb;
	    delete hm;
	    delete hp;

	    
	    // do the RMS
	    min = h44->GetMinimum();
	    max = h44->GetMaximum();
	    h55 = new TH1F("h55","h55",int(max-min+1),min,max+1);
	    // std::cout << max << " " << min << " " << int(max-min+1) << " " << nbin << std::endl;
	    for (int j=0;j!=nbin;j++){
	      h55->Fill(h44->GetBinContent(j+1));
	    }
	    mean = h55->GetMean();
	    rms = h55->GetRMS();
	    for (int j=0;j!=h55->GetNbinsX();j++){
	      int bin_center = h55->GetBinCenter(j+1);
	      if (bin_center < mean - 4.5*rms || bin_center > mean + 4.5*rms){
	     	h55->SetBinContent(j+1,0);
	      }
	    }
	    
	    mean = h55->GetMean();
	    if (h55->GetRMS()>0)
	      rms = h55->GetRMS();
	    
	    delete h55;
	    // remove +- 3sigma one
	    for (int j=0;j!=nbin;j++){
	      float content = h44->GetBinContent(j+1);
	      if (fabs(content-mean)>std::min(std::max(protection_factor*rms,upper_adc_limit),min_adc_limit)){
		//	    	h44->SetBinContent(j+1,0);
	    	//signals.push_back(j);
	    	signalsBool.at(j) = 1;
	    	// add the front and back padding
	    	for (int k=0;k!=pad_window_ub;k++){
	    	  int bin = j+k+1;
	    	  if (bin > nbin-1) bin = nbin-1;
	    	  signalsBool.at(bin) = 1;
	    	  //auto it = find(signals.begin(),signals.end(),bin);
	    	  //if (it == signals.end())
	    	  //signals.push_back(bin);
	    	}
	    	for (int k=0;k!=pad_window_uf;k++){
	    	  int bin = j-k-1;
	    	  if (bin <0) bin = 0;
	    	  signalsBool.at(bin) = 1;
	    	  //it = find(signals.begin(),signals.end(),bin);
	    	  //if (it == signals.end())
	    	  //signals.push_back(bin);
	    	}
	      }
	    }
	    
	    std::vector< std::vector<int> > rois_u;
	    {
	      // partition waveform indices into consecutive regions with
	      // signalsBool true.
	      
	      bool inside = false;
	      for (int ind=0; ind<nbin; ++ind) {
		if (inside) {
		  if (signalsBool[ind]) { // still inside
                    rois_u.back().push_back(ind);
		  }else{
		    inside = false;
		  }
		}
		else {                  // outside the Rio
		  if (signalsBool[ind]) { // just entered ROI
                    std::vector<int> roi;
                    roi.push_back(ind);
                    rois_u.push_back(roi);
		    inside = true;
		  }
		}
	      }
	      std::map<int, bool> flag_replace;
	      for (auto roi: rois_u){
		flag_replace[roi.front()] = false;
	      }

	      // // new deconvolution ...
	      // {
	      // 	// use ROI to get a new waveform ...
	      // 	TH1F *h44_temp = (TH1F*)h44->Clone("h44_temp");
	      // 	h44_temp->Reset();
	      // 	for (auto roi: rois){
	      // 	  const int bin0 = std::max(roi.front()-1, 0);
	      // 	  const int binf = std::min(roi.back()+1, nbin-1);
	      // 	  const double m0 = h44->GetBinContent(bin0+1);
	      // 	  const double mf = h44->GetBinContent(binf+1);
	      // 	  const double roi_run = binf - bin0;
	      // 	  const double roi_rise = mf - m0;
	      // 	  for (auto bin : roi) {
	      // 	    const double m = m0 + (bin - bin0)/roi_run*roi_rise;
	      // 	    h44_temp->SetBinContent(bin+1,h44->GetBinContent(bin+1)- m);
	      // 	  }
	      // 	}
	      // 	// do the deconvolution with a very loose low-frequency filter
	      // 	TH1 *hm = h44_temp->FFT(0,"MAG");
	      // 	TH1 *hp = h44_temp->FFT(0,"PH");
	      // 	for (Int_t j=0;j!=nbin;j++){
	      // 	  double freq=0;
	      // 	  if (j < nbin/2.){
	      // 	    freq = j/(1.*nbin)*2.;
	      // 	  }else{
	      // 	    freq = (nbin - j)/(1.*nbin)*2.;
	      // 	  }
	      // 	  double rho = hm->GetBinContent(j+1)/hmr_u->GetBinContent(j+1) *filter_time->Eval(freq)*filter_low_loose->Eval(freq);
	      // 	  double phi = hp->GetBinContent(j+1) - hpr_u->GetBinContent(j+1);
	      // 	  value_re[j] = rho*cos(phi)/nbin;
	      // 	  value_im[j] = rho*sin(phi)/nbin;
	      // 	}
	      // 	Int_t n = nbin;
	      // 	TVirtualFFT *ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
	      // 	ifft->SetPointsComplex(value_re,value_im);
	      // 	ifft->Transform();
	      // 	TH1 *fb = TH1::TransformHisto(ifft,0,"Re");
		
	      // 	for (auto roi: rois){
	      // 	  const int bin0 = std::max(roi.front()-1, 0);
	      // 	  const int binf = std::min(roi.back()+1, nbin-1);
	      // 	  flag_replace[roi.front()] = false;
		  
	      // 	  double max_val=0;
		  
		  
	      // 	  // to be modified ... 
	      // 	  for (int i=bin0; i<=binf; i++){
	      // 	    int time_bin = i-uplane_time_shift;
	      // 	    if (time_bin <0) time_bin += nbin;
	      // 	    if (time_bin >=nbin) time_bin -= nbin;
		    
	      // 	    if (i==bin0){
	      // 	      max_val = fb->GetBinContent(time_bin+1);
	      // 	      // max_adc_val = medians.at(i);
	      // 	      // min_adc_val = medians.at(i);
	      // 	    }else{
	      // 	      if (fb->GetBinContent(time_bin+1) > max_val) max_val = fb->GetBinContent(time_bin+1);
	      // 	      // if (medians.at(i) > max_adc_val) max_adc_val = medians.at(i);
	      // 	      // if (medians.at(i) < min_adc_val) min_adc_val = medians.at(i);
	      // 	    }
	      // 	  }
		  
	      // 	  //std::cout << "Xin: " << upper_decon_limit1 << " " << max_val << std::endl;
		  
	      // 	  if ( max_val > upper_decon_limit1)
	      // 	    flag_replace[roi.front()] = true;
		  
	      // 	}
		
	      // 	//judge if a roi is good or not, shrink things back properly ...
	      // 	delete h44_temp;
	      // 	delete hm;
	      // 	delete hp;
	      // 	delete fb;
	      // 	delete ifft;
	      // }
	      
	      // Replace medians for above regions with interpolation on values
	      // just outside each region.
	      for (auto roi : rois_u) {
		// original code used the bins just outside the ROI
		const int bin0 = std::max(roi.front()-1, 0);
		const int binf = std::min(roi.back()+1, nbin-1);
		if (flag_replace[roi.front()]){
		  const double m0 = h44->GetBinContent(bin0+1);//medians[bin0];
		  const double mf = h44->GetBinContent(binf+1);//medians[binf];
		  const double roi_run = binf - bin0;
		  const double roi_rise = mf - m0;
		  for (auto bin : roi) {
		    const double m = m0 + (bin - bin0)/roi_run*roi_rise;
		    h44->SetBinContent(bin+1,m);
		    //medians.at(bin) = m;
		  }
		}
	      }
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
	  

      	  
	  for (int k=0;k!=uplane_all.at(i).size();k++){
	    double scaling=0;
	    if (ave_coef!=0)
	      scaling = coef_all.at(k)/ave_coef;
	    if (scaling < 0) scaling = 0;
	    if (scaling > 1.5) scaling = 1.5;

	    
	    // new protection ... 
	    {
	      TH1F *h44_temp = (TH1F*)h44->Clone("h44_temp");
	      h44_temp->Reset();
	      
	      for (auto roi: rois_u){
		const int bin0 = std::max(roi.front()-1, 0);
		const int binf = std::min(roi.back()+1, nbin-1);
	      	  const double m0 = hu[uplane_all.at(i).at(k)]->GetBinContent(bin0+1);
	      	  const double mf = hu[uplane_all.at(i).at(k)]->GetBinContent(binf+1);
	      	  const double roi_run = binf - bin0;
	      	  const double roi_rise = mf - m0;
	      	  for (auto bin : roi) {
	      	    const double m = m0 + (bin - bin0)/roi_run*roi_rise;
	      	    h44_temp->SetBinContent(bin+1,hu[uplane_all.at(i).at(k)]->GetBinContent(bin+1) - m);
	      	  }
	      }
	      // do the deconvolution ....
	      TH1 *hm = h44_temp->FFT(0,"MAG");
	      TH1 *hp = h44_temp->FFT(0,"PH");
	      for (Int_t j=0;j!=nbin;j++){
		double freq=0;
		if (j < nbin/2.){
		  freq = j/(1.*nbin)*2.;
		}else{
		  freq = (nbin - j)/(1.*nbin)*2.;
		}
		double rho = hm->GetBinContent(j+1)/hmr_u->GetBinContent(j+1) *filter_time->Eval(freq)*filter_low_loose->Eval(freq);
		double phi = hp->GetBinContent(j+1) - hpr_u->GetBinContent(j+1);
		value_re[j] = rho*cos(phi)/nbin;
		value_im[j] = rho*sin(phi)/nbin;
	      }
	      Int_t n = nbin;
	      TVirtualFFT *ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
	      ifft->SetPointsComplex(value_re,value_im);
	      ifft->Transform();
	      TH1 *fb = TH1::TransformHisto(ifft,0,"Re");

	      std::map<int, bool> flag_replace;
	      for (auto roi: rois_u){
		flag_replace[roi.front()] = false;
	      }

	      for (auto roi: rois_u){
		const int bin0 = std::max(roi.front()-1, 0);
		const int binf = std::min(roi.back()+1, nbin-1);
	      		  
		double max_val=0;
		double min_val=0;
		
		for (int i=bin0; i<=binf; i++){
		  int time_bin = i-uplane_time_shift;
		  if (time_bin <0) time_bin += nbin;
		  if (time_bin >=nbin) time_bin -= nbin;
		    
		  if (i==bin0){
		    max_val = fb->GetBinContent(time_bin+1);
		    min_val = fb->GetBinContent(time_bin+1);
		  }else{
		    if (fb->GetBinContent(time_bin+1) > max_val) max_val = fb->GetBinContent(time_bin+1);
		    if (fb->GetBinContent(time_bin+1) < min_val) min_val = fb->GetBinContent(time_bin+1);
		  }
		}
		  
		if ( max_val > upper_decon_limit1 && fabs(min_val) < max_val*0.8)
		  flag_replace[roi.front()] = true;
	      }


	      TH1F *h44_temp1 = (TH1F*)h44->Clone("h44_temp1");
	      for (auto roi : rois_u) {
		// original code used the bins just outside the ROI
		const int bin0 = std::max(roi.front()-1, 0);
		const int binf = std::min(roi.back()+1, nbin-1);
		if (flag_replace[roi.front()]){
		  const double m0 = h44_temp1->GetBinContent(bin0+1);
		  const double mf = h44_temp1->GetBinContent(binf+1);
		  const double roi_run = binf - bin0;
		  const double roi_rise = mf - m0;
		  for (auto bin : roi) {
		    const double m = m0 + (bin - bin0)/roi_run*roi_rise;
		    h44_temp1->SetBinContent(bin+1,m);
		  }
		}
	      }
	      

	      for (int j=0;j!=nbin;j++){
		if (fabs(hu[uplane_all.at(i).at(k)]->GetBinContent(j+1))>0.001)
		  hu[uplane_all.at(i).at(k)]->SetBinContent(j+1,hu[uplane_all.at(i).at(k)]->GetBinContent(j+1)-h44_temp1->GetBinContent(j+1) * scaling);
	      }
	      
	      
	      delete h44_temp;
	      delete h44_temp1;
	      delete hm;
	      delete hp;
	      delete fb;
	      delete ifft;
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
	  
	  // add the code to 
	  std::vector<int> signals;
	  std::vector<bool> signalsBool;
	  for (int j=0;j!=nbin;j++)
	    signalsBool.push_back(0);
	  
	  TH1 *hm = h44->FFT(0,"MAG");
	  TH1 *hp = h44->FFT(0,"PH");
	  
	  for (Int_t j=0;j!=nbin;j++){
	    double freq=0;
	    if (j < nbin/2.){
	      freq = j/(1.*nbin)*2.;
	    }else{
	      freq = (nbin - j)/(1.*nbin)*2.;
	    }
	    
	    double rho = hm->GetBinContent(j+1)/hmr_v->GetBinContent(j+1) *filter_time->Eval(freq)*filter_low->Eval(freq);
	    double phi = hp->GetBinContent(j+1) - hpr_v->GetBinContent(j+1);
	    
	    //   if(freq < 0.03) rho = 0;
	    if (j==0) rho = 0;
	    value_re[j] = rho*cos(phi)/nbin;
	    value_im[j] = rho*sin(phi)/nbin;
	  }
	  Int_t n = nbin;
	  TVirtualFFT *ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
	  ifft->SetPointsComplex(value_re,value_im);
	  ifft->Transform();
	  TH1 *fb = TH1::TransformHisto(ifft,0,"Re");
	    
	    
	  // do the RMS
	  float min = fb->GetMinimum();
	  float max = fb->GetMaximum();
	  TH1F *h55 = new TH1F("h55","h55",int((max-min)*100+1),min,max+1e-5);
	  // std::cout << max << " " << min << " " << int(max-min+1) << " " << nbin << std::endl;
	  for (int j=0;j!=nbin;j++){
	    h55->Fill(fb->GetBinContent(j+1));
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
	  if (h55->GetRMS()>0)
	    rms = h55->GetRMS();
	  //  std::cout << 0 << " " << i << " " << rms << std::endl;
	  delete h55;
	    
	  for (int j=0;j!=nbin;j++){
	    float content = fb->GetBinContent(j+1);
	    if ((content-mean)>std::max(protection_factor*rms,upper_decon_limit)){
	      int time_bin = j + vplane_time_shift;
	      if (time_bin >= nbin) time_bin -= nbin;
	      //h44->SetBinContent(time_bin+1,0);
	      signalsBool.at(time_bin) = 1;
	      // add the front and back padding
	      for (int k=0;k!=pad_window_vb;k++){
		int bin = time_bin+k+1;
		if (bin > nbin-1) bin = nbin-1;
		signalsBool.at(bin) = 1;
	      }
	      for (int k=0;k!=pad_window_vf;k++){
		int bin = time_bin-k-1;
		if (bin <0) bin = 0;
		signalsBool.at(bin) = 1;
	      }
	    }
	  }
	  
	  delete ifft;
	  delete fb;
	  delete hm;
	  delete hp;
	  
      	   // do the RMS
      	  min = h44->GetMinimum();
      	  max = h44->GetMaximum();
      	  h55 = new TH1F("h55","h55",int(max-min+1),min,max+1);
      	  //std::cout << max << " " << min << " " << int(max-min+1) << std::endl;
      	  for (int j=0;j!=nbin;j++){
      	    h55->Fill(h44->GetBinContent(j+1));
      	  }
	  mean = h55->GetMean();
	  rms = h55->GetRMS();
      	  for (int j=0;j!=h55->GetNbinsX();j++){
      	    int bin_center = h55->GetBinCenter(j+1);
      	    if (bin_center < mean - 4.5*rms || bin_center > mean + 4.5*rms){
      	      h55->SetBinContent(j+1,0);
      	    }
      	  }
	  
      	  mean = h55->GetMean();
	  if (h55->GetRMS()>0)
	    rms = h55->GetRMS();
	  
	  // std::cout << 1 << " " << i << " " << rms << std::endl;

      	  delete h55;
      	  // remove +- 3sigma one
      	  
      	  for (int j=0;j!=nbin;j++){
      	    float content = h44->GetBinContent(j+1);
      	    if (fabs(content-mean)>std::min(std::max(protection_factor*rms,upper_adc_limit),min_adc_limit)){
	      //      	      h44->SetBinContent(j+1,0);
      	      //signals.push_back(j);
      	      signalsBool.at(j) = 1;
	    
      	      // add the front and back padding
      	      for (int k=0;k!=pad_window_vb;k++){
      	        int bin = j+k+1;
      	        if (bin > nbin-1) bin = nbin-1;
      	        signalsBool.at(bin) = 1;
      	        //auto it = find(signals.begin(),signals.end(),bin);
      	        //if (it == signals.end())
      	  	  //signals.push_back(bin);
	      }
	      for (int k=0;k!=pad_window_vf;k++){
      	        int bin = j-k-1;
      	        if (bin <0) bin = 0;
      	        signalsBool.at(bin) = 1;
      	        //it = find(signals.begin(),signals.end(),bin);
      	        //if (it == signals.end())
      	  	//signals.push_back(bin);
      	      }
      	    }
      	  }



	  

	  
	  std::vector< std::vector<int> > rois_v;
	  {
	    // partition waveform indices into consecutive regions with
	    // signalsBool true.
	    
	    bool inside = false;
	    for (int ind=0; ind<nbin; ++ind) {
	      if (inside) {
		if (signalsBool[ind]) { // still inside
		  rois_v.back().push_back(ind);
		}else{
		  inside = false;
		}
	      }
	      else {                  // outside the Rio
		if (signalsBool[ind]) { // just entered ROI
		  std::vector<int> roi;
                    roi.push_back(ind);
                    rois_v.push_back(roi);
		    inside = true;
		}
	      }
	    }

	    std::map<int, bool> flag_replace;
	    
	    for (auto roi: rois_v){
	      flag_replace[roi.front()] = false;
	    }
	    
	  
	    
	    // Replace medians for above regions with interpolation on values
	    // just outside each region.
	    for (auto roi : rois_v) {
	      // original code used the bins just outside the ROI
	      const int bin0 = std::max(roi.front()-1, 0);
	      const int binf = std::min(roi.back()+1, nbin-1);
	      if (flag_replace[roi.front()]){
		const double m0 = h44->GetBinContent(bin0+1);//medians[bin0];
		const double mf = h44->GetBinContent(binf+1);//medians[binf];
		const double roi_run = binf - bin0;
		const double roi_rise = mf - m0;
		for (auto bin : roi) {
		  const double m = m0 + (bin - bin0)/roi_run*roi_rise;
		  h44->SetBinContent(bin+1,m);
		  //medians.at(bin) = m;
		}
	      }
	    }
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

      	 
	  for (int k=0;k!=vplane_all.at(i).size();k++){
	    double scaling=0;
	    if (ave_coef!=0)
	      scaling = coef_all.at(k)/ave_coef;
	    if (scaling < 0) scaling = 0;
	    if (scaling > 1.5) scaling = 1.5;



	    // new protection ... 
	    {
	      TH1F *h44_temp = (TH1F*)h44->Clone("h44_temp");
	      h44_temp->Reset();
	      
	      for (auto roi: rois_v){
		const int bin0 = std::max(roi.front()-1, 0);
		const int binf = std::min(roi.back()+1, nbin-1);
		const double m0 = hv[vplane_all.at(i).at(k)]->GetBinContent(bin0+1);
		const double mf = hv[vplane_all.at(i).at(k)]->GetBinContent(binf+1);
		const double roi_run = binf - bin0;
		const double roi_rise = mf - m0;
		for (auto bin : roi) {
		  const double m = m0 + (bin - bin0)/roi_run*roi_rise;
		  h44_temp->SetBinContent(bin+1,hv[vplane_all.at(i).at(k)]->GetBinContent(bin+1) - m);
		}
	      }
	      // do the deconvolution ....
	      TH1 *hm = h44_temp->FFT(0,"MAG");
	      TH1 *hp = h44_temp->FFT(0,"PH");
	      for (Int_t j=0;j!=nbin;j++){
		double freq=0;
		if (j < nbin/2.){
		  freq = j/(1.*nbin)*2.;
		}else{
		  freq = (nbin - j)/(1.*nbin)*2.;
		}
		double rho = hm->GetBinContent(j+1)/hmr_v->GetBinContent(j+1) *filter_time->Eval(freq)*filter_low_loose->Eval(freq);
		double phi = hp->GetBinContent(j+1) - hpr_v->GetBinContent(j+1);
		value_re[j] = rho*cos(phi)/nbin;
		value_im[j] = rho*sin(phi)/nbin;
	      }
	      Int_t n = nbin;
	      TVirtualFFT *ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
	      ifft->SetPointsComplex(value_re,value_im);
	      ifft->Transform();
	      TH1 *fb = TH1::TransformHisto(ifft,0,"Re");

	      std::map<int, bool> flag_replace;
	      for (auto roi: rois_v){
		flag_replace[roi.front()] = false;
	      }

	      for (auto roi: rois_v){
		const int bin0 = std::max(roi.front()-1, 0);
		const int binf = std::min(roi.back()+1, nbin-1);
	      		  
		double max_val=0;
		double min_val = 0;
		for (int i=bin0; i<=binf; i++){
		  int time_bin = i-vplane_time_shift;
		  if (time_bin <0) time_bin += nbin;
		  if (time_bin >=nbin) time_bin -= nbin;
		  
		  if (i==bin0){
		    max_val = fb->GetBinContent(time_bin+1);
		    min_val = fb->GetBinContent(time_bin+1);
		  }else{
		    if (fb->GetBinContent(time_bin+1) > max_val) max_val = fb->GetBinContent(time_bin+1);
		    if (fb->GetBinContent(time_bin+1) < min_val) min_val = fb->GetBinContent(time_bin+1);
		    
		  }
		}
		
		if ( max_val > upper_decon_limit1 && fabs(min_val) < 0.8*max_val)
		  flag_replace[roi.front()] = true;
	      }

	      
	      TH1F *h44_temp1 = (TH1F*)h44->Clone("h44_temp1");
	      for (auto roi : rois_v) {
		// original code used the bins just outside the ROI
		const int bin0 = std::max(roi.front()-1, 0);
		const int binf = std::min(roi.back()+1, nbin-1);
		if (flag_replace[roi.front()]){
		  const double m0 = h44_temp1->GetBinContent(bin0+1);
		  const double mf = h44_temp1->GetBinContent(binf+1);
		  const double roi_run = binf - bin0;
		  const double roi_rise = mf - m0;
		  for (auto bin : roi) {
		    const double m = m0 + (bin - bin0)/roi_run*roi_rise;
		    h44_temp1->SetBinContent(bin+1,m);
		  }
		}
	      }
	      
	      for (int j=0;j!=nbin;j++){
		if (fabs(hv[vplane_all.at(i).at(k)]->GetBinContent(j+1))>0.001)
		  hv[vplane_all.at(i).at(k)]->SetBinContent(j+1,hv[vplane_all.at(i).at(k)]->GetBinContent(j+1)-h44_temp1->GetBinContent(j+1)* scaling);
	      }
	      
	      
	      
	      delete h44_temp;
	      delete h44_temp1;
	      delete hm;
	      delete hp;
	      delete fb;
	      delete ifft;
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
	
      	//std::cout << "W " << i << " " << wplane_all.at(i).at(0) << " ";
	
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
	  if (h55->GetRMS()>0)
	    rms = h55->GetRMS();

	  //std::cout << mean << " " << rms << " " << std::endl;
	  
	  //std::cout << 2 << " " << i << " " << rms << std::endl;

      	  delete h55;
      	  // remove +- 3sigma one
      	  std::vector<int> signals;
      	  std::vector<bool> signalsBool;
      	  for (int j=0;j!=nbin;j++)
      		signalsBool.push_back(0);

      	  for (int j=0;j!=nbin;j++){
      	    float content = h44->GetBinContent(j+1);
      	    if (fabs(content-mean)>std::min(protection_factor*rms,min_adc_limit)){
	      //      	      h44->SetBinContent(j+1,0);
      	      //signals.push_back(j);
      	      signalsBool.at(j) = 1;
	    
      	      // add the front and back padding
      	      for (int k=0;k!=pad_window_w;k++){
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
      	  // for (int j=0;j!=nbin;j++)
      	  // 	if( signalsBool.at(j) == 1 )
      	  // 		signals.push_back(j);
	  
      	  // // adaptive baseline 
      	  // for (int j=0;j!=signals.size();j++){
      	  //   int bin = signals.at(j);
      	  //   int prev_bin=bin;
      	  //   int next_bin=bin;
	    
      	  //   int flag = 1;
      	  //   while(flag){
      	  //     prev_bin--;
      	  //     if (find(signals.begin(),signals.end(),prev_bin)==signals.end() || prev_bin <=0){
      	  // 	flag = 0;
      	  //     }
      	  //   }

      	  //   // prev_bin = prev_bin - pad_window;
      	  //   // if (prev_bin <0) prev_bin = 0;

      	  //   flag =1;
      	  //   while(flag){
      	  //     next_bin++;
      	  //     if (find(signals.begin(),signals.end(),next_bin)==signals.end() || next_bin >=nbin-1){
      	  // 	flag = 0;
      	  //     }
      	  //   }

      	  //   // next_bin = next_bin + pad_window;
      	  //   // if (next_bin > nbin-1) next_bin = nbin-1; 


      	  //   // if (prev_bin>0) prev_bin --;
      	  //   // if (next_bin<nbin-1) next_bin++; 

      	  //   float prev_content, next_content;
      	  //   // if (prev_bin >=4){
      	  //   //   prev_content = (h44->GetBinContent(prev_bin+1) + h44->GetBinContent(prev_bin) + h44->GetBinContent(prev_bin-1) + 
      	  //   // 		      h44->GetBinContent(prev_bin-2) + h44->GetBinContent(prev_bin-3))/5.;
      	  //   // }else{
      	  //     prev_content = h44->GetBinContent(prev_bin+1);
      	  //   // }
      	  //   // if (next_bin <= nbin-5){
      	  //   //   next_content = (h44->GetBinContent(next_bin+1) + h44->GetBinContent(next_bin+2) + h44->GetBinContent(next_bin+3)+
      	  //   // 		      h44->GetBinContent(next_bin+4) + h44->GetBinContent(next_bin+5))/5.;
      	  //   // }else{
      	  //     next_content = h44->GetBinContent(next_bin+1);
      	  //   // }
	    

      	  //   //std::cout << prev_bin << " " << bin << " " << next_bin << " " << signals.size() << std::endl;
      	  //   float content = prev_content + (bin - prev_bin)/ (next_bin - prev_bin*1.0) 
      	  //     * (next_content - prev_content);

      	  //   h44->SetBinContent(bin+1,content);
      	  // }
	  {
	      // partition waveform indices into consecutive regions with
	      // signalsBool true.
	      std::vector< std::vector<int> > rois;
	      bool inside = false;
	      for (int ind=0; ind<nbin; ++ind) {
		if (inside) {
		  if (signalsBool[ind]) { // still inside
                    rois.back().push_back(ind);
		  }else{
		    inside = false;
		  }
		}
		else {                  // outside the Rio
		  if (signalsBool[ind]) { // just entered ROI
                    std::vector<int> roi;
                    roi.push_back(ind);
                    rois.push_back(roi);
		    inside = true;
		  }
		}
	      }
	      // Replace medians for above regions with interpolation on values
	      // just outside each region.
	      for (auto roi : rois) {
		// original code used the bins just outside the ROI
		const int bin0 = std::max(roi.front()-1, 0);
		const int binf = std::min(roi.back()+1, nbin-1);
		const double m0 = h44->GetBinContent(bin0+1);//medians[bin0];
		const double mf = h44->GetBinContent(binf+1);//medians[binf];
		const double roi_run = binf - bin0;
		const double roi_rise = mf - m0;
		for (auto bin : roi) {
		  const double m = m0 + (bin - bin0)/roi_run*roi_rise;
		  h44->SetBinContent(bin+1,m);
		  //medians.at(bin) = m;
		}
	      }
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
	   
      	  
	  for (int k=0;k!=wplane_all.at(i).size();k++){
	    double scaling=0;
	    if (ave_coef!=0)
	      scaling = coef_all.at(k)/ave_coef;
	    if (scaling < 0) scaling = 0;
	    if (scaling > 1.5) scaling = 1.5;
	    for (int j=0;j!=nbin;j++){
      	      if (fabs(hw[wplane_all.at(i).at(k)]->GetBinContent(j+1))>0.001)
      		hw[wplane_all.at(i).at(k)]->SetBinContent(j+1,hw[wplane_all.at(i).at(k)]->GetBinContent(j+1)-h44->GetBinContent(j+1)* scaling);
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


      }


     


      // deal with the w plane to remove the PMT signal (negative pulse ...)
      for (int i=0;i!=nwire_w;i++){
	//std::cout << i << std::endl;
      	IDPMTSignalCollection(hw[i],wrms_map[i],i);
      }
      for (int i=0;i!=nwire_u;i++){
      	IDPMTSignalInduction(hu[i],urms_map[i],0,i);
      }
      for (int i=0;i!=nwire_v;i++){
      	IDPMTSignalInduction(hv[i],vrms_map[i],1,i);
      }
      
      for (int i=0;i!=PMT_ROIs.size();i++){
       	PMT_ROIs.at(i)->sort_wires(4);
	int flag_qx = 0;
       	if (PMT_ROIs.at(i)->get_sorted_uwires().size() > 0 && PMT_ROIs.at(i)->get_sorted_vwires().size() > 0){
	  for (int j=0;j!=PMT_ROIs.at(i)->get_sorted_uwires().size();j++){
	    if (PMT_ROIs.at(i)->get_average_uwires_peak_height(j) < 2.0* PMT_ROIs.at(i)->get_average_wwires_peak_height() &&
		PMT_ROIs.at(i)->get_max_uwires_peak_height(j) < 0.5 * PMT_ROIs.at(i)->get_max_wwires_peak_height() && 
		PMT_ROIs.at(i)->get_sorted_uwires().at(j).size() <= PMT_ROIs.at(i)->get_sorted_wwires().size()
		){
	      flag_qx = 1;
	      for (int k=0;k!=PMT_ROIs.at(i)->get_sorted_uwires().at(j).size();k++){
		RemovePMTSignal(hu[PMT_ROIs.at(i)->get_sorted_uwires().at(j).at(k)],PMT_ROIs.at(i)->get_start_bin(),PMT_ROIs.at(i)->get_end_bin());
	      }

	      // std::cout << i << " " <<  PMT_ROIs.at(i)->get_peaks().at(0) << " " << PMT_ROIs.at(i)->get_peaks().size() << " " << PMT_ROIs.at(i)->get_sorted_uwires().at(j).size() << " " << j << " " << PMT_ROIs.at(i)->get_average_uwires_peak_height(j) << " " <<   PMT_ROIs.at(i)->get_average_wwires_peak_height() << " " << PMT_ROIs.at(i)->get_max_uwires_peak_height(j) << " " <<  PMT_ROIs.at(i)->get_max_wwires_peak_height() << " " << PMT_ROIs.at(i)->get_sorted_uwires().at(j).size() << " " <<  PMT_ROIs.at(i)->get_sorted_wwires().size() << std::endl;
		
		

	    }
	  }

	  for (int j=0;j!=PMT_ROIs.at(i)->get_sorted_vwires().size();j++){
	    if (PMT_ROIs.at(i)->get_average_vwires_peak_height(j) < 2.0 * PMT_ROIs.at(i)->get_average_wwires_peak_height() && 
		PMT_ROIs.at(i)->get_max_vwires_peak_height(j) < 0.5 * PMT_ROIs.at(i)->get_max_wwires_peak_height() &&
		PMT_ROIs.at(i)->get_sorted_vwires().at(j).size() <= PMT_ROIs.at(i)->get_sorted_wwires().size() 
		){
	      flag_qx = 1;
	      for (int k=0;k!=PMT_ROIs.at(i)->get_sorted_vwires().at(j).size();k++){
		RemovePMTSignal(hv[PMT_ROIs.at(i)->get_sorted_vwires().at(j).at(k)],PMT_ROIs.at(i)->get_start_bin(),PMT_ROIs.at(i)->get_end_bin());
	      }

	      //  std::cout << i << " " <<  PMT_ROIs.at(i)->get_peaks().at(0) << " " << PMT_ROIs.at(i)->get_peaks().size() << " " << PMT_ROIs.at(i)->get_sorted_vwires().at(j).size() << " " << j << " " << PMT_ROIs.at(i)->get_average_vwires_peak_height(j) << " " <<   PMT_ROIs.at(i)->get_average_wwires_peak_height() << " " << PMT_ROIs.at(i)->get_max_vwires_peak_height(j) << " " <<  PMT_ROIs.at(i)->get_max_wwires_peak_height() << " " << PMT_ROIs.at(i)->get_sorted_vwires().at(j).size() << " " <<  PMT_ROIs.at(i)->get_sorted_wwires().size() << std::endl;
	    }
	  }

	  if (PMT_ROIs.at(i)->get_sorted_wwires().size() >=6){
	    for (int j=0;j!=PMT_ROIs.at(i)->get_sorted_wwires().size();j++){
	      //	    std::cout << PMT_ROIs.at(i)->get_sorted_wwires().size() << std::endl;
	      RemovePMTSignal(hw[PMT_ROIs.at(i)->get_sorted_wwires().at(j)],PMT_ROIs.at(i)->get_start_bin(),PMT_ROIs.at(i)->get_end_bin(),1);
	    }
	  }
	  
     
       	  //if (flag_qx == 1) std::cout << i << " " <<  PMT_ROIs.at(i)->get_peaks().at(0) << " " <<  PMT_ROIs.at(i)->get_peaks().size() << " " << PMT_ROIs.at(i)->get_uwires().size() << " " << PMT_ROIs.at(i)->get_vwires().size() << " " << PMT_ROIs.at(i)->get_sorted_uwires().size() << " " << PMT_ROIs.at(i)->get_sorted_vwires().size() << " " << PMT_ROIs.at(i)->get_sorted_wwires().size() << " " << PMT_ROIs.at(i)->get_average_uwires_peak_height() << " " << PMT_ROIs.at(i)->get_average_vwires_peak_height() << " " << PMT_ROIs.at(i)->get_average_wwires_peak_height() << " " << PMT_ROIs.at(i)->get_max_uwires_peak_height() << " " << PMT_ROIs.at(i)->get_max_vwires_peak_height() << " " << PMT_ROIs.at(i)->get_max_wwires_peak_height() << std::endl;
      	}
       	delete PMT_ROIs.at(i);
      }
      PMT_ROIs.clear();

      


      
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
      
      
      // loop through the chirping channels ... 
      // save chirping ones into the lf_noisy_channels
      for (auto it = uchirp_map.begin(); it!= uchirp_map.end(); it++){
	int channel = it->first;
	if (it->second.first!=0 || it->second.second<bins_per_frame-1){
	  lf_noisy_channels.insert(channel);
	}
      }
      for (auto it = vchirp_map.begin(); it!= vchirp_map.end(); it++){
	int channel = it->first+nwire_u;
	if (it->second.first!=0 || it->second.second<bins_per_frame-1){
	  lf_noisy_channels.insert(channel);
	}
      }
      
      for (auto it = ided_rc_uplane.begin(); it!= ided_rc_uplane.end(); it++){
	int channel = *it;
	lf_noisy_channels.insert(channel);
      }
      for (auto it = ided_rc_vplane.begin(); it!= ided_rc_vplane.end(); it++){
	int channel = *it+nwire_u;
	lf_noisy_channels.insert(channel);
      }
     

      // Try to ID  lf noisy channels
      for (int i=0;i!=nwire_u;i++){
	//std::cout << i << std::endl;
	if (ID_lf_noisy(hu[i])) lf_noisy_channels.insert(i);
	if (ID_lf_noisy(hv[i])) lf_noisy_channels.insert(nwire_u+i);
      }

      
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

      delete hmr_u;
      delete hmr_v;
      delete hpr_u;
      delete hpr_v;
      delete hu_resp;
      delete hv_resp;
      delete filter_time;
      delete filter_low;
      
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


void WireCellSst::DatauBooNEFrameDataSource::fix_ADC_shift(int chid, TH1F *h1){
  
  // judge if there is a shift
  //  std::cout << chid << std::endl;
  std::vector<int> counter(12,0);
  std::vector<int> s,s1(12,0);
  int nbin = 500;
  for (int i=0;i!=nbin;i++){
    int x = h1->GetBinContent(i+1);
    
    //if (chid ==1152) std::cout << i << " " << x << std::endl;

    s.clear();
    do
      {
	s.push_back( (x & 1));
      } while (x >>= 1);
    s.resize(12);
    
    for (int j=0;j!=12;j++){
      counter.at(j) += abs(s.at(j) - s1.at(j));
    }
    s1=s;
  }
  int nshift = 0;
  for (int i=0;i!=12;i++){
    if (counter.at(i) < nbin/2. - nbin/2. *sqrt(1./nbin) * 7.5){
      nshift ++;
    }else{
      break;
    }
  }
  // if (chid ==1152) std::cout << chid << " " << nshift << std::endl;
  
  if (nshift !=0 && nshift<11){ // difficult to do it if shift by 11 bins ... 
    std::cout << chid << " " << "ADC shifted " << nshift << " bits" << std::endl;
    

    // copy the data, and do the shift back
    const int nl = bins_per_frame;//h1->GetNbinsX();
    int x[nl], x_orig[nl];
    for (int i=0;i!=nl;i++){
      x_orig[i] = h1->GetBinContent(i+1);
      x[i] = h1->GetBinContent(i+1);
    }

    std::set<int> fillings;
    int exp_diff = pow(2,12-nshift)*0.8;
      
    int filling;// = lowest_bits(x[0], nshift);
    int mean = 0;
    for (int i=0;i!=nl;i++){
      filling = lowest_bits(x_orig[i], nshift);
      int y = shift_right(x_orig[i], nshift, filling, 12);
      fillings.insert(filling);
      // cout << y << " ";
      // filling = lowest_bits(x[i], nshift);
      // filling = x[i] & int(pow(2, nshift)-1);
      x[i] = y;
      mean += x[i];
    }
    mean = mean/nl;
    
    // if (chid == 1281){
    //   std::cout << nl << " " << nshift << " " << mean << " " << exp_diff << " " << fillings.size() << std::endl;
    //   for (auto it = fillings.begin(); it!=fillings.end();it++){
    // 	std::cout << (*it) << std::endl;
    //   }
    // }

    // examine the results ... 
    int prev1_bin_content = mean;
    int prev_bin_content = mean;
    int next_bin_content = mean;
    int next1_bin_content = mean;
    int curr_bin_content;
    
    for (Int_t i=0;i<nl;i=i+1){
      curr_bin_content = x[i];
      // when to judge if one is likely bad ... 
      if (abs(curr_bin_content-mean)>exp_diff && 
	  (abs(curr_bin_content - prev_bin_content) > exp_diff ||
	   abs(curr_bin_content - next_bin_content) > exp_diff)
	  ){
	Int_t exp_value = ( (2*prev_bin_content - prev1_bin_content) +
			    (prev_bin_content + next_bin_content)/2. + 
			    (prev_bin_content * 2./3. + next1_bin_content/3.))/3.;
	for (auto it = fillings.begin(); it!=fillings.end(); it++){
	  int y = shift_right(x_orig[i], nshift, (*it), 12);
	  // when to switch ... 
	  if (fabs(y-exp_value) < fabs(x[i] - exp_value)){
	    x[i] = y;//hs->SetBinContent(i+1,y);
	  }
	}
      }
      // if (chid == 1281 && i < 10){
      // 	std::cout <<x_orig[i] << " " << x[i] << std::endl;
      // }
      //hs->SetBinContent(i+1,mean);
      
      prev1_bin_content = prev_bin_content;
      prev_bin_content = x[i];
      if (i+2 < nl){
	next_bin_content = x[i+2];
      }else{
	next_bin_content = mean;
      }
      if (i+3 < nl){
	next1_bin_content = x[i+3];
      }else{
	next1_bin_content = mean;
      }
    }


    // change the histogram ...
    for (int i=0;i!=nl;i++){
      h1->SetBinContent(i+1,x[i]);
    }
  }
}

int WireCellSst::DatauBooNEFrameDataSource::shift_right(int value, int n, int filling, int totalBit)
{
  return (value >> n) | (filling << (totalBit-n));
}

int  WireCellSst::DatauBooNEFrameDataSource::lowest_bits(int value, int n)
{
  return ( value & ((1 << n) - 1) );
}

int WireCellSst::DatauBooNEFrameDataSource::jump_no_noise(int frame_number)
{
  
  //  frame.clear();
  //return frame_number;

  
  if (load_results_from_file) return frame_number;

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
    
    // default is non-misconfigued
    flag_mis_config = 0; 

    if (siz > 0 && frame_number < siz) {
      
      TH1F **hu=0;
      TH1F **hv=0;
      TH1F **hw=0;
      
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
      
      
      // initialize the response function
      Double_t hu_res_array[120]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.364382, 0.387949, 0.411053, 0.433979, 0.456863, 0.479746, 0.502641, 0.52554, 0.548441, 0.57134, 0.591765, 0.609448, 0.626848, 0.644094, 0.661364, 0.678859, 0.695231, 0.710462, 0.726147, 0.742373, 0.761332, 0.783313, 0.806325, 0.830412, 0.857676, 0.888412, 0.920705, 0.954624, 0.990242, 1.02766, 1.06121, 1.09027, 1.12037, 1.15157, 1.18392, 1.21748, 1.25229, 1.28824, 1.32509, 1.36256, 1.40051, 1.43907, 1.47857, 1.51933, 1.56134, 1.60404, 1.72665, 1.94005, 2.16994, 2.42041, 2.69475, 3.07222, 3.67375, 4.60766, 5.91864, 7.30178, 8.3715, 8.94736, 8.93705, 8.40339, 7.2212, 5.76382, 3.8931, 1.07893, -3.52481, -11.4593, -20.4011, -29.1259, -34.9544, -36.9358, -35.3303, -31.2068, -25.8614, -20.3613, -15.3794, -11.2266, -7.96091, -5.50138, -3.71143, -2.44637, -1.57662, -0.99733, -0.62554, -0.393562, -0.249715, -0.15914, -0.100771, -0.062443, -0.037283, -0.0211508, -0.0112448, -0.00552085, -0.00245133, -0.000957821, -0.000316912, -8.51679e-05, -2.21299e-05, -1.37496e-05, -1.49806e-05, -1.36935e-05, -9.66758e-06, -5.20773e-06, -7.4787e-07, 3.71199e-06, 8.17184e-06, 1.26317e-05, 1.70916e-05, 2.15514e-05, 2.60113e-05, 3.04711e-05};
      Double_t hv_res_array[120]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0865303, 0.0925559, 0.0983619, 0.104068, 0.109739, 0.115403, 0.121068, 0.126735, 0.132403, 0.138072, 0.143739, 0.149408, 0.155085, 0.160791, 0.166565, 0.172454, 0.178514, 0.184795, 0.191341, 0.198192, 0.205382, 0.212944, 0.220905, 0.229292, 0.238129, 0.247441, 0.257256, 0.267601, 0.278502, 0.28999, 0.298745, 0.304378, 0.310105, 0.315921, 0.321818, 0.327796, 0.333852, 0.339967, 0.346098, 0.352169, 0.358103, 0.363859, 0.36945, 0.374915, 0.380261, 0.385401, 0.39016, 0.394378, 0.39804, 0.401394, 0.405145, 0.410714, 0.4205, 0.437951, 0.467841, 0.516042, 0.587738, 0.694157, 0.840763, 1.01966, 1.22894, 1.5612, 2.12348, 3.31455, 5.59355, 9.10709, 14.1756, 18.4603, 19.9517, 17.4166, 10.6683, 1.40656, -10.0638, -19.034, -23.654, -24.0558, -21.4418, -17.3229, -12.9485, -9.08912, -6.05941, -3.86946, -2.38669, -1.43678, -0.853335, -0.503951, -0.296551, -0.173029, -0.0990099, -0.0547172, -0.0287882, -0.0142758, -0.00661815, -0.00284757, -0.00115702, -0.000456456, -0.000183439, -8.04214e-05, -4.20533e-05, -2.62903e-05, -1.64098e-05, -6.68039e-06, 3.04903e-06, 1.27784e-05, 2.25079e-05, 3.22373e-05, 4.19667e-05, 5.16961e-05, 6.14255e-05, 7.11549e-05};
      
      TH1F *hu_resp = new TH1F("hu_resp","hu_resp",bins_per_frame,0,bins_per_frame);
      TH1F *hv_resp = new TH1F("hv_resp","hv_resp",bins_per_frame,0,bins_per_frame);
      for (Int_t i=0;i!=120;i++){
	hu_resp->SetBinContent(i+1,hu_res_array[i]);
	hv_resp->SetBinContent(i+1,hv_res_array[i]);
      }

     
      TH1 *hmr_u = hu_resp->FFT(0,"MAG");
      TH1 *hpr_u = hu_resp->FFT(0,"PH");
      
      TH1 *hmr_v = hv_resp->FFT(0,"MAG");
      TH1 *hpr_v = hv_resp->FFT(0,"PH");
      double value_re[9600]={0.0};
      double value_im[9600]={0.0};
      
      TF1 *filter_time = new TF1("filter_time","(x>0.0)*exp(-0.5*pow(x/[0],[1]))");
      double par[2]={1.43555e+01/200.*2.,4.95096e+00};
      filter_time->SetParameters(par);
      TF1 *filter_low = new TF1("filter_low","1-exp(-pow(x/0.06,2))");
      
      

      std::cout << "Load Data " << std::endl;
      // load into frame
      int nchannels = channelid->size();
      
      for (size_t ind=0; ind < nchannels; ++ind) {
	TH1F* signal = dynamic_cast<TH1F*>(esignal->At(ind));
	if (!signal) continue;
	
	
	WireCell::Trace trace;
	trace.chid = channelid->at(ind);
	
	//	trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
	//trace.charge.resize(bins_per_frame, 0.0);

	TH1F *htemp=0;
	float threshold=0;

	
	if (trace.chid < nwire_u){
	  htemp = hu[trace.chid];
	}else if (trace.chid < nwire_u + nwire_v){
	  htemp = hv[trace.chid - nwire_u];
	}else{
	  htemp = hw[trace.chid - nwire_u - nwire_v];
	}
	
	for (int ibin=0; ibin != bins_per_frame; ibin++) {
	  htemp->SetBinContent(ibin+1,signal->GetBinContent(ibin+1));
	}
	
	//	fix_ADC_shift(htemp);

      }
      
      int nu = 2400, nv = 2400, nw = 3456;
      int ntotal = nu + nv + nw;
      
      if (0){

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


      std::cout << "Adaptive Baseline " << uchirp_map.size() << " " << 
      	vchirp_map.size() << " " << wchirp_map.size() << std::endl;
      // do the adaptive baseline ... 


      // do the adaptive baseline for the bad RC channels ... 
      
      
      std::cout << "Noisy Channel " << std::endl;
      // deal with the noisy signal, and put them into chirping map 

      }

      
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

      delete hmr_u;
      delete hmr_v;
      delete hpr_u;
      delete hpr_v;
      delete hu_resp;
      delete hv_resp;
      delete filter_time;
      delete filter_low;
      
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

void WireCellSst::DatauBooNEFrameDataSource::IDPMTSignalInduction(TH1F* hist, float rms, int plane, int channel){
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
    for (int i=0;i!=PMT_ROIs.size();i++){
      WireCell::PMTNoiseROI* ROI= PMT_ROIs.at(i);
      for (int j=0;j!= ROI->get_peaks().size();j++){
	int peak = ROI->get_peaks().at(j);
	int peak_m1 = peak - 1; if (peak_m1 <0) peak_m1 = 0;
	int peak_m2 = peak - 2; if (peak_m2 <0) peak_m2 = 0;
	int peak_m3 = peak - 3; if (peak_m3 <0) peak_m3 = 0;
	int peak_p1 = peak + 1; if (peak_p1 >= hist->GetNbinsX()) peak_p1 = hist->GetNbinsX()-1;
	int peak_p2 = peak + 2; if (peak_p2 >= hist->GetNbinsX()) peak_p2 = hist->GetNbinsX()-1;
	int peak_p3 = peak + 3; if (peak_p3 >= hist->GetNbinsX()) peak_p3 = hist->GetNbinsX()-1;
	if (fabs(hist->GetBinContent(peak+1))> 5 * rms1 && 
	    fabs(hist->GetBinContent(peak+1)) + fabs(hist->GetBinContent(peak_m1+1)) + fabs(hist->GetBinContent(peak_p1+1)) > fabs(hist->GetBinContent(peak_m1+1)) + fabs(hist->GetBinContent(peak_m2+1)) + fabs(hist->GetBinContent(peak+1)) &&
	    fabs(hist->GetBinContent(peak+1)) + fabs(hist->GetBinContent(peak_m1+1)) + fabs(hist->GetBinContent(peak_p1+1)) > fabs(hist->GetBinContent(peak_p1+1)) + fabs(hist->GetBinContent(peak_p2+1)) + fabs(hist->GetBinContent(peak+1)) && 
	    fabs(hist->GetBinContent(peak+1)) + fabs(hist->GetBinContent(peak_m1+1)) + fabs(hist->GetBinContent(peak_p1+1)) > fabs(hist->GetBinContent(peak_m2+1)) + fabs(hist->GetBinContent(peak_m3+1)) + fabs(hist->GetBinContent(peak_m1+1)) && 
	    fabs(hist->GetBinContent(peak+1)) + fabs(hist->GetBinContent(peak_m1+1)) + fabs(hist->GetBinContent(peak_p1+1)) > fabs(hist->GetBinContent(peak_p2+1)) + fabs(hist->GetBinContent(peak_p3+1)) + fabs(hist->GetBinContent(peak_p1+1)) ){

	  //	    fabs(hist->GetBinContent(peak+1))>= fabs(hist->GetBinContent(peak_m2+1)) && 
	  // fabs(hist->GetBinContent(peak+1))>= fabs(hist->GetBinContent(peak_p1+1)) && 
	  // fabs(hist->GetBinContent(peak+1))>= fabs(hist->GetBinContent(peak_p2+1)) && 
	  // fabs(hist->GetBinContent(peak+1))>= fabs(hist->GetBinContent(peak_p3+1)) && 
	  // fabs(hist->GetBinContent(peak+1))>= fabs(hist->GetBinContent(peak_m3+1)) 
	  // ){
	  if (plane == 0 ){
	    ROI->insert_uwires(channel,fabs(hist->GetBinContent(peak+1)));
	    break;
	  }else{
	    ROI->insert_vwires(channel,fabs(hist->GetBinContent(peak+1)));
	    break;
	  }
	}
      }
    }
  }
}

void WireCellSst::DatauBooNEFrameDataSource::RemovePMTSignal(TH1F* hist, int start_bin, int end_bin, int flag){

  //  return;

  
  int pad_window = 5;
  
  int flag_start = 0;
	    
  //adaptive baseline
  float start_content = hist->GetBinContent(start_bin+1);
  for (int j=start_bin;j>=start_bin - pad_window;j--){
    if (j<0) continue;
    if (fabs(hist->GetBinContent(j+1)) < fabs(start_content)){
      start_bin = j;
      start_content = hist->GetBinContent(start_bin+1);
    }
  }
  float end_content = hist->GetBinContent(end_bin+1);
  for (int j=end_bin; j<=end_bin+pad_window;j++){
    if (j>=hist->GetNbinsX()) continue;
    if (fabs(hist->GetBinContent(j+1)) < fabs(end_content)){
      end_bin = j;
      end_content = hist->GetBinContent(end_bin+1);
    }
  }
  
  for (int j=start_bin;j<=end_bin;j++){
    float content = start_content + (end_content - start_content) * (j - start_bin) / (end_bin - start_bin*1.0);
    if (flag==0){
      hist->SetBinContent(j+1,content);
    }else{
      if (hist->GetBinContent(j+1)<0)
	hist->SetBinContent(j+1,content);
    }
  }
}


void WireCellSst::DatauBooNEFrameDataSource::IDPMTSignalCollection(TH1F* hist, float rms, int channel){

  int pad_window = 5;
  int min_window_length = 4;
  
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
    int flag_start = 0;
    int start_bin=0;
    int end_bin=0;
    int peak_bin=0;

    //std::cout << channel << " " << rms1 << std::endl;
    
    for (int i=0;i!=hist->GetNbinsX();i++){
      float content = hist->GetBinContent(i+1);

      if (flag_start == 0){
	if (content < -5 *rms1){
	  start_bin = i;
	  flag_start = 1;
	}	
      }else{
	if (content >= -5 * rms1){
	  end_bin = i-1;
	  if (end_bin > start_bin+min_window_length){
	    float min = hist->GetBinContent(start_bin+1);
	    peak_bin = start_bin;
	    for (int j=start_bin+1;j!=end_bin;j++){
	      if (hist->GetBinContent(j+1) < min)
		peak_bin = j;
	    }
	  	 
	    WireCell::PMTNoiseROI *ROI = new WireCell::PMTNoiseROI(start_bin,end_bin,peak_bin,channel,hist->GetBinContent(peak_bin+1));
	    
	    //std::cout << start_bin << " " << end_bin << std::endl;
	    if (PMT_ROIs.size()==0){
	      PMT_ROIs.push_back(ROI);
	    }else{
	      bool flag_merge = false;
	      for (int i=0;i!=PMT_ROIs.size();i++){
		flag_merge = PMT_ROIs.at(i)->merge_ROI(*ROI);
		if (flag_merge){
		  delete ROI;
		  break;
		}
	      }
	      if (!flag_merge){
		PMT_ROIs.push_back(ROI);
	      }
	    }
	    // std::cout << "h " << PMT_ROIs.size() << std::endl;
	    
	    flag_start = 0;
	    
	    // //adaptive baseline
	    // float start_content = hist->GetBinContent(start_bin+1);
	    // for (int j=start_bin;j>=start_bin - pad_window;j--){
	    //   if (j<0) continue;
	    //   if (fabs(hist->GetBinContent(j+1)) < fabs(start_content)){
	    // 	start_bin = j;
	    // 	start_content = hist->GetBinContent(start_bin+1);
	    //   }
	    // }
	    // float end_content = hist->GetBinContent(end_bin+1);
	    // for (int j=end_bin; j<=end_bin+pad_window;j++){
	    //   if (j>=hist->GetNbinsX()) continue;
	    //   if (fabs(hist->GetBinContent(j+1)) < fabs(end_content)){
	    // 	end_bin = j;
	    // 	end_content = hist->GetBinContent(end_bin+1);
	    //   }
	    // }
	    // for (int j=start_bin;j<=end_bin;j++){
	    //   float content = start_content + (end_content - start_content) * (j - start_bin) / (end_bin - start_bin*1.0);
	    //   hist->SetBinContent(j+1,content);
	    // }
	    
	  }
	}
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
