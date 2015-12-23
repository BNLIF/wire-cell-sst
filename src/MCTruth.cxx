#include "WireCellSst/MCTruth.h"

WireCellSst::MCTruth::MCTruth(std::string rootfile)
{
  TFile *file = new TFile(rootfile.c_str(), "read");
  mcTree = (TTree*)file->Get("/Event/Sim");
  if (!mcTree ) std::cout<<"cannot find /Event/Sim tree."<<std::endl;
  mc_daughters= new std::vector<std::vector<int> >;
  mcTree->SetBranchAddress("eventNo",&eventNo);
  mcTree->SetBranchAddress("runNo",&runNo);
  mcTree->SetBranchAddress("subRunNo",&subrunNo);
  mcTree->SetBranchAddress("mc_Ntrack", &mc_Ntrack);  // number of tracks in MC
  mcTree->SetBranchAddress("mc_id", &mc_id);  // track id; size == mc_Ntrack
  mcTree->SetBranchAddress("mc_pdg", &mc_pdg);  // track particle pdg; size == mc_Ntrack
  mcTree->SetBranchAddress("mc_process", &mc_process);  // track generation process code; size == mc_Ntrack
  mcTree->SetBranchAddress("mc_mother", &mc_mother);  // mother id of this track; size == mc_Ntrack
  mcTree->SetBranchAddress("mc_daughters", &mc_daughters);  // daughters id of this track; vector
  mcTree->SetBranchAddress("mc_startXYZT", &mc_startXYZT);  // start position of this track; size == mc_Ntrack
  mcTree->SetBranchAddress("mc_endXYZT", &mc_endXYZT);  // start position of this track; size == mc_Ntrack
  mcTree->SetBranchAddress("mc_startMomentum", &mc_startMomentum);  // start momentum of this track; size == mc_Ntrack
  mcTree->SetBranchAddress("mc_endMomentum", &mc_endMomentum);  // start momentum of this track; size == mc_Ntrack
  if (mcTree->GetBranch("mc_trackPosition")) {
    mc_trackPosition = new TObjArray();
    mcTree->SetBranchAddress("mc_trackPosition",mc_trackPosition);
  }
}

WireCellSst::MCTruth::~MCTruth()
{
  if (mc_daughters != NULL)
    delete mc_daughters;
  if (mc_trackPosition != NULL)
    delete mc_trackPosition;
}

void WireCellSst::MCTruth::GetEntry(int i)
{
  if (i < mcTree->GetEntries()) {
    mcTree->GetEntry(i);
  } else {
    std::cout<<i<<" exceeds the number of entries, which is "<<mcTree->GetEntries()<<std::endl;
  }
}
