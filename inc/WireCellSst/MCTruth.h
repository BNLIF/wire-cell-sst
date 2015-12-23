#ifndef WIRECELLSST_MCTRUTH_H
#define WIRECELLSST_MCTRUTH_H

#include "TClonesArray.h"
#include "TObjArray.h"
#include <vector>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TTree.h"

#define MAX_TRACKS 10000

namespace WireCellSst {

  class MCTruth {
  public:
    MCTruth(std::string rootfile);
    ~MCTruth();
    void GetEntry(int i);

  private:
    TTree *mcTree;
    
  public:
    int runNo;
    int subrunNo;
    int eventNo;
    int mc_Ntrack;  // number of tracks in MC
    int mc_id[MAX_TRACKS];  // track id; size == mc_Ntrack
    int mc_pdg[MAX_TRACKS];  // track particle pdg; size == mc_Ntrack
    int mc_process[MAX_TRACKS];  // track generation process code; size == mc_Ntrack
    int mc_mother[MAX_TRACKS];  // mother id of this track; size == mc_Ntrack
    float mc_startXYZT[MAX_TRACKS][4];  // start position of this track; size == mc_Ntrack
    float mc_endXYZT[MAX_TRACKS][4];  // end position of this track; size == mc_Ntrack
    float mc_startMomentum[MAX_TRACKS][4];  // start momentum of this track; size == mc_Ntrack
    float mc_endMomentum[MAX_TRACKS][4];  // end momentum of this track; size == mc_Ntrack
    std::vector<std::vector<int> > *mc_daughters;  // daughters id of this track; vector
    TObjArray* mc_trackPosition;    
  };
  
}

#endif
