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
    void Rotate_Shift(float x_center=0, float y_center=0, float z_center=0, float rotate_angle=0, float x_shift=0, float y_shift=0, float z_shift=0);
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

    int mc_isnu; // is neutrino interaction
    int mc_nGeniePrimaries; // number of Genie primaries
    int mc_nu_pdg; // pdg code of neutrino
    int mc_nu_ccnc; // cc or nc
    int mc_nu_mode; // mode: http://nusoft.fnal.gov/larsoft/doxsvn/html/MCNeutrino_8h_source.html
    int mc_nu_intType; // interaction type
    int mc_nu_target; // target interaction
    int mc_hitnuc; // hit nucleon
    int mc_hitquark; // hit quark
    double mc_nu_Q2; // Q^2
    double mc_nu_W; // W
    double mc_nu_X; // X
    double mc_nu_Y; // Y
    double mc_nu_Pt; // Pt
    double mc_nu_Theta; // angle relative to lepton
    float mc_nu_pos[4];  // interaction position of nu
    float mc_nu_mom[4];  // interaction momentum of nu


  };
  
}

#endif
