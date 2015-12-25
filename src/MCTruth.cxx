#include "WireCellSst/MCTruth.h"
#include "TLorentzVector.h"

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

  mcTree->SetBranchAddress("mc_isnu", &mc_isnu);
  mcTree->SetBranchAddress("mc_nGeniePrimaries", &mc_nGeniePrimaries);
  mcTree->SetBranchAddress("mc_nu_pdg", &mc_nu_pdg);
  mcTree->SetBranchAddress("mc_nu_ccnc", &mc_nu_ccnc);
  mcTree->SetBranchAddress("mc_nu_mode", &mc_nu_mode);
  mcTree->SetBranchAddress("mc_nu_intType", &mc_nu_intType);
  mcTree->SetBranchAddress("mc_nu_target", &mc_nu_target);
  mcTree->SetBranchAddress("mc_hitnuc", &mc_hitnuc);
  mcTree->SetBranchAddress("mc_hitquark", &mc_hitquark);
  mcTree->SetBranchAddress("mc_nu_Q2", &mc_nu_Q2);
  mcTree->SetBranchAddress("mc_nu_W", &mc_nu_W);
  mcTree->SetBranchAddress("mc_nu_X", &mc_nu_X);
  mcTree->SetBranchAddress("mc_nu_Y", &mc_nu_Y);
  mcTree->SetBranchAddress("mc_nu_Pt", &mc_nu_Pt);
  mcTree->SetBranchAddress("mc_nu_Theta", &mc_nu_Theta);
  mcTree->SetBranchAddress("mc_nu_pos", mc_nu_pos);
  mcTree->SetBranchAddress("mc_nu_mom", mc_nu_mom);




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

void WireCellSst::MCTruth::Rotate_Shift(float x_center, float y_center, float z_center, float rotate_angle, float x_shift, float y_shift, float z_shift){

  Double_t temp_x, temp_y, temp_z;
  for (int i=0;i!=mc_Ntrack;i++){
    //mc_startXYZT
    temp_x = mc_startXYZT[i][0];
    temp_y = mc_startXYZT[i][1];
    temp_z = mc_startXYZT[i][2];
    mc_startXYZT[i][0] = (temp_x - x_center)*cos(rotate_angle) - (temp_z - z_center)*sin(rotate_angle) + x_center + x_shift;
    mc_startXYZT[i][1] = temp_y + y_shift;
    mc_startXYZT[i][2] = (temp_x - x_center)*sin(rotate_angle) + (temp_z - z_center)*cos(rotate_angle) + z_center + z_shift;
    //mc_endXYZT
    temp_x = mc_endXYZT[i][0];
    temp_y = mc_endXYZT[i][1];
    temp_z = mc_endXYZT[i][2];
    mc_endXYZT[i][0] = (temp_x - x_center)*cos(rotate_angle) - (temp_z - z_center)*sin(rotate_angle) + x_center + x_shift;
    mc_endXYZT[i][1] = temp_y + y_shift;
    mc_endXYZT[i][2] = (temp_x - x_center)*sin(rotate_angle) + (temp_z- z_center)*cos(rotate_angle) + z_center + z_shift;
    //mc_startMomentum
    temp_x = mc_startMomentum[i][0];
    temp_y = mc_startMomentum[i][1];
    temp_z = mc_startMomentum[i][2];
    mc_startMomentum[i][0] = temp_x*cos(rotate_angle) - temp_z*sin(rotate_angle);
    mc_startMomentum[i][1] = temp_y;
    mc_startMomentum[i][2] = temp_x*cos(rotate_angle) + temp_z*sin(rotate_angle);

    //mc_endMomentum
    temp_x = mc_endMomentum[i][0];
    temp_y = mc_endMomentum[i][1];
    temp_z = mc_endMomentum[i][2];
    mc_endMomentum[i][0] = temp_x*cos(rotate_angle) - temp_z*sin(rotate_angle);
    mc_endMomentum[i][1] = temp_y;
    mc_endMomentum[i][2] = temp_x*cos(rotate_angle) + temp_z*sin(rotate_angle);

    //mc_trackPosition
    if (mcTree->GetBranch("mc_trackPosition")) {
      TClonesArray *trackPoints = (TClonesArray*)(*mc_trackPosition)[i];
      int nPoints = trackPoints->GetEntries();
      for (int j=0;j!=nPoints;j++){
	TLorentzVector *p = (TLorentzVector*)(*trackPoints)[j];
	//cout << p->X() << " " << p->Y() << " " << p->Z() << std::endl;
	Double_t temp_x, temp_y, temp_z;
	temp_x = (p->X() - x_center)*cos(rotate_angle) - (p->Z() - z_center)*sin(rotate_angle) + x_center + x_shift;
	temp_y = p->Y() + y_shift;
	temp_z = (p->X() - x_center)*sin(rotate_angle) + (p->Z() - z_center)*cos(rotate_angle) + z_center + z_shift;
	p->SetXYZT(temp_x,temp_y,temp_z,p->T());
      }
    }
    
  }
  
}
