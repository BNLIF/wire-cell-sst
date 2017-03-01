#include "WireCellSst/MCTruth.h"
#include "TLorentzVector.h"

WireCellSst::MCTruth::MCTruth(std::string rootfile)
  : mcTree(0)
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
    mcTree->SetBranchAddress("mc_trackPosition",&mc_trackPosition);
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


  mc_oldVertex[0] = -1;
  mc_oldVertex[0] = -1;
  mc_oldVertex[0] = -1;
}

float WireCellSst::MCTruth::find_neutrino_true_energy(int event_no){
  mcTree->GetEntry(event_no);
  return mc_nu_mom[3];
}

float WireCellSst::MCTruth::find_neutrino_visible_energy(int event_no){
  mcTree->GetEntry(event_no);
  // Hack for now ... 
  float visE =0;
  
   for (int i=0;i!=mc_Ntrack;i++){
     if (mc_mother[i] == 0){
       if (mc_pdg[i]==2212){
	 visE += sqrt(pow(mc_startMomentum[i][3],2)-pow(0.939,2));
       }
     }
   }

   return find_neutrino_true_energy(event_no) - visE;
}


WireCell::Point WireCellSst::MCTruth::find_neutrino_vertex(int event_no){
  mcTree->GetEntry(event_no);
  
  WireCell::Point vertex(0,0,0);
  
  vertex[0] = mc_nu_pos[0];
  vertex[1] = mc_nu_pos[1];
  vertex[2] = mc_nu_pos[2];


  return vertex;
}

WireCell::MCParticleSelection WireCellSst::MCTruth::find_primary_photons(int event_no){
  // mcTree->GetEntry(event_no);
  WireCell::MCParticleSelection photons;
  WireCell::Point primary_vertex = find_primary_vertex(event_no);
  
  for (int i=0;i!=mc_Ntrack;i++){
    if (mc_pdg[i] != 22) continue;

    float dis;
    dis = sqrt(pow(mc_startXYZT[i][0] - primary_vertex[0],2) 
	       + pow(mc_startXYZT[i][1] - primary_vertex[1],2) 
	       + pow(mc_startXYZT[i][2] - primary_vertex[2],2));
    //std::cout << dis << std::endl;
    if (dis < 0.5 ){ // cm
      WireCell::MCParticle *photon = new WireCell::MCParticle();
      photon->pdg = mc_pdg[i];
      
      // temporary fix 
      mc_startMomentum[i][2] = sqrt(pow(mc_startMomentum[i][3],2) - pow(mc_startMomentum[i][0],2) - pow(mc_startMomentum[i][1],2));
      for (int j=0;j!=4;j++){
	photon->startXYZT[j] = mc_startXYZT[i][j];
	photon->endXYZT[j] = mc_endXYZT[i][j];
	photon->startMomentum[j] = mc_startMomentum[i][j];
      }
      
      photons.push_back(photon);
    }
  }

  return photons;
}

WireCell::MCParticle* WireCellSst::MCTruth::find_primary_electron(int event_no){
  mcTree->GetEntry(event_no);
  
  for (int i=0;i!=mc_Ntrack;i++){
    if (mc_mother[i] == 0 && fabs(mc_pdg[i])==11){
      WireCell::MCParticle *electron = new WireCell::MCParticle();
      electron->pdg = mc_pdg[i];
      
      // a fix to start_Momentum[3] due to a previous bug ...
      mc_startMomentum[i][2] = sqrt(pow(mc_startMomentum[i][3],2) - pow(0.511e-3,2) - pow(mc_startMomentum[i][0],2) - pow(mc_startMomentum[i][1],2));
      //

      for (int j=0;j!=4;j++){
	electron->startXYZT[j] = mc_startXYZT[i][j];
	electron->endXYZT[j] = mc_endXYZT[i][j];
	electron->startMomentum[j] = mc_startMomentum[i][j];
      }

      
      

      if (mcTree->GetBranch("mc_trackPosition")) {
      	TClonesArray *trackPoints = (TClonesArray*)(*mc_trackPosition)[i];
      	int nPoints = trackPoints->GetEntries();
      	for (int j=0;j!=nPoints;j++){
      	  TLorentzVector *p = (TLorentzVector*)(*trackPoints)[j];
      	  WireCell::Point point(p->X(),p->Y(),p->Z());
      	  electron->trajectory.push_back(point);
      	  //   cout << p->X() << " " << p->Y() << " " << p->Z() << std::endl;
      	}
      }


      return electron;
    }
  }
  
  return 0;
  
}

WireCell::Point WireCellSst::MCTruth::find_primary_vertex(int event_no){
  
  mcTree->GetEntry(event_no);
  
  WireCell::Point vertex(0,0,0);
  
  for (int i = 0; i!= mc_Ntrack;i++){
    if (mc_mother[i] == 0) { // mother is primary particle ...
      vertex[0] = mc_startXYZT[i][0];
      vertex[1] = mc_startXYZT[i][1];
      vertex[2] = mc_startXYZT[i][2];
      break;
    }
  }


  return vertex;
}


WireCellSst::MCTruth::MCTruth(TTree *TMC)
{
  mcTree = TMC;
  if (!mcTree ) 
    std::cout<<"cannot find /Event/Sim tree."<<std::endl;
  mc_daughters= new std::vector<std::vector<int> >;
  // mcTree->SetBranchAddress("eventNo",&eventNo);
  // mcTree->SetBranchAddress("runNo",&runNo);
  // mcTree->SetBranchAddress("subRunNo",&subrunNo);
  eventNo = -1;
  runNo = -1;
  subrunNo = -1;

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
    mcTree->SetBranchAddress("mc_trackPosition",&mc_trackPosition);
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
  
  mcTree->SetBranchAddress("mc_oldVertex", &mc_oldVertex);
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

  Double_t temp_x=0, temp_y=0, temp_z=0;
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
    mc_startMomentum[i][2] = temp_x*sin(rotate_angle) + temp_z*cos(rotate_angle);

    //mc_endMomentum
    temp_x = mc_endMomentum[i][0];
    temp_y = mc_endMomentum[i][1];
    temp_z = mc_endMomentum[i][2];
    mc_endMomentum[i][0] = temp_x*cos(rotate_angle) - temp_z*sin(rotate_angle);
    mc_endMomentum[i][1] = temp_y;
    mc_endMomentum[i][2] = temp_x*sin(rotate_angle) + temp_z*cos(rotate_angle);

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
