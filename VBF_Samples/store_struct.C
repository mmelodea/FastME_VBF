#include "fme_particle_class.h"

void store_struct(TString file_name){
  
  TFile *f = TFile::Open(file_name);
  TTree *t = (TTree*)f->Get("VBF");
  Int_t ParticleID[6];
  Double_t RecoParticle[6][3][2], Zon_mass, Zoff_mass, ZZ_mass;
  t->SetBranchAddress("RecoParticle",&RecoParticle);
  t->SetBranchAddress("ParticleID",&ParticleID);
  t->SetBranchAddress("Zon_mass",&Zon_mass);
  t->SetBranchAddress("Zoff_mass",&Zoff_mass);
  t->SetBranchAddress("H_reco_mass",&ZZ_mass);
  
  TTree *nt = new TTree("VBF","New VBF samples with matrix inserted on struct");
  nt->SetDirectory(0);
  particle Particle;
  nt->Branch("Particle",&Particle);
  nt->Branch("Zon_mass",&Zon_mass,"Zon_mass/D");
  nt->Branch("Zoff_mass",&Zoff_mass,"Zoff_mass/D");
  nt->Branch("ZZ_mass",&ZZ_mass,"ZZ_mass/D");
  
  particle tmp_particle;
  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    
    for(int j=0; j<6; j++){
      
      Particle.id[j]      = ParticleID[j];
      
      Particle.pt[j]      = RecoParticle[j][0][0];
      Particle.pt_res[j]  = RecoParticle[j][0][1];
      
      Particle.eta[j]     = RecoParticle[j][1][0];
      Particle.eta_res[j] = RecoParticle[j][1][1];
      
      Particle.phi[j]     = RecoParticle[j][2][0];
      Particle.phi_res[j] = RecoParticle[j][2][1];
      
    }
    nt->Fill();
  }
  
  TFile *fn = new TFile("New"+file_name,"recreate");
  nt->Write();
  fn->Close();
  
}
