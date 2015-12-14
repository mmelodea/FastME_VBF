#include "fme_particle_class.h"

void test(void){
  TFile *f = TFile::Open("StructToTProcPool/NewVBF_QED4_QCD0_SIG_FastME_2.root");
  TTree *t = (TTree*)f->Get("VBF");
  
  particle Particle;
  t->SetBranchAddress("Particle",&Particle.id);
  
  t->GetEntry(0);
  for(int i=0; i<6; i++)
  cout<<"Particles: "<<Particle.id[i]<<endl;
}
