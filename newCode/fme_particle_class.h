///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::///
///:: Library to class of particles properties used by FastME ::///
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::///


#ifndef fme_particle_class_h
#define fme_particle_class_h

#define nObjs 6

typedef struct{
  Int_t    id[nObjs];
  Double_t pt[nObjs];
  Double_t pt_res[nObjs];
  Double_t eta[nObjs];
  Double_t eta_res[nObjs];
  Double_t phi[nObjs];
  Double_t phi_res[nObjs];
} particle;

#endif
