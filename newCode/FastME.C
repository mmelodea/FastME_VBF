///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::	    FAST MATRIX ELEMENT TO VBF TOPOLOGY		::::::::::::::::::::::::::::::::::
///::::::::::::::::::::::::::::::::::	    Code Author: Miquéias M. de Almeida		::::::::::::::::::::::::::::::::::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
///::															::
///:: To run the analysis make in Root (compilation mode is the fastest):						::
///:: root -l -b -q FastME.C+												::
///:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

///NOTE The first version of my FastME didn't take care about the lepton flavors when comparing them. Even so, the results
///NOTE found were similar to the new ones with separation of flavors. This is a interisting behavior to discuss later.

///To print informations about DR functions running
#define debug		false

///The header to refer to particles properties to be used in analysis
#include "fme_particle_class.h"

#include <iostream>
#include <string>
#include <ctime>
#include <exception>
#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TCanvas.h>

#define pi		3.14159265358979312
#define pedestal 	-99					///Reset Value to Variables


///Scale Factors to normalize Deltas
#define scale_dPt	50.
#define scale_dEta	5.
#define scale_dPhi	pi

///Flag to control dPhi use on event distance
#define use_dPhi	false

///Flag to control number of objects in the final state
#define nObjs		6

using namespace std;



///========== Compute Distance Between Events - Minimum Distance Method ===============
///Takes the combination that gives minimum distance when equal final states
Double_t ComputeDR_MinDist(particle *DtParticle, particle *McParticle){
  
  ///Checks the final states
  int dt_el= 0, dt_mu= 0, mc_el= 0, mc_mu=0;
  for(int u=0; u<nObjs; u++){
    if(abs(DtParticle->id[u]) == 11) dt_el++;
    if(abs(DtParticle->id[u]) == 13) dt_mu++;
    if(abs(McParticle->id[u]) == 11) mc_el++;
    if(abs(McParticle->id[u]) == 13) mc_mu++;
  }
  
  ///To distiguish between right and wrong framework
  if( dt_el != mc_el || dt_mu != mc_mu ) return -1;
  
  ///If Data and MC final state are equals compute the distance between events
  if( debug ){
    cout<<"\n\n(#e,#u)>>>> Data ("<<dt_el<<", "<<dt_mu<<")\t\tMc ("<<mc_el<<", "<<mc_mu<<")"<<endl;
    for(int f=0; f<6; f++){
      cout<<"Data("<<f<<") >> Pt: "<<DtParticle->pt[f]<<"\t\tEta: "<<DtParticle->eta[f]<<endl;
      cout<<"Mc("<<f<<") >> Pt: "<<McParticle->pt[f]<<"\t\tEta: "<<McParticle->eta[f]<<endl;
    }
  }
  bool mc_approved;
  Double_t particles_distance, min_particles_distance;
  Double_t dPt=0, dEta=0, dPhi=0;
  Double_t sum_dPt2=0, sum_dEta2=0, sum_dPhi2=0, event_distance=-1;
  
  int min_imc;
  std::vector<int> vmin_imc;
  vmin_imc.clear();
  for(int idt=0; idt<nObjs; idt++){
    
    min_imc = -1, particles_distance = -1; min_particles_distance = 1.E15;
    for(int imc=0; imc<nObjs; imc++){
      mc_approved = true;
      if(int(vmin_imc.size()) > 0)
	for(int g=0; g<int(vmin_imc.size()); g++)
	  if(imc == vmin_imc[g]) mc_approved = false;
	
      if(mc_approved == true){
	///Avoid different Data-MC particles comparison
	if(DtParticle->id[idt] != McParticle->id[imc]) continue;
	if( debug ) cout<<"DataPos: "<<idt<<"  ID: "<<DtParticle->id[idt]<<"  MCPos: "<<imc<<"   ID: "<<McParticle->id[imc]<<" >>> ";

	dPt  = (DtParticle->pt[idt]-McParticle->pt[imc])/(scale_dPt*DtParticle->pt_res[idt]);
	dEta = (DtParticle->eta[idt]-McParticle->eta[imc])/(scale_dEta*DtParticle->eta_res[idt]);
	if(use_dPhi == true){
	  dPhi = (DtParticle->phi[idt]-McParticle->phi[imc]);
	  if(fabs(dPhi) > pi)
	    dPhi = (2*pi-fabs(dPhi))/(scale_dPhi*DtParticle->phi[idt]);
	  else
	    dPhi = dPhi/(scale_dPhi*DtParticle->phi[idt]);
	}
	
	particles_distance = sqrt(dPt*dPt + dEta*dEta);
	if( debug ) cout<<"particles_distance = "<<particles_distance<<endl;
	if(use_dPhi == true) particles_distance = sqrt(dPt*dPt + dEta*dEta + dPhi*dPhi);
	if(particles_distance < min_particles_distance){
	  min_imc = imc;
	  min_particles_distance = particles_distance;
	}
      }
    }
      
    ///Monitor of chosen MCs to avoid object recounting
    vmin_imc.push_back(min_imc);
    if( debug ) cout<<"Chosen->>  MCPos: "<<min_imc<<"   ID: "<<McParticle->id[min_imc]<<endl;
      
    dPt  = (DtParticle->pt[idt]-McParticle->pt[min_imc])/(scale_dPt*DtParticle->pt_res[idt]);
    dEta = (DtParticle->eta[idt]-McParticle->eta[min_imc])/(scale_dEta*DtParticle->eta_res[idt]);
    if(use_dPhi == true){
      dPhi = (DtParticle->phi[idt]-McParticle->phi[min_imc]);
      if(fabs(dPhi) > pi)
	dPhi = (2*pi-fabs(dPhi))/(scale_dPhi*DtParticle->phi_res[idt]);
      else
	dPhi = dPhi/(scale_dPhi*DtParticle->phi_res[idt]);
    }
    sum_dPt2  += dPt*dPt;
    sum_dEta2 += dEta*dEta;
    if(use_dPhi == true) sum_dPhi2 += dPhi*dPhi;
  }
    
  if(use_dPhi == false) event_distance = sqrt(sum_dPt2 + sum_dEta2);
  else event_distance = sqrt(sum_dPt2 + sum_dEta2 + sum_dPhi2);
  
  return event_distance;
}
///================================================================================================
///================================================================================================


/*
///============ Compute Distance Between Data-MC Events - Media Method ============================
///Takes the media of the possible combinations when equal final states
Double_t ComputeDR_Media(Int_t FS, Int_t DataParticleID[nObjs], Float_t Data[nObjs][nProp][nVals],
			 Int_t McParticleID[nObjs], Float_t MC[nObjs][nProp][nVals]){
  
  Double_t dPt=0, dEta=0, dPhi=0;
  Double_t sum_dPt2=0, sum_dEta2=0, sum_dPhi2=0, event_distance=-1;
  
  if( debug ) cout<<"\nFS: "<<FS<<endl;
  for(int idt=0; idt<nObjs; idt++)
  for(int imc=0; imc<nObjs; imc++){
    ///Avoid different Data-MC particles comparison
    if(DataParticleID[idt] != McParticleID[imc]) continue;
    if( debug ) cout<<"DataPos: "<<idt<<"  ID: "<<DataParticleID[idt]<<"  MCPos: "<<imc<<"   ID: "<<McParticleID[imc]<<endl;
  
    if(FS == 2){
      
      dPt  = (Data[idt][0][0]-MC[imc][0][0])/(scale_dPt*Data[idt][0][1]);
      dEta = (Data[idt][1][0]-MC[imc][1][0])/(scale_dEta*Data[idt][1][1]);
      if(use_dPhi == true){
	dPhi = (Data[idt][2][0]-MC[imc][2][0]);
	if(fabs(dPhi) > pi)
	  dPhi = (2*pi-fabs(dPhi))/(scale_dPhi*Data[idt][2][1]);
	else
	  dPhi = dPhi/(scale_dPhi*Data[idt][2][1]);
      }
      
    }
    
    ///Takes the media of combinations in 4e and 4mu final states
    else if(FS == 0 || FS == 1){
      
      dPt  = (Data[idt][0][0]-MC[imc][0][0])/(2*scale_dPt*Data[idt][0][1]);
      dEta = (Data[idt][1][0]-MC[imc][1][0])/(2*scale_dEta*Data[idt][1][1]);
      if(use_dPhi == true){
	dPhi = (Data[idt][2][0]-MC[imc][2][0]);
	if(fabs(dPhi) > pi)
	  dPhi = (2*pi-fabs(dPhi))/(2*scale_dPhi*Data[idt][2][1]);
	else
	  dPhi = dPhi/(2*scale_dPhi*Data[idt][2][1]);
      }
      
    }

    else{
      cout<<"[Error] Irregular Final State passed to ComputeDR function!"<<endl;
      throw exception();
    }
    
    sum_dPt2  += dPt*dPt;
    sum_dEta2 += dEta*dEta;
    if(use_dPhi == true)
      sum_dPhi2 += dPhi*dPhi;
  }      
  
  if(use_dPhi == false) event_distance = sqrt(sum_dPt2 + sum_dEta2);
  else event_distance = sqrt(sum_dPt2 + sum_dEta2 + sum_dPhi2);

  if(event_distance == -1) throw exception();
  else return event_distance;
}
///================================================================================================
*/


///========= Compute Discriminant Values ==========================================================
///Based on Distance
Double_t PsbD(Double_t min_dr_sig, Double_t min_dr_bkg){
  Double_t DD = min_dr_bkg/(min_dr_sig + min_dr_bkg);
  return DD;
}

///Based on Weight
Double_t PsbW(Double_t sig_event_weight, Double_t sigXS, Double_t bkg_event_weight, Double_t bkgXS){
  ///Splits the discriminant calculus according to final state (the cross section are different)
  Double_t DW = (sig_event_weight/sigXS)/(sig_event_weight/sigXS + bkg_event_weight/bkgXS);
  return DW;
}
///================================================================================================

int FastME(){
   
  ///--------------	Preparing Inputs	-----------------------------------------
  TString SamplesPath = "../VBF_Samples/Struct/";
  TString Out_Name    = "BkgLikeData_EWbkg";
  TString Data_Path   = SamplesPath + "VBF_QED4_QCD0_SIG_FastME_1.root";
  //TString Data_Path   = SamplesPath + "VBF_QED4_QCD0_BKG_FastME_1.root";
  
  ///MCs
  TString MC_Sig_Path = SamplesPath + "VBF_QED4_QCD0_SIG_FastME_2.root";
  TString MC_Bkg_Path = SamplesPath + "VBF_QED4_QCD0_BKG_FastME_2.root";
  
  TString Tree_Name   = "VBF";
  TString Objs_Branch = "Particle";
  TString EW_Branch   = "EventWeight";
  TString Num_Trials  = "Ntrials";
  ///------------------------------------------------------------------------------------
  
  ///Flag to control what discriminant to save
  ///(0: based on distance, 1: based on weight, 2: both)
  int Compute = 0;
  
  ///Flag to control the Data-MC comparisons to compute DR
  ///(0: minimum distance, 1: media of combinations)
  int DrMethod = 0;
  
  ///Flags to control final state usage
  ///(to use only one channel - combinations not implemented!)
  //bool fs4e 	= false;
  //bool fs4u 	= false;
  //bool fs2e2u 	= false;  
  
  
  ///Opening the input files
  TFile *fData = TFile::Open(Data_Path);
  TFile *fMC_Sig = TFile::Open(MC_Sig_Path);
  TFile *fMC_Bkg = TFile::Open(MC_Bkg_Path);
  TTree *Data_Tree = (TTree*)fData->Get(Tree_Name);
  TTree *MC_Sig_Tree = (TTree*)fMC_Sig->Get(Tree_Name);
  TTree *MC_Bkg_Tree = (TTree*)fMC_Bkg->Get(Tree_Name);
  
  ///Creating the Data Tree
  particle *DtParticle = nullptr;
  Data_Tree->SetBranchAddress(Objs_Branch,&DtParticle);
  int ndata = Data_Tree->GetEntries();
  
  ///Creating the Signal Tree
  Int_t MC_SIG_NTRIALS;
  Double_t MC_SIG_WEIGHT;
  particle *McSigParticle = nullptr;
  MC_Sig_Tree->SetBranchAddress(Objs_Branch,&McSigParticle);
  if(Compute == 1 || Compute == 2){
    MC_Sig_Tree->SetBranchAddress(EW_Branch,&MC_SIG_WEIGHT);
    MC_Sig_Tree->SetBranchAddress(Num_Trials,&MC_SIG_NTRIALS);
  }
  int nsig = MC_Sig_Tree->GetEntries();
  
  ///Creating the Background Tree
  Int_t MC_BKG_NTRIALS;
  Double_t MC_BKG_WEIGHT;
  particle *McBkgParticle = nullptr;
  MC_Bkg_Tree->SetBranchAddress(Objs_Branch,&McBkgParticle);
  if(Compute == 1 || Compute == 2){
    MC_Bkg_Tree->SetBranchAddress(EW_Branch,&MC_BKG_WEIGHT);
    MC_Bkg_Tree->SetBranchAddress(Num_Trials,&MC_BKG_NTRIALS);
  }
  int nbkg = MC_Bkg_Tree->GetEntries();
  
  
  cout<<"\n::::::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<":::::: Starting FastME Processing ::::::"<<endl;
  cout<<"::::::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<":: DrComputation Mode: ";
       if(DrMethod == 0) cout<<"Minimum Distance"<<endl;
  else if(DrMethod == 1) cout<<"Media"<<endl;
  if(use_dPhi == false)
    cout<<":: No Using dPhi!"<<endl;
  cout<<"::--------------------------------------"<<endl;
  cout<<":: #Data Events:   "<<ndata<<endl;
  cout<<":: #MC Sig Events: "<<nsig<<endl;
  cout<<":: #MC BKG Events: "<<nbkg<<endl;
  cout<<"::--------------------------------------"<<endl;
/*  
  ///Compute the final states yields in the ntuples
  ///And the Cross Section to compute weights
  int d4e=0, d4u=0, d2e2u=0;
  for(int fs=0; fs<ndata; fs++){
    Data_Tree->GetEntry(fs);
    if(Dat==0) d4e    += 1;
    if(DataFS==1) d4u    += 1;
    if(DataFS==2) d2e2u  += 1;
  }
  int s4e=0, s4u=0, s2e2u=0;
  Double_t sig_4eSum_weight=0, sig_4uSum_weight=0, sig_2e2uSum_weight=0;
  Double_t sig_4eSum_trials=0, sig_4uSum_trials=0, sig_2e2uSum_trials=0;
  for(int fs=0; fs<nsig; fs++){
    MC_Sig_Tree->GetEntry(fs);
    if(MC_SIG_FS==0){
      s4e    += 1;
      if(Compute == 1 || Compute == 2){
	sig_4eSum_weight += MC_SIG_WEIGHT;
	sig_4eSum_trials += MC_SIG_NTRIALS;
      }
    }
    if(MC_SIG_FS==1){
      s4u    += 1;
      if(Compute == 1 || Compute == 2){
	sig_4uSum_weight += MC_SIG_WEIGHT;
	sig_4uSum_trials += MC_SIG_NTRIALS;
      }
    }
    if(MC_SIG_FS==2){
      s2e2u  += 1;
      if(Compute == 1 || Compute == 2){
	sig_2e2uSum_weight += MC_SIG_WEIGHT;
	sig_2e2uSum_trials += MC_SIG_NTRIALS;
      }
    }      
  }
  int b4e=0, b4u=0, b2e2u=0;
  Double_t bkg_4eSum_weight=0, bkg_4uSum_weight=0, bkg_2e2uSum_weight=0;
  Double_t bkg_4eSum_trials=0, bkg_4uSum_trials=0, bkg_2e2uSum_trials=0;
  for(int fs=0; fs<nbkg; fs++){
    MC_Bkg_Tree->GetEntry(fs);
    if(MC_BKG_FS==0){
      b4e    += 1;
      if(Compute == 1 || Compute == 2){
	bkg_4eSum_weight += MC_BKG_WEIGHT;
	bkg_4eSum_trials += MC_BKG_NTRIALS;
      }
    }
    if(MC_BKG_FS==1){
      b4u    += 1;
      if(Compute == 1 || Compute == 2){
	bkg_4uSum_weight += MC_BKG_WEIGHT;
	bkg_4uSum_trials += MC_BKG_NTRIALS;
      }
    }
    if(MC_BKG_FS==2){
      b2e2u  += 1;
      if(Compute == 1 || Compute == 2){
	bkg_2e2uSum_weight += MC_BKG_WEIGHT;
	bkg_2e2uSum_trials += MC_BKG_NTRIALS;
      }
    }
  }
  Double_t SigTrials[3], BkgTrials[3], SigXS[3], BkgXS[3];
  if(Compute == 1 || Compute == 2){
    SigTrials[0] = sig_4eSum_trials;			BkgTrials[0] = bkg_4eSum_trials;
    SigTrials[1] = sig_4uSum_trials;			BkgTrials[1] = bkg_4uSum_trials;
    SigTrials[2] = sig_2e2uSum_trials;			BkgTrials[2] = bkg_2e2uSum_trials; 

    ///Computing and storing XS
    SigXS[0] = sig_4eSum_weight/SigTrials[0];		BkgXS[0] = bkg_4eSum_weight/BkgTrials[0];
    SigXS[1] = sig_4uSum_weight/SigTrials[1];		BkgXS[1] = bkg_4uSum_weight/BkgTrials[1];
    SigXS[2] = sig_2e2uSum_weight/SigTrials[2];		BkgXS[2] = bkg_2e2uSum_weight/BkgTrials[2];
  }
    
  ///Plot on PC Screen the XS
  if(Compute == 1 || Compute == 2){
    cout<<":: XS Computed from the Ntuples"<<endl;
    cout<<":: Sig   "<<"4e: "<<SigXS[0]<<"\t4mu: "<<SigXS[1]<<"\t2e2mu: "<<SigXS[2]<<endl;
    cout<<":: Bkg   "<<"4e: "<<BkgXS[0]<<"\t4mu: "<<BkgXS[1]<<"\t2e2mu: "<<BkgXS[2]<<endl;
    cout<<"::--------------------------------------"<<endl;
  }
  
  cout<<":: Final State Yields in the Ntuples"<<endl;
  cout<<":: Data  "<<"4e: "<<d4e<<"\t4mu: "<<d4u<<"\t2e2mu: "<<d2e2u<<endl;
  cout<<":: Sig   "<<"4e: "<<s4e<<"\t4mu: "<<s4u<<"\t2e2mu: "<<s2e2u<<endl;
  cout<<":: Bkg   "<<"4e: "<<b4e<<"\t4mu: "<<b4u<<"\t2e2mu: "<<b2e2u<<endl;
  cout<<"::______________________________________"<<endl;
  cout<<"::--------------------------------------"<<endl;
*/  
  
  ///Creating Tree to store Fast Matrix Element results
  Double_t dr_test, minDR_toSig, minDR_toBkg, psb_distance, sig_event_weight;
  Double_t min_dr_sig, min_dr_bkg, bkg_event_weight, psb_weight;
  TTree *FME_out = new TTree("FastME_Results","Fast Matrix Element Results");
  FME_out->SetDirectory(0);
  FME_out->Branch("minDR_toSig",&minDR_toSig);
  FME_out->Branch("minDR_toBkg",&minDR_toBkg);
  FME_out->Branch("PSB_Distance",&psb_distance);
  if(Compute == 1 || Compute == 2){
    FME_out->Branch("WSig_ToEvent",&sig_event_weight);
    FME_out->Branch("WBkg_ToEvent",&bkg_event_weight);
    FME_out->Branch("PSB_Weight",&psb_weight);
  }

  TH1D *h = new TH1D("Distances", "Distances between Data and MC events", 200, 0, 20);
  ///For timming the process
  time_t start, stop, delaied;
  double seconds, elapsed_time;
  string unity = "s";
  time(&start);
  ///--------------------------
  
  ///Starts the analysis
  ndata = 100;
  for(int i=0; i<ndata; i++){
    time(&delaied);
    if(i % (ndata/10) == 0 && i != 0) cout<<":: Remaining Data: "<<ndata-i<<"\t||Elapsed: "<<difftime(delaied,start)<<"s"<<endl;
    Data_Tree->GetEntry(i);
        
    if(i < 2) cout<<"Dt("<<i<<") >> Pt: "<<DtParticle->pt[0]<<"\t\tEta: "<<DtParticle->eta[0]<<endl;
    ///Reseting Variables
    minDR_toSig   	= pedestal;
    minDR_toBkg   	= pedestal;
    psb_distance  	= pedestal;
    sig_event_weight 	= pedestal;
    bkg_event_weight 	= pedestal;
    psb_weight		= pedestal;
    min_dr_sig    	= 1.E15;
    min_dr_bkg    	= 1.E15;
    
    ///:::::::::::::  Tests MC Signal  ::::::::::::::::::::::::::::::::::::::::::::::::::
    nsig = 100;
    for(int s=0; s<nsig; s++){
      MC_Sig_Tree->GetEntry(s);
      if(i<2 && s < 2) cout<<"Mc1("<<s<<") >> Pt: "<<McSigParticle->pt[0]<<"\t\tEta: "<<McSigParticle->eta[0]<<endl;
      
	   if(DrMethod == 0) dr_test = ComputeDR_MinDist(DtParticle, McSigParticle);
      //else if(DrMethod == 1) dr_test = ComputeDR_Media(DtParticle, McSigParticle);
      if(dr_test >= 0) h->Fill(dr_test);  
      if(dr_test >= 0 && dr_test < min_dr_sig){
	min_dr_sig = dr_test;
	sig_event_weight = MC_SIG_WEIGHT;
      }
    }
    ///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    

    ///::::::::::  Using MC Background  :::::::::::::::::::::::::::::::::::::::::::::::::
    nbkg = 100;
    for(int b=0; b<nbkg; b++){
      MC_Bkg_Tree->GetEntry(b);
      if(i < 2 && b < 2) cout<<"Mc2("<<b<<") >> Pt: "<<McBkgParticle->pt[0]<<"\t\tEta: "<<McBkgParticle->eta[0]<<endl;

	   if(DrMethod == 0) dr_test = ComputeDR_MinDist(DtParticle, McBkgParticle);
      //else if(DrMethod == 1) dr_test = ComputeDR_Media(DtParticle, McBkgParticle);
      if(dr_test >= 0) h->Fill(dr_test);  
      if(dr_test >= 0 && dr_test < min_dr_bkg){
	min_dr_bkg = dr_test;
	bkg_event_weight = MC_BKG_WEIGHT;
      }
    }
    ///::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    
    ///Getting the discriminant value based on event distance
    if(min_dr_sig == 0 && min_dr_bkg == 0){ 
      min_dr_sig = 1.;
      min_dr_bkg = 1.;
    }
    psb_distance = PsbD(min_dr_sig, min_dr_bkg);
    if(psb_distance < 0 || psb_distance > 1){
      cout<<"[Warning] PSB_Distance: "<<psb_distance<<"    SigD: "<<min_dr_sig<<"  BkgD: "<<min_dr_bkg<<endl;
    }
    minDR_toSig  = min_dr_sig;
    minDR_toBkg  = min_dr_bkg;
    
/*    
    ///Getting the discriminant value based on most close MC event weight
    if(Compute == 1 || Compute == 2){
      if(sig_event_weight == 0 && bkg_event_weight == 0){
	sig_event_weight = 1.;
	bkg_event_weight = 1.;
      }
      psb_weight = PsbW(sig_event_weight/SigTrials[DataFS], SigXS[DataFS],
			bkg_event_weight/BkgTrials[DataFS], BkgXS[DataFS]);
  }
    
    FinalState = DataFS;
*/
    ///Stores in the Tree the results
    //FME_out->Fill();
  }
  
  h->Draw();
  //TFile *FME_Results = new TFile(Out_Name+".root","recreate");
  //FME_out->Write();
  //FME_Results->Close();
  
  
  ///Stoping timming
  time(&stop);
  cout<<":::::::::: Process Finished ::::::::::::"<<endl;    
  seconds = difftime(stop,start);
  if(seconds < 60) elapsed_time = seconds;
  if(seconds >= 60){
    elapsed_time = seconds/60.;
    unity = "min.";
  }
  if(seconds >= 3600){
    elapsed_time = seconds/3600.;
    unity = "h";
  }
  cout<<":: Time Consumed: "<<elapsed_time<<unity<<endl;
  cout<<"::--------------------------------------\n"<<endl;
  
  
  return 0; //if process well finished
}
///====================================================================================================
