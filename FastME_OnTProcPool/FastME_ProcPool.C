///To print informations about DR functions running
#define debug		false

#include <iostream>
#include <string>
#include <exception>
#include <vector>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStopwatch.h>

///Headers to TProcPool
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TProcPool.h>
#include <PoolUtils.h>
#include <TSystem.h>


#define pi		3.14159265358979312
#define pedestal 	-99					///Reset Value to Variables


///Scale Factors to normalize Deltas
#define scale_dPt	50.
#define scale_dEta	5.
#define scale_dPhi	pi

using namespace std;


///========= Compute Discriminant Value ==========================================================
///Based on Distance
Double_t PsbD(Double_t min_dr_sig, Double_t min_dr_bkg){
  Double_t DD = min_dr_bkg/(min_dr_sig + min_dr_bkg);
  return DD;
}
///================================================================================================


///Parameters (from left to right):
///1. Address to Data sample;
///2. Vector with address of MC samples;
///3. Number of Data events;
///4. Number of MC types (if all have different classification that is just vector size);
///5. Number of final state particles;
///6. Name of the output file to store FastME analysis results.
void FastME_ProcPool(string Data_Path, vector<string> MCs, const Int_t N_DT, const Int_t N_TMC, const Int_t N_FSParticles, TString out_name){
  
  ///TProcPool declaration to objects to be analised
  auto workItem = [Data_Path,N_DT,N_TMC,N_FSParticles](TTreeReader &tread) -> TObject* {
    TStopwatch t2;
    
    ///Defines 2D histogram to stores minimum distances
    TH2D *mdists = new TH2D("mdists","Minimum Data-MC distances found",N_DT,0,N_DT,N_TMC,0,N_TMC);
    mdists->GetXaxis()->SetTitle("Data Events");
    mdists->GetYaxis()->SetTitle("MC Types (0-1:Sig, 1-n:Bkg)");
    mdists->SetDirectory(0);
    
    ///Addresses the MC branches to be used
    TTreeReaderValue<Int_t>    McType(tread, "McType"); ///McType for Signal=0 and Background >0
    TTreeReaderArray<Int_t>    McId(tread, "ParticleId");
    TTreeReaderArray<Double_t> McPt(tread, "ParticlePt");
    TTreeReaderArray<Double_t> McEta(tread, "ParticleEta");
    
    ///Addresses the Data branches to be used
    TFile *fData = TFile::Open((TString)Data_Path);
    TTreeReader refReader("VBF",fData);
    TTreeReaderArray<Int_t>    DataId(refReader, "ParticleId");
    TTreeReaderArray<Double_t> DataPt(refReader, "ParticlePt");
    TTreeReaderArray<Double_t> DataEta(refReader, "ParticleEta");


    ///Loop on Data events
    for(Int_t dt=0; dt<N_DT; dt++){
    
      if(dt!= 0 && dt%(N_DT/10) == 0){ 
	cout<<":: [Remaining Data]: "<<N_DT-dt<<"\t\t[Elapsed Time]: ";
	t2.Stop();
	t2.Print();
	t2.Continue();
      }
      refReader.SetEntry(dt); ///Move on Data loop
      Double_t min_distance = 1.E15, event_distance=-1;
      Int_t f_type=-1;
      Int_t nMonteCarlo = tread.GetEntries(true);
      
      for(Int_t mc=0; mc<nMonteCarlo; mc++){
	tread.SetEntry(mc); ///Move on MC loop
        
	
      ///============================================================================================================
      ///::::::::::::::::::::::::: Fast Matrix Element method to compute Data - MC distance :::::::::::::::::::::::::
      ///============================================================================================================

	///Checks the final states
	int dt_el= 0, dt_mu= 0, mc_el= 0, mc_mu= 0;
	for(int u=0; u<N_FSParticles; u++){
	  if(abs(DataId[u]) == 11) dt_el++;
	  if(abs(DataId[u]) == 13) dt_mu++;
	  if(abs(McId[u])   == 11) mc_el++;
	  if(abs(McId[u])   == 13) mc_mu++;
	}
	///Avoid wrong framework (no works with different final states)
	if( dt_el != mc_el || dt_mu != mc_mu ) continue;

	bool mc_approved;
	Double_t particles_distance, min_particles_distance;
	Double_t dPt=0, dEta=0, dPhi=0;
	Double_t sum_dPt2=0, sum_dEta2=0, sum_dPhi2=0, event_distance=-1;
	int min_imc;
	vector<int> vmin_imc;
	vmin_imc.clear();
	
	for(int idt=0; idt<N_FSParticles; idt++){
	  min_particles_distance = 1.E15;
	  particles_distance = -1; 
	  min_imc = -1; 
    
	  for(int imc=0; imc<N_FSParticles; imc++){
	    mc_approved = true;
	    if(int(vmin_imc.size()) > 0)
	    for(int g=0; g<int(vmin_imc.size()); g++)
	    if(imc == vmin_imc[g]) mc_approved = false;
	
	    if(mc_approved == true){
	      ///Avoid different Data-MC particles comparison
	      if(DataId[idt] != McId[imc]) continue;
	      if( debug ) cout<<"DataPos: "<<idt<<"  ID: "<<DataId[idt]<<"  MCPos: "<<imc<<"   ID: "<<McId[imc]<<" >>> ";

	      dPt  = (DataPt[idt]-McPt[imc])/(scale_dPt);//DataPtRes[idt]);
	      dEta = (DataEta[idt]-McEta[imc])/(scale_dEta);//DataEtaRes[idt]);
	
	      particles_distance = sqrt(dPt*dPt + dEta*dEta);
	      if( debug ) cout<<"particles_distance = "<<particles_distance<<endl;
	      if(particles_distance < min_particles_distance){
		min_imc = imc;
		min_particles_distance = particles_distance;
	      }
	    }
	  }///Ends MC event loop
      
	  ///Monitor of chosen MCs to avoid object recounting
	  vmin_imc.push_back(min_imc);
	  if( debug ) cout<<"Chosen->>  MCPos: "<<min_imc<<"   ID: "<<McId[min_imc]<<endl;
      
	  dPt  = (DataPt[idt]-McPt[min_imc])/(scale_dPt);//DataPtRes[idt]);
	  dEta = (DataEta[idt]-McEta[min_imc])/(scale_dEta);//DataEtaRes[idt]);
	  sum_dPt2  += dPt*dPt;
	  sum_dEta2 += dEta*dEta;
	}///Ends Data event loop  
  
	///Compute final Data-MC events distance
	event_distance = sqrt(sum_dPt2 + sum_dEta2);
      ///============================================================================================================

    
	///Finds minimum distance computed and closet MC event type (signal/background)
	if(event_distance < min_distance){
	  min_distance = event_distance;
	  f_type = *McType;
	}
      }///End MC sample loop

      
      ///Stores the minimum distances found
      mdists->Fill(dt,f_type,min_distance);
      if( debug ) cout<<"dt: "<<dt<<"\tf_type: "<<f_type<<"\tmin_distance: "<<min_distance<<endl;
    }///End Data sample loop
    
    t2.Stop();
    delete fData;
    return mdists;
  };
  
  ///_____________________________ For timming the process ________________________________________________________
  TStopwatch t;
  t.Start();
  cout<<"\n\n";
  cout<<"=================================================================================================="<<endl;
  cout<<"::::::::::::::::::::::::[ Starting Fast Matrix Element Analysis ]:::::::::::::::::::::::::::::::::"<<endl;
  cout<<"=================================================================================================="<<endl;
  ///--------------------------------------------------------------------------------------------------------------

  ///Calls analysis through TProcPool
  ///TObject returned of TTree type (but not all member classes are accessibles)
  UInt_t ncores = MCs.size();
  TProcPool workers(ncores);
  auto f_hist = workers.ProcTree(MCs, workItem);
  
  ///This allows one to get all TH2 member classes
  TFile *tmp = TFile::Open(out_name+".root","recreate");
  f_hist->Write();
  if( debug ) f_hist->Draw();
  TH2D *mdists2 = (TH2D*)tmp->Get("mdists");

  ///_______________________ Compute discriminant from MDMCED _____________________________________________________
  cout<<":: [Time to Make Analysis]: "; t.Stop(); t.Print(); t.Continue();
  cout<<"\n::::::::::::::::::::::::::::::::[ Computing discriminant ]::::::::::::::::::::::::::::::::::::::::"<<endl;
  ///--------------------------------------------------------------------------------------------------------------
  TH1D *hdisc = new TH1D("hdisc","",50,0,1);
  for(Int_t data=0; data<N_DT; data++){
    if(data%(N_DT/10) == 0){
      cout<<":: [Remaining Data]: "<<N_DT-data<<"\t\t[Elapsed Time]: ";
      t.Stop(); t.Print(); t.Continue();
    }

    Double_t min_dr_sig = mdists2->GetBinContent(data+1,1);
    Double_t min_dr_bkg = 1.E15;
    ///Needs to find minimum distance to MC
    for(Int_t mcs=1; mcs<N_TMC; mcs++)
      if( mdists2->GetBinContent(data+1,mcs+1) < min_dr_bkg )
	min_dr_bkg = mdists2->GetBinContent(data+1,mcs+1);

    hdisc->Fill( PsbD(min_dr_sig, min_dr_bkg) );
    if( debug ) cout<<"SigMin: "<< min_dr_sig <<"\t\tBkgMin: "<< min_dr_bkg << endl;
  }
  hdisc->Write();
  tmp->Close();
  
  ///________________________________ Stoping timming ________________________________________________________
  cout<<"\n::::::::::::::::::::::::::::::::::::[ Process Finished ]::::::::::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<":: [Analysis Total Time]: "; t.Stop(); t.Print();
  cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"<<endl;
  ///---------------------------------------------------------------------------------------------------------
}
///====================================================================================================
