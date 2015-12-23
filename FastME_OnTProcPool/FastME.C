///####################################################################################################################
///===================================== [Code Author: Miqueias M. de Almeida] ========================================
///####################################################################################################################


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


///========== Compute Discriminant Value ===============
///Based on Distance
Double_t PsbD(Double_t min_dr_sig, Double_t min_dr_bkg){
  Double_t DD = min_dr_bkg/(min_dr_sig + min_dr_bkg);
  return DD;
}
///=====================================================






///********************************************************************************************************************
///Parameters (from left to right):
///1. Address to Data sample;
///2. Vector with address of MC samples;
///3. Name of tree containig the events;
///4. Vector with name of each MC;
///5. Number of final state particles;
///6. Number of cores to be used (be carefull, this option must respect machine core number available);
///7. Name of the output file to store FastME analysis results.

void FastME(string Data_Path, vector<string> MCs, TString TreeName, vector<string> MC_Names,
		     Int_t N_FSParticles, UInt_t N_Cores, TString out_name){
  
  ///Getting some numbers
  TFile *fData = TFile::Open((TString)Data_Path);
  TTreeReader tmpReader(TreeName,fData);
  Int_t nData = tmpReader.GetEntries(true);

  Int_t N_MCT = MC_Names.size();
  Int_t N_MC = MCs.size();
  Int_t NMCEV[N_MC];
  for(Int_t ne=0; ne<N_MC; ne++){
    TFile *tmpf = TFile::Open((TString)MCs[ne]);
    TTreeReader tmpReader(TreeName,tmpf);
    NMCEV[ne] = tmpReader.GetEntries(true);
    tmpf->Close();
  }

  ///TProcPool declaration to objects to be analised  
  auto workItem = [fData, TreeName, nData, N_MCT, N_FSParticles](TTreeReader &tread) -> TObject* {
    TStopwatch t2;
        
    ///Defines 2D histogram to stores minimum distances
    TH2D *mdists = new TH2D("mdists","Minimum Data-MC distances found",nData,0,nData,N_MCT,0,N_MCT);
    mdists->SetDirectory(0);    
  
    ///Addresses the MC branches to be used
    TTreeReaderValue<Int_t>    McType(tread, "McType"); ///McType for Signal=0 and Background >0
    TTreeReaderArray<Int_t>    McId(tread, "ParticleId");
    TTreeReaderArray<Double_t> McPt(tread, "ParticlePt");
    TTreeReaderArray<Double_t> McEta(tread, "ParticleEta");
    
    ///Addresses the Data branches to be used
    TTreeReader refReader(TreeName,fData);
    TTreeReaderArray<Int_t>    DataId(refReader, "ParticleId");
    TTreeReaderArray<Double_t> DataPt(refReader, "ParticlePt");
    TTreeReaderArray<Double_t> DataEta(refReader, "ParticleEta");


    ///Loop on Data events
    for(Int_t dt=0; dt<nData; dt++){
    
      if( (dt!= 0 && dt%(nData/10) == 0) || (nData-dt) == 1 ){ 
	cout<< Form(":: [Remaining Data]:  %i\t\t[Elapsed Time]:  ",nData-dt);
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
        if( debug ) cout<<"dt_el: "<<dt_el<<"\tmc_el: "<<mc_el<<"\tdt_mu: "<<dt_mu<<"\tmc_mu: "<<mc_mu<<endl;
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
  cout<<"=============================================================================================="<<endl;
  cout<<"::::::::::::::::::::::::::[ Fast Matrix Element Analysis Started]:::::::::::::::::::::::::::::"<<endl;
  cout<<"=============================================================================================="<<endl;
  cout<< Form(":: Data Events:      %i",nData) <<endl;
  cout<< Form(":: MC Samples:       %i\t[",N_MC); 
            for(Int_t ne=0; ne<N_MC; ne++){
              if( ne<(N_MC-1) ) cout<< ne << "= " << NMCEV[ne] << ",  ";
              if( ne==(N_MC-1) ) cout<< ne << "= " << NMCEV[ne] << "]" <<endl;
            }
  cout<< Form(":: Final State:      %i  Objects",N_FSParticles) <<endl;
  cout<< Form(":: Cores to Use:     %i  Cores",N_Cores) <<endl;
  cout<<"--------------------------------------------------------------------------------------------------"<<endl;
  cout<<" [Analysing events...]"<<endl;  
  ///--------------------------------------------------------------------------------------------------------------

  ///Calls analysis through TProcPool
  TProcPool workers(N_Cores);
  auto f_hist = (TH2D*)workers.ProcTree(MCs, workItem);
  f_hist->GetXaxis()->SetTitle("Data Events");
  for(int mcn=0; mcn<int(MC_Names.size()); mcn++)
    f_hist->GetYaxis()->SetBinLabel(mcn+1,(TString)MC_Names[mcn]);
  

  ///_______________________ Compute discriminant from MDMCED _____________________________________________________
  cout<<":: [Distance Computing Time]: "; t.Stop(); t.Print(); t.Continue();
  cout<<"\n::::::::::::::::::::::::::::::::[ Computing discriminant ]::::::::::::::::::::::::::::::::::::::::"<<endl;
  ///--------------------------------------------------------------------------------------------------------------
  Int_t Event;
  Double_t G_PsbD_MinDist, G_PsbD_Media;
  vector<Double_t> PsbD_MinDist, PsbD_Media;
  TTree *tree = new TTree("FastME","Fast Matrix Element Analysis Results");
  tree->SetDirectory(0);
  tree->Branch("Event",&Event,"Event/I");
  tree->Branch("G_PsbD_MinDist",&G_PsbD_MinDist,"G_PsbD_MinDist/D");
  //tree->Branch("G_PsbD_Media",&G_PsbD_Media,"G_PsbD_Media/D");
  tree->Branch("PsbD_MinDist","vector<Double_t> PsbD_MinDist",&PsbD_MinDist);
  //tree->Branch("PsbD_Media","vector<Double_t> PsbD_Media",&PsbD_Media);
  for(Int_t data=0; data<nData; data++){
    Event = data;
    if(data%(nData/10) == 0)
      cout<< Form(":: [Remaining Data]:   %i",nData-data) <<endl;

    Double_t min_dr_sig = f_hist->GetBinContent(data+1,1);
    Double_t min_dr_bkg = 1.E15;
    PsbD_MinDist.clear();
    PsbD_Media.clear();
    ///Finds closet MC
    for(Int_t mcs=1; mcs<N_MCT; mcs++){
      PsbD_MinDist.push_back( PsbD(min_dr_sig, f_hist->GetBinContent(data+1,mcs+1)) );
      //PsbD_Media.push_back( PsbD(min_dr_sig, f_hist->GetBinContent(data+1,mcs+1)) );
      if( f_hist->GetBinContent(data+1,mcs+1) < min_dr_bkg )
	min_dr_bkg = f_hist->GetBinContent(data+1,mcs+1);
    }

    G_PsbD_MinDist = PsbD(min_dr_sig, min_dr_bkg);
    //G_PsbD_Media   = PsbD(min_dr_sig, min_dr_bkg);
    if( debug ) cout<<"GSigMin: "<< min_dr_sig <<"\t\tGBkgMin: "<< min_dr_bkg << endl;
    tree->Fill();
  }
  
  ///________________________________ Stoping timming ________________________________________________________
  cout<<"\n::::::::::::::::::::::::::::::::::::[ Process Finished ]::::::::::::::::::::::::::::::::::::::::::::"<<endl;
  cout<<":: [Analysis Total Time]: "; t.Stop(); t.Print();
  cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n"<<endl;
  ///---------------------------------------------------------------------------------------------------------

  ///Saving FastME results
  TFile *tmp = TFile::Open(out_name+".root","recreate");
  f_hist->Write();
  tree->Write();
  tmp->Close();
  
}
///********************************************************************************************************************