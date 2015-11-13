{
  ///To plot discriminant distributions
  
 
  TFile *fsig = TFile::Open("SigLikeData_QCDbkg.root");
  TFile *fbkg = TFile::Open("BkgLikeData_QCDbkg.root");

  TTree *tsig = (TTree*)fsig->Get("FastME_Results");
  TTree *tbkg = (TTree*)fbkg->Get("FastME_Results");
  
  TH1D *sig_psbD = new TH1D("sig_psbD","",50,0,1);
  sig_psbD->SetLineWidth(2);
  //TH1D *sig_psbW = new TH1D("sig_psbW","",50,0,1);
  //sig_psbW->SetLineWidth(2);
  TH1D *bkg_psbD = new TH1D("bkg_psbD","",50,0,1); 
  bkg_psbD->SetLineWidth(2);
  bkg_psbD->SetLineColor(kRed);
  bkg_psbD->GetXaxis()->SetTitle("P_{SB}(Distance)");
  bkg_psbD->GetYaxis()->SetTitle("Events/0.02 (Normalized)");
  //TH1D *bkg_psbW = new TH1D("bkg_psbW","",50,0,1);
  //bkg_psbW->SetLineWidth(2);
  //bkg_psbW->SetLineColor(kRed);
  //bkg_psbW->GetXaxis()->SetTitle("P_{SB}(Weight)");
  //bkg_psbW->GetYaxis()->SetTitle("Events/0.02");

  tsig->Project("sig_psbD","PSB_Distance");
  //tsig->Project("sig_psbW","PSB_Weight");
  tbkg->Project("bkg_psbD","PSB_Distance");
  //tbkg->Project("bkg_psbW","PSB_Weight");

  Double_t SigPsbDistance, BkgPsbDistance;
  tsig->SetBranchAddress("PSB_Distance",&SigPsbDistance);
  tbkg->SetBranchAddress("PSB_Distance",&BkgPsbDistance);
  if( tsig->GetEntries() != tbkg->GetEntries() ) cout<<"!!Different number of events!!"<<endl;
  int nevents = tsig->GetEntries();
  
  Double_t cutoff;
  const int discret = 100;
  float TPR[discret], FPR[discret], TP, FP, TN, FN;
  for(int j=0; j<discret; j++){
    cutoff = j/float(discret);
    TP = FP = TN = FN = 0;
    for(int i=0; i<nevents; i++){
      tsig->GetEntry(i);
      if( SigPsbDistance > cutoff ) TP++;
      if( SigPsbDistance < cutoff ) FN++;
      tbkg->GetEntry(i);
      if( BkgPsbDistance > cutoff ) FP++;
      if( BkgPsbDistance < cutoff ) TN++;
    }
    TPR[j] = TP/float(TP + FN);
    FPR[j] = FP/float(FP + TN);
    
    //if(TPR[j] == 0.95 || TPR[j] == 0.90 || TPR[j] == 0.85 || TPR[j] == 0.80)
      //cout<<"Cutoff= "<<cutoff<<"   TPR= "<<TPR[j]<<"   FPR= "<<FPR[j]<<endl;
  }
  TGraph *roc = new TGraph(discret,FPR,TPR);
  roc->SetMarkerStyle(4);
  roc->SetMarkerSize(0.9);
  roc->GetXaxis()->SetTitle("False Positive Rate");
  roc->GetYaxis()->SetTitle("True Positive Rate");
  
  TLine *l1 = new TLine(0,0,1,1);
  l1->SetLineColor(kBlack);
  l1->SetLineWidth(2);
  l1->SetLineStyle(2);
  
  TLegend *leg = new TLegend(0.2,0.75,0.5,0.85);
  leg->AddEntry(sig_psbD,"VBF Sig","l");
  leg->AddEntry(bkg_psbD,"VBF QCD Bkg","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);

  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas();
  cv->Divide(2,1);
  cv->cd(1);
  bkg_psbD->DrawNormalized();
  sig_psbD->DrawNormalized("same");
  leg->Draw();
  cv->cd(2);
  roc->Draw("AP");
  l1->Draw();

 
///To plot Efficiency vs Purity
/*
  Int_t SigFinalState, BkgFinalState;
  Double_t psbD_Sig, psbW_Sig, psbD_Bkg, psbW_Bkg;
  TFile *fme_sig = TFile::Open("ggH4l_14TeV_powheg_FMEe4MC_Scale_dEta5_dPt50.root");
  //TFile *fme_sig = TFile::Open("ggH4l_14TeV_madgraph_FMEe4MC_Scale_dEta5_dPt50.root");
  TTree *tsig = (TTree*)fme_sig->Get("FastME_Results");
  tsig->SetBranchAddress("PSB_Distance",&psbD_Sig);
  tsig->SetBranchAddress("PSB_Weight",&psbW_Sig);
  int nsig = tsig->GetEntries();
  
  TFile *fme_bkg = TFile::Open("ZZ4l_14TeV_powheg_FMEe4MC_Scale_dEta5_dPt50.root");
  //TFile *fme_bkg = TFile::Open("ZZ4l_14TeV_madgraph_FMEe4MC_Scale_dEta5_dPt50.root");
  TTree *tbkg = (TTree*)fme_bkg->Get("FastME_Results");
  tbkg->SetBranchAddress("PSB_Distance",&psbD_Bkg);
  tbkg->SetBranchAddress("PSB_Weight",&psbW_Bkg);
  int nbkg = tbkg->GetEntries();
  
  TH2D *eff_pur_D = new TH2D("eff_pur_D","",1000,0,1,1000,0,1);
  eff_pur_D->SetMarkerColor(kBlack);
  TH2D *eff_pur_W = new TH2D("eff_pur_W","",1000,0,1,1000,0,1);
  eff_pur_W->SetMarkerColor(kRed);

  int npoints = 100, n1=0, n2=0;
  float purityD[npoints], efficiencyD[npoints], purityW[npoints], efficiencyW[npoints], x1[npoints], x2[npoints];
  if(nsig != nbkg) cout<<"[Warning]: Different event number coming from two sources..."<<endl;
  for(int cut=0; cut<npoints; cut++){
    float vcut = cut/float(npoints);

    float sig_acp_psbD = 0, sig_acp_psbW = 0;
    float bkg_acp_psbD = 0, bkg_acp_psbW = 0;
    for(int i=0; i<nsig; i++){
      tsig->GetEntry(i);
      tbkg->GetEntry(i);
      
      ///Less than cut is signal-like
      ///(Efficiency)
      if(psbD_Sig > vcut){
	if(SigFinalState == 0) sig_acp_psbD += 5.05519E+01/3350.;
	if(SigFinalState == 1) sig_acp_psbD += 5.05519E+01/3350.;
	if(SigFinalState == 2) sig_acp_psbD += 5.05519E+01/3350.;
      }
      if(psbD_Bkg > vcut){
	if(BkgFinalState == 0) bkg_acp_psbD += 1.88348E-02/3350.;
	if(BkgFinalState == 1) bkg_acp_psbD += 1.88348E-02/3350.;
	if(BkgFinalState == 2) bkg_acp_psbD += 4.27071E-02/3350.;
      }
      if(psbW_Sig > vcut){
	if(SigFinalState == 0) sig_acp_psbW += 5.05519E+01/3350.;
	if(SigFinalState == 1) sig_acp_psbW += 5.05519E+01/3350.;
	if(SigFinalState == 2) sig_acp_psbW += 5.05519E+01/3350.;
      }
      if(psbW_Bkg > vcut){
	if(BkgFinalState == 0) bkg_acp_psbW += 1.88348E-02/3350.;
	if(BkgFinalState == 1) bkg_acp_psbW += 1.88348E-02/3350.;
	if(BkgFinalState == 2) bkg_acp_psbW += 4.27071E-02/3350.;
      }
    }
    
    ///Purity
    if(sig_acp_psbD != 0){
      x1[n1] = vcut;
      efficiencyD[n1] = (sig_acp_psbD + bkg_acp_psbD)/float(nsig + nbkg);
      purityD[n1] = sig_acp_psbD/float(sig_acp_psbD + bkg_acp_psbD);
      if(cut % 2 == 0 && vcut >= 0.34 && vcut <= 0.68)
	eff_pur_D->Fill(purityD[n1],efficiencyD[n1],vcut);
      n1 += 1;
    }
    if(sig_acp_psbW != 0){
      x2[n2] = vcut;
      efficiencyW[n2] = (sig_acp_psbW + bkg_acp_psbW)/float(nsig + nbkg);
      purityW[n2] = sig_acp_psbW/float(sig_acp_psbW + bkg_acp_psbW);
      if(cut % 4 == 0)
	eff_pur_W->Fill(purityW[n2],efficiencyW[n2],vcut);
      n2 += 1;
    }
  }

  ///Efficiency and Purity vs Discriminant Cut

  TGraph *effD_cut = new TGraph(n1,x1,efficiencyD);
  effD_cut->SetMarkerColor(kBlack);
  effD_cut->SetMarkerStyle(20);
  TGraph *purD_cut = new TGraph(n1,x1,purityD);
  purD_cut->SetMarkerColor(kRed);
  purD_cut->SetMarkerStyle(20);
  
  TGraph *effW_cut = new TGraph(n2,x2,efficiencyW);
  effW_cut->SetMarkerColor(kViolet);
  effW_cut->SetMarkerStyle(20);
  TGraph *purW_cut = new TGraph(n2,x2,purityW);
  purW_cut->SetMarkerColor(kBlue);
  purW_cut->SetMarkerStyle(20);

  TMultiGraph DistAnalysis("DistAnalysis","");
  DistAnalysis.Add(effD_cut);
  DistAnalysis.Add(purD_cut);
  DistAnalysis.Add(effW_cut);
  DistAnalysis.Add(purW_cut);
  DistAnalysis.Draw("AP");
  
  TLegend *leg = new TLegend(0.6,0.85,0.9,0.95);
  leg->AddEntry(effD_cut,"P_{SB}(Distance) Efficiency","p");
  leg->AddEntry(purD_cut,"P_{SB}(Distance) Purity","p");
  leg->AddEntry(effW_cut,"P_{SB}(Weight) Efficiency","p");
  leg->AddEntry(purW_cut,"P_{SB}(Weight) Purity","p");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->Draw();
*/

 ///Efficiency vs Purity without show cuts
 /*
 TGraph *effD_vs_purD = new TGraph(n1,purityD,efficiencyD);
 effD_vs_purD->SetMarkerColor(kBlack);
 TGraph *effW_vs_purW = new TGraph(n2,purityW,efficiencyW);
 effW_vs_purW->SetMarkerColor(kRed);
 
 TMultiGraph DistAnalysis("DistAnalysis","");
 DistAnalysis.Add(effD_vs_purD);
 DistAnalysis.Add(effW_vs_purW);
 DistAnalysis.Draw("AP");
 */ 

 ///With cuts
 //eff_pur_D->Draw("Text");
 //eff_pur_W->Draw("Text,same");
 
 
 ///Plot Discriminant vs ZZ mass
/* 
 Int_t sig_fs, bkg_fs;
 Double_t sig_ZZ_mass, sig_PsbD, sig_PsbW, bkg_ZZ_mass, bkg_PsbD, bkg_PsbW;
 //TFile *fsig = TFile::Open("NtuplesPowhegPythia/ggH4l_14TeV_powheg-pythia8.root");
 TFile *fsig = TFile::Open("NtuplesMadGraph/gg4l_madgraph.root");
 TTree *tsig = (TTree*)fsig->Get("LHE_Tree");
 tsig->SetBranchAddress("ZZ_mass",&sig_ZZ_mass);
 int nsig = tsig->GetEntries();

 //TFile *sig_fme = TFile::Open("ggH4l_14TeV_powheg_FMEe4MC_Scale_dEta5_dPt50.root");
 TFile *sig_fme = TFile::Open("ggH4l_14TeV_madgraph_FMEe4MC_Scale_dEta5_dPt50.root");
 TTree *tsig_fme = (TTree*)sig_fme->Get("FastME_Results");
 tsig_fme->SetBranchAddress("PSB_Distance",&sig_PsbD);
 tsig_fme->SetBranchAddress("PSB_Weight",&sig_PsbW);

 
 //TFile *fbkg = TFile::Open("NtuplesPowhegPythia/ZZTo4l_14TeV_powheg.root");
 TFile *fbkg = TFile::Open("NtuplesMadGraph/qq4l_madgraph.root");
 TTree *tbkg = (TTree*)fbkg->Get("LHE_Tree");
 tbkg->SetBranchAddress("ZZ_mass",&bkg_ZZ_mass);
 int nbkg = tbkg->GetEntries();

 //TFile *bkg_fme = TFile::Open("ZZ4l_14TeV_powheg_FMEe4MC_Scale_dEta5_dPt50.root");
 TFile *bkg_fme = TFile::Open("ZZ4l_14TeV_madgraph_FMEe4MC_Scale_dEta5_dPt50.root");
 TTree *tbkg_fme = (TTree*)bkg_fme->Get("FastME_Results");
 tbkg_fme->SetBranchAddress("PSB_Distance",&bkg_PsbD);
 tbkg_fme->SetBranchAddress("PSB_Weight",&bkg_PsbW);
 
  
  TH2D *sig_psbD_ZZmass = new TH2D("sig_psbD_ZZmass","",75,50,800,50,0,1);
  sig_psbD_ZZmass->SetMarkerColor(kBlue);
  sig_psbD_ZZmass->GetXaxis()->SetTitle("ZZ Mass (GeV)");
  sig_psbD_ZZmass->GetYaxis()->SetTitle("P_{SB}(Distance)");
  sig_psbD_ZZmass->SetTitle("gg #rightarrow ZZ #rightarrow 4l");
  
  TH2D *sig_psbW_ZZmass = new TH2D("sig_psbW_ZZmass","",75,50,800,50,0,1);
  sig_psbW_ZZmass->SetMarkerColor(kBlue);
  sig_psbW_ZZmass->GetXaxis()->SetTitle("ZZ Mass (GeV)");
  sig_psbW_ZZmass->GetYaxis()->SetTitle("P_{SB}(Weight)");
  sig_psbW_ZZmass->SetTitle("gg #rightarrow ZZ #rightarrow 4l");
  
  TH2D *bkg_psbD_ZZmass = new TH2D("bkg_psbD_ZZmass","",75,50,800,50,0,1);
  bkg_psbD_ZZmass->SetMarkerColor(kRed);
  bkg_psbD_ZZmass->GetXaxis()->SetTitle("ZZ Mass (GeV)");
  bkg_psbD_ZZmass->GetYaxis()->SetTitle("P_{SB}(Distance)");
  bkg_psbD_ZZmass->SetTitle("qq #rightarrow ZZ #rightarrow 4l");

  TH2D *bkg_psbW_ZZmass = new TH2D("bkg_psbW_ZZmass","",75,50,800,50,0,1);
  bkg_psbW_ZZmass->SetMarkerColor(kRed);
  bkg_psbW_ZZmass->GetXaxis()->SetTitle("ZZ Mass (GeV)");
  bkg_psbW_ZZmass->GetYaxis()->SetTitle("P_{SB}(Weight)");
  bkg_psbW_ZZmass->SetTitle("qq #rightarrow ZZ #rightarrow 4l");

  
  for(int i=0; i<nsig; i++){
    tsig->GetEntry(i);
    tsig_fme->GetEntry(i);
    sig_psbD_ZZmass->Fill(sig_ZZ_mass,sig_PsbD);
    sig_psbW_ZZmass->Fill(sig_ZZ_mass,sig_PsbW);
  }
  for(int i=0; i<nbkg; i++){
    tbkg->GetEntry(i);
    tbkg_fme->GetEntry(i);
    bkg_psbD_ZZmass->Fill(bkg_ZZ_mass,bkg_PsbD);
    bkg_psbW_ZZmass->Fill(bkg_ZZ_mass,bkg_PsbW);    
  }
  
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas();
  cv->Divide(2,2);
  cv->cd(1);
  sig_psbD_ZZmass->Draw("Colz");
  cv->cd(2);
  sig_psbW_ZZmass->Draw("Colz");
  cv->cd(3);
  bkg_psbD_ZZmass->Draw("Colz");
  cv->cd(4);
  bkg_psbW_ZZmass->Draw("Colz");
*/
}
