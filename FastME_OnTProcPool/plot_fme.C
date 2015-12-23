{
  ///To plot discriminant distributions
  
 
  TFile *fsig = TFile::Open("fme_results_sig.root");
  TFile *fbkg = TFile::Open("fme_results_bkg.root");
  
  TH2D *h2sig = (TH2D*)fsig->Get("mdists");
  TH2D *h2bkg = (TH2D*)fbkg->Get("mdists");
  
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,1);
  c1->cd(1);
  h2sig->Draw("Colz");
  c1->cd(2);
  h2bkg->Draw("Colz");
  
  Double_t GSigPsbDistance, GBkgPsbDistance;
  vector<Double_t> *SigPsbDistance, *BkgPsbDistance;
  TTree *tsig = (TTree*)fsig->Get("FastME");
  TTree *tbkg = (TTree*)fbkg->Get("FastME");
  tsig->SetBranchAddress("G_PsbD_MinDist",&GSigPsbDistance);
  tsig->SetBranchAddress("PsbD_MinDist",&SigPsbDistance);
  tbkg->SetBranchAddress("G_PsbD_MinDist",&GBkgPsbDistance);
  tbkg->SetBranchAddress("PsbD_MinDist",&BkgPsbDistance);
  if( tsig->GetEntries() != tbkg->GetEntries() ) cout<<"!!Different number of events!!"<<endl;
  int nevents = tsig->GetEntries();
  
  TH1D *sig_psbD = new TH1D("sig_psbD","",50,0,1);
  TH1D *bkg_psbD = new TH1D("bkg_psbD","",50,0,1); 
  TH1D *bkg1_psbD = new TH1D("bkg1_psbD","",50,0,1); 
  TH1D *bkg2_psbD = new TH1D("bkg2_psbD","",50,0,1); 
  bkg_psbD->SetLineColor(kRed);
  bkg1_psbD->SetLineColor(kBlue);
  bkg2_psbD->SetLineColor(kGreen);
  bkg_psbD->GetXaxis()->SetTitle("P_{SB}(Distance)");
  bkg_psbD->GetYaxis()->SetTitle("Events/0.02 (Normalized)");
  tsig->Project("sig_psbD","G_PsbD_MinDist");
  tsig->Project("bkg1_psbD","PsbD_MinDist[0]");
  tsig->Project("bkg2_psbD","PsbD_MinDist[1]");
  tbkg->Project("bkg_psbD","G_PsbD_MinDist");
  
  Double_t cutoff, integral=0;
  const int discret = 10000;
  float TPR[discret], FPR[discret], TP, FP, TN, FN;
  for(int j=0; j<discret; j++){
    cutoff = j/float(discret);
    TP = FP = TN = FN = 0;
    for(int i=0; i<nevents; i++){
      tsig->GetEntry(i);
      if( GSigPsbDistance > cutoff ) TP++;
      if( GSigPsbDistance < cutoff ) FN++;
      tbkg->GetEntry(i);
      if( GBkgPsbDistance > cutoff ) FP++;
      if( GBkgPsbDistance < cutoff ) TN++;
    }
    TPR[j] = TP/float(TP + FN);
    FPR[j] = FP/float(FP + TN);
    if(j>0)
     integral += fabs(FPR[j-1]-FPR[j])*TPR[j];
  }
  cout<<"Area in the ROC plot: "<<integral<<endl;
  
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
  leg->AddEntry(bkg_psbD,"VBF Bkg (EW + QCD)","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);

  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas();
  cv->Divide(2,1);
  cv->cd(1);
  sig_psbD->DrawNormalized();
  bkg1_psbD->DrawNormalized("same");
  bkg2_psbD->DrawNormalized("same");
  bkg_psbD->DrawNormalized("same");
  leg->Draw();
  cv->cd(2);
  roc->Draw("AP");
  l1->Draw();
}