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
  
  
  TH1D *h1sig = (TH1D*)fsig->Get("hdisc");
  TH1D *h1bkg = (TH1D*)fbkg->Get("hdisc");
  h1bkg->SetLineColor(kRed);
  int n = h1sig->GetNbinsX();
  float tpr[50], fpr[50];
  for(int i=0; i<n; i++){
    float fn=0, tn=0, fp=0, tp=0;
    for(int j=0; j<n; j++){
      if(j<i){
	fn += h1sig->GetBinContent(j);
	tn += h1bkg->GetBinContent(j);
      }
      if(j>i){
	fp += h1bkg->GetBinContent(j);
	tp += h1sig->GetBinContent(j);
      }
    }
    tpr[i] = tp/(tp+fn);
    fpr[i] = fp/(fp+tn);
  }
  
  TGraph *roc = new TGraph(n,fpr,tpr);
  roc->SetMarkerStyle(4);
  roc->SetMarkerSize(0.9);
  roc->GetXaxis()->SetTitle("False Positive Rate");
  roc->GetYaxis()->SetTitle("True Positive Rate");
  
  TLine *l1 = new TLine(0,0,1,1);
  l1->SetLineColor(kBlack);
  l1->SetLineWidth(2);
  l1->SetLineStyle(2);
  
  TLegend *leg = new TLegend(0.2,0.75,0.5,0.85);
  leg->AddEntry(h1sig,"VBF Sig","l");
  leg->AddEntry(h1bkg,"VBF Bkg (EW + QCD)","l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);

  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas();
  cv->Divide(2,1);
  cv->cd(1);
  h1sig->DrawNormalized();
  h1bkg->DrawNormalized("same");
  leg->Draw();
  cv->cd(2);
  roc->Draw("AP");
  l1->Draw();
}
