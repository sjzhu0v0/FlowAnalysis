#include "MHead.h"

void Sigma_Plot() {
  TFile* f = new TFile("/home/sjzhu/Documents/Image/QC/sigma.root");
  TH1F* h_22o = (TH1F*)f->Get("sigma_22o");
  TH1F* h_23u = (TH1F*)f->Get("sigma_23u");

  TCanvas* c = new TCanvas("c", "", 800, 600);
  h_22o->SetMarkerColor(kRed);
  h_22o->SetLineColor(kRed);
  h_23u->Draw();  
  h_22o->Draw("same");

  TCanvas* c2 = new TCanvas("c2", "", 800, 600);
  TH1F* h_ratio = (TH1F*)h_22o->Clone("h_ratio");
  h_ratio->GetXaxis()->SetRangeUser(0,15);
  h_ratio->GetYaxis()->SetRangeUser(0,3.);
  h_ratio->Divide(h_23u);
  h_ratio->Draw();

}