#include "MHead.h"

TFile *file_input = nullptr;
TFile *file_output = nullptr;

void GetNonuniformEffect(string input_file = "input/non_uniform.root",
                         string output_file = "output/nonuniform_effect.root") {
                          SetPolttingStyle();
  file_input = new TFile(input_file.c_str(), "READ");

  if (output_file != "") {
    file_output = new TFile(output_file.c_str(), "RECREATE");
  }

  vector<string> vec_str_name_obj;
  vector<TObject *> vec_obj = GetObjectRecursive(file_input, vec_str_name_obj);

  vector<string> vec_str_name_obj_output;
  vector<TH1 *> vec_h1 = GetObjectVector<TH1>(
      vec_obj, vec_str_name_obj, vec_str_name_obj_output, (TString) "Phi_Eta");

  for (auto &h1 : vec_h1)
    h1->Print();

  TH2F *h2_input = (TH2F *)vec_h1[2];
  cout << h2_input->GetXaxis()->GetNbins() << endl;
  h2_input->RebinX(100);
// h2_input->GetYaxis()->SetRangeUser(0., 2.*TMath::Pi());
  // h3_input->RebinY();
  h2_input->Draw();
  TH2F *h2_output = (TH2F *)h2_input->Clone("h2_nonuniform_effect");
  h2_output->GetXaxis()->SetTitle("#eta");
  h2_output->GetYaxis()->SetTitle("#phi");

  TH2F* h2_2draw = new TH2F("h2_2draw","",100,-0.9,0.9,100,0,2.*TMath::Pi());

  // h2_output->GetZaxis()->SetTitle("Multiplicity");

  for (int iBinEta = 1; iBinEta <= h2_input->GetXaxis()->GetNbins();
       iBinEta++) {
    TH1D *h1_input = (TH1D *)h2_input->ProjectionY(
        Form("h1_input_%d", iBinEta), iBinEta, iBinEta);
    h1_input->Scale(h1_input->GetNbinsX() / h1_input->Integral()/2.);
    for (int iBinPhi = 1; iBinPhi <= h2_input->GetYaxis()->GetNbins();
         iBinPhi++) {
      double content = h1_input->GetBinContent(iBinPhi);
      h2_output->SetBinContent(iBinEta, iBinPhi, content);
      h2_output->SetBinError(iBinEta, iBinPhi, 0.);
    }
  }
  TCanvas* c0 = new TCanvas("c0","",800,600);
  h2_output->GetYaxis()->SetRangeUser(0,2.*TMath::Pi());
  h2_output->GetXaxis()->SetRangeUser(-0.9,0.9);
  h2_output->Draw("colz");

  TH1D* h1_left = (TH1D*)h2_output->ProjectionY("h1_left",1,1);
  TH1D* h1_right = (TH1D*)h2_output->ProjectionY("h1_right",2,2);
  h1_left->SetLineColor(kRed);
  TCanvas* c1 = new TCanvas("c1","",800,600);
  h1_left->GetYaxis()->SetRangeUser(0,1.5);
  h1_left->Draw();
  h1_right->Draw("same");
  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->SetMargin(0.9);
  leg->SetLineColor(0);
  leg->AddEntry(h1_left,"#eta<0","l");
  leg->AddEntry(h1_right,"#eta>0","l");
  leg->Draw("same");
  // if (file_output != nullptr) {
  //   file_output->cd();
  //   h2_input->Write();
  //   h2_output->Write();
  //   file_output->Close();
  // }
}
