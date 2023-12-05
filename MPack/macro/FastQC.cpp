#include "MHead.h"

TFile* file_input = nullptr;
TFile* file_output = nullptr;

void FastQC( TString input_file="/data/work/work_alice/FlowReference/AnalysisResults.root", TString tag="",TString output_file="") {
  file_input = new TFile(input_file, "READ");

  if (output_file != "") {
    file_output = new TFile(output_file, "RECREATE");
  }  

  vector<string> vec_str_name_obj;
  vector<TObject*> vec_obj = GetObjectRecursive(file_input, vec_str_name_obj);

  vector<string> vec_str_name_obj_output;
  vector<TH1*> vec_h1 = GetObjectVector<TH1>(vec_obj, vec_str_name_obj, vec_str_name_obj_output, tag);

  // Draw all th1 in vec_h1 in separated canvas and save it
  for (auto name : vec_str_name_obj)
    cout << name << endl;

  for (auto &h1 : vec_h1) {
    TCanvas *c = new TCanvas();
    if (h1->GetEntries() == 0) {
      cout << "Skip empty histogram: " << h1->GetName() << endl;
      continue;
    }
    gPad->SetLogz();
    h1->Draw("COLZ");
    c->SaveAs(Form("FastQC/%s_%d.pdf", h1->GetName(),GenerateUID()));
    delete c;
  }
  
  if (file_output != nullptr) {
    file_output->cd();
    for (auto &obj : vec_obj) {
      obj->Write();
    }
    file_output->Close();
  }

}
