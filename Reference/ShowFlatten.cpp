#include "MHead.h"

#define verbose 2
#define roothttp

void ShowFlatten(string path_file_input =
                     "/data/work/work_alice/Analysis/JpsiFlow/input/"
                     "AnalysisResults_new2.root") {
  SetPolttingStyle();

  TFile *file_input = new TFile(path_file_input.c_str(), "READ");
  vector<string> vec_str_name_obj;
  vector<TObject *> vec_obj = GetObjectRecursive(file_input, vec_str_name_obj);
  vector<TH1 *> vec_h1_reference =
      GetObjectVector<TH1>(vec_obj, vec_str_name_obj, vector<string>{},
                           (TString) "fh3EtaPhiMultReferenceCorrected");
  TH3F* h3_reference = (TH3F*)vec_h1_reference[0];
  cout << h3_reference->GetXaxis()->GetNbins() << endl;
  cout << h3_reference->GetYaxis()->GetNbins() << endl;
  cout << h3_reference->GetZaxis()->GetNbins() << endl;
  h3_reference->RebinX(90);
  TH1D* h1_reference1 = (TH1D*)h3_reference->ProjectionY("h1_reference1",1,1);
  TH1D* h1_reference2 = (TH1D*)h3_reference->ProjectionY("h1_reference2",2,2);

  h1_reference1->Draw();
  h1_reference1->SetLineColor(kRed);
  h1_reference2->Draw("same");
  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(h1_reference1,"#eta < 0","l");
  leg->AddEntry(h1_reference2,"#eta > 0","l");
  leg->Draw();
}