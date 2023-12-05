#include "MHead.h"

TFile *file_input = nullptr;
TFile *file_output = nullptr;

void httpserver(
    TString input_file =
        "/data/work/work_alice/FlowReference/AnalysisResults.root") {
  file_input = new TFile(input_file, "READ");

  vector<string> vec_str_name_obj;
  vector<TObject *> vec_obj = GetObjectRecursive(file_input, vec_str_name_obj);

  // Draw all th1 in vec_h1 in separated canvas and save it
  for (auto name : vec_str_name_obj)
    cout << name << endl;

  THttpServer *serv = new THttpServer(Form("http:8090?top=%s", "job"));
  serv->SetReadOnly(kFALSE);
  gBenchmark->Start("job");

  // register the histogram in the server
  for (int i = 0; i < (int)vec_obj.size(); i++) {
    if (vec_obj[i]->InheritsFrom("TH1")) {
      TH1 *h1 = (TH1 *)vec_obj[i];
      serv->Register(vec_str_name_obj[i].c_str(), h1);
    }
  }
}