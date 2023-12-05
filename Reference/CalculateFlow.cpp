#include "MHead.h"

TFile *file_input = nullptr;
TFile *file_output = nullptr;

void CalculateFlow(
    TString input_file =
        "/data/work/work_alice/FlowReference/AnalysisResults.root",
    TString tag = "Combined",
    TString output_file = "output/flow_reference.root") {
  file_input = new TFile(input_file, "READ");

  if (output_file != "") {
    file_output = new TFile(output_file, "RECREATE");
  }

  vector<string> vec_str_name_obj;
  vector<TObject *> vec_obj = GetObjectRecursive(file_input, vec_str_name_obj);

  vector<string> vec_str_name_obj_output;
  vector<TH1 *> vec_h1 = GetObjectVector<TH1>(vec_obj, vec_str_name_obj,
                                              vec_str_name_obj_output, tag);

  for (auto &h1 : vec_h1)
    h1->Print();

  vector<TList *> vec_list_output;
  for (int i_bin_mult = 1; i_bin_mult <= vec_h1[0]->GetYaxis()->GetNbins();
       i_bin_mult++) {
    TList *list_output = new TList();
    list_output->SetName(Form("mult_%d", i_bin_mult));
    list_output->SetOwner();
    vec_list_output.push_back(list_output);
  }

  for (int i_bin_mult = 1; i_bin_mult <= vec_h1[0]->GetYaxis()->GetNbins();
       i_bin_mult++) {
    for (int i_hist = 0; i_hist < vec_h1.size(); i_hist += 2) {
      TGraphErrors *gr = new TGraphErrors();
      TGraphErrors *gr_flow = new TGraphErrors();
      gr->SetName(Form("%s_%d", vec_h1[i_hist]->GetName(), i_bin_mult));
      gr->SetTitle(Form("%s_%d;%s;%s", vec_h1[i_hist]->GetTitle(), i_bin_mult,
                        vec_h1[i_hist]->GetXaxis()->GetTitle(),
                        vec_h1[i_hist]->GetZaxis()->GetTitle()));
      gr_flow->SetName(Form("%s_%d", vec_h1[i_hist]->GetName(), i_bin_mult));
      gr_flow->SetTitle(Form("%s_%s_%d;%s;%s", vec_h1[i_hist]->GetTitle(), "v2",
                             i_bin_mult, vec_h1[i_hist]->GetXaxis()->GetTitle(),
                             "v_{2}"));
      TProfile2D *prof2d_im = (TProfile2D *)vec_h1[i_hist];
      TProfile2D *prof2d_re = (TProfile2D *)vec_h1[i_hist + 1];
      for (int i_bin_x_prof2d = 1;
           i_bin_x_prof2d <= prof2d_im->GetXaxis()->GetNbins();
           i_bin_x_prof2d++) {
        int i_bin_prof2d = prof2d_im->GetBin(i_bin_x_prof2d, i_bin_mult);
        double im = prof2d_im->GetBinContent(i_bin_prof2d);
        double re = prof2d_re->GetBinContent(i_bin_prof2d);
        double err_im = prof2d_im->GetBinError(i_bin_prof2d);
        double err_re = prof2d_re->GetBinError(i_bin_prof2d);
        double *result = ProessComplex(im, re, err_im, err_re);
        double *result_flow =
            ProcessPowerFunction(result[0], 1.0 / 2.0, result[1]);
        gr->SetPoint(i_bin_x_prof2d - 1,
                     prof2d_im->GetXaxis()->GetBinLowEdge(i_bin_x_prof2d),
                     result[0]);
        gr->SetPointError(i_bin_x_prof2d - 1, 0, result[1]);
        gr_flow->SetPoint(i_bin_x_prof2d - 1,
                          prof2d_im->GetXaxis()->GetBinLowEdge(i_bin_x_prof2d),
                          result_flow[0]);
        gr_flow->SetPointError(i_bin_x_prof2d - 1, 0, result_flow[1]);
      }
      vec_list_output[i_bin_mult - 1]->Add(gr);
      vec_list_output[i_bin_mult - 1]->Add(gr_flow);
    }
  }

  if (file_output != nullptr) {
    file_output->cd();
    for (auto &list_output : vec_list_output) {
      list_output->Write();
    }
    file_output->Close();
  }
}