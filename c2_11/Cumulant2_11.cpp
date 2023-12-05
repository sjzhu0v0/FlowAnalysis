#include "MHead.h"

#define verbose 2
// #define roothttp

void Cumulant2_11(
    string path_file_input = "/data/work/work_alice/Analysis/JpsiFlow/input/"
                             "AnalysisResults_new.root",
    string path_file_jpsiReco = "/data/work/work_alice/Analysis/JpsiFlow/c2_11/"
                                "output/JpsiReco.root",
    string output_file =
        "/data/work/work_alice/Analysis/JpsiFlow/c2_11/output/JpsiFlow.root") {
#ifdef roothttp
  THttpServer *server = new THttpServer("http:8090?top=job");
  server->SetReadOnly(kFALSE);
  gBenchmark->Start("job");
#endif

  TFile *file_input = new TFile(path_file_input.c_str(), "READ");
  TFile *file_jpsiReco = new TFile(path_file_jpsiReco.c_str(), "READ");
  TFile *file_output = new TFile(output_file.c_str(), "RECREATE");

  vector<string> vec_str_name_obj;
  vector<TObject *> vec_obj = GetObjectRecursive(file_input, vec_str_name_obj);
  vector<TH1 *> vec_h1_reference =
      GetObjectVector<TH1>(vec_obj, vec_str_name_obj, vector<string>{},
                           (TString) "ProfileFlowReference");
  vector<TH1 *> vec_h1_target = GetObjectVector<TH1>(
      vec_obj, vec_str_name_obj, vector<string>{}, (TString) "Target");

  vector<string> vec_str_name_obj_jpsiReco;
  vector<TObject *> vec_obj_jpsiReco =
      GetObjectRecursive(file_jpsiReco, vec_str_name_obj_jpsiReco);
  vector<TH1F *> vec_h1_jpsiReco =
      GetObjectVector<TH1F>(vec_obj_jpsiReco, vec_str_name_obj_jpsiReco,
                            vector<string>{}, (TString) "");

  cout << endl;
  for (int i = 0; i < vec_h1_reference.size(); i++) {
    cout << "index:\t" << i << "\t\t";
    vec_h1_reference[i]->Print();
  }

  cout << endl;
  for (int i = 0; i < vec_h1_target.size(); i++) {
    cout << "index:\t" << i << "\t\t";
    vec_h1_target[i]->Print();
  }

  cout << endl;
  for (auto &h1 : vec_h1_jpsiReco)
    h1->Print();

  TH1F *h1_signal2Total_Jpsi =
      (TH1F *)vec_h1_jpsiReco[0]->Clone("h1_signal2Total_Jpsi");
  h1_signal2Total_Jpsi->SetTitle(
      "h1_signal2Total_Jpsi;#it{p}_{T} (GeV/#it{c});#it{S}/#it{N}");

  for (int i_bin_pt = 1; i_bin_pt <= vec_h1_jpsiReco[0]->GetXaxis()->GetNbins();
       i_bin_pt++) {
    double signal = vec_h1_jpsiReco[0]->GetBinContent(i_bin_pt);
    double total = vec_h1_jpsiReco[1]->GetBinContent(i_bin_pt);
    double signal2Total = signal / total;
    double err_signal = vec_h1_jpsiReco[0]->GetBinError(i_bin_pt);
    double err_total = vec_h1_jpsiReco[1]->GetBinError(i_bin_pt);
    // (t-2s)/t^3*deltaS^2+s^2/t^4*deltaT^2
    double *result = ProcessDivision(signal, total, err_signal, err_total);
    h1_signal2Total_Jpsi->SetBinContent(i_bin_pt, result[0]);
    h1_signal2Total_Jpsi->SetBinError(i_bin_pt, result[1]);
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  h1_signal2Total_Jpsi->Draw();

  // profile: {Im, Re} {reference_c11, target_c2_11, target_bg_c2_11}
  // Im : reference_c11
  TProfile2D *prof2d_reference_c11_im = (TProfile2D *)vec_h1_reference[4];
  // Re : reference_c11
  TProfile2D *prof2d_reference_c11_re = (TProfile2D *)vec_h1_reference[5];
  // Im : target_c2_11
  TProfile2D *prof2d_target_c2_11_im = (TProfile2D *)vec_h1_target[0];
  // Re : target_c2_11
  TProfile2D *prof2d_target_c2_11_re = (TProfile2D *)vec_h1_target[1];
  // Im : target_bg_c2_11
  TProfile2D *prof2d_target_bg_c2_11_im = (TProfile2D *)vec_h1_target[24];
  // Re : target_bg_c2_11
  TProfile2D *prof2d_target_bg_c2_11_re = (TProfile2D *)vec_h1_target[25];

  int n_bins_deltaEta = prof2d_target_c2_11_im->GetXaxis()->GetNbins();
  int n_bins_pt = prof2d_target_c2_11_im->GetYaxis()->GetNbins();
  double *bins_deltaEta =
      prof2d_target_c2_11_im->GetXaxis()->GetXbins()->fArray;
  double *bins_pt = prof2d_target_c2_11_im->GetYaxis()->GetXbins()->fArray;
  // target h2{delta_eta;pt}: c2_11 c2 v2
  // target_bg h2{delta_eta;pt}: c2_11 c2 v2
  // reference h1{delta_eta}: c11 v1
  TH2F *h2_target_c2_11 =
      new TH2F("h2_target_c2_11", ";#eta width;p_{T};c2_11", n_bins_deltaEta,
               bins_deltaEta, n_bins_pt, bins_pt);
  TH2F *h2_target_v2 =
      new TH2F("h2_target_v2", ";#eta width;p_{T}", n_bins_deltaEta,
               bins_deltaEta, n_bins_pt, bins_pt);

  TH2F *h2_target_c2_11_bg =
      new TH2F("h2_target_c2_11_bg", ";#eta width;p_{T};c2_11", n_bins_deltaEta,
               bins_deltaEta, n_bins_pt, bins_pt);
  TH2F *h2_target_v2_bg =
      new TH2F("h2_target_v2_bg", ";#eta width;p_{T}", n_bins_deltaEta,
               bins_deltaEta, n_bins_pt, bins_pt);

  TH1F *h1_v2_jpsi =
      new TH1F("h1_v2_jpsi", ";p_{T};v2 J/#psi", n_bins_pt, bins_pt);
  TH1F *h1_v2_ls =
      new TH1F("h1_v2_ls", ";p_{T};v2 like sign", n_bins_pt, bins_pt);
  TH1F *h1_v2_us =
      new TH1F("h1_v2_us", ";p_{T};v2 unlike sign", n_bins_pt, bins_pt);
  TH1F *h1_c211_us =
      new TH1F("h1_c211_us", ";p_{T};c{2,11} unlike sign", n_bins_pt, bins_pt);
  TH1F *h1_c211_ls =
      new TH1F("h1_c211_ls", ";p_{T};c{2,11} like sign", n_bins_pt, bins_pt);

  TH1F *h1_reference_c11 = new TH1F("h1_reference_c11", ";#Delta#eta;c{11}",
                                    n_bins_deltaEta, bins_deltaEta);
  TH1F *h1_reference_v1 = new TH1F("h1_reference_v1", ";#Delta#eta;v1",
                                   n_bins_deltaEta, bins_deltaEta);

  // c2_11 c11->c2_11/c11=c2->v2;c11->v1^2

  for (int i_bin_etaWidth = 1; i_bin_etaWidth <= n_bins_deltaEta;
       i_bin_etaWidth++) {
    for (int i_bin_pt = 1; i_bin_pt <= n_bins_pt; i_bin_pt++) {
      int i_bin_prof2d =
          prof2d_target_c2_11_im->GetBin(i_bin_etaWidth, i_bin_pt);
      double im_target_c2_11 =
          prof2d_target_c2_11_im->GetBinContent(i_bin_prof2d);
      double re_target_c2_11 =
          prof2d_target_c2_11_re->GetBinContent(i_bin_prof2d);
      double im_target_c2_err =
          prof2d_target_c2_11_im->GetBinError(i_bin_prof2d);
      double re_target_c2_err =
          prof2d_target_c2_11_re->GetBinError(i_bin_prof2d);
      double *result_target_c2_11 = ProessComplex(
          im_target_c2_11, re_target_c2_11, im_target_c2_err, re_target_c2_err);
      h2_target_c2_11->SetBinContent(i_bin_etaWidth, i_bin_pt,
                                     result_target_c2_11[0]);
      h2_target_c2_11->SetBinError(i_bin_etaWidth, i_bin_pt,
                                   result_target_c2_11[1]);
    }
  }

  // for bg
  for (int i_bin_etaWidth = 1; i_bin_etaWidth <= n_bins_deltaEta;
       i_bin_etaWidth++) {
    for (int i_bin_pt = 1; i_bin_pt <= n_bins_pt; i_bin_pt++) {
      int i_bin_prof2d =
          prof2d_target_bg_c2_11_im->GetBin(i_bin_etaWidth, i_bin_pt);
      double im_target_c2_11 =
          prof2d_target_bg_c2_11_im->GetBinContent(i_bin_prof2d);
      double re_target_c2_11 =
          prof2d_target_bg_c2_11_re->GetBinContent(i_bin_prof2d);
      double im_target_c2_err =
          prof2d_target_bg_c2_11_im->GetBinError(i_bin_prof2d);
      double re_target_c2_err =
          prof2d_target_bg_c2_11_re->GetBinError(i_bin_prof2d);
      double *result_target_c2_11 = ProessComplex(
          im_target_c2_11, re_target_c2_11, im_target_c2_err, re_target_c2_err);
      h2_target_c2_11_bg->SetBinContent(i_bin_etaWidth, i_bin_pt,
                                        result_target_c2_11[0]);
      h2_target_c2_11_bg->SetBinError(i_bin_etaWidth, i_bin_pt,
                                      result_target_c2_11[1]);
    }
  }

  for (int i_bin_deltaEta = 1; i_bin_deltaEta <= n_bins_deltaEta;
       i_bin_deltaEta++) {
    int i_bin_prof2d = prof2d_reference_c11_im->GetBin(i_bin_deltaEta, 1);
    double im_reference_c11 =
        prof2d_reference_c11_im->GetBinContent(i_bin_prof2d);
    double re_reference_c11 =
        prof2d_reference_c11_re->GetBinContent(i_bin_prof2d);
    double im_reference_c11_err =
        prof2d_reference_c11_im->GetBinError(i_bin_prof2d);
    double re_reference_c11_err =
        prof2d_reference_c11_re->GetBinError(i_bin_prof2d);
    double *result_reference_c11 =
        ProessComplex(im_reference_c11, re_reference_c11, im_reference_c11_err,
                      re_reference_c11_err);
    h1_reference_c11->SetBinContent(i_bin_deltaEta, result_reference_c11[0]);
    h1_reference_c11->SetBinError(i_bin_deltaEta, result_reference_c11[1]);
    double *result_reference_v1 = ProcessPowerFunction(
        result_reference_c11[0], 0.5, result_reference_c11[1]);
    h1_reference_v1->SetBinContent(i_bin_deltaEta, result_reference_v1[0]);
    h1_reference_v1->SetBinError(i_bin_deltaEta, result_reference_v1[1]);
  }

  for (int i_bin_etaWidth = 1; i_bin_etaWidth <= n_bins_deltaEta;
       i_bin_etaWidth++) {
    for (int i_bin_pt = 1; i_bin_pt <= n_bins_pt; i_bin_pt++) {
      double c2_11 = h2_target_c2_11->GetBinContent(i_bin_etaWidth, i_bin_pt);
      double c11 = h1_reference_c11->GetBinContent(i_bin_etaWidth);
      double c2_11_err = h2_target_c2_11->GetBinError(i_bin_etaWidth, i_bin_pt);
      double c11_err = h1_reference_c11->GetBinError(i_bin_etaWidth);
      double *result_v2 = ProcessDivision(c2_11, c11, c2_11_err, c11_err);
      h2_target_v2->SetBinContent(i_bin_etaWidth, i_bin_pt, result_v2[0]);
      h2_target_v2->SetBinError(i_bin_etaWidth, i_bin_pt, result_v2[1]);

      double c2_11_bg =
          h2_target_c2_11_bg->GetBinContent(i_bin_etaWidth, i_bin_pt);
      double c2_11_bg_err =
          h2_target_c2_11_bg->GetBinError(i_bin_etaWidth, i_bin_pt);
      double *result_v2_bg =
          ProcessDivision(c2_11_bg, c11, c2_11_bg_err, c11_err);
      h2_target_v2_bg->SetBinContent(i_bin_etaWidth, i_bin_pt, result_v2_bg[0]);
      h2_target_v2_bg->SetBinError(i_bin_etaWidth, i_bin_pt, result_v2_bg[1]);
    }
  }

  int i_bin_eta_required = 3;
  for (int i_bin_pt = 1; i_bin_pt <= n_bins_pt; i_bin_pt++) {
    double v2 = h2_target_v2->GetBinContent(i_bin_eta_required, i_bin_pt);
    double v2_err = h2_target_v2->GetBinError(i_bin_eta_required, i_bin_pt);
    h1_v2_us->SetBinContent(i_bin_pt, v2);
    h1_v2_us->SetBinError(i_bin_pt, v2_err);
    double v2_bg = h2_target_v2_bg->GetBinContent(i_bin_eta_required, i_bin_pt);
    double v2_bg_err =
        h2_target_v2_bg->GetBinError(i_bin_eta_required, i_bin_pt);
    h1_v2_ls->SetBinContent(i_bin_pt, v2_bg);
    h1_v2_ls->SetBinError(i_bin_pt, v2_bg_err);
    double signal2Total = h1_signal2Total_Jpsi->GetBinContent(i_bin_pt);
    double signal2Total_err = h1_signal2Total_Jpsi->GetBinError(i_bin_pt);
    double v2_jpsi = signal2Total * v2 + (1 - signal2Total) * v2_bg;
    double v2_jpsi_err =
        sqrt(pow(signal2Total * v2_err, 2) + pow(v2 * signal2Total_err, 2) +
             pow((1 - signal2Total) * v2_bg_err, 2) +
             pow(v2_bg * signal2Total_err, 2));
    cout << v2_jpsi << endl;
    cout << v2_jpsi_err << endl;
    h1_v2_jpsi->SetBinContent(i_bin_pt, v2_jpsi);
    h1_v2_jpsi->SetBinError(i_bin_pt, v2_jpsi_err);

    double c2_11 = h2_target_c2_11->GetBinContent(i_bin_eta_required, i_bin_pt);
    double c2_11_err = h2_target_c2_11->GetBinError(i_bin_eta_required, i_bin_pt);
    double c2_11_bg =
        h2_target_c2_11_bg->GetBinContent(i_bin_eta_required, i_bin_pt);
    double c2_11_bg_err =
        h2_target_c2_11_bg->GetBinError(i_bin_eta_required, i_bin_pt);
    h1_c211_us->SetBinContent(i_bin_pt, c2_11);
    h1_c211_us->SetBinError(i_bin_pt, c2_11_err);
    h1_c211_ls->SetBinContent(i_bin_pt, c2_11_bg);
    h1_c211_ls->SetBinError(i_bin_pt, c2_11_bg_err);
  }

  TCanvas *c_result = new TCanvas("c_result", "", 800, 800);
  gStyle->SetOptStat(0);
  c_result->Divide(2, 2);
  c_result->cd(1);
  h1_v2_jpsi->Draw();
  // h1_c211_ls->Draw();

  //   gPad->SetLogy();
  //   h2_target_c2_11->Draw("colz");
  //   vec_h1_jpsiReco[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  //   vec_h1_jpsiReco[0]->GetYaxis()->SetTitle("Signal");
  //   vec_h1_jpsiReco[0]->Draw("colz");

  c_result->cd(2);
  // h1_v2_ls->Draw();
  h1_c211_us->Draw("E");

  //   gPad->SetLogy();
  //   h2_target_c2_11_bg->Draw("colz");
  //   vec_h1_jpsiReco[1]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  //   vec_h1_jpsiReco[1]->GetYaxis()->SetTitle("Total");
  //   vec_h1_jpsiReco[1]->Draw("colz");
  c_result->cd(3);
  h1_v2_us->Draw();
  //   h1_signal2Total_Jpsi->GetYaxis()->SetRangeUser(0, 0.2);
  //   h1_signal2Total_Jpsi->Draw();
  //   h1_reference_c11->Draw();
  c_result->cd(4);
  // h1_v2_jpsi->Draw();
  h1_reference_v1->Draw();
}