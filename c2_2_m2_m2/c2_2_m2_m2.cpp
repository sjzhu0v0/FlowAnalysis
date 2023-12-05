#include "MHead.h"

#define verbose 2
#define roothttp

void c2_2_m2_m2(
    string path_file_input = "/data/work/work_alice/Analysis/JpsiFlow/input/"
                             "AnalysisResults_new2.root",
    string path_file_jpsiReco = "/data/work/work_alice/Analysis/JpsiFlow/c2_m2/"
                                "output/JpsiReco.root",
    string path_file_output =
        "/data/work/work_alice/Analysis/JpsiFlow/c2_11/output/JpsiFlow.root") {
#ifdef roothttp
  THttpServer *server = new THttpServer("http:8090?top=job");
  server->SetReadOnly(kFALSE);
  gBenchmark->Start("job");
#endif

  SetPolttingStyle();

  TFile *file_input = new TFile(path_file_input.c_str(), "READ");
  TFile *file_jpsiReco = new TFile(path_file_jpsiReco.c_str(), "READ");
  TFile *file_output = new TFile(path_file_output.c_str(), "RECREATE");

  vector<string> vec_str_name_obj;
  vector<TObject *> vec_obj = GetObjectRecursive(file_input, vec_str_name_obj);

  vector<TH1 *> vec_h1_reference =
      GetObjectVector<TH1>(vec_obj, vec_str_name_obj, vector<string>{},
                           (TString) "ProfileFlowReference");
  vector<TH1 *> vec_h1_target =
      GetObjectVector<TH1>(vec_obj, vec_str_name_obj, vector<string>{},
                           (TString) "fProfileTarget_EtaWidth_Pt_Mass");
  vector<TH1 *> vec_h1_targetbg =
      GetObjectVector<TH1>(vec_obj, vec_str_name_obj, vector<string>{},
                           (TString) "fProfileTargetBg_EtaWidth_Pt_Mass");
  vector<TH1 *> vec_h1_mass = GetObjectVector<TH1>(
      vec_obj, vec_str_name_obj, vector<string>{}, (TString) "Eta_Pt_Mass");

  vector<string> vec_str_name_obj_jpsiReco;
  vector<TObject *> vec_obj_jpsiReco =
      GetObjectRecursive(file_jpsiReco, vec_str_name_obj_jpsiReco);
  vector<TObject *> vec_h1_jpsiReco =
      GetObjectVector<TObject>(vec_obj_jpsiReco, vec_str_name_obj_jpsiReco,
                               vector<string>{}, (TString) "");
  vector<TH1 *> vec_h1_mass_fit =
      GetObjectVector<TH1>(vec_obj_jpsiReco, vec_str_name_obj_jpsiReco,
                           vector<string>{}, (TString) "Rawhisto");
  vector<TH1 *> vec_h1_mass_uls =
      GetObjectVector<TH1>(vec_obj_jpsiReco, vec_str_name_obj_jpsiReco,
                           vector<string>{}, (TString) "h1_uls");
  vector<TF1 *> vec_f1_resBkg =
      GetObjectVector<TF1>(vec_obj_jpsiReco, vec_str_name_obj_jpsiReco,
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
  for (int i = 0; i < vec_h1_targetbg.size(); i++) {
    cout << "index:\t" << i << "\t\t";
    vec_h1_targetbg[i]->Print();
  }

  cout << endl;
  for (auto &h1 : vec_h1_jpsiReco)
    h1->Print();

  cout << endl;
  for (auto &h1 : vec_h1_mass)
    h1->Print();

  // reference flow
  TProfile2D *fProfileFlowReferenceSubEvent_DeltaEta_Mult_c2_m2_im =
      (TProfile2D *)vec_h1_reference[6];
  TProfile2D *fProfileFlowReferenceSubEvent_DeltaEta_Mult_c2_m2_re =
      (TProfile2D *)vec_h1_reference[7];
  TProfile2D *fProfileFlowReferenceSubEvent_DeltaEta_Mult_c22_m2m2_im =
      (TProfile2D *)vec_h1_reference[8];
  TProfile2D *fProfileFlowReferenceSubEvent_DeltaEta_Mult_c22_m2m2_re =
      (TProfile2D *)vec_h1_reference[9];

  double *bins_deltaEta =
      fProfileFlowReferenceSubEvent_DeltaEta_Mult_c22_m2m2_im->GetXaxis()
          ->GetXbins()
          ->fArray;
  int n_bins_deltaEta =
      fProfileFlowReferenceSubEvent_DeltaEta_Mult_c22_m2m2_im->GetXaxis()
          ->GetNbins();

  TH1F *h1_v2_2_reference =
      new TH1F("h1_v2_2_reference", ";#Delta #eta;v_{2}{2}_{2 sub events}",
               n_bins_deltaEta, bins_deltaEta);
  TH1F *h1_v2_4_reference =
      new TH1F("h1_v2_4_reference", ";#Delta #eta;v_{2}{4}_{2 sub events}",
               n_bins_deltaEta, bins_deltaEta);
  TH1F *h1_c2_4_pow_reference =
      new TH1F("h1_c2_4_pow_reference", ";#Delta #eta;v_{2}{4}_{2 sub events}",
               n_bins_deltaEta, bins_deltaEta);
  for (int i = 1; i <= n_bins_deltaEta; i++) {
    double c2_m2_im =
        fProfileFlowReferenceSubEvent_DeltaEta_Mult_c2_m2_im->GetBinContent(i,
                                                                            1);
    double c2_m2_re =
        fProfileFlowReferenceSubEvent_DeltaEta_Mult_c2_m2_re->GetBinContent(i,
                                                                            1);
    double c22_m2m2_im =
        fProfileFlowReferenceSubEvent_DeltaEta_Mult_c22_m2m2_im->GetBinContent(
            i, 1);
    double c22_m2m2_re =
        fProfileFlowReferenceSubEvent_DeltaEta_Mult_c22_m2m2_re->GetBinContent(
            i, 1);

    double c2_m2_im_err =
        fProfileFlowReferenceSubEvent_DeltaEta_Mult_c2_m2_im->GetBinError(i, 1);
    double c2_m2_re_err =
        fProfileFlowReferenceSubEvent_DeltaEta_Mult_c2_m2_re->GetBinError(i, 1);
    double c22_m2m2_im_err =
        fProfileFlowReferenceSubEvent_DeltaEta_Mult_c22_m2m2_im->GetBinError(i,
                                                                             1);
    double c22_m2m2_re_err =
        fProfileFlowReferenceSubEvent_DeltaEta_Mult_c22_m2m2_re->GetBinError(i,
                                                                             1);

    TComplex c22_m2m2(c22_m2m2_re, c22_m2m2_im);
    TComplex c2_m2(c2_m2_re, c2_m2_im);
    TComplex c22_m2m2_err(c22_m2m2_re_err, c22_m2m2_im_err);
    TComplex c2_m2_err(c2_m2_re_err, c2_m2_im_err);
    TComplex c_4_2 = c22_m2m2 - 2. * c2_m2 * c2_m2;

    TComplex *c2_m2_square =
        ProcessComplexProduct(c2_m2, c2_m2, c2_m2_err, c2_m2_err);
    double err_c_4_2_re =
        sqrt(pow(c22_m2m2_err.Re(), 2) + pow(2. * c2_m2_square[1].Re(), 2));
    double err_c_4_2_im =
        sqrt(pow(c22_m2m2_err.Im(), 2) + pow(2. * c2_m2_square[1].Im(), 2));
    TComplex c_4_2_err(err_c_4_2_re, err_c_4_2_im);
    double *result_c_4_2 =
        ProcessComplex(c_4_2.Im(), c_4_2.Re(), c_4_2_err.Im(), c_4_2_err.Re());
    double *result_c_2_2 =
        ProcessComplex(c2_m2.Im(), c2_m2.Re(), c2_m2_err.Im(), c2_m2_err.Re());

    double *result_v2_4 =
        ProcessPowerFunction(result_c_4_2[0], 0.25, result_c_4_2[1]);
    double *result_v2_2 =
        ProcessPowerFunction(result_c_2_2[0], 0.5, result_c_2_2[1]);
    double *result_c2_4_pow =
        ProcessPowerFunction(result_c_2_2[0], 0.75, result_c_2_2[1]);

    h1_v2_2_reference->SetBinContent(i, result_v2_2[0]);
    h1_v2_2_reference->SetBinError(i, result_v2_2[1]);
    h1_v2_4_reference->SetBinContent(i, result_v2_4[0]);
    h1_v2_4_reference->SetBinError(i, result_v2_4[1]);
    h1_c2_4_pow_reference->SetBinContent(i, result_c2_4_pow[0]);
    h1_c2_4_pow_reference->SetBinError(i, result_c2_4_pow[1]);
  }

  TCanvas *c1_reference = new TCanvas("c1_reference", "", 400, 800);
  c1_reference->Divide(1, 2);
  c1_reference->cd(1);
  h1_v2_2_reference->GetYaxis()->SetRangeUser(0., 0.2);
  h1_v2_2_reference->Draw();
  c1_reference->cd(2);
  h1_v2_4_reference->GetYaxis()->SetRangeUser(0., 0.2);
  h1_v2_4_reference->Draw();

  //////////////////////////////////////////////////////////////////
  TH3F *h3Target_Eta_Pt_Mass = (TH3F *)vec_h1_mass[0];
  TH3F *h3TargetBg_Eta_Pt_Mass = (TH3F *)vec_h1_mass[1];

  double *bins_eta_invmass =
      h3Target_Eta_Pt_Mass->GetXaxis()->GetXbins()->fArray;
  double *bins_pt_invmass =
      h3Target_Eta_Pt_Mass->GetYaxis()->GetXbins()->fArray;
  double *bins_mass_invmass =
      h3Target_Eta_Pt_Mass->GetZaxis()->GetXbins()->fArray;
  int n_bins_eta_invmass = h3Target_Eta_Pt_Mass->GetXaxis()->GetNbins();
  int n_bins_pt_invmass = h3Target_Eta_Pt_Mass->GetYaxis()->GetNbins();
  int n_bins_mass_invmass = h3Target_Eta_Pt_Mass->GetZaxis()->GetNbins();

  // target flow
  // tprofile3d for target: c2_2_m2_m2, c2_0_m2_m0, c2_0_m0_m2, c0_2_m2_m0,
  // c0_2_m0_m2
  TProfile3D *fProfileTarget_EtaWidth_Pt_Mass_c2_2_m2_m2_im =
      (TProfile3D *)vec_h1_target[2];
  TProfile3D *fProfileTarget_EtaWidth_Pt_Mass_c2_2_m2_m2_re =
      (TProfile3D *)vec_h1_target[3];
  TProfile3D *fProfileTarget_EtaWidth_Pt_Mass_c2_0_m2_m0_im =
      (TProfile3D *)vec_h1_target[4];
  TProfile3D *fProfileTarget_EtaWidth_Pt_Mass_c2_0_m2_m0_re =
      (TProfile3D *)vec_h1_target[5];
  TProfile3D *fProfileTarget_EtaWidth_Pt_Mass_c2_0_m0_m2_im =
      (TProfile3D *)vec_h1_target[8];
  TProfile3D *fProfileTarget_EtaWidth_Pt_Mass_c2_0_m0_m2_re =
      (TProfile3D *)vec_h1_target[9];
  TProfile3D *fProfileTarget_EtaWidth_Pt_Mass_c0_2_m2_m0_im =
      (TProfile3D *)vec_h1_target[10];
  TProfile3D *fProfileTarget_EtaWidth_Pt_Mass_c0_2_m2_m0_re =
      (TProfile3D *)vec_h1_target[11];
  TProfile3D *fProfileTarget_EtaWidth_Pt_Mass_c0_2_m0_m2_im =
      (TProfile3D *)vec_h1_target[6];
  TProfile3D *fProfileTarget_EtaWidth_Pt_Mass_c0_2_m0_m2_re =
      (TProfile3D *)vec_h1_target[7];
  // target flow bg
  // tprofile3d for target: c2_2_m2_m2, c2_0_m2_m0, c2_0_m0_m2, c0_2_m2_m0,
  // c0_2_m0_m2

  TProfile3D *fProfileTargetbg_EtaWidth_Pt_Mass_c2_2_m2_m2_im =
      (TProfile3D *)vec_h1_targetbg[2];
  TProfile3D *fProfileTargetbg_EtaWidth_Pt_Mass_c2_2_m2_m2_re =
      (TProfile3D *)vec_h1_targetbg[3];
  TProfile3D *fProfileTargetbg_EtaWidth_Pt_Mass_c2_0_m2_m0_im =
      (TProfile3D *)vec_h1_targetbg[4];
  TProfile3D *fProfileTargetbg_EtaWidth_Pt_Mass_c2_0_m2_m0_re =
      (TProfile3D *)vec_h1_targetbg[5];
  TProfile3D *fProfileTargetbg_EtaWidth_Pt_Mass_c2_0_m0_m2_im =
      (TProfile3D *)vec_h1_targetbg[8];
  TProfile3D *fProfileTargetbg_EtaWidth_Pt_Mass_c2_0_m0_m2_re =
      (TProfile3D *)vec_h1_targetbg[9];
  TProfile3D *fProfileTargetbg_EtaWidth_Pt_Mass_c0_2_m2_m0_im =
      (TProfile3D *)vec_h1_targetbg[10];
  TProfile3D *fProfileTargetbg_EtaWidth_Pt_Mass_c0_2_m2_m0_re =
      (TProfile3D *)vec_h1_targetbg[11];
  TProfile3D *fProfileTargetbg_EtaWidth_Pt_Mass_c0_2_m0_m2_im =
      (TProfile3D *)vec_h1_targetbg[6];
  TProfile3D *fProfileTargetbg_EtaWidth_Pt_Mass_c0_2_m0_m2_re =
      (TProfile3D *)vec_h1_targetbg[7];

  double *bins_widthEta =
      fProfileTarget_EtaWidth_Pt_Mass_c2_2_m2_m2_im->GetXaxis()
          ->GetXbins()
          ->fArray;
  double *bins_pt = fProfileTarget_EtaWidth_Pt_Mass_c2_2_m2_m2_im->GetYaxis()
                        ->GetXbins()
                        ->fArray;
  double *bins_mass = fProfileTarget_EtaWidth_Pt_Mass_c2_2_m2_m2_im->GetZaxis()
                          ->GetXbins()
                          ->fArray;
  int n_bins_eta_width =
      fProfileTarget_EtaWidth_Pt_Mass_c2_2_m2_m2_im->GetXaxis()->GetNbins();
  int n_bins_pt =
      fProfileTarget_EtaWidth_Pt_Mass_c2_2_m2_m2_im->GetYaxis()->GetNbins();
  int n_bins_mass =
      fProfileTarget_EtaWidth_Pt_Mass_c2_2_m2_m2_im->GetZaxis()->GetNbins();

  TH3D *h3_widthEta_pt_mass_target =
      new TH3D("h3_widthEta_pt_mass_target",
               ";#eta width;#it{p}_{T};#it{m}_{e^{+}e^{-}}", n_bins_eta_width,
               bins_widthEta, n_bins_pt, bins_pt, n_bins_mass, bins_mass);
  TH3D *h3_widthEta_pt_mass_targetbg =
      new TH3D("h3_widthEta_pt_mass_targetbg",
               ";#eta width;#it{p}_{T};#it{m}_{e^{+}e^{-}}", n_bins_eta_width,
               bins_widthEta, n_bins_pt, bins_pt, n_bins_mass, bins_mass);

  int index_bin_mass_low =
      h3Target_Eta_Pt_Mass->GetZaxis()->FindBin(bins_mass[0]);
  int index_bin_mass_high =
      h3Target_Eta_Pt_Mass->GetZaxis()->FindBin(bins_mass[n_bins_mass - 1]);

  for (int i_bin_widthEta = 1; i_bin_widthEta <= n_bins_eta_width;
       i_bin_widthEta++) {
    for (int i_bin_pt = 1; i_bin_pt <= n_bins_pt; i_bin_pt++) {
      double widthEta = bins_widthEta[i_bin_widthEta - 1];
      double pt = bins_pt[i_bin_pt - 1];
      double edge_eta_low = -widthEta / 2.;
      double edge_eta_high = widthEta / 2.;
      int index_bin_eta_invMass_low =
          h3Target_Eta_Pt_Mass->GetXaxis()->FindBin(edge_eta_low);
      int index_bin_eta_invMass_high =
          h3Target_Eta_Pt_Mass->GetXaxis()->FindBin(edge_eta_high);
      for (int i_bin_mass = index_bin_mass_low;
           i_bin_mass <= index_bin_mass_high; i_bin_mass++) {
        double count_uls = 0.;
        double count_uls_err = 0.;
        double count_ls = 0.;
        double count_ls_err = 0.;
        for (int i_bin_eta = index_bin_eta_invMass_low;
             i_bin_eta < index_bin_eta_invMass_high; i_bin_eta++) {
          count_uls += h3Target_Eta_Pt_Mass->GetBinContent(i_bin_eta, i_bin_pt,
                                                           i_bin_mass);
          count_uls_err += pow(h3Target_Eta_Pt_Mass->GetBinError(
                                   i_bin_eta, i_bin_pt, i_bin_mass),
                               2);
          count_ls += h3TargetBg_Eta_Pt_Mass->GetBinContent(i_bin_eta, i_bin_pt,
                                                            i_bin_mass);
          count_ls_err += pow(h3TargetBg_Eta_Pt_Mass->GetBinError(
                                  i_bin_eta, i_bin_pt, i_bin_mass),
                              2);
        }
        count_uls_err = sqrt(count_uls_err);
        count_ls_err = sqrt(count_ls_err);
        int index_bin_mass_target =
            h3_widthEta_pt_mass_target->GetZaxis()->FindBin(
                h3Target_Eta_Pt_Mass->GetZaxis()->GetBinCenter(i_bin_mass));
        h3_widthEta_pt_mass_target->SetBinContent(
            i_bin_widthEta, i_bin_pt, index_bin_mass_target, count_uls);
        h3_widthEta_pt_mass_target->SetBinError(
            i_bin_widthEta, i_bin_pt, index_bin_mass_target, count_uls_err);
        h3_widthEta_pt_mass_targetbg->SetBinContent(
            i_bin_widthEta, i_bin_pt, index_bin_mass_target, count_ls);
        h3_widthEta_pt_mass_targetbg->SetBinError(
            i_bin_widthEta, i_bin_pt, index_bin_mass_target, count_ls_err);
      }
    }
  }

  TH3D *h3target_widthEta_pt_mass_C2_4 =
      new TH3D("h3target_widthEta_pt_mass_C2_4",
               "h3target_widthEta_pt_mass_C2_4;#eta width;p_{T} "
               "[GeV/c];Mass_{e^{+}e^{-}}[GeV/c]",
               n_bins_eta_width, bins_widthEta, n_bins_pt, bins_pt, n_bins_mass,
               bins_mass); // c2_2_m2_m2
  TH3D *h3targetbg_widthEta_pt_mass_C2_4 =
      new TH3D("h3targetbg_widthEta_pt_mass_C2_4",
               "h3targetbg_widthEta_pt_mass_C2_4;#eta width;p_{T} "
               "[GeV/c];Mass_{e^{+}e^{-}}[GeV/c]",
               n_bins_eta_width, bins_widthEta, n_bins_pt, bins_pt, n_bins_mass,
               bins_mass);
  TH3D *h3target_widthEta_pt_mass_v2_4 =
      new TH3D("h3target_widthEta_pt_mass_v2_4",
               "h3target_widthEta_pt_mass_v2_4;#eta width;p_{T} "
               "[GeV/c];Mass_{e^{+}e^{-}}[GeV/c]",
               n_bins_eta_width, bins_widthEta, n_bins_pt, bins_pt, n_bins_mass,
               bins_mass); // c2_2_m2_m2
  TH3D *h3targetbg_widthEta_pt_mass_v2_4 =
      new TH3D("h3targetbg_widthEta_pt_mass_v2_4",
               "h3targetbg_widthEta_pt_mass_v2_4;#eta width;p_{T} "
               "[GeV/c];Mass_{e^{+}e^{-}}[GeV/c]",
               n_bins_eta_width, bins_widthEta, n_bins_pt, bins_pt, n_bins_mass,
               bins_mass);
  TH3D *h3pre_widthEta_pt_mass_C2_4 =
      new TH3D("h3pre_widthEta_pt_mass_C2_4",
               "h3pre_widthEta_pt_mass_C2_4;#eta width;p_{T} "
               "[GeV/c];Mass_{e^{+}e^{-}}[GeV/c]",
               n_bins_eta_width, bins_widthEta, n_bins_pt, bins_pt, n_bins_mass,
               bins_mass);
  TH3D *h3pre_widthEta_pt_mass_NC2_4 =
      new TH3D("h3pre_widthEta_pt_mass_NC2_4",
               "h3pre_widthEta_pt_mass_NC2_4;#eta width;p_{T} "
               "[GeV/c];Mass_{e^{+}e^{-}}[GeV/c]",
               n_bins_eta_width, bins_widthEta, n_bins_pt, bins_pt, n_bins_mass,
               bins_mass);
  TH3D *h3pre_widthEta_pt_mass_v2_4 =
      new TH3D("h3pre_widthEta_pt_mass_v2_4", "", n_bins_eta_width,
               bins_widthEta, n_bins_pt, bins_pt, n_bins_mass, bins_mass);
  TH3D *h3pre_widthEta_pt_mass_Nv2_4 =
      new TH3D("h3pre_widthEta_pt_mass_Nv2_4", "", n_bins_eta_width,
               bins_widthEta, n_bins_pt, bins_pt, n_bins_mass, bins_mass);

  TH3D *h3_widthEta_pt_mass_Nv2_4_NmNbg =
      new TH3D("h3_widthEta_pt_mass_Nv2_4_NmNbg", "", n_bins_eta_width,
               bins_widthEta, n_bins_pt, bins_pt, n_bins_mass, bins_mass);

  for (int i_bin_widthEta = 1; i_bin_widthEta <= n_bins_eta_width;
       i_bin_widthEta++) {
    for (int i_bin_pt = 1; i_bin_pt <= n_bins_pt; i_bin_pt++) {
      for (int i_bin_mass = 1; i_bin_mass <= n_bins_mass; i_bin_mass++) {
        double etaWidth = bins_widthEta[i_bin_widthEta - 1];
        int index_deltaEta_reference =
            h1_v2_4_reference->GetXaxis()->FindBin(etaWidth);
        double c2_4_pow_reference =
            h1_c2_4_pow_reference->GetBinContent(index_deltaEta_reference);
        double c2_4_pow_reference_err =
            h1_c2_4_pow_reference->GetBinError(index_deltaEta_reference);
        // cout << v2_reference << endl;
        MComplex mc_c2_4_pow_reference(c2_4_pow_reference, 0.,
                                       c2_4_pow_reference_err, 0.);
        double c2_2_m2_m2_im_target =
            fProfileTarget_EtaWidth_Pt_Mass_c2_2_m2_m2_im->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_2_m2_m2_re_target =
            fProfileTarget_EtaWidth_Pt_Mass_c2_2_m2_m2_re->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m2_m0_im_target =
            fProfileTarget_EtaWidth_Pt_Mass_c2_0_m2_m0_im->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m2_m0_re_target =
            fProfileTarget_EtaWidth_Pt_Mass_c2_0_m2_m0_re->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m0_m2_im_target =
            fProfileTarget_EtaWidth_Pt_Mass_c2_0_m0_m2_im->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m0_m2_re_target =
            fProfileTarget_EtaWidth_Pt_Mass_c2_0_m0_m2_re->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m2_m0_im_target =
            fProfileTarget_EtaWidth_Pt_Mass_c0_2_m2_m0_im->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m2_m0_re_target =
            fProfileTarget_EtaWidth_Pt_Mass_c0_2_m2_m0_re->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m0_m2_im_target =
            fProfileTarget_EtaWidth_Pt_Mass_c0_2_m0_m2_im->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m0_m2_re_target =
            fProfileTarget_EtaWidth_Pt_Mass_c0_2_m0_m2_re->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_2_m2_m2_im_target_err =
            fProfileTarget_EtaWidth_Pt_Mass_c2_2_m2_m2_im->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_2_m2_m2_re_target_err =
            fProfileTarget_EtaWidth_Pt_Mass_c2_2_m2_m2_re->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m2_m0_im_target_err =
            fProfileTarget_EtaWidth_Pt_Mass_c2_0_m2_m0_im->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m2_m0_re_target_err =
            fProfileTarget_EtaWidth_Pt_Mass_c2_0_m2_m0_re->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m0_m2_im_target_err =
            fProfileTarget_EtaWidth_Pt_Mass_c2_0_m0_m2_im->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m0_m2_re_target_err =
            fProfileTarget_EtaWidth_Pt_Mass_c2_0_m0_m2_re->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m2_m0_im_target_err =
            fProfileTarget_EtaWidth_Pt_Mass_c0_2_m2_m0_im->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m2_m0_re_target_err =
            fProfileTarget_EtaWidth_Pt_Mass_c0_2_m2_m0_re->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m0_m2_im_target_err =
            fProfileTarget_EtaWidth_Pt_Mass_c0_2_m0_m2_im->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m0_m2_re_target_err =
            fProfileTarget_EtaWidth_Pt_Mass_c0_2_m0_m2_re->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);

        double N_target = h3_widthEta_pt_mass_target->GetBinContent(
            i_bin_widthEta, i_bin_pt, i_bin_mass);
        double N_target_err = h3_widthEta_pt_mass_target->GetBinError(
            i_bin_widthEta, i_bin_pt, i_bin_mass);
        MComplex c2_2_m2_m2_target(c2_2_m2_m2_re_target, c2_2_m2_m2_im_target,
                                   c2_2_m2_m2_re_target_err,
                                   c2_2_m2_m2_im_target_err);
        MComplex c2_0_m2_m0_target(c2_0_m2_m0_re_target, c2_0_m2_m0_im_target,
                                   c2_0_m2_m0_re_target_err,
                                   c2_0_m2_m0_im_target_err);
        MComplex c2_0_m0_m2_target(c2_0_m0_m2_re_target, c2_0_m0_m2_im_target,
                                   c2_0_m0_m2_re_target_err,
                                   c2_0_m0_m2_im_target_err);
        MComplex c0_2_m2_m0_target(c0_2_m2_m0_re_target, c0_2_m2_m0_im_target,
                                   c0_2_m2_m0_re_target_err,
                                   c0_2_m2_m0_im_target_err);
        MComplex c0_2_m0_m2_target(c0_2_m0_m2_re_target, c0_2_m0_m2_im_target,
                                   c0_2_m0_m2_re_target_err,
                                   c0_2_m0_m2_im_target_err);
        MComplex C2_4_target = c2_2_m2_m2_target -
                               (c2_0_m2_m0_target * c0_2_m0_m2_target) -
                               (c2_0_m0_m2_target * c0_2_m2_m0_target);
        h3target_widthEta_pt_mass_C2_4->SetBinContent(
            i_bin_widthEta, i_bin_pt, i_bin_mass, C2_4_target.Mag());
        h3target_widthEta_pt_mass_C2_4->SetBinError(
            i_bin_widthEta, i_bin_pt, i_bin_mass, C2_4_target.MagErr());

        double v2_4_target = C2_4_target.Mag() / mc_c2_4_pow_reference.Re();
        double v2_4_target_err =
            sqrt(pow(C2_4_target.MagErr() / mc_c2_4_pow_reference.Re(), 2) +
                 pow(C2_4_target.Mag() * mc_c2_4_pow_reference.ReErr() /
                         pow(mc_c2_4_pow_reference.Re(), 2),
                     2));
        h3target_widthEta_pt_mass_v2_4->SetBinContent(i_bin_widthEta, i_bin_pt,
                                                      i_bin_mass, v2_4_target);
        h3target_widthEta_pt_mass_v2_4->SetBinError(
            i_bin_widthEta, i_bin_pt, i_bin_mass, v2_4_target_err);

        //////////////// target bg ////////////////////////////////////
        double c2_2_m2_m2_im_targetbg =
            fProfileTargetbg_EtaWidth_Pt_Mass_c2_2_m2_m2_im->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_2_m2_m2_re_targetbg =
            fProfileTargetbg_EtaWidth_Pt_Mass_c2_2_m2_m2_re->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m2_m0_im_targetbg =
            fProfileTargetbg_EtaWidth_Pt_Mass_c2_0_m2_m0_im->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m2_m0_re_targetbg =
            fProfileTargetbg_EtaWidth_Pt_Mass_c2_0_m2_m0_re->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m0_m2_im_targetbg =
            fProfileTargetbg_EtaWidth_Pt_Mass_c2_0_m0_m2_im->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m0_m2_re_targetbg =
            fProfileTargetbg_EtaWidth_Pt_Mass_c2_0_m0_m2_re->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m2_m0_im_targetbg =
            fProfileTargetbg_EtaWidth_Pt_Mass_c0_2_m2_m0_im->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m2_m0_re_targetbg =
            fProfileTargetbg_EtaWidth_Pt_Mass_c0_2_m2_m0_re->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m0_m2_im_targetbg =
            fProfileTargetbg_EtaWidth_Pt_Mass_c0_2_m0_m2_im->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m0_m2_re_targetbg =
            fProfileTargetbg_EtaWidth_Pt_Mass_c0_2_m0_m2_re->GetBinContent(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_2_m2_m2_im_targetbg_err =
            fProfileTargetbg_EtaWidth_Pt_Mass_c2_2_m2_m2_im->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_2_m2_m2_re_targetbg_err =
            fProfileTargetbg_EtaWidth_Pt_Mass_c2_2_m2_m2_re->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m2_m0_im_targetbg_err =
            fProfileTargetbg_EtaWidth_Pt_Mass_c2_0_m2_m0_im->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m2_m0_re_targetbg_err =
            fProfileTargetbg_EtaWidth_Pt_Mass_c2_0_m2_m0_re->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m0_m2_im_targetbg_err =
            fProfileTargetbg_EtaWidth_Pt_Mass_c2_0_m0_m2_im->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c2_0_m0_m2_re_targetbg_err =
            fProfileTargetbg_EtaWidth_Pt_Mass_c2_0_m0_m2_re->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m2_m0_im_targetbg_err =
            fProfileTargetbg_EtaWidth_Pt_Mass_c0_2_m2_m0_im->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m2_m0_re_targetbg_err =
            fProfileTargetbg_EtaWidth_Pt_Mass_c0_2_m2_m0_re->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m0_m2_im_targetbg_err =
            fProfileTargetbg_EtaWidth_Pt_Mass_c0_2_m0_m2_im->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double c0_2_m0_m2_re_targetbg_err =
            fProfileTargetbg_EtaWidth_Pt_Mass_c0_2_m0_m2_re->GetBinError(
                i_bin_widthEta, i_bin_pt, i_bin_mass);
        double N_targetbg = h3_widthEta_pt_mass_targetbg->GetBinContent(
            i_bin_widthEta, i_bin_pt, i_bin_mass);
        double N_targetbg_err = h3_widthEta_pt_mass_targetbg->GetBinError(
            i_bin_widthEta, i_bin_pt, i_bin_mass);
        MComplex c2_2_m2_m2_targetbg(
            c2_2_m2_m2_re_targetbg, c2_2_m2_m2_im_targetbg,
            c2_2_m2_m2_re_targetbg_err, c2_2_m2_m2_im_targetbg_err);
        MComplex c2_0_m2_m0_targetbg(
            c2_0_m2_m0_re_targetbg, c2_0_m2_m0_im_targetbg,
            c2_0_m2_m0_re_targetbg_err, c2_0_m2_m0_im_targetbg_err);
        MComplex c2_0_m0_m2_targetbg(
            c2_0_m0_m2_re_targetbg, c2_0_m0_m2_im_targetbg,
            c2_0_m0_m2_re_targetbg_err, c2_0_m0_m2_im_targetbg_err);
        MComplex c0_2_m2_m0_targetbg(
            c0_2_m2_m0_re_targetbg, c0_2_m2_m0_im_targetbg,
            c0_2_m2_m0_re_targetbg_err, c0_2_m2_m0_im_targetbg_err);
        MComplex c0_2_m0_m2_targetbg(
            c0_2_m0_m2_re_targetbg, c0_2_m0_m2_im_targetbg,
            c0_2_m0_m2_re_targetbg_err, c0_2_m0_m2_im_targetbg_err);
        MComplex C2_4_targetbg = c2_2_m2_m2_targetbg -
                                 (c2_0_m2_m0_targetbg * c0_2_m0_m2_targetbg) -
                                 (c2_0_m0_m2_targetbg * c0_2_m2_m0_targetbg);
        h3targetbg_widthEta_pt_mass_C2_4->SetBinContent(
            i_bin_widthEta, i_bin_pt, i_bin_mass, C2_4_targetbg.Mag());
        h3targetbg_widthEta_pt_mass_C2_4->SetBinError(
            i_bin_widthEta, i_bin_pt, i_bin_mass, C2_4_targetbg.MagErr());

        double v2_4_targetbg = C2_4_targetbg.Mag() / mc_c2_4_pow_reference.Re();
        double v2_4_targetbg_err =
            sqrt(pow(C2_4_targetbg.MagErr() / mc_c2_4_pow_reference.Re(), 2) +
                 pow(C2_4_targetbg.Mag() * mc_c2_4_pow_reference.ReErr() /
                         pow(mc_c2_4_pow_reference.Re(), 2),
                     2));
        h3targetbg_widthEta_pt_mass_v2_4->SetBinContent(
            i_bin_widthEta, i_bin_pt, i_bin_mass, v2_4_targetbg);
        h3targetbg_widthEta_pt_mass_v2_4->SetBinError(
            i_bin_widthEta, i_bin_pt, i_bin_mass, v2_4_targetbg_err);

        ///////////////////////////// pre flow ////////////////////////////////
        MComplex mc_target(N_target, 0., N_target_err, 0.);
        MComplex mc_targetbg(N_targetbg, 0., N_targetbg_err, 0.);

        h3_widthEta_pt_mass_Nv2_4_NmNbg->SetBinContent(
            i_bin_widthEta, i_bin_pt, i_bin_mass,
            (mc_target - mc_targetbg).Re());
        h3_widthEta_pt_mass_Nv2_4_NmNbg->SetBinError(
            i_bin_widthEta, i_bin_pt, i_bin_mass,
            (mc_target - mc_targetbg).ReErr());

        MComplex NC4_2_pre =
            (mc_target * C2_4_target) - (mc_targetbg * C2_4_targetbg);
        h3pre_widthEta_pt_mass_NC2_4->SetBinContent(
            i_bin_widthEta, i_bin_pt, i_bin_mass, NC4_2_pre.Mag());
        h3pre_widthEta_pt_mass_NC2_4->SetBinError(
            i_bin_widthEta, i_bin_pt, i_bin_mass, NC4_2_pre.MagErr());

        MComplex C4_2_pre = NC4_2_pre / (mc_target - mc_targetbg);
        // cout << (mc_target - mc_targetbg).Re() << endl;
        double err1_C4_2_pre =
            N_target / (N_target - N_targetbg) * C2_4_target.ReErr();
        double err2_C4_2_pre =
            N_targetbg / (N_target - N_targetbg) * C2_4_targetbg.ReErr();
        double err3_C4_2_pre = N_target / pow(N_target - N_targetbg, 2) *
                               C2_4_targetbg.Re() * N_targetbg_err;
        double err4_C4_2_pre = N_targetbg / pow(N_target - N_targetbg, 2) *
                               C2_4_target.Re() * N_target_err;
        double err_C4_2_pre =
            sqrt(pow(err1_C4_2_pre, 2) + pow(err2_C4_2_pre, 2) +
                 pow(err3_C4_2_pre, 2) + pow(err4_C4_2_pre, 2));
        C4_2_pre.SetReErr(err_C4_2_pre);
        err1_C4_2_pre =
            N_target / (N_target - N_targetbg) * C2_4_target.ImErr();
        err2_C4_2_pre =
            N_targetbg / (N_target - N_targetbg) * C2_4_targetbg.ImErr();
        err3_C4_2_pre = N_target / pow(N_target - N_targetbg, 2) *
                        C2_4_targetbg.Im() * N_targetbg_err;
        err4_C4_2_pre = N_targetbg / pow(N_target - N_targetbg, 2) *
                        C2_4_target.Im() * N_target_err;
        err_C4_2_pre = sqrt(pow(err1_C4_2_pre, 2) + pow(err2_C4_2_pre, 2) +
                            pow(err3_C4_2_pre, 2) + pow(err4_C4_2_pre, 2));
        C4_2_pre.SetImErr(err_C4_2_pre);
        h3pre_widthEta_pt_mass_C2_4->SetBinContent(i_bin_widthEta, i_bin_pt,
                                                   i_bin_mass, C4_2_pre.Mag());
        h3pre_widthEta_pt_mass_C2_4->SetBinError(i_bin_widthEta, i_bin_pt,
                                                 i_bin_mass, C4_2_pre.Mag());

        MComplex Nv2_4_pre = NC4_2_pre / mc_c2_4_pow_reference;
        MComplex mag_C4_2_pre(C4_2_pre.Mag(), 0, err_C4_2_pre, 0);
        MComplex v2_4_pre = mag_C4_2_pre / mc_c2_4_pow_reference;
        h3pre_widthEta_pt_mass_Nv2_4->SetBinContent(
            i_bin_widthEta, i_bin_pt, i_bin_mass, Nv2_4_pre.Mag());
        h3pre_widthEta_pt_mass_v2_4->SetBinContent(i_bin_widthEta, i_bin_pt,
                                                   i_bin_mass, v2_4_pre.Mag());
        h3pre_widthEta_pt_mass_Nv2_4->SetBinError(
            i_bin_widthEta, i_bin_pt, i_bin_mass, Nv2_4_pre.MagErr());
        h3pre_widthEta_pt_mass_v2_4->SetBinError(i_bin_widthEta, i_bin_pt,
                                                 i_bin_mass, v2_4_pre.MagErr());
      }
    }
  }

  TCanvas *c1_C2_4 = new TCanvas("c1_C2_4", "c1_C2_4", 600, 600);
  c1_C2_4->Divide(2, 2);
  c1_C2_4->cd(1);
  h3target_widthEta_pt_mass_C2_4
      ->ProjectionZ("h3target_widthEta_pt_mass_C2_4_z", 4, 4, 1, 1)
      ->Draw();
  c1_C2_4->cd(2);
  h3targetbg_widthEta_pt_mass_C2_4
      ->ProjectionZ("h3targetbg_widthEta_pt_mass_C2_4_z", 4, 4, 1, 1)
      ->Draw();
  c1_C2_4->cd(3);
  h3pre_widthEta_pt_mass_NC2_4
      ->ProjectionZ("h3pre_widthEta_pt_mass_NC2_4_z", 4, 4, 1, 1)
      ->Draw();

  c1_C2_4->cd(4);
  h3pre_widthEta_pt_mass_C2_4
      ->ProjectionZ("h3pre_widthEta_pt_mass_C2_4_z", 4, 4, 1, 1)
      ->Draw();
  server->Register("c1_C2_4", c1_C2_4);

  double binWidth_requested = 0.8;
  int index_binWidth_requested =
      h3target_widthEta_pt_mass_v2_4->GetXaxis()->FindBin(binWidth_requested);

  TH2F *h2_v2_4_sig_mass_pt = new TH2F(
      "h2_v2_4_target", ";Mass_{e^{+}e^{-}}[GeV/c^{2}];#it{p}_{T} [GeV/c]",
      n_bins_mass, bins_mass, n_bins_pt, bins_pt);

  vector<TH1D *> vec_h1_mass_fit_flow;

  for (int i_bin_pt = 1; i_bin_pt <= n_bins_pt; i_bin_pt++) {
    TH1D *h1_targetbg_mass =
        (TH1D *)h3targetbg_widthEta_pt_mass_v2_4->ProjectionZ(
            Form("h1_targetbg_mass_%d", i_bin_pt), index_binWidth_requested,
            index_binWidth_requested, i_bin_pt, i_bin_pt);
    TF1 f1_Bkg_flow("f1_Bkg_flow", "pol1", bins_mass[0],
                    bins_mass[n_bins_mass]);
    h1_targetbg_mass->Fit(&f1_Bkg_flow, "R");
    vec_h1_mass_fit_flow.push_back(h1_targetbg_mass);
    for (int i_bin_mass = 1; i_bin_mass <= n_bins_mass; i_bin_mass++) {
      double value_binCenter_mass =
          h3target_widthEta_pt_mass_v2_4->GetZaxis()->GetBinCenter(i_bin_mass);
      int index_bin_mass_fit =
          vec_h1_mass_fit[i_bin_pt - 1]->GetXaxis()->FindBin(
              value_binCenter_mass);
      double count_signal =
          vec_h1_mass_fit[i_bin_pt - 1]->GetBinContent(index_bin_mass_fit);
      double count_signal_err =
          vec_h1_mass_fit[i_bin_pt - 1]->GetBinError(index_bin_mass_fit);
      count_signal -= vec_f1_resBkg[i_bin_pt - 1]->Eval(value_binCenter_mass);
      double count_uls =
          vec_h1_mass_uls[i_bin_pt - 1]->GetBinContent(index_bin_mass_fit);
      double count_uls_err =
          vec_h1_mass_uls[i_bin_pt - 1]->GetBinError(index_bin_mass_fit);

      double v2_target = h3target_widthEta_pt_mass_v2_4->GetBinContent(
          index_binWidth_requested, i_bin_pt, i_bin_mass);
      double v2_target_err = h3target_widthEta_pt_mass_v2_4->GetBinError(
          index_binWidth_requested, i_bin_pt, i_bin_mass);
      double v2_targetbg = f1_Bkg_flow.Eval(value_binCenter_mass);
      // double v2_targetbg = h3targetbg_widthEta_pt_mass_v2_4->GetBinContent(
      //     index_binWidth_requested, i_bin_pt, i_bin_mass);
      // double v2_targetbg_err = h3targetbg_widthEta_pt_mass_v2_4->GetBinError(
      //     index_binWidth_requested, i_bin_pt, i_bin_mass);

      double v2_sig =
          (v2_target * count_uls - v2_targetbg * (count_uls - count_signal)) /
          (count_signal);
      double ntot2nsig = count_uls / count_signal;
      double ntot2nsig_err = sqrt(
          pow(count_uls_err / count_signal, 2) +
          pow(count_signal_err * count_uls / count_signal / count_signal, 2));
      double v2_sig_err1 = ntot2nsig * v2_target_err;
      double v2_sig_err2 = ntot2nsig_err * v2_target;
      // double v2_sig_err3 = (ntot2nsig - 1) * v2_targetbg_err;
      double v2_sig_err4 = v2_targetbg * ntot2nsig_err;
      double v2_sig_err =
          sqrt(pow(v2_sig_err1, 2) + pow(v2_sig_err2, 2) + pow(v2_sig_err4, 2));
      h2_v2_4_sig_mass_pt->SetBinContent(i_bin_mass, i_bin_pt, v2_sig);
      h2_v2_4_sig_mass_pt->SetBinError(i_bin_mass, i_bin_pt, v2_sig_err);
    }
  }

  int index_bin_mass_flow_low = h2_v2_4_sig_mass_pt->GetXaxis()->FindBin(2.72);
  int index_bin_mass_flow_high = h2_v2_4_sig_mass_pt->GetXaxis()->FindBin(3.2);

  TH1F *h1_v2_4_sig_pt =
      new TH1F("h1_v2_4_sig_pt", ";#it{p}_{T} [GeV/c];v_{2}{4}_{3 sub events}",
               n_bins_pt, bins_pt);
  for (int i_bin_pt = 1; i_bin_pt <= n_bins_pt; i_bin_pt++) {
    double value = 0;
    double value_err = 0;
    for (int i_bin_mass = index_bin_mass_flow_low;
         i_bin_mass <= index_bin_mass_flow_high; i_bin_mass++) {
      double bin_err = h2_v2_4_sig_mass_pt->GetBinError(i_bin_mass, i_bin_pt);
      double bin_content =
          h2_v2_4_sig_mass_pt->GetBinContent(i_bin_mass, i_bin_pt);
      value += bin_content / bin_err / bin_err;
      value_err += 1. / bin_err / bin_err;
    }
    value /= value_err;
    value_err = sqrt(1. / value_err);
    cout << "pt:\t" << bins_pt[i_bin_pt - 1] << "\t" << value << "\t"
         << value_err << endl;
    h1_v2_4_sig_pt->SetBinContent(i_bin_pt, value);
    h1_v2_4_sig_pt->SetBinError(i_bin_pt, value_err);
  }
  TCanvas *c_sig_flow = new TCanvas("c_sig_flow", "c_sig_flow", 600, 600);
  c_sig_flow->Divide(2, 2);
  c_sig_flow->cd(1);
  h2_v2_4_sig_mass_pt->Draw("colz");
  c_sig_flow->cd(2);
  h1_v2_4_sig_pt->Draw();
  c_sig_flow->cd(3);
  TH1D *h1_h2_v2_4_sig_mass_pt =
      (TH1D *)h2_v2_4_sig_mass_pt->ProjectionX("h2_v2_4_sig_mass_x", 2, 2);
  h1_h2_v2_4_sig_mass_pt->GetYaxis()->SetRangeUser(-1, 1);
  h1_h2_v2_4_sig_mass_pt->GetYaxis()->SetTitle("v_{2}{4}_{Signal}");
  TF1 *f1_h1_h2_v2_4_sig_mass_pt =
      new TF1("f1_h1_h2_v2_4_sig_mass_pt", "pol0", 2.72, 3.2);
  h1_h2_v2_4_sig_mass_pt->Fit(f1_h1_h2_v2_4_sig_mass_pt, "R");
  h1_h2_v2_4_sig_mass_pt->Draw();

  c_sig_flow->cd(4);
  vec_h1_mass_fit_flow[0]->Draw();
  for (int i = 1; i < vec_h1_mass_fit_flow.size(); i++) {
    vec_h1_mass_fit_flow[i]->Draw("same");
  }
  server->Register("c_sig_flow", c_sig_flow);

  TCanvas *c1_C2_4_QA = new TCanvas("c1_C2_4_QA", "c1_C2_4_QA", 600, 600);
  c1_C2_4_QA->Divide(2, 2);
  c1_C2_4_QA->cd(1);
  h3_widthEta_pt_mass_Nv2_4_NmNbg
      ->ProjectionZ("h3_widthEta_pt_mass_Nv2_4_NmNbg_z", 4, 4, 1, 1)
      ->Draw();
  c1_C2_4_QA->cd(2);
  h3pre_widthEta_pt_mass_v2_4
      ->ProjectionZ("h3pre_widthEta_pt_mass_v2_4_z", 4, 4, 1, 1)
      ->Draw();
  c1_C2_4_QA->cd(3);
  TH1D *h1_z_h3target_widthEta_pt_mass_v2_4 =
      (TH1D *)h3target_widthEta_pt_mass_v2_4->ProjectionZ(
          "h3target_widthEta_pt_mass_v2_4_z", 5, 5, 2, 2);
  h1_z_h3target_widthEta_pt_mass_v2_4->GetXaxis()->SetRangeUser(-1, 1);
  h1_z_h3target_widthEta_pt_mass_v2_4->Draw();

  c1_C2_4_QA->cd(4);
  TH1D *h1_h3targetbg_widthEta_pt_mass_v2_2_z =
      h3targetbg_widthEta_pt_mass_v2_4->ProjectionZ(
          "h3targetbg_widthEta_pt_mass_v2_4_z", 5, 5, 2, 2);
  TF1 *f1 = new TF1("f1", "pol1", 2.4, 3.6);
  h1_h3targetbg_widthEta_pt_mass_v2_2_z->Fit(f1, "R");
  h1_h3targetbg_widthEta_pt_mass_v2_2_z->Draw();
  server->Register("c1_C2_4_QA", c1_C2_4_QA);
}