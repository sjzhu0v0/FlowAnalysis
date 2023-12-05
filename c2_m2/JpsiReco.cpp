#include "MHead.h"

#include "../src/Fitter.cpp"
#include "../src/libmytool.cpp"
#include "Fitter.h"
#include "libmytool.h"

#define verbose 2
#define roothttp

double width_eta_jpsi = 0.8;

mysignal Jpsi;
fitter *myfunc = new fitter();
TFile *resultfile = nullptr;
TF1 *fsignal = nullptr;
TF1 *fbkg = nullptr;
TF1 *fglobal = nullptr;

TH1F *hSE_unlike = nullptr;
TH1F *hSE_like = nullptr;

typedef struct {
  TCanvas *canvas;
  mysignal signal;
  TH1 *Rawhisto;
  TH1D* h1_unlikeSign;
  TH1D* h1_likeSign;
  TF1 *fResidualBkg;
} ContainerRecoJpsi;

void Init() {
  // initiate the fit function
  std::vector<double> signal_lowlimit = {1.0, 3.0, 0.02, 0.2, 0.2};
  std::vector<double> signal_uplimit = {1e8, 3.2, 0.065, 10, 100};
  myfunc->set_massrange(1.8, 4.2);
  myfunc->set_binwidth(0.04);
  myfunc->set_sideband(2.6, 3.2);
  myfunc->set_fitmethod("unlike_like");
  fsignal = myfunc->get_func("CB1");
  fbkg = myfunc->get_func("pol2");
  std::function<double(double *, double *)> mfglobal = [&](double *x,
                                                           double *p) {
    int npar = fsignal->GetNpar();
    return fsignal->EvalPar(x, &p[0]) + fbkg->EvalPar(x, &p[npar]);
  };
  fglobal =
      new TF1("global", mfglobal, 0, 10, fsignal->GetNpar() + fbkg->GetNpar());
  myfunc->set_fitfunction(fsignal, fbkg, fglobal, fbkg);
  // set parlimits
  for (int ipar = 0; ipar < fsignal->GetNpar(); ipar++) {
    myfunc->fGlobal->SetParLimits(ipar, signal_lowlimit[ipar],
                                  signal_uplimit[ipar]);
  }
  // set signal
  Jpsi.kMasswindow[0] = 2.72;
  Jpsi.kMasswindow[1] = 3.2;
  Jpsi.kWindowcolor = 1;
  myfunc->fSignal->SetParameters(1e6, 3.06, 0.04, 0.3, 1);
}

ContainerRecoJpsi FitJpsi(TH1D *hSE_unlike, TH1D *hSE_like) {
  ContainerRecoJpsi container;
  Init();
  myfunc->fit(hSE_unlike, hSE_like);
  Jpsi = myfunc->calculate(Jpsi);
  container.signal = Jpsi;
  TCanvas *c1 = myfunc->save_fitcanvas();
  c1->cd(1);
  tool::SetLatex(0.21, 0.85, Form("ALICE Run3 pp #sqrt{s}=13.6 TeV"), 22, 0.055,
                 1);
  tool::SetLatex(0.21, 0.77, Form("%s", "QA"), 22, 0.055, 2);
  tool::SetLatex(0.21, 0.69, Form("p_{T} integrated, |#eta| < 0.4"), 22, 0.055,
                 2);
  c1->cd(2);
  tool::SetLatex(
      0.21, 0.92,
      Form("N_{J/#psi} = %.0f #pm %.0f", Jpsi.kRawsignal, Jpsi.kRawsignal_err),
      22, 0.045, 1);
  tool::SetLatex(
      0.21, 0.84,
      Form("S/B = %.2f #pm %.2f", Jpsi.kS_over_B, Jpsi.kS_over_B_err), 22,
      0.045, 1);
  tool::SetLatex(0.21, 0.76,
                 Form("#frac{S}{#sqrt{S + 2B}} = %.2f", Jpsi.kSignificance), 22,
                 0.045, 1);
  tool::SetLine(Jpsi.kMasswindow[0], 0, Jpsi.kMasswindow[0],
                myfunc->Rawhisto->GetMaximum() * 0.6, Jpsi.kWindowcolor, 2, 2);
  tool::SetLine(Jpsi.kMasswindow[1], 0, Jpsi.kMasswindow[1],
                myfunc->Rawhisto->GetMaximum() * 0.6, Jpsi.kWindowcolor, 2, 2);

  tool::SetLatex(0.71, 0.65,
                 Form("#mu = %.3f #pm %.3f", myfunc->fGlobal->GetParameter(1),
                      myfunc->fGlobal->GetParError(1)),
                 22, 0.045, 1);
  tool::SetLatex(0.71, 0.57,
                 Form("#sigma = %.3f #pm %.3f",
                      myfunc->fGlobal->GetParameter(2),
                      myfunc->fGlobal->GetParError(2)),
                 22, 0.045, 1);
  container.canvas = (TCanvas *)c1->Clone();
  container.Rawhisto =
      (TH1 *)myfunc->Rawhisto->Clone(Form("Rawhisto_%d", GenerateUID()));
  container.fResidualBkg = (TF1 *)myfunc->fResidualBkg->Clone(
      Form("fResidualBkg_%d", GenerateUID()));
  return container;
}

void JpsiReco(string input_file =
                  "/data/work/work_alice/Analysis/JpsiFlow/input/"
                  "AnalysisResults_new2.root",
              string output_file = "/data/work/work_alice/Analysis/JpsiFlow/"
                                   "c2_m2/output/JpsiReco.root") {
#ifdef roothttp
  THttpServer *server = new THttpServer("http:8090?top=job");
  server->SetReadOnly(kFALSE);
  gBenchmark->Start("job");
#endif

  TFile *file_input = new TFile(input_file.c_str(), "READ");
  TFile *file_output = new TFile(output_file.c_str(), "RECREATE");
  vector<string> vec_str_name_obj;
  vector<TObject *> vec_obj = GetObjectRecursive(file_input, vec_str_name_obj);
  vector<TH1 *> vec_h1 = GetObjectVector<TH1>(
      vec_obj, vec_str_name_obj, vector<string>{}, (TString) "JpsiFlow");
  vector<TH3F *> vec_h1_mass = GetObjectVector<TH3F>(
      vec_obj, vec_str_name_obj, vector<string>{}, (TString) "Eta_Pt_Mass");

#if verbose == 2
  for (auto &h1 : vec_h1)
    h1->Print();

  cout << endl;

  for (auto &h1 : vec_h1_mass)
    h1->Print();
#endif

  vector<TH1D *> vec_h1_mass_pt_ls;
  vector<TH1D *> vec_h1_mass_pt_uls;
  vector<ContainerRecoJpsi> vec_container_reco_jpsi;

  // vec_h1_mass[0]->GetXaxis()->SetRangeUser(-0.4, 0.4);
  // vec_h1_mass[1]->GetXaxis()->SetRangeUser(-0.4, 0.4);
  vec_h1_mass[0]->ProjectionZ("h1", 0, -1, 0, -1)->Draw();
  cout << vec_h1_mass[0]->ProjectionZ("h1", 0, -1, 0, -1)->Integral() << endl;
  cout << vec_h1_mass[1]->ProjectionZ("h1_c", 0, -1, 0, -1)->Integral() << endl;
  // vec_h1_mass[1]->ProjectionZ("h1", 0, -1, 0, -1)->Draw();
  // ContainerRecoJpsi container_reco_jpsi =
  //     FitJpsi(vec_h1_mass[0]->ProjectionZ("h1_a", 0, -1, 0, -1),
  //             vec_h1_mass[1]->ProjectionZ("h1_b", 0, -1, 0, -1));

  for (int i_pt = 1; i_pt <= vec_h1_mass[0]->GetYaxis()->GetNbins(); i_pt++) {
    TH1D *h1_ls =
        vec_h1_mass[1]->ProjectionZ(Form("h1_ls_%d", i_pt), 6, 13, i_pt, i_pt);
    TH1D *h1_uls =
        vec_h1_mass[0]->ProjectionZ(Form("h1_uls_%d", i_pt), 6, 13, i_pt, i_pt);
    vec_h1_mass_pt_ls.push_back(h1_ls);
    vec_h1_mass_pt_uls.push_back(h1_uls);
    ContainerRecoJpsi container_reco_jpsi = FitJpsi(h1_uls, h1_ls);
    container_reco_jpsi.Rawhisto->SetName(Form("Rawhisto_Pt_%d", i_pt));
    container_reco_jpsi.fResidualBkg->SetName(Form("fResidualBkg_Pt_%d", i_pt));
    container_reco_jpsi.h1_unlikeSign = h1_uls;
    container_reco_jpsi.h1_likeSign = h1_ls;
    vec_container_reco_jpsi.push_back(container_reco_jpsi);
  }

#ifdef roothttp
  for (int i_container = 0; i_container < vec_container_reco_jpsi.size();
       i_container++) {
    server->Register(Form("jpsi%d", i_container),
                     vec_container_reco_jpsi[i_container].canvas);
  }
#endif
  TH1F *h1_N_jpsi =
      new TH1F("h1_N_jpsi", "", vec_h1_mass[0]->GetYaxis()->GetNbins(),
               vec_h1_mass[0]->GetYaxis()->GetXbins()->fArray);
  TH1F *h1_N_total =
      new TH1F("h1_N_total", "", vec_h1_mass[0]->GetYaxis()->GetNbins(),
               vec_h1_mass[0]->GetYaxis()->GetXbins()->fArray);
  for (int i_pt = 1; i_pt <= vec_h1_mass[0]->GetYaxis()->GetNbins(); i_pt++) {
    h1_N_jpsi->SetBinContent(
        i_pt, vec_container_reco_jpsi[i_pt - 1].signal.kRawsignal);
    h1_N_jpsi->SetBinError(
        i_pt, vec_container_reco_jpsi[i_pt - 1].signal.kRawsignal_err);
    h1_N_total->SetBinContent(i_pt, vec_h1_mass_pt_uls[i_pt - 1]->Integral());
    h1_N_total->SetBinError(i_pt,
                            sqrt(vec_h1_mass_pt_uls[i_pt - 1]->Integral()));
  }

  server->Register("N_jpsi", h1_N_jpsi);
  server->Register("N_total", h1_N_total);

#ifdef roothttp
  gBenchmark->Show("job");
#endif

  if (file_output != nullptr) {
    for (auto &container_reco_jpsi : vec_container_reco_jpsi) {
      container_reco_jpsi.canvas->Write();
      container_reco_jpsi.Rawhisto->Write();
      container_reco_jpsi.fResidualBkg->Write();
      container_reco_jpsi.h1_unlikeSign->Write();
      container_reco_jpsi.h1_likeSign->Write();
    }
    h1_N_jpsi->Write();
    h1_N_total->Write();
    file_output->Close();
  }
}