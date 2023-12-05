#include "MHead.h"

void FlowTest() {
  TFile* f_input = new TFile("../../FlowTest/AnalysisResults.root"); 
  THashList* hash_list = (THashList*)f_input->Get("analysis-same-event-pairing/output");
  THashList* hash_list_flow = (THashList*)f_input->Get("analysis-dilepton-hadron/output");
  hash_list->Print();
  TH1* h_unlike = (TH1*)((TList*)hash_list->FindObject("PairsBarrelSEPM_Jpsi_TPCPost_calib_debug4"))->FindObject("Mass");
  h_unlike->GetXaxis()->SetRangeUser(2, 4);
  h_unlike->GetYaxis()->SetRangeUser(0., 140);
  h_unlike->Draw();
  TH1* h_like1 = (TH1*)((TList*)hash_list->FindObject("PairsBarrelSEPP_Jpsi_TPCPost_calib_debug4"))->FindObject("Mass");
  // h_like1->Draw("same");
  TH1* h_like2 = (TH1*)((TList*)hash_list->FindObject("PairsBarrelSEMM_Jpsi_TPCPost_calib_debug4"))->FindObject("Mass");
  h_like2->Add(h_like1);
  h_like2->Draw("same");
  double number_total = h_unlike->Integral(h_unlike->GetXaxis()->FindBin(2.8), h_unlike->GetXaxis()->FindBin(3.2));
  double number_bkg = h_like2->Integral(h_unlike->GetXaxis()->FindBin(2.8), h_unlike->GetXaxis()->FindBin(3.2));

  vector<TProfile*> vec_pro;
  vector<string> vec_str_name_profile = {"profile_QqbarLikeSignRe", "profile_QqbarUnlikeSignRe", "profile_QLikeSignRe", "profile_QUnlikeSignRe", "profile_qbarLikeSignRe", "profile_qbarUnlikeSignRe", "profile_QqbarLikeSignIm", "profile_QqbarUnlikeSignIm", "profile_QLikeSignIm", "profile_QUnlikeSignIm", "profile_qbarLikeSignIm", "profile_qbarUnlikeSignIm"};
  for (int i = 0; i < vec_str_name_profile.size(); i++) {
    TProfile* pro = (TProfile*)((TList*)hash_list_flow->FindObject("DileptonHadronFlow"))->FindObject(vec_str_name_profile[i].c_str());
    vec_pro.push_back(pro);
  }
  // QqbarRe - QRe*qbarRe + QIm*qbarIm
  // err = sqrt(err(QqbarRe)^2 + err(QRe)^2*qbarRe^2 + err(qbarRe)^2*QRe^2 + err(QIm)^2*qbarIm^2 + err(qbarIm)^2*QIm^2
  // for likesign
  TProfile* pro_flow2_likesign = (TProfile*)vec_pro[0]->Clone("pro_flow2_likesign");
  TProfile* pro_flow2_unlikesign = (TProfile*)vec_pro[1]->Clone("pro_flow2_unlikesign");
  for (int i_bin = 1; i_bin <= vec_pro[0]->GetNbinsX();i_bin++)
  {
    double QqbarRe = vec_pro[0]->GetBinContent(i_bin);
    double QRe = vec_pro[2]->GetBinContent(i_bin); 
    double qbarRe = vec_pro[4]->GetBinContent(i_bin);
    double QIm = vec_pro[8]->GetBinContent(i_bin);
    double qbarIm = vec_pro[10]->GetBinContent(i_bin);
    double flow2 = QqbarRe - QRe*qbarRe + QIm*qbarIm;
    pro_flow2_likesign->SetBinContent(i_bin, flow2);
    double err = sqrt(pow(vec_pro[0]->GetBinError(i_bin), 2) + pow(vec_pro[2]->GetBinError(i_bin)*qbarRe, 2) + pow(vec_pro[4]->GetBinError(i_bin)*QRe, 2) + pow(vec_pro[8]->GetBinError(i_bin)*qbarIm, 2) + pow(vec_pro[10]->GetBinError(i_bin)*QIm, 2));
    pro_flow2_likesign->SetBinError(i_bin, err);
  }
  // for unlikesign
  for (int i_bin = 1; i_bin <= vec_pro[1]->GetNbinsX();i_bin++)
  {
    double QqbarRe = vec_pro[1]->GetBinContent(i_bin);
    double QRe = vec_pro[3]->GetBinContent(i_bin); 
    double qbarRe = vec_pro[5]->GetBinContent(i_bin);
    double QIm = vec_pro[9]->GetBinContent(i_bin);
    double qbarIm = vec_pro[11]->GetBinContent(i_bin);
    double flow2 = QqbarRe - QRe*qbarRe + QIm*qbarIm;
    pro_flow2_unlikesign->SetBinContent(i_bin, flow2);
    double err = sqrt(pow(vec_pro[1]->GetBinError(i_bin), 2) + pow(vec_pro[3]->GetBinError(i_bin)*qbarRe, 2) + pow(vec_pro[5]->GetBinError(i_bin)*QRe, 2) + pow(vec_pro[9]->GetBinError(i_bin)*qbarIm, 2) + pow(vec_pro[11]->GetBinError(i_bin)*QIm, 2));
    pro_flow2_unlikesign->SetBinError(i_bin, err);
  }
  

  // profile all 
  TProfile* pro_flow2 = (TProfile*)pro_flow2_likesign->Clone("pro_flow2");
  for (int i_bin = 1; i_bin <= vec_pro[1]->GetNbinsX();i_bin++)
  {
    double flow2 = (number_total-number_bkg)/number_total*pro_flow2_unlikesign->GetBinContent(i_bin) - number_bkg/number_total*pro_flow2_likesign->GetBinContent(i_bin);
    double err1 = sqrt(pow((number_total-number_bkg)/number_total*pro_flow2_unlikesign->GetBinError(i_bin), 2) + pow(number_bkg/number_total*pro_flow2_likesign->GetBinError(i_bin), 2));
    double err2 = (number_bkg/number_total)*sqrt(1./number_bkg+1./number_total)*sqrt(pow(pro_flow2_unlikesign->GetBinContent(i_bin),2)+pow(pro_flow2_likesign->GetBinContent(i_bin),2));
    double err = sqrt(pow(err1, 2) + pow(err2, 2));
    pro_flow2->SetBinContent(i_bin, flow2);
    pro_flow2->SetBinError(i_bin, err);
  }
  // for all
  pro_flow2->Draw();
  pro_flow2->GetYaxis()->SetRangeUser(-0.00005, 0.00005);

}