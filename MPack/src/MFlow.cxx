#include "MFlow.h"

MFlow::MFlow(int maxHarmonic, int maxPower)
{
  TComplex zero(0., 0.);
  for (int i = 0; i < maxHarmonic; i++) {
    std::vector<TComplex> temp;
    for (int j = 0; j < maxPower; j++) {
      temp.push_back(zero);
    }
    fQvectors.push_back(temp);
  }
}

MFlow::~MFlow()
{
}

void MFlow::Fill(double phi, double weight)
{
  for (int i = 0; i < (int)fQvectors.size(); i++) {
    for (int j = 0; j < (int)fQvectors[i].size(); j++) {
      fQvectors[i][j] += TComplex::Exp(i * phi * TComplex::I()) * pow(weight, j);
    }
  }
}

TComplex MFlow::Q(int harmonic, int power)
{
  if (abs(harmonic) >= (int)fQvectors.size() || abs(power) >= (int)fQvectors[harmonic].size()) {
    return TComplex(0., 0.);
  }
  if (harmonic < 0) {
    return TComplex::Conjugate(fQvectors[abs(harmonic)][power]);
  }
  return fQvectors[harmonic][power];
}

void MFlow::Reset()
{
  for (int i = 0; i < (int)fQvectors.size(); i++) {
    for (int j = 0; j < (int)fQvectors[i].size(); j++) {
      fQvectors[i][j] = TComplex(0., 0.);
    }
  }
}

void MFlow::Print()
{
  cout << "MFlow::Print(): " << this << endl;
  for (int i = 0; i < (int)fQvectors.size(); i++) {
    for (int j = 0; j < (int)fQvectors[i].size(); j++) {
      cout << "Q(" << i << "," << j << ") = " << fQvectors[i][j] << endl;
    }
  }
}

MFlow& MFlow::operator+=(const MFlow& other)
{
  if (fQvectors.size() != other.fQvectors.size() ||
      fQvectors[0].size() != other.fQvectors[0].size()) {
    return *this;
  }

  for (int i = 0; i < (int)fQvectors.size(); i++) {
    for (int j = 0; j < (int)fQvectors[i].size(); j++) {
      fQvectors[i][j] += other.fQvectors[i][j];
    }
  }

  return *this;
}

MFlow& MFlow::operator-=(const MFlow& other)
{
  if (fQvectors.size() != other.fQvectors.size() ||
      fQvectors[0].size() != other.fQvectors[0].size()) {
    return *this;
  }

  for (int i = 0; i < (int)fQvectors.size(); i++) {
    for (int j = 0; j < (int)fQvectors[i].size(); j++) {
      fQvectors[i][j] -= other.fQvectors[i][j];
    }
  }
  return *this;
}

MFlow MFlow::operator+(const MFlow& other)
{
  MFlow temp = *this;
  temp += other;
  return temp;
}

MFlow MFlow::operator-(const MFlow& other)
{
  MFlow temp = *this;
  temp -= other;
  return temp;
}

MFlow MFlow::operator=(const MFlow& other)
{
  for (int i = 0; i < (int)fQvectors.size(); i++) {
    for (int j = 0; j < (int)fQvectors[i].size(); j++) {
      fQvectors[i][j] = other.fQvectors[i][j];
    }
  }
  return *this;
}

TString MFGF::fgVariableUnits[] = {"#it{p}_{T}[GeV/c]", "#it{N}_{ch}", "#eta", "Weight", "Im(Cumulant)", "Re(Cumulant)"};
TString MFGF::fgVariableNames[] = {"Pt", "Mult", "Eta", "Weight", "CumulantIm", "CumulantRe"};

std::array<MFlow, MFGV::NBinsEta> MFGF::fArrayFlowReferenceEta;
std::array<std::array<MFlow, MFGV::NBinsPt>, MFGV::NBinsEta> MFGF::fArrayFlowTargetEtaPt;
std::array<std::array<MFlow, MFGV::NBinsPt>, MFGV::NBinsEta> MFGF::fArrayFlowTargetBgEtaPt;
std::array<std::array<std::array<MFlow, MFGV::NBinsPt>, MFGV::NBinsEta>, MFGV::NBinsMass> MFGF::fArrayFlowTargetEtaPtMass;
std::array<std::array<std::array<MFlow, MFGV::NBinsPt>, MFGV::NBinsEta>, MFGV::NBinsMass> MFGF::fArrayFlowTargetBgEtaPtMass;

double MFGF::fMinAbsEtaReference = 0;
double MFGF::fMaxAbsEtaReference = 0.4;
double MFGF::fMinAbsEtaTarget = 0.4;
double MFGF::fMaxAbsEtaTarget = 0.9;
double MFGF::fMultiplicity = 30;
double MFGF::fMass = 0;

TList* MFGF::fListOutput;
TList* MFGF::fListQA;
TList* MFGF::fListReference;
TList* MFGF::fListTarget;
TList* MFGF::fListCorrelation;

TFile* MFGF::fFileInput;
// input histograms
TH2F* MFGF::fh2NonUniformEtaPhiMultReference;
TH3F* MFGF::fh3EtaPhiMultReference;
TH3F* MFGF::fh3EtaPhiMultReferenceCorrected;
TH3F* MFGF::fh3EtaPhiMultBeforeCut;
TH1F* MFGF::fh1CutTrack;

int MFGF::NBinsDeltaEtaReferenceSubEvent = 0;
double* MFGF::BinDeltaEtaReferenceSubEvent;
int MFGF::NBinsEtaWidthReferenceTotal = 0;
double* MFGF::BinEtaWidthReferenceTotal;

////////////////////////////////////////////////////////////
// Reference Flow
TProfile2D* MFGF::fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[MFGV::N_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant][2];
TProfile2D* MFGF::fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[MFGV::N_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant][2];
TH2F* MFGF::fh2FlowReferenceTotal_EtaWidth_Mult_WeightSum[MFGV::N_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant][2];
TH2F* MFGF::fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum[MFGV::N_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant][2];

// Target Flow
TProfile2D* MFGF::fProfileTarget_EtaWidth_Pt_cumulant[MFGV::N_fProfileTarget_EtaWidth_Pt_cumulant][2];
TProfile2D* MFGF::fProfileTargetBg_EtaWidth_Pt_cumulant[MFGV::N_fProfileTarget_EtaWidth_Pt_Mass_cumulant][2];
TH2F* MFGF::fh2Target_EtaWidth_Pt_WeightSum[MFGV::N_fProfileTarget_EtaWidth_Pt_cumulant][2];
TH2F* MFGF::fh2TargetBg_EtaWidth_Pt_WeightSum[MFGV::N_fProfileTarget_EtaWidth_Pt_Mass_cumulant][2];

TProfile3D* MFGF::fProfileTarget_EtaWidth_Pt_Mass_cumulant[MFGV::N_fProfileTargetBg_EtaWidth_Pt_cumulant][2];
TProfile3D* MFGF::fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[MFGV::N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant][2];
TH3F* MFGF::fh3Target_EtaWidth_Pt_Mass_WeightSum[MFGV::N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant][2];
TH3F* MFGF::fh3TargetBg_EtaWidth_Pt_Mass_WeightSum[MFGV::N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant][2];
///////////////////////////////////////////////////////////////
TH3F* MFGF::fh3Target_Eta_Pt_Mass;
TH3F* MFGF::fh3TargetBg_Eta_Pt_Mass;
TH1F* MFGF::fh1DileptonMass_NoCut;
TH1F* MFGF::fh1DileptonMass_MultCut;

// Correlation
TH1F* MFGF::fh1DeltaPhiReference;
TH1F* MFGF::fh1DeltaPhiTargetReference;
TH1F* MFGF::fh1DeltaPhiTargetBgReference;

void MFGF::HistogramInit(bool doNeedNonUniformEtaPhiMultReference)
{
  if (doNeedNonUniformEtaPhiMultReference) {
    fFileInput = new TFile("nonuniform_effect.root", "READ");
    fh2NonUniformEtaPhiMultReference = (TH2F*)fFileInput->Get("h2_nonuniform_effect");
  }
  cout << "Histogram Initialing" << endl;

  fListOutput = new TList();
  fListOutput->SetName("JpsiFlow");
  fListOutput->SetOwner(kTRUE);

  fListQA = new TList();
  fListQA->SetName("QA");
  // fListQA->SetOwner(kTRUE);
  fListOutput->Add(fListQA);

  fListReference = new TList();
  fListReference->SetName("Reference");
  // fListReference->SetOwner(kTRUE);
  fListOutput->Add(fListReference);

  fListTarget = new TList();
  fListTarget->SetName("Target");
  // fListTarget->SetOwner(kTRUE);
  fListOutput->Add(fListTarget);

  fListCorrelation = new TList();
  fListCorrelation->SetName("Correlation");
  // fListCorrelation->SetOwner(kTRUE);
  fListOutput->Add(fListCorrelation);

  // bin 180,-0.9,0.9
  // bin 360,0,2pi
  // bin 100,-1,1
  // bin 100,2.8,3.2
  double bin_eta_temp[180 + 1];
  for (int i = 0; i < 180 + 1; i++)
    bin_eta_temp[i] = -0.9 + 0.01 * i;

  double bin_phi_temp[360 + 1];
  for (int i = 0; i < 360 + 1; i++)
    bin_phi_temp[i] = 0 + i * 2 * TMath::Pi() / 360.;

  // double bin_cumulant_temp[100 + 1];
  // for (int i = 0; i < 100 + 1; i++)
  //   bin_cumulant_temp[i] = -1 + 0.02 * i;

  double bin_mass[50 + 1];
  for (int i = 0; i < 50 + 1; i++)
    bin_mass[i] = 2. + 0.04 * i;

  double bin_mass_for_flow[MFGV::NBinsMass + 1];
  for (int i = 0; i < MFGV::NBinsMass + 1; i++)
    bin_mass_for_flow[i] = MFGV::EdgesMass[0] + i * (MFGV::EdgesMass[1] - MFGV::EdgesMass[0]) / (double)MFGV::NBinsMass;

  ///////////////////////////////////////////////////////////////////
  fh3EtaPhiMultReference = new TH3F("fh3EtaPhiMultReference", "fh3EtaPhiMultReference;#eta;#phi;Multiplicity", 180, bin_eta_temp, 360, bin_phi_temp, MFGV::NBinsMult, MFGV::BinMult);
  fh3EtaPhiMultReferenceCorrected = new TH3F("fh3EtaPhiMultReferenceCorrected", "fh3EtaPhiMultReferenceCorrected;#eta;#phi;Multiplicity", 180, bin_eta_temp, 360, bin_phi_temp, MFGV::NBinsMult, MFGV::BinMult);

  fh3EtaPhiMultBeforeCut = new TH3F("fh3EtaPhiMultBeforeCut", "fh3EtaPhiMultBeforeCut;#eta;#phi;Multiplicity", 180, bin_eta_temp, 360, bin_phi_temp, MFGV::NBinsMult, MFGV::BinMult);

  fh1CutTrack = new TH1F("fh1CutTrack", "Cut Track", 11, -1, 10);
  /////////////////////////////////////////////////////////////////////

  if (MFGV::NBinsEta % 2 == 0) {
    NBinsDeltaEtaReferenceSubEvent = MFGV::NBinsEta / 2;
    NBinsEtaWidthReferenceTotal = MFGV::NBinsEta / 2;
  } else {
    NBinsDeltaEtaReferenceSubEvent = (MFGV::NBinsEta - 1) / 2;
    NBinsEtaWidthReferenceTotal = (MFGV::NBinsEta + 1) / 2;
  }

  BinDeltaEtaReferenceSubEvent = new double[NBinsDeltaEtaReferenceSubEvent + 1];
  BinEtaWidthReferenceTotal = new double[NBinsEtaWidthReferenceTotal + 1];

  for (int i = 0; i < NBinsDeltaEtaReferenceSubEvent; i++) {
    BinDeltaEtaReferenceSubEvent[NBinsDeltaEtaReferenceSubEvent - i - 1] = MFGV::BinEta[MFGV::NBinsEta - i - 1] - MFGV::BinEta[i + 1];
  }
  BinDeltaEtaReferenceSubEvent[NBinsDeltaEtaReferenceSubEvent] = 2 * BinDeltaEtaReferenceSubEvent[NBinsDeltaEtaReferenceSubEvent - 1] - BinDeltaEtaReferenceSubEvent[NBinsDeltaEtaReferenceSubEvent - 2];

  for (int i = 0; i < NBinsEtaWidthReferenceTotal; i++) {
    int i_left = i;
    int i_right = MFGV::NBinsEta - i;
    BinEtaWidthReferenceTotal[NBinsEtaWidthReferenceTotal - i - 1] = MFGV::BinEta[i_right] - MFGV::BinEta[i_left];
  }
  BinEtaWidthReferenceTotal[NBinsEtaWidthReferenceTotal] = 2 * BinEtaWidthReferenceTotal[NBinsEtaWidthReferenceTotal - 1] - BinEtaWidthReferenceTotal[NBinsEtaWidthReferenceTotal - 2];
  ////////////////////// Profile Definition For Flow //////////////////////////

  for (int iprofile = 0; iprofile < MFGV::N_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant; iprofile++) {
    fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[iprofile][0] = new TProfile2D(Form("fProfileFlowReferenceTotal_EtaWidth_Mult_%s_im", MFGV::Name_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[iprofile].c_str()), Form("fProfileFlowReferenceTotal_EtaWidth_Mult_%s_im;#eta;Multiplicity", MFGV::Name_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsMult, MFGV::BinMult);

    fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[iprofile][1] = new TProfile2D(Form("fProfileFlowReferenceTotal_EtaWidth_Mult_%s_re", MFGV::Name_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[iprofile].c_str()), Form("fProfileFlowReferenceTotal_EtaWidth_Mult_%s_re;#eta;Multiplicity", MFGV::Name_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsMult, MFGV::BinMult);

    fh2FlowReferenceTotal_EtaWidth_Mult_WeightSum[iprofile][0] = new TH2F(Form("fh2FlowReferenceTotal_EtaWidth_Mult_WeightSum_%s", MFGV::Name_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[iprofile].c_str()), Form("fh2FlowReferenceTotal_EtaWidth_Mult_WeightSum_%s;#eta;Multiplicity", MFGV::Name_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsMult, MFGV::BinMult);

    fh2FlowReferenceTotal_EtaWidth_Mult_WeightSum[iprofile][1] = new TH2F(Form("fh2FlowReferenceTotal_EtaWidth_Mult_EventSum_%s", MFGV::Name_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[iprofile].c_str()), Form("fh2FlowReferenceTotal_EtaWidth_Mult_EventSum_%s;#eta;Multiplicity", MFGV::Name_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsMult, MFGV::BinMult);
  }

  for (int iprofile = 0; iprofile < MFGV::N_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant; iprofile++) {
    fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[iprofile][0] = new TProfile2D(Form("fProfileFlowReferenceSubEvent_DeltaEta_Mult_%s_im", MFGV::Name_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[iprofile].c_str()), Form("fProfileFlowReferenceSubEvent_DeltaEta_Mult_%s_im;#Delta#eta;Multiplicity", MFGV::Name_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[iprofile].c_str()), NBinsDeltaEtaReferenceSubEvent, BinDeltaEtaReferenceSubEvent, MFGV::NBinsMult, MFGV::BinMult);

    fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[iprofile][1] = new TProfile2D(Form("fProfileFlowReferenceSubEvent_DeltaEta_Mult_%s_re", MFGV::Name_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[iprofile].c_str()), Form("fProfileFlowReferenceSubEvent_DeltaEta_Mult_%s_re;#Delta#eta;Multiplicity", MFGV::Name_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[iprofile].c_str()), NBinsDeltaEtaReferenceSubEvent, BinDeltaEtaReferenceSubEvent, MFGV::NBinsMult, MFGV::BinMult);

    fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum[iprofile][0] = new TH2F(Form("fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum_%s", MFGV::Name_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[iprofile].c_str()), Form("fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum_%s;#Delta#eta;Multiplicity", MFGV::Name_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[iprofile].c_str()), NBinsDeltaEtaReferenceSubEvent, BinDeltaEtaReferenceSubEvent, MFGV::NBinsMult, MFGV::BinMult);

    fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum[iprofile][1] = new TH2F(Form("fh2FlowReferenceSubEvent_DeltaEta_Mult_EventSum_%s", MFGV::Name_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[iprofile].c_str()), Form("fh2FlowReferenceSubEvent_DeltaEta_Mult_EventSum_%s;#Delta#eta;Multiplicity", MFGV::Name_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[iprofile].c_str()), NBinsDeltaEtaReferenceSubEvent, BinDeltaEtaReferenceSubEvent, MFGV::NBinsMult, MFGV::BinMult);
  }

  for (int iprofile = 0; iprofile < MFGV::N_fProfileTarget_EtaWidth_Pt_cumulant; iprofile++) {
    fProfileTarget_EtaWidth_Pt_cumulant[iprofile][0] = new TProfile2D(Form("fProfileTarget_EtaWidth_Pt_%s_im", MFGV::Name_fProfileTarget_EtaWidth_Pt_cumulant[iprofile].c_str()), Form("fProfileTarget_EtaWidth_Pt_%s_im;#eta;#it{p}_{T}[GeV/c]", MFGV::Name_fProfileTarget_EtaWidth_Pt_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt);

    fProfileTarget_EtaWidth_Pt_cumulant[iprofile][1] = new TProfile2D(Form("fProfileTarget_EtaWidth_Pt_%s_re", MFGV::Name_fProfileTarget_EtaWidth_Pt_cumulant[iprofile].c_str()), Form("fProfileTarget_EtaWidth_Pt_%s_re;#eta;#it{p}_{T}[GeV/c]", MFGV::Name_fProfileTarget_EtaWidth_Pt_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt);

    fh2Target_EtaWidth_Pt_WeightSum[iprofile][0] = new TH2F(Form("fh2Target_EtaWidth_Pt_WeightSum_%s", MFGV::Name_fProfileTarget_EtaWidth_Pt_cumulant[iprofile].c_str()), Form("fh2Target_EtaWidth_Pt_WeightSum_%s;#eta;#it{p}_{T}[GeV/c]", MFGV::Name_fProfileTarget_EtaWidth_Pt_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt);

    fh2Target_EtaWidth_Pt_WeightSum[iprofile][1] = new TH2F(Form("fh2Target_EtaWidth_Pt_EventSum_%s", MFGV::Name_fProfileTarget_EtaWidth_Pt_cumulant[iprofile].c_str()), Form("fh2Target_EtaWidth_Pt_EventSum_%s;#eta;#it{p}_{T}[GeV/c]", MFGV::Name_fProfileTarget_EtaWidth_Pt_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt);
  }

  for (int iprofile = 0; iprofile < MFGV::N_fProfileTarget_EtaWidth_Pt_Mass_cumulant; iprofile++) {
    fProfileTarget_EtaWidth_Pt_Mass_cumulant[iprofile][0] = new TProfile3D(Form("fProfileTarget_EtaWidth_Pt_Mass_%s_im", MFGV::Name_fProfileTarget_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), Form("fProfileTarget_EtaWidth_Pt_Mass_%s_im;#eta;#it{p}_{T}[GeV/c];Mass_{e^{+}e^{-}}[GeV/c^{2}]", MFGV::Name_fProfileTarget_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt, MFGV::NBinsMass, bin_mass_for_flow);

    fProfileTarget_EtaWidth_Pt_Mass_cumulant[iprofile][1] = new TProfile3D(Form("fProfileTarget_EtaWidth_Pt_Mass_%s_re", MFGV::Name_fProfileTarget_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), Form("fProfileTarget_EtaWidth_Pt_Mass_%s_re;#eta;#it{p}_{T}[GeV/c];Mass_{e^{+}e^{-}}[GeV/c^{2}]", MFGV::Name_fProfileTarget_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt, MFGV::NBinsMass, bin_mass_for_flow);

    fh3Target_EtaWidth_Pt_Mass_WeightSum[iprofile][0] = new TH3F(Form("fh3Target_EtaWidth_Pt_Mass_WeightSum_%s", MFGV::Name_fProfileTarget_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), Form("fh3Target_EtaWidth_Pt_Mass_WeightSum_%s;#eta;#it{p}_{T}[GeV/c];Mass_{e^{+}e^{-}}[GeV/c^{2}]", MFGV::Name_fProfileTarget_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt, MFGV::NBinsMass, bin_mass_for_flow);

    fh3Target_EtaWidth_Pt_Mass_WeightSum[iprofile][1] = new TH3F(Form("fh3Target_EtaWidth_Pt_Mass_WeightSum_%s", MFGV::Name_fProfileTarget_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), Form("fh3Target_EtaWidth_Pt_Mass_EventSum_%s;#eta;#it{p}_{T}[GeV/c];Mass_{e^{+}e^{-}}[GeV/c^{2}]", MFGV::Name_fProfileTarget_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt, MFGV::NBinsMass, bin_mass_for_flow);
  }

  for (int iprofile = 0; iprofile < MFGV::N_fProfileTargetBg_EtaWidth_Pt_cumulant; iprofile++) {
    fProfileTargetBg_EtaWidth_Pt_cumulant[iprofile][0] = new TProfile2D(Form("fProfileTargetBg_EtaWidth_Pt_%s_im", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_cumulant[iprofile].c_str()), Form("fProfileTargetBg_EtaWidth_Pt_%s_im;#eta;#it{p}_{T}[GeV/c]", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt);

    fProfileTargetBg_EtaWidth_Pt_cumulant[iprofile][1] = new TProfile2D(Form("fProfileTargetBg_EtaWidth_Pt_%s_re", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_cumulant[iprofile].c_str()), Form("fProfileTargetBg_EtaWidth_Pt_%s_re;#eta;#it{p}_{T}[GeV/c]", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt);

    fh2TargetBg_EtaWidth_Pt_WeightSum[iprofile][0] = new TH2F(Form("fh2TargetBg_EtaWidth_Pt_WeightSum_%s", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_cumulant[iprofile].c_str()), Form("fh2TargetBg_EtaWidth_Pt_WeightSum_%s;#eta;#it{p}_{T}[GeV/c]", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt);
    
    fh2TargetBg_EtaWidth_Pt_WeightSum[iprofile][1] = new TH2F(Form("fh2TargetBg_EtaWidth_Pt_EventSum_%s", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_cumulant[iprofile].c_str()), Form("fh2TargetBg_EtaWidth_Pt_EventSum_%s;#eta;#it{p}_{T}[GeV/c]", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt);
  }

  for (int iprofile = 0; iprofile < MFGV::N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant; iprofile++) {
    fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[iprofile][0] = new TProfile3D(Form("fProfileTargetBg_EtaWidth_Pt_Mass_%s_im", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), Form("fProfileTargetBg_EtaWidth_Pt_Mass_%s_im;#eta;#it{p}_{T}[GeV/c];Mass_{e^{+}e^{-}}[GeV/c^{2}]", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt, MFGV::NBinsMass, bin_mass_for_flow);

    fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[iprofile][1] = new TProfile3D(Form("fProfileTargetBg_EtaWidth_Pt_Mass_%s_re", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), Form("fProfileTargetBg_EtaWidth_Pt_Mass_%s_re;#eta;#it{p}_{T}[GeV/c];Mass_{e^{+}e^{-}}[GeV/c^{2}]", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt, MFGV::NBinsMass, bin_mass_for_flow);

    fh3TargetBg_EtaWidth_Pt_Mass_WeightSum[iprofile][0] = new TH3F(Form("fh3TargetBg_EtaWidth_Pt_Mass_WeightSum_%s", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), Form("fh3TargetBg_EtaWidth_Pt_Mass_WeightSum_%s;#eta;#it{p}_{T}[GeV/c];Mass_{e^{+}e^{-}}[GeV/c^{2}]", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt, MFGV::NBinsMass, bin_mass_for_flow);

    fh3TargetBg_EtaWidth_Pt_Mass_WeightSum[iprofile][0] = new TH3F(Form("fh3TargetBg_EtaWidth_Pt_Mass_EventSum_%s", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), Form("fh3TargetBg_EtaWidth_Pt_Mass_EventSum_%s;#eta;#it{p}_{T}[GeV/c];Mass_{e^{+}e^{-}}[GeV/c^{2}]", MFGV::Name_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[iprofile].c_str()), NBinsEtaWidthReferenceTotal, BinEtaWidthReferenceTotal, MFGV::NBinsPt, MFGV::BinPt, MFGV::NBinsMass, bin_mass_for_flow);
  }

  /////////////////////////////////////////////////////////////////////////////
  fh3Target_Eta_Pt_Mass = new TH3F("fh3Target_Eta_Pt_Mass", "fh3Target_Eta_Pt_Mass", MFGV::NBinsEta, MFGV::BinEta, MFGV::NBinsPt, MFGV::BinPt, 50, bin_mass);
  fh3Target_Eta_Pt_Mass->GetXaxis()->SetTitle("#eta");
  fh3Target_Eta_Pt_Mass->GetYaxis()->SetTitle("#it{p}_{T}[GeV/c]");
  fh3Target_Eta_Pt_Mass->GetZaxis()->SetTitle("Mass [GeV/c^{2}]");

  fh3TargetBg_Eta_Pt_Mass = new TH3F("fh3TargetBg_Eta_Pt_Mass", "fh3TargetBg_Eta_Pt_Mass", MFGV::NBinsEta, MFGV::BinEta, MFGV::NBinsPt, MFGV::BinPt, 50, bin_mass);
  fh3TargetBg_Eta_Pt_Mass->GetXaxis()->SetTitle("#eta");
  fh3TargetBg_Eta_Pt_Mass->GetYaxis()->SetTitle("#it{p}_{T}[GeV/c]");
  fh3TargetBg_Eta_Pt_Mass->GetZaxis()->SetTitle("Mass [GeV/c^{2}]");

  fh1DileptonMass_NoCut = new TH1F("fh1DileptonMass_NoCut", "fh1DileptonMass_NoCut", 200, 2, 4);
  fh1DileptonMass_MultCut = new TH1F("fh1DileptonMass_MultCut", "fh1DileptonMass_MultCut", 200, 2, 4);

  fh1DeltaPhiReference = new TH1F("fh1DeltaPhiReference", "fh1DeltaPhiReference", 100, -2. * TMath::Pi(), 2. * TMath::Pi());
  fh1DeltaPhiTargetReference = new TH1F("fh1DeltaPhiTargetReference", "fh1DeltaPhiTargetReference", 800, -4. * TMath::Pi(), 4. * TMath::Pi());
  fh1DeltaPhiTargetBgReference = new TH1F("fh1DeltaPhiTargetBgReference", "fh1DeltaPhiTargetBgReference", 800, -4. * TMath::Pi(), 4. * TMath::Pi());

  ////////////////////////////////////////////////////////////////

  fListQA->Add(fh1DileptonMass_NoCut);
  fListQA->Add(fh1DileptonMass_MultCut);
  fListQA->Add(fh3EtaPhiMultReference);
  fListQA->Add(fh3EtaPhiMultReferenceCorrected);
  fListQA->Add(fh3EtaPhiMultBeforeCut);
  fListQA->Add(fh1CutTrack);

  ////////////////////////////////////////////////////////////////
  for (int iprofile = 0; iprofile < MFGV::N_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant; iprofile++) {
    fListReference->Add(fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[iprofile][0]);
    fListReference->Add(fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[iprofile][1]);
    fListReference->Add(fh2FlowReferenceTotal_EtaWidth_Mult_WeightSum[iprofile][0]);
    fListReference->Add(fh2FlowReferenceTotal_EtaWidth_Mult_WeightSum[iprofile][1]);
  }

  for (int iprofile = 0; iprofile < MFGV::N_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant; iprofile++) {
    fListReference->Add(fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[iprofile][0]);
    fListReference->Add(fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[iprofile][1]);
    fListReference->Add(fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum[iprofile][0]);
    fListReference->Add(fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum[iprofile][1]);
  }

  ////////////////////////////////////////////////////////////
  for (int iprofile = 0; iprofile < MFGV::N_fProfileTarget_EtaWidth_Pt_cumulant; iprofile++) {
    fListTarget->Add(fProfileTarget_EtaWidth_Pt_cumulant[iprofile][0]);
    fListTarget->Add(fProfileTarget_EtaWidth_Pt_cumulant[iprofile][1]);
    fListTarget->Add(fh2Target_EtaWidth_Pt_WeightSum[iprofile][0]);
    fListTarget->Add(fh2Target_EtaWidth_Pt_WeightSum[iprofile][1]);
  }

  for (int iprofile = 0; iprofile < MFGV::N_fProfileTarget_EtaWidth_Pt_Mass_cumulant; iprofile++) {
    fListTarget->Add(fProfileTarget_EtaWidth_Pt_Mass_cumulant[iprofile][0]);
    fListTarget->Add(fProfileTarget_EtaWidth_Pt_Mass_cumulant[iprofile][1]);
    fListTarget->Add(fh3Target_EtaWidth_Pt_Mass_WeightSum[iprofile][0]);
    fListTarget->Add(fh3Target_EtaWidth_Pt_Mass_WeightSum[iprofile][1]);
  }

  for (int iprofile = 0; iprofile < MFGV::N_fProfileTargetBg_EtaWidth_Pt_cumulant; iprofile++) {
    fListTarget->Add(fProfileTargetBg_EtaWidth_Pt_cumulant[iprofile][0]);
    fListTarget->Add(fProfileTargetBg_EtaWidth_Pt_cumulant[iprofile][1]);
    fListTarget->Add(fh2TargetBg_EtaWidth_Pt_WeightSum[iprofile][0]);
    fListTarget->Add(fh2TargetBg_EtaWidth_Pt_WeightSum[iprofile][1]);
  }

  for (int iprofile = 0; iprofile < MFGV::N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant; iprofile++) {
    fListTarget->Add(fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[iprofile][0]);
    fListTarget->Add(fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[iprofile][1]);
    fListTarget->Add(fh3TargetBg_EtaWidth_Pt_Mass_WeightSum[iprofile][0]);
    fListTarget->Add(fh3TargetBg_EtaWidth_Pt_Mass_WeightSum[iprofile][1]);
  }

  fListTarget->Add(fh3Target_Eta_Pt_Mass);
  fListTarget->Add(fh3TargetBg_Eta_Pt_Mass);

  /////////////////////////////////////////////////////////////

  fListCorrelation->Add(fh1DeltaPhiReference);
  fListCorrelation->Add(fh1DeltaPhiTargetReference);
  fListCorrelation->Add(fh1DeltaPhiTargetBgReference);
}

void MFGF::ArrayFlowInit()
{
  cout << "MFlow Initialing" << endl;
  // fFlowReference init
  for (int i = 0; i < MFGV::NBinsEta; i++) {
    fArrayFlowReferenceEta[i] = MFlow();
  }

  // fFlowTarget
  for (int i = 0; i < MFGV::NBinsEta; i++) {
    for (int k = 0; k < MFGV::NBinsPt; k++) {
      fArrayFlowTargetEtaPt[i][k] = MFlow();
      fArrayFlowTargetBgEtaPt[i][k] = MFlow();
    }
  }

  for (int i = 0; i < MFGV::NBinsEta; i++) {
    for (int k = 0; k < MFGV::NBinsPt; k++) {
      for (int j = 0; j < MFGV::NBinsMass; j++) {
        fArrayFlowTargetEtaPtMass[j][i][k] = MFlow();
        fArrayFlowTargetBgEtaPtMass[j][i][k] = MFlow();
      }
    }
  }

  cout << "MFlow Initialing Done" << endl;
}

void MFGF::ArrayFlowReferenceFill(double eta, double phi, double weight)
{
  fh3EtaPhiMultReference->Fill(eta, phi, fMultiplicity);
  fh3EtaPhiMultReferenceCorrected->Fill(eta, phi, fMultiplicity, weight);
  // cout << "eta: " << eta << " phi: " << phi << " weight: " << weight << endl;
  if (eta < MFGV::BinEta[0] || eta > MFGV::BinEta[MFGV::NBinsEta]) {
    return;
  }
  int etaBin = 0;
  for (int i = 0; i < MFGV::NBinsEta; i++) {
    if (eta < MFGV::BinEta[i + 1]) {
      etaBin = i;
      break;
    }
  }
  fArrayFlowReferenceEta[etaBin].Fill(phi, weight);
}

void MFGF::ArrayFlowTargetFill(double eta, double pt, double phi, double weight)
{
  // fh3Target_Eta_Pt_Mass->Fill(eta, pt, fMass);
  if (eta < MFGV::BinEta[0] || eta > MFGV::BinEta[MFGV::NBinsEta]) {
    return;
  }
  if (pt < MFGV::BinPt[0] || pt > MFGV::BinPt[MFGV::NBinsPt]) {
    return;
  }
  if (MFGF::fMass < MFGV::EdgesMass[0] || MFGF::fMass > MFGV::EdgesMass[1]) {
    return;
  }
  int etaBin = 0;
  int ptBin = 0;
  int massBin = 0;
  for (int i = 0; i < MFGV::NBinsEta; i++) {
    if (eta < MFGV::BinEta[i + 1]) {
      etaBin = i;
      break;
    }
  }
  for (int i = 0; i < MFGV::NBinsPt; i++) {
    if (pt < MFGV::BinPt[i + 1]) {
      ptBin = i;
      break;
    }
  }
  for (int i = 0; i < MFGV::NBinsMass; i++) {
    double edge_right = MFGV::EdgesMass[0] + (i + 1) * (MFGV::EdgesMass[1] - MFGV::EdgesMass[0]) / (double)MFGV::NBinsMass;
    if (fMass < edge_right) {
      massBin = i;
      break;
    }
  }
  fArrayFlowTargetEtaPt[etaBin][ptBin].Fill(phi, weight);
  fArrayFlowTargetEtaPtMass[massBin][etaBin][ptBin].Fill(phi, weight);
}

void MFGF::ArrayFlowTargetBgFill(double eta, double pt, double phi, double weight)
{
  // fh3Target_Eta_Pt_Mass->Fill(eta, pt, fMass);
  if (eta < MFGV::BinEta[0] || eta > MFGV::BinEta[MFGV::NBinsEta]) {
    return;
  }
  if (pt < MFGV::BinPt[0] || pt > MFGV::BinPt[MFGV::NBinsPt]) {
    return;
  }
  if (MFGF::fMass < MFGV::EdgesMass[0] || MFGF::fMass > MFGV::EdgesMass[1]) {
    return;
  }
  int etaBin = 0;
  int ptBin = 0;
  int massBin = 0;
  for (int i = 0; i < MFGV::NBinsEta; i++) {
    if (eta < MFGV::BinEta[i + 1]) {
      etaBin = i;
      break;
    }
  }
  for (int i = 0; i < MFGV::NBinsPt; i++) {
    if (pt < MFGV::BinPt[i + 1]) {
      ptBin = i;
      break;
    }
  }
  for (int i = 0; i < MFGV::NBinsMass; i++) {
    double edge_right = MFGV::EdgesMass[0] + (i + 1) * (MFGV::EdgesMass[1] - MFGV::EdgesMass[0]) / (double)MFGV::NBinsMass;
    if (fMass < edge_right) {
      massBin = i;
      break;
    }
  }
  fArrayFlowTargetBgEtaPt[etaBin][ptBin].Fill(phi, weight);
  fArrayFlowTargetBgEtaPtMass[massBin][etaBin][ptBin].Fill(phi, weight);
}

void MFGF::ArrayFlowReferenceReset()
{
  for (int i = 0; i < MFGV::NBinsEta; i++) {
    fArrayFlowReferenceEta[i].Reset();
  }
}

void MFGF::ArrayFlowTargetReset()
{
  for (int i = 0; i < MFGV::NBinsEta; i++) {
    for (int k = 0; k < MFGV::NBinsPt; k++) {
      fArrayFlowTargetEtaPt[i][k].Reset();
    }
  }

  for (int i = 0; i < MFGV::NBinsEta; i++) {
    for (int k = 0; k < MFGV::NBinsPt; k++) {
      for (int j = 0; j < MFGV::NBinsMass; j++) {
        fArrayFlowTargetEtaPtMass[j][i][k].Reset();
      }
    }
  }
}

void MFGF::ArrayFlowTargetBgReset()
{
  for (int i = 0; i < MFGV::NBinsEta; i++) {
    for (int k = 0; k < MFGV::NBinsPt; k++) {
      fArrayFlowTargetBgEtaPt[i][k].Reset();
    }
  }

  for (int i = 0; i < MFGV::NBinsEta; i++) {
    for (int k = 0; k < MFGV::NBinsPt; k++) {
      for (int j = 0; j < MFGV::NBinsMass; j++) {
        fArrayFlowTargetBgEtaPtMass[j][i][k].Reset();
      }
    }
  }
}

void MFGF::CalculateFlowReference()
{
  // Total calculate
  for (int i_bin_eta_width = 0; i_bin_eta_width < NBinsEtaWidthReferenceTotal; i_bin_eta_width++) {
    double eta_width = BinEtaWidthReferenceTotal[i_bin_eta_width] / 2. + BinEtaWidthReferenceTotal[i_bin_eta_width + 1] / 2.;
    MFlow flow_total;
    flow_total.Reset();
    for (int i_bin_eta = NBinsEtaWidthReferenceTotal - i_bin_eta_width - 1; i_bin_eta <= NBinsEtaWidthReferenceTotal + i_bin_eta_width - NBinsEtaWidthReferenceTotal % 2; i_bin_eta++) {
      flow_total += fArrayFlowReferenceEta[i_bin_eta];
    }

    double weight_c1m1 = (flow_total.Q(0, 1) * flow_total.Q(0, 1) - flow_total.Q(0, 2)).Re();
    double weight2_c1m1 = (flow_total.Q(0, 2) * flow_total.Q(0, 2) - flow_total.Q(0, 4)).Re();

    double c1m1_im = (flow_total.Q(1, 1) * flow_total.Q(-1, 1) - flow_total.Q(0, 2)).Im();
    double c1m1_re = (flow_total.Q(1, 1) * flow_total.Q(-1, 1) - flow_total.Q(0, 2)).Re();
    double c2m2_im = (flow_total.Q(2, 1) * flow_total.Q(-2, 1) - flow_total.Q(0, 2)).Im();
    double c2m2_re = (flow_total.Q(2, 1) * flow_total.Q(-2, 1) - flow_total.Q(0, 2)).Re();

    if (weight_c1m1 > 0.0001) {
      c1m1_im /= weight_c1m1;
      c1m1_re /= weight_c1m1;
      c2m2_im /= weight_c1m1;
      c2m2_re /= weight_c1m1;
      fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[0][0]->Fill(eta_width, fMultiplicity, c1m1_im, weight_c1m1);
      fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[0][1]->Fill(eta_width, fMultiplicity, c1m1_re, weight_c1m1);
      fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[1][0]->Fill(eta_width, fMultiplicity, c2m2_im, weight_c1m1);
      fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[1][1]->Fill(eta_width, fMultiplicity, c2m2_re, weight_c1m1);
      fh2FlowReferenceTotal_EtaWidth_Mult_WeightSum[0][0]->Fill(eta_width, fMultiplicity, weight_c1m1);
      fh2FlowReferenceTotal_EtaWidth_Mult_WeightSum[0][1]->Fill(eta_width, fMultiplicity, weight2_c1m1/weight_c1m1);
      fh2FlowReferenceTotal_EtaWidth_Mult_WeightSum[1][0]->Fill(eta_width, fMultiplicity, weight_c1m1);
      fh2FlowReferenceTotal_EtaWidth_Mult_WeightSum[1][1]->Fill(eta_width, fMultiplicity, weight2_c1m1/weight_c1m1);
    }
  }

  // subevent calculate
  for (int i_bin_delta_eta = 0; i_bin_delta_eta < NBinsDeltaEtaReferenceSubEvent; i_bin_delta_eta++) {
    MFlow flow_subevent_left;
    MFlow flow_subevent_right;
    flow_subevent_left.Reset();
    flow_subevent_right.Reset();
    int index_left_max = NBinsDeltaEtaReferenceSubEvent - i_bin_delta_eta;
    int index_right_min = MFGV::NBinsEta + i_bin_delta_eta - NBinsDeltaEtaReferenceSubEvent;
    for (int i_bin_eta = 0; i_bin_eta < index_left_max; i_bin_eta++) {
      flow_subevent_left += fArrayFlowReferenceEta[i_bin_eta];
    }
    for (int i_bin_eta = index_right_min; i_bin_eta < MFGV::NBinsEta; i_bin_eta++) {
      flow_subevent_right += fArrayFlowReferenceEta[i_bin_eta];
    }

    double delta_eta = BinDeltaEtaReferenceSubEvent[i_bin_delta_eta] / 2. + BinDeltaEtaReferenceSubEvent[i_bin_delta_eta + 1] / 2.;

    double weight_c1_m1 = (flow_subevent_left.Q(0, 1) * flow_subevent_right.Q(0, 1)).Re();
    double weight2_c1_m1 = (flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 2)).Re();
    double c1_m1_im = (flow_subevent_left.Q(1, 1) * flow_subevent_right.Q(-1, 1)).Im();
    double c1_m1_re = (flow_subevent_left.Q(1, 1) * flow_subevent_right.Q(-1, 1)).Re();

    double c2_m2_im = (flow_subevent_left.Q(2, 1) * flow_subevent_right.Q(-2, 1)).Im();
    double c2_m2_re = (flow_subevent_left.Q(2, 1) * flow_subevent_right.Q(-2, 1)).Re();

    double weight_c11_m1m1 = (flow_subevent_left.Q(0, 1) * flow_subevent_right.Q(0, 1) * flow_subevent_left.Q(0, 1) * flow_subevent_right.Q(0, 1) - flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 1) * flow_subevent_right.Q(0, 1) - flow_subevent_left.Q(0, 1) * flow_subevent_left.Q(0, 1) * flow_subevent_right.Q(0, 2) + flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 2)).Re();
    double weight2_c11_m1m1 = (flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 2) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 2) - flow_subevent_left.Q(0, 4) * flow_subevent_right.Q(0, 2) * flow_subevent_right.Q(0, 2) - flow_subevent_left.Q(0, 2) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 4) + flow_subevent_left.Q(0, 4) * flow_subevent_right.Q(0, 4)).Re();
    double c22_m2m2_im = (flow_subevent_left.Q(2, 1) * flow_subevent_right.Q(-2, 1) * flow_subevent_left.Q(2, 1) * flow_subevent_right.Q(-2, 1) - flow_subevent_left.Q(4, 2) * flow_subevent_right.Q(-2, 1) * flow_subevent_right.Q(-2, 1) - flow_subevent_left.Q(2, 1) * flow_subevent_left.Q(2, 1) * flow_subevent_right.Q(-4, 2) + flow_subevent_left.Q(4, 2) * flow_subevent_right.Q(-4, 2)).Im();
    double c22_m2m2_re = (flow_subevent_left.Q(2, 1) * flow_subevent_right.Q(-2, 1) * flow_subevent_left.Q(2, 1) * flow_subevent_right.Q(-2, 1) - flow_subevent_left.Q(4, 2) * flow_subevent_right.Q(-2, 1) * flow_subevent_right.Q(-2, 1) - flow_subevent_left.Q(2, 1) * flow_subevent_left.Q(2, 1) * flow_subevent_right.Q(-4, 2) + flow_subevent_left.Q(4, 2) * flow_subevent_right.Q(-4, 2)).Re();

    if (weight_c1_m1 > 0.0001) {
      c1_m1_im /= weight_c1_m1;
      c1_m1_re /= weight_c1_m1;
      c2_m2_im /= weight_c1_m1;
      c2_m2_re /= weight_c1_m1;
      fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[0][0]->Fill(delta_eta, fMultiplicity, c1_m1_im, weight_c1_m1);
      fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[0][1]->Fill(delta_eta, fMultiplicity, c1_m1_re, weight_c1_m1);
      fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[1][0]->Fill(delta_eta, fMultiplicity, c2_m2_im, weight_c1_m1);
      fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[1][1]->Fill(delta_eta, fMultiplicity, c2_m2_re, weight_c1_m1);
      fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum[0][0]->Fill(delta_eta, fMultiplicity, weight_c1_m1);
      fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum[0][1]->Fill(delta_eta, fMultiplicity, weight2_c1_m1/weight_c1_m1);
      fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum[1][0]->Fill(delta_eta, fMultiplicity, weight_c1_m1);
      fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum[1][1]->Fill(delta_eta, fMultiplicity, weight2_c1_m1/weight_c1_m1);
    }

    if (weight_c11_m1m1 > 0.0001) {
      c22_m2m2_im /= weight_c11_m1m1;
      c22_m2m2_re /= weight_c11_m1m1;
      fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[2][0]->Fill(delta_eta, fMultiplicity, c22_m2m2_im, weight_c11_m1m1);
      fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[2][1]->Fill(delta_eta, fMultiplicity, c22_m2m2_re, weight_c11_m1m1);
      fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum[2][0]->Fill(delta_eta, fMultiplicity, weight_c11_m1m1);
      fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum[2][1]->Fill(delta_eta, fMultiplicity, weight2_c11_m1m1/weight_c11_m1m1);
    }
  }
};

void MFGF::CalculateFlowTarget()
{
  for (int i_bin_delta_eta = 0; i_bin_delta_eta < NBinsDeltaEtaReferenceSubEvent; i_bin_delta_eta++) {
    MFlow flow_subevent_left;
    MFlow flow_subevent_right;
    flow_subevent_left.Reset();
    flow_subevent_right.Reset();
    int index_left_max = NBinsDeltaEtaReferenceSubEvent - i_bin_delta_eta;
    int index_right_min = MFGV::NBinsEta + i_bin_delta_eta - NBinsDeltaEtaReferenceSubEvent;
    for (int i_bin_eta = 0; i_bin_eta < index_left_max; i_bin_eta++) {
      flow_subevent_left += fArrayFlowReferenceEta[i_bin_eta];
    }
    for (int i_bin_eta = index_right_min; i_bin_eta < MFGV::NBinsEta; i_bin_eta++) {
      flow_subevent_right += fArrayFlowReferenceEta[i_bin_eta];
    }

    double delta_eta = BinDeltaEtaReferenceSubEvent[i_bin_delta_eta] / 2. + BinDeltaEtaReferenceSubEvent[i_bin_delta_eta + 1] / 2.;

    // Target Flow
    if (index_left_max == index_right_min)
      continue;

    MFlow flow_Target[MFGV::NBinsPt] = {MFlow()};
    MFlow flow_TargetBg[MFGV::NBinsPt] = {MFlow()};

    // reset
    for (int i_bin_pt = 0; i_bin_pt < MFGV::NBinsPt; i_bin_pt++) {
      flow_Target[i_bin_pt].Reset();
      flow_TargetBg[i_bin_pt].Reset();
    }

    for (int i_bin_eta = index_left_max; i_bin_eta < index_right_min; i_bin_eta++) {
      for (int i_bin_pt = 0; i_bin_pt < MFGV::NBinsPt; i_bin_pt++) {
        flow_Target[i_bin_pt] += fArrayFlowTargetEtaPt[i_bin_eta][i_bin_pt];
        flow_TargetBg[i_bin_pt] += fArrayFlowTargetBgEtaPt[i_bin_eta][i_bin_pt];
      }
    }

    for (int i_bin_pt = 0; i_bin_pt < MFGV::NBinsPt; i_bin_pt++) {
      double weight_c2_m1_m1 = (flow_Target[i_bin_pt].Q(0, 1) * flow_subevent_left.Q(0, 1) * flow_subevent_right.Q(0, 1)).Re();
      double weight2_c2_m1_m1 = (flow_Target[i_bin_pt].Q(0, 2) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 2)).Re();
      double c2_m1_m1_im = (flow_Target[i_bin_pt].Q(2, 1) * flow_subevent_left.Q(-1, 1) * flow_subevent_right.Q(-1, 1)).Im();
      double c2_m1_m1_re = (flow_Target[i_bin_pt].Q(2, 1) * flow_subevent_left.Q(-1, 1) * flow_subevent_right.Q(-1, 1)).Re();

      double weight_c1_1_m1_m1 = (flow_Target[i_bin_pt].Q(0, 1) * (flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1)) * flow_subevent_left.Q(0, 1) * flow_subevent_right.Q(0, 1) - flow_Target[i_bin_pt].Q(0, 1) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 1) - flow_Target[i_bin_pt].Q(0, 1) * flow_subevent_right.Q(0, 2) * flow_subevent_left.Q(0, 1)).Re();
      double weight2_c1_1_m1_m1 = (flow_Target[i_bin_pt].Q(0, 2) * (flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2)) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 2) - flow_Target[i_bin_pt].Q(0, 2) * flow_subevent_left.Q(0, 4) * flow_subevent_right.Q(0, 2) - flow_Target[i_bin_pt].Q(0, 2) * flow_subevent_right.Q(0, 4) * flow_subevent_left.Q(0, 2)).Re();
      double c2_2_m2_m2_im = (flow_Target[i_bin_pt].Q(2, 1) * (flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_left.Q(-2, 1) * flow_subevent_right.Q(-2, 1) - flow_Target[i_bin_pt].Q(2, 1) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(-2, 1) - flow_Target[i_bin_pt].Q(2, 1) * flow_subevent_right.Q(0, 2) * flow_subevent_left.Q(-2, 1)).Im();
      double c2_2_m2_m2_re = (flow_Target[i_bin_pt].Q(2, 1) * (flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_left.Q(-2, 1) * flow_subevent_right.Q(-2, 1) - flow_Target[i_bin_pt].Q(2, 1) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(-2, 1) - flow_Target[i_bin_pt].Q(2, 1) * flow_subevent_right.Q(0, 2) * flow_subevent_left.Q(-2, 1)).Re();

      double weight_c1_0_m1_0 = (flow_Target[i_bin_pt].Q(0, 1) * flow_subevent_left.Q(0, 1)).Re();
      double weight2_c1_0_m1_0 = (flow_Target[i_bin_pt].Q(0, 2) * flow_subevent_left.Q(0, 2)).Re();
      double c2_0_m2_0_im = (flow_Target[i_bin_pt].Q(2, 1) * flow_subevent_left.Q(-2, 1)).Im();
      double c2_0_m2_0_re = (flow_Target[i_bin_pt].Q(2, 1) * flow_subevent_left.Q(-2, 1)).Re();

      double weight_c0_1_0_m1 = ((flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1)) * flow_subevent_right.Q(0, 1) - flow_subevent_right.Q(0, 2)).Re();
      double weight2_c0_1_0_m1 = ((flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2)) * flow_subevent_right.Q(0, 2) - flow_subevent_right.Q(0, 4)).Re();
      double c0_2_0_m2_im = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_right.Q(-2, 1) - flow_subevent_right.Q(0, 2)).Im();
      double c0_2_0_m2_re = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_right.Q(-2, 1) - flow_subevent_right.Q(0, 2)).Re();

      double weight_c1_0_0_m1 = (flow_Target[i_bin_pt].Q(0, 1) * flow_subevent_right.Q(0, 1)).Re();
      double weight2_c1_0_0_m1 = (flow_Target[i_bin_pt].Q(0, 2) * flow_subevent_right.Q(0, 2)).Re();
      double c2_0_0_m2_im = (flow_Target[i_bin_pt].Q(2, 1) * flow_subevent_right.Q(-2, 1)).Im();
      double c2_0_0_m2_re = (flow_Target[i_bin_pt].Q(2, 1) * flow_subevent_right.Q(-2, 1)).Re();

      double weight_c0_1_m1_0 = ((flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1)) * flow_subevent_left.Q(0, 1) - flow_subevent_left.Q(0, 2)).Re();
      double weight2_c0_1_m1_0 = ((flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2)) * flow_subevent_left.Q(0, 2) - flow_subevent_left.Q(0, 4)).Re();
      double c0_2_m2_0_im = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_right.Q(-2, 1) - flow_subevent_right.Q(0, 2)).Im();
      double c0_2_m2_0_re = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_right.Q(-2, 1) - flow_subevent_right.Q(0, 2)).Re();

      double weight_c1_m1 = (flow_Target[i_bin_pt].Q(0, 1) * (flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1))).Re();
      double weight2_c1_m1 = (flow_Target[i_bin_pt].Q(0, 2) * (flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2))).Re();
      double c2_m2_im = (flow_Target[i_bin_pt].Q(2, 1) * (flow_subevent_left.Q(-2, 1) + flow_subevent_right.Q(-2, 1))).Im();
      double c2_m2_re = (flow_Target[i_bin_pt].Q(2, 1) * (flow_subevent_left.Q(-2, 1) + flow_subevent_right.Q(-2, 1))).Re();

      double weights[MFGV::N_fProfileTarget_EtaWidth_Pt_cumulant] = {weight_c2_m1_m1, weight_c1_1_m1_m1, weight_c1_0_m1_0, weight_c0_1_0_m1, weight_c1_0_0_m1, weight_c0_1_m1_0, weight_c1_m1};
      double weight2s[MFGV::N_fProfileTarget_EtaWidth_Pt_cumulant] = {weight2_c2_m1_m1, weight2_c1_1_m1_m1, weight2_c1_0_m1_0, weight2_c0_1_0_m1, weight2_c1_0_0_m1, weight2_c0_1_m1_0, weight2_c1_m1};
      double ims[MFGV::N_fProfileTarget_EtaWidth_Pt_cumulant] = {c2_m1_m1_im, c2_2_m2_m2_im, c2_0_m2_0_im, c0_2_0_m2_im, c2_0_0_m2_im, c0_2_m2_0_im, c2_m2_im};
      double res[MFGV::N_fProfileTarget_EtaWidth_Pt_cumulant] = {c2_m1_m1_re, c2_2_m2_m2_re, c2_0_m2_0_re, c0_2_0_m2_re, c2_0_0_m2_re, c0_2_m2_0_re, c2_m2_re};

      for (int i = 0; i < MFGV::N_fProfileTarget_EtaWidth_Pt_cumulant; i++) {
        if (weights[i] > 0.0001) {
          ims[i] /= weights[i];
          res[i] /= weights[i];
          fProfileTarget_EtaWidth_Pt_cumulant[i][0]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], ims[i], weights[i]);
          fProfileTarget_EtaWidth_Pt_cumulant[i][1]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], res[i], weights[i]);
          fh2Target_EtaWidth_Pt_WeightSum[i][0]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], weights[i]);
          fh2Target_EtaWidth_Pt_WeightSum[i][1]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], weight2s[i]/weights[i]);
        }
      }

      for (int imass = 0; imass < MFGV::NBinsMass; imass++) {
        MFlow flow_Target_mass = MFlow();
        double mass = MFGV::EdgesMass[0] + imass * (MFGV::EdgesMass[1] - MFGV::EdgesMass[0]) / (double)MFGV::NBinsMass;
        flow_Target_mass.Reset();

        for (int j_bin_eta = index_left_max; j_bin_eta < index_right_min; j_bin_eta++) {
          flow_Target_mass += fArrayFlowTargetEtaPtMass[imass][j_bin_eta][i_bin_pt];
        }

        double weight_c2_m1_m1_mass = (flow_Target_mass.Q(0, 1) * flow_subevent_left.Q(0, 1) * flow_subevent_right.Q(0, 1)).Re();
        double weight2_c2_m1_m1_mass = (flow_Target_mass.Q(0, 2) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 2)).Re();
        double c2_m1_m1_im_mass = (flow_Target_mass.Q(2, 1) * flow_subevent_left.Q(-1, 1) * flow_subevent_right.Q(-1, 1)).Im();
        double c2_m1_m1_re_mass = (flow_Target_mass.Q(2, 1) * flow_subevent_left.Q(-1, 1) * flow_subevent_right.Q(-1, 1)).Re();

        double weight_c1_1_m1_m1_mass = (flow_Target_mass.Q(0, 1) * (flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1)) * flow_subevent_left.Q(0, 1) * flow_subevent_right.Q(0, 1) - flow_Target_mass.Q(0, 1) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 1) - flow_Target_mass.Q(0, 1) * flow_subevent_right.Q(0, 2) * flow_subevent_left.Q(0, 1)).Re();
        double weight2_c1_1_m1_m1_mass = (flow_Target_mass.Q(0, 2) * (flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2)) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 2) - flow_Target_mass.Q(0, 2) * flow_subevent_left.Q(0, 4) * flow_subevent_right.Q(0, 2) - flow_Target_mass.Q(0, 2) * flow_subevent_right.Q(0, 4) * flow_subevent_left.Q(0, 2)).Re();
        double c2_2_m2_m2_im_mass = (flow_Target_mass.Q(2, 1) * (flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_left.Q(-2, 1) * flow_subevent_right.Q(-2, 1) - flow_Target_mass.Q(2, 1) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(-2, 1) - flow_Target_mass.Q(2, 1) * flow_subevent_right.Q(0, 2) * flow_subevent_left.Q(-2, 1)).Im();
        double c2_2_m2_m2_re_mass = (flow_Target_mass.Q(2, 1) * (flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_left.Q(-2, 1) * flow_subevent_right.Q(-2, 1) - flow_Target_mass.Q(2, 1) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(-2, 1) - flow_Target_mass.Q(2, 1) * flow_subevent_right.Q(0, 2) * flow_subevent_left.Q(-2, 1)).Re();

        double weight_c1_0_m1_0_mass = (flow_Target_mass.Q(0, 1) * flow_subevent_left.Q(0, 1)).Re();
        double weight2_c1_0_m1_0_mass = (flow_Target_mass.Q(0, 2) * flow_subevent_left.Q(0, 2)).Re();
        double c2_0_m2_0_im_mass = (flow_Target_mass.Q(2, 1) * flow_subevent_left.Q(-2, 1)).Im();
        double c2_0_m2_0_re_mass = (flow_Target_mass.Q(2, 1) * flow_subevent_left.Q(-2, 1)).Re();

        double weight_c0_1_0_m1_mass = ((flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1)) * flow_subevent_right.Q(0, 1) - flow_subevent_right.Q(0, 2)).Re();
        double weight2_c0_1_0_m1_mass = ((flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2)) * flow_subevent_right.Q(0, 2) - flow_subevent_right.Q(0, 4)).Re();
        double c0_2_0_m2_im_mass = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_right.Q(-2, 1) - flow_subevent_right.Q(0, 2)).Im();
        double c0_2_0_m2_re_mass = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_right.Q(-2, 1) - flow_subevent_right.Q(0, 2)).Re();

        double weight_c1_0_0_m1_mass = (flow_Target_mass.Q(0, 1) * flow_subevent_right.Q(0, 1)).Re();
        double weight2_c1_0_0_m1_mass = (flow_Target_mass.Q(0, 2) * flow_subevent_right.Q(0, 2)).Re();
        double c2_0_0_m2_im_mass = (flow_Target_mass.Q(2, 1) * flow_subevent_right.Q(-2, 1)).Im();
        double c2_0_0_m2_re_mass = (flow_Target_mass.Q(2, 1) * flow_subevent_right.Q(-2, 1)).Re();

        double weight_c0_1_m1_0_mass = ((flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1)) * flow_subevent_left.Q(0, 1) - flow_subevent_left.Q(0, 2)).Re();
        double weight2_c0_1_m1_0_mass = ((flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2)) * flow_subevent_left.Q(0, 2) - flow_subevent_left.Q(0, 4)).Re();
        double c0_2_m2_0_im_mass = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_left.Q(-2, 1) - flow_subevent_left.Q(0, 2)).Im();
        double c0_2_m2_0_re_mass = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_left.Q(-2, 1) - flow_subevent_left.Q(0, 2)).Re();

        double weight_c1_m1_mass = (flow_Target_mass.Q(0, 1) * (flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1))).Re();
        double weight2_c1_m1_mass = (flow_Target_mass.Q(0, 2) * (flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2))).Re();
        double c2_m2_im_mass = (flow_Target_mass.Q(2, 1) * (flow_subevent_left.Q(-2, 1) + flow_subevent_right.Q(-2, 1))).Im();
        double c2_m2_re_mass = (flow_Target_mass.Q(2, 1) * (flow_subevent_left.Q(-2, 1) + flow_subevent_right.Q(-2, 1))).Re();

        double weights_mass[MFGV::N_fProfileTarget_EtaWidth_Pt_Mass_cumulant] = {weight_c2_m1_m1_mass, weight_c1_1_m1_m1_mass, weight_c1_0_m1_0_mass, weight_c0_1_0_m1_mass, weight_c1_0_0_m1_mass, weight_c0_1_m1_0_mass, weight_c1_m1_mass};
        double weight2s_mass[MFGV::N_fProfileTarget_EtaWidth_Pt_Mass_cumulant] = {weight2_c2_m1_m1_mass, weight2_c1_1_m1_m1_mass, weight2_c1_0_m1_0_mass, weight2_c0_1_0_m1_mass, weight2_c1_0_0_m1_mass, weight2_c0_1_m1_0_mass, weight2_c1_m1_mass};
        double ims_mass[MFGV::N_fProfileTarget_EtaWidth_Pt_Mass_cumulant] = {c2_m1_m1_im_mass, c2_2_m2_m2_im_mass, c2_0_m2_0_im_mass, c0_2_0_m2_im_mass, c2_0_0_m2_im_mass, c0_2_m2_0_im_mass, c2_m2_im_mass};
        double res_mass[MFGV::N_fProfileTarget_EtaWidth_Pt_Mass_cumulant] = {c2_m1_m1_re_mass, c2_2_m2_m2_re_mass, c2_0_m2_0_re_mass, c0_2_0_m2_re_mass, c2_0_0_m2_re_mass, c0_2_m2_0_re_mass, c2_m2_re_mass};

        for (int i = 0; i < MFGV::N_fProfileTarget_EtaWidth_Pt_Mass_cumulant; i++) {
          if (weights_mass[i] > 0.0001) {
            ims_mass[i] /= weights_mass[i];
            res_mass[i] /= weights_mass[i];
            fProfileTarget_EtaWidth_Pt_Mass_cumulant[i][0]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], mass, ims_mass[i], weights_mass[i]);
            fProfileTarget_EtaWidth_Pt_Mass_cumulant[i][1]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], mass, res_mass[i], weights_mass[i]);
            fh3Target_EtaWidth_Pt_Mass_WeightSum[i][0]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], mass, weights_mass[i]);
            fh3Target_EtaWidth_Pt_Mass_WeightSum[i][1]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], mass, weight2s_mass[i]/weights_mass[i]);
          }
        }
      }
    }

    for (int i_bin_pt = 0; i_bin_pt < MFGV::NBinsPt; i_bin_pt++) {
      double weight_c2_m1_m1 = (flow_TargetBg[i_bin_pt].Q(0, 1) * flow_subevent_left.Q(0, 1) * flow_subevent_right.Q(0, 1)).Re();
      double weight2_c2_m1_m1 = (flow_TargetBg[i_bin_pt].Q(0, 2) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 2)).Re();
      double c2_m1_m1_im = (flow_TargetBg[i_bin_pt].Q(2, 1) * flow_subevent_left.Q(-1, 1) * flow_subevent_right.Q(-1, 1)).Im();
      double c2_m1_m1_re = (flow_TargetBg[i_bin_pt].Q(2, 1) * flow_subevent_left.Q(-1, 1) * flow_subevent_right.Q(-1, 1)).Re();

      double weight_c1_1_m1_m1 = (flow_TargetBg[i_bin_pt].Q(0, 1) * (flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1)) * flow_subevent_left.Q(0, 1) * flow_subevent_right.Q(0, 1) - flow_TargetBg[i_bin_pt].Q(0, 1) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 1) - flow_TargetBg[i_bin_pt].Q(0, 1) * flow_subevent_right.Q(0, 2) * flow_subevent_left.Q(0, 1)).Re();
      double weight2_c1_1_m1_m1 = (flow_TargetBg[i_bin_pt].Q(0, 2) * (flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2)) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 2) - flow_TargetBg[i_bin_pt].Q(0, 2) * flow_subevent_left.Q(0, 4) * flow_subevent_right.Q(0, 2) - flow_TargetBg[i_bin_pt].Q(0, 2) * flow_subevent_right.Q(0, 4) * flow_subevent_left.Q(0, 2)).Re();
      double c2_2_m2_m2_im = (flow_TargetBg[i_bin_pt].Q(2, 1) * (flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_left.Q(-2, 1) * flow_subevent_right.Q(-2, 1) - flow_TargetBg[i_bin_pt].Q(2, 1) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(-2, 1) - flow_TargetBg[i_bin_pt].Q(2, 1) * flow_subevent_right.Q(0, 2) * flow_subevent_left.Q(-2, 1)).Im();
      double c2_2_m2_m2_re = (flow_TargetBg[i_bin_pt].Q(2, 1) * (flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_left.Q(-2, 1) * flow_subevent_right.Q(-2, 1) - flow_TargetBg[i_bin_pt].Q(2, 1) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(-2, 1) - flow_TargetBg[i_bin_pt].Q(2, 1) * flow_subevent_right.Q(0, 2) * flow_subevent_left.Q(-2, 1)).Re();

      double weight_c1_0_m1_0 = (flow_TargetBg[i_bin_pt].Q(0, 1) * flow_subevent_left.Q(0, 1)).Re();
      double weight2_c1_0_m1_0 = (flow_TargetBg[i_bin_pt].Q(0, 2) * flow_subevent_left.Q(0, 2)).Re();
      double c2_0_m2_0_im = (flow_TargetBg[i_bin_pt].Q(2, 1) * flow_subevent_left.Q(-2, 1)).Im();
      double c2_0_m2_0_re = (flow_TargetBg[i_bin_pt].Q(2, 1) * flow_subevent_left.Q(-2, 1)).Re();

      double weight_c0_1_0_m1 = ((flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1)) * flow_subevent_right.Q(0, 1) - flow_subevent_right.Q(0, 2)).Re();
      double weight2_c0_1_0_m1 = ((flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2)) * flow_subevent_right.Q(0, 2) - flow_subevent_right.Q(0, 4)).Re();
      double c0_2_0_m2_im = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_right.Q(-2, 1) - flow_subevent_right.Q(0, 2)).Im();
      double c0_2_0_m2_re = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_right.Q(-2, 1) - flow_subevent_right.Q(0, 2)).Re();

      double weight_c1_0_0_m1 = (flow_TargetBg[i_bin_pt].Q(0, 1) * flow_subevent_right.Q(0, 1)).Re();
      double weight2_c1_0_0_m1 = (flow_TargetBg[i_bin_pt].Q(0, 2) * flow_subevent_right.Q(0, 2)).Re();
      double c2_0_0_m2_im = (flow_TargetBg[i_bin_pt].Q(2, 1) * flow_subevent_right.Q(-2, 1)).Im();
      double c2_0_0_m2_re = (flow_TargetBg[i_bin_pt].Q(2, 1) * flow_subevent_right.Q(-2, 1)).Re();

      double weight_c0_1_m1_0 = ((flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1)) * flow_subevent_left.Q(0, 1) - flow_subevent_left.Q(0, 2)).Re();
      double weight2_c0_1_m1_0 = ((flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2)) * flow_subevent_left.Q(0, 2) - flow_subevent_left.Q(0, 4)).Re();
      double c0_2_m2_0_im = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_right.Q(-2, 1) - flow_subevent_right.Q(0, 2)).Im();
      double c0_2_m2_0_re = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_right.Q(-2, 1) - flow_subevent_right.Q(0, 2)).Re();

      double weight_c1_m1 = (flow_Target[i_bin_pt].Q(0, 1) * (flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1))).Re();
      double weight2_c1_m1 = (flow_Target[i_bin_pt].Q(0, 2) * (flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2))).Re();
      double c2_m2_im = (flow_TargetBg[i_bin_pt].Q(2, 1) * (flow_subevent_left.Q(-2, 1) + flow_subevent_right.Q(-2, 1))).Im();
      double c2_m2_re = (flow_TargetBg[i_bin_pt].Q(2, 1) * (flow_subevent_left.Q(-2, 1) + flow_subevent_right.Q(-2, 1))).Re();

      double weights[MFGV::N_fProfileTargetBg_EtaWidth_Pt_cumulant] = {weight_c2_m1_m1, weight_c1_1_m1_m1, weight_c1_0_m1_0, weight_c0_1_0_m1, weight_c1_0_0_m1, weight_c0_1_m1_0, weight_c1_m1};
      double weight2s[MFGV::N_fProfileTargetBg_EtaWidth_Pt_cumulant] = {weight2_c2_m1_m1, weight2_c1_1_m1_m1, weight2_c1_0_m1_0, weight2_c0_1_0_m1, weight2_c1_0_0_m1, weight2_c0_1_m1_0, weight2_c1_m1};
      double ims[MFGV::N_fProfileTargetBg_EtaWidth_Pt_cumulant] = {c2_m1_m1_im, c2_2_m2_m2_im, c2_0_m2_0_im, c0_2_0_m2_im, c2_0_0_m2_im, c0_2_m2_0_im, c2_m2_im};
      double res[MFGV::N_fProfileTargetBg_EtaWidth_Pt_cumulant] = {c2_m1_m1_re, c2_2_m2_m2_re, c2_0_m2_0_re, c0_2_0_m2_re, c2_0_0_m2_re, c0_2_m2_0_re, c2_m2_re};

      for (int i = 0; i < MFGV::N_fProfileTargetBg_EtaWidth_Pt_cumulant; i++) {
        if (weights[i] > 0.0001) {
          ims[i] /= weights[i];
          res[i] /= weights[i];
          fProfileTargetBg_EtaWidth_Pt_cumulant[i][0]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], ims[i], weights[i]);
          fProfileTargetBg_EtaWidth_Pt_cumulant[i][1]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], res[i], weights[i]);
          fh2TargetBg_EtaWidth_Pt_WeightSum[i][0]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], weights[i]);
          fh2TargetBg_EtaWidth_Pt_WeightSum[i][1]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], weight2s[i]/weights[i]);
        }
      }

      for (int imass = 0; imass < MFGV::NBinsMass; imass++) {
        MFlow flow_TargetBg_mass = MFlow();
        double mass = MFGV::EdgesMass[0] + imass * (MFGV::EdgesMass[1] - MFGV::EdgesMass[0]) / (double)MFGV::NBinsMass;
        flow_TargetBg_mass.Reset();

        for (int j_bin_eta = index_left_max; j_bin_eta < index_right_min; j_bin_eta++) {
          flow_TargetBg_mass += fArrayFlowTargetBgEtaPtMass[imass][j_bin_eta][i_bin_pt];
        }

        double weight_c2_m1_m1_mass = (flow_TargetBg_mass.Q(0, 1) * flow_subevent_left.Q(0, 1) * flow_subevent_right.Q(0, 1)).Re();
        double weight2_c2_m1_m1_mass = (flow_TargetBg_mass.Q(0, 2) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 2)).Re();
        double c2_m1_m1_im_mass = (flow_TargetBg_mass.Q(2, 1) * flow_subevent_left.Q(-1, 1) * flow_subevent_right.Q(-1, 1)).Im();
        double c2_m1_m1_re_mass = (flow_TargetBg_mass.Q(2, 1) * flow_subevent_left.Q(-1, 1) * flow_subevent_right.Q(-1, 1)).Re();

        double weight_c1_1_m1_m1_mass = (flow_TargetBg_mass.Q(0, 1) * (flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1)) * flow_subevent_left.Q(0, 1) * flow_subevent_right.Q(0, 1) - flow_TargetBg_mass.Q(0, 1) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 1) - flow_TargetBg_mass.Q(0, 1) * flow_subevent_right.Q(0, 2) * flow_subevent_left.Q(0, 1)).Re();
        double weight2_c1_1_m1_m1_mass = (flow_TargetBg_mass.Q(0, 2) * (flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2)) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(0, 2) - flow_TargetBg_mass.Q(0, 2) * flow_subevent_left.Q(0, 4) * flow_subevent_right.Q(0, 2) - flow_TargetBg_mass.Q(0, 2) * flow_subevent_right.Q(0, 4) * flow_subevent_left.Q(0, 2)).Re();
        double c2_2_m2_m2_im_mass = (flow_TargetBg_mass.Q(2, 1) * (flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_left.Q(-2, 1) * flow_subevent_right.Q(-2, 1) - flow_TargetBg_mass.Q(2, 1) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(-2, 1) - flow_TargetBg_mass.Q(2, 1) * flow_subevent_right.Q(0, 2) * flow_subevent_left.Q(-2, 1)).Im();
        double c2_2_m2_m2_re_mass = (flow_TargetBg_mass.Q(2, 1) * (flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_left.Q(-2, 1) * flow_subevent_right.Q(-2, 1) - flow_TargetBg_mass.Q(2, 1) * flow_subevent_left.Q(0, 2) * flow_subevent_right.Q(-2, 1) - flow_TargetBg_mass.Q(2, 1) * flow_subevent_right.Q(0, 2) * flow_subevent_left.Q(-2, 1)).Re();

        double weight_c1_0_m1_0_mass = (flow_TargetBg_mass.Q(0, 1) * flow_subevent_left.Q(0, 1)).Re();
        double weight2_c1_0_m1_0_mass = (flow_TargetBg_mass.Q(0, 2) * flow_subevent_left.Q(0, 2)).Re();
        double c2_0_m2_0_im_mass = (flow_TargetBg_mass.Q(2, 1) * flow_subevent_left.Q(-2, 1)).Im();
        double c2_0_m2_0_re_mass = (flow_TargetBg_mass.Q(2, 1) * flow_subevent_left.Q(-2, 1)).Re();

        double weight_c0_1_0_m1_mass = ((flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1)) * flow_subevent_right.Q(0, 1) - flow_subevent_right.Q(0, 2)).Re();
        double weight2_c0_1_0_m1_mass = ((flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2)) * flow_subevent_right.Q(0, 2) - flow_subevent_right.Q(0, 4)).Re();
        double c0_2_0_m2_im_mass = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_right.Q(-2, 1) - flow_subevent_right.Q(0, 2)).Im();
        double c0_2_0_m2_re_mass = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_right.Q(-2, 1) - flow_subevent_right.Q(0, 2)).Re();

        double weight_c1_0_0_m1_mass = (flow_TargetBg_mass.Q(0, 1) * flow_subevent_right.Q(0, 1)).Re();
        double weight2_c1_0_0_m1_mass = (flow_TargetBg_mass.Q(0, 2) * flow_subevent_right.Q(0, 2)).Re();
        double c2_0_0_m2_im_mass = (flow_TargetBg_mass.Q(2, 1) * flow_subevent_right.Q(-2, 1)).Im();
        double c2_0_0_m2_re_mass = (flow_TargetBg_mass.Q(2, 1) * flow_subevent_right.Q(-2, 1)).Re();

        double weight_c0_1_m1_0_mass = ((flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1)) * flow_subevent_left.Q(0, 1) - flow_subevent_left.Q(0, 2)).Re();
        double weight2_c0_1_m1_0_mass = ((flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2)) * flow_subevent_left.Q(0, 2) - flow_subevent_left.Q(0, 4)).Re();
        double c0_2_m2_0_im_mass = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_left.Q(-2, 1) - flow_subevent_left.Q(0, 2)).Im();
        double c0_2_m2_0_re_mass = ((flow_subevent_left.Q(2, 1) + flow_subevent_right.Q(2, 1)) * flow_subevent_left.Q(-2, 1) - flow_subevent_left.Q(0, 2)).Re();

        double weight_c1_m1_mass = (flow_TargetBg_mass.Q(0, 1) * (flow_subevent_left.Q(0, 1) + flow_subevent_right.Q(0, 1))).Re();
        double weight2_c1_m1_mass = (flow_TargetBg_mass.Q(0, 2) * (flow_subevent_left.Q(0, 2) + flow_subevent_right.Q(0, 2))).Re();
        double c2_m2_im_mass = (flow_TargetBg_mass.Q(2, 1) * (flow_subevent_left.Q(-2, 1) + flow_subevent_right.Q(-2, 1))).Im();
        double c2_m2_re_mass = (flow_TargetBg_mass.Q(2, 1) * (flow_subevent_left.Q(-2, 1) + flow_subevent_right.Q(-2, 1))).Re();

        double weights_mass[MFGV::N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant] = {weight_c2_m1_m1_mass, weight_c1_1_m1_m1_mass, weight_c1_0_m1_0_mass, weight_c0_1_0_m1_mass, weight_c1_0_0_m1_mass, weight_c0_1_m1_0_mass, weight_c1_m1_mass};
        double weight2s_mass[MFGV::N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant] = {weight2_c2_m1_m1_mass, weight2_c1_1_m1_m1_mass, weight2_c1_0_m1_0_mass, weight2_c0_1_0_m1_mass, weight2_c1_0_0_m1_mass, weight2_c0_1_m1_0_mass, weight2_c1_m1_mass};
        double ims_mass[MFGV::N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant] = {c2_m1_m1_im_mass, c2_2_m2_m2_im_mass, c2_0_m2_0_im_mass, c0_2_0_m2_im_mass, c2_0_0_m2_im_mass, c0_2_m2_0_im_mass, c2_m2_im_mass};
        double res_mass[MFGV::N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant] = {c2_m1_m1_re_mass, c2_2_m2_m2_re_mass, c2_0_m2_0_re_mass, c0_2_0_m2_re_mass, c2_0_0_m2_re_mass, c0_2_m2_0_re_mass, c2_m2_re_mass};

        for (int i = 0; i < MFGV::N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant; i++) {
          if (weights_mass[i] > 0.0001) {
            ims_mass[i] /= weights_mass[i];
            res_mass[i] /= weights_mass[i];
            fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[i][0]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], mass, ims_mass[i], weights_mass[i]);
            fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[i][1]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], mass, res_mass[i], weights_mass[i]);
            fh3TargetBg_EtaWidth_Pt_Mass_WeightSum[i][0]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], mass, weights_mass[i]);
            fh3TargetBg_EtaWidth_Pt_Mass_WeightSum[i][1]->Fill(delta_eta, MFGV::BinPt[i_bin_pt], mass, weight2s_mass[i]/weights_mass[i]);
          }
        }
      }
    }
  }
};
