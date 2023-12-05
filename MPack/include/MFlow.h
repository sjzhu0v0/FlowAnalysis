#ifndef MFLOW_H
#define MFLOW_H

#include "TComplex.h"
#include "iostream"
#include "TMath.h"
#include "TString.h"
#include "TList.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "string"

using namespace std;

namespace MFlowGlobalVariables
{
const int MaxHarmonic = 3;
const int MaxPower = 3;
const int NBinsPt = 2;
const int NBinsMult = 1;
const int NBinsEta = 18;
const int NBinsMass = 30;

const double BinPt[NBinsPt + 1] = {0., 5.  , 40.};
const double BinMult[NBinsMult + 1] = {30., 10000.};
const double BinEta[NBinsEta + 1] = {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
                                     0.0,
                                     0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
const double EdgesMass[2] = {2.4, 3.6};

/////////////////////////////////////////////////////////////////////////

const int N_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant = 2;
const int N_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant = 3;
const string Name_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[N_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant] = {"c1m1", "c2m2"};
const string Name_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[N_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant] = {"c1_m1", "c2_m2", "c22_m2m2"};

const int N_fProfileTarget_EtaWidth_Pt_cumulant = 7;
const int N_fProfileTarget_EtaWidth_Pt_Mass_cumulant = 7;
const string Name_fProfileTarget_EtaWidth_Pt_cumulant[N_fProfileTarget_EtaWidth_Pt_cumulant] = {"c2_m1_m1", "c2_2_m2_m2", "c2_0_m2_0", "c0_2_0_m2", "c2_0_0_m2", "c0_2_m2_0","c2_m2"};
const string Name_fProfileTarget_EtaWidth_Pt_Mass_cumulant[N_fProfileTarget_EtaWidth_Pt_Mass_cumulant] = {"c2_m1_m1", "c2_2_m2_m2", "c2_0_m2_0", "c0_2_0_m2", "c2_0_0_m2", "c0_2_m2_0","c2_m2"};

const int N_fProfileTargetBg_EtaWidth_Pt_cumulant = 7;
const int N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant = 7;
const string Name_fProfileTargetBg_EtaWidth_Pt_cumulant[N_fProfileTargetBg_EtaWidth_Pt_cumulant] = {"c2_m1_m1", "c2_2_m2_m2", "c2_0_m2_0", "c0_2_0_m2", "c2_0_0_m2", "c0_2_m2_0","c2_m2"};
const string Name_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant] = {"c2_m1_m1", "c2_2_m2_m2", "c2_0_m2_0", "c0_2_0_m2", "c2_0_0_m2", "c0_2_m2_0","c2_m2"};
}; // namespace MFlowGlobalVariables

namespace MFGV = MFlowGlobalVariables;

class MFlow
{
 private:
  std::vector<std::vector<TComplex>> fQvectors;

 public:
  MFlow(int maxHarmonic = MFGV::MaxHarmonic, int maxPower = MFGV::MaxPower);
  void Fill(double phi, double weight);
  TComplex Q(int harmonic, int power);
  void Reset();
  void Print();

  MFlow& operator+=(const MFlow& other);
  MFlow& operator-=(const MFlow& other);
  MFlow operator+(const MFlow& other);
  MFlow operator-(const MFlow& other);
  MFlow operator=(const MFlow& other);

  ~MFlow();
};

class MFGF
{
 public:
  enum TypeVariable {
    kNothing = -1,
    kPt = 0,
    kMult,
    kEta,
    kWeight,
    kCumulantIm,
    kCumulantRe,
    // four particle correlations
    kNVariables
  };

  static TString fgVariableNames[kNVariables];
  static TString fgVariableUnits[kNVariables];

  static std::array<MFlow, MFGV::NBinsEta> fArrayFlowReferenceEta;
  static std::array<std::array<MFlow, MFGV::NBinsPt>, MFGV::NBinsEta> fArrayFlowTargetEtaPt;
  static std::array<std::array<MFlow, MFGV::NBinsPt>, MFGV::NBinsEta> fArrayFlowTargetBgEtaPt;
  static std::array<std::array<std::array<MFlow, MFGV::NBinsPt>, MFGV::NBinsEta>, MFGV::NBinsMass> fArrayFlowTargetEtaPtMass;
  static std::array<std::array<std::array<MFlow, MFGV::NBinsPt>, MFGV::NBinsEta>, MFGV::NBinsMass> fArrayFlowTargetBgEtaPtMass;

  static double fMinAbsEtaReference;
  static double fMaxAbsEtaReference;
  static double fMinAbsEtaTarget;
  static double fMaxAbsEtaTarget;
  static double fMultiplicity;
  static double fMass;

  static TList* fListOutput;
  static TList* fListQA;
  static TList* fListReference;
  static TList* fListTarget;
  static TList* fListCorrelation;

  static TFile* fFileInput;

  // input histograms
  static TH2F* fh2NonUniformEtaPhiMultReference;

  static TH3F* fh3EtaPhiMultReference;
  static TH3F* fh3EtaPhiMultReferenceCorrected;
  static TH3F* fh3EtaPhiMultBeforeCut;
  static TH1F* fh1CutTrack;

  static int NBinsDeltaEtaReferenceSubEvent;
  static double* BinDeltaEtaReferenceSubEvent;
  static int NBinsEtaWidthReferenceTotal;
  static double* BinEtaWidthReferenceTotal;

  // static TProfile2D* fProfileFlowReferenceTotal_Eta_Mult[6];
  // Total_c1m1; Total_c2m2;
  // SubEvent_c1_m1; SubEvent_c2_m2; SubEvent_c22_m2m2
  static TProfile2D* fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant[MFGV::N_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant][2];
  static TProfile2D* fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant[MFGV::N_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant][2];
  static TH2F* fh2FlowReferenceTotal_EtaWidth_Mult_WeightSum[MFGV::N_fProfileFlowReferenceTotal_EtaWidth_Mult_cumulant][2];
  static TH2F* fh2FlowReferenceSubEvent_DeltaEta_Mult_WeightSum[MFGV::N_fProfileFlowReferenceSubEvent_DeltaEta_Mult_cumulant][2];

  // Target Flow
  static TProfile2D* fProfileTarget_EtaWidth_Pt_cumulant[MFGV::N_fProfileTarget_EtaWidth_Pt_cumulant][2];
  static TProfile2D* fProfileTargetBg_EtaWidth_Pt_cumulant[MFGV::N_fProfileTarget_EtaWidth_Pt_Mass_cumulant][2];
  static TH2F* fh2Target_EtaWidth_Pt_WeightSum[MFGV::N_fProfileTarget_EtaWidth_Pt_cumulant][2];
  static TH2F* fh2TargetBg_EtaWidth_Pt_WeightSum[MFGV::N_fProfileTarget_EtaWidth_Pt_Mass_cumulant][2];

  static TProfile3D* fProfileTarget_EtaWidth_Pt_Mass_cumulant[MFGV::N_fProfileTargetBg_EtaWidth_Pt_cumulant][2];
  static TProfile3D* fProfileTargetBg_EtaWidth_Pt_Mass_cumulant[MFGV::N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant][2];
  static TH3F* fh3Target_EtaWidth_Pt_Mass_WeightSum[MFGV::N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant][2];
  static TH3F* fh3TargetBg_EtaWidth_Pt_Mass_WeightSum[MFGV::N_fProfileTargetBg_EtaWidth_Pt_Mass_cumulant][2];

  static TH3F* fh3Target_Eta_Pt_Mass;
  static TH3F* fh3TargetBg_Eta_Pt_Mass;
  static TH1F* fh1DileptonMass_NoCut;
  static TH1F* fh1DileptonMass_MultCut;

  // Correlation
  static TH1F* fh1DeltaPhiReference;
  static TH1F* fh1DeltaPhiTargetReference;
  static TH1F* fh1DeltaPhiTargetBgReference;

  static void ArrayFlowInit();
  static void HistogramInit(bool doNeedNonUniformEtaPhiMultReference = false);
  static void ArrayFlowReferenceFill(double eta, double phi, double weight);
  static void ArrayFlowTargetFill(double eta, double pt, double phi, double weight);
  static void ArrayFlowTargetBgFill(double eta, double pt, double phi, double weight);
  static void ArrayFlowReferenceReset();
  static void ArrayFlowTargetReset();
  static void ArrayFlowTargetBgReset();
  static void CalculateFlowReference();
  static void CalculateFlowTarget();
}; // namespace MFlowGlobaFunctions

#endif // end of #ifndef MFLOW_H
