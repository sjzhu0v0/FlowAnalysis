//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 22 14:17:53 2023 by ROOT version 6.24/00
// from TTree Event/Event Information
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef MEvent_h
#define MEvent_h

#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TROOT.h>

// Header file for the classes stored in the TTree if any.

class MTrack {
private:
public:
  int fid;
  int fstatus;
  int fmother1;
  int fmother2;
  int fdaughter1;
  int fdaughter2;
  int fcol;
  int facol;
  float fpx;
  float fpy;
  float fpz;
  float fe;
  float fm;
  float fscale;
  float fxProd;
  float fyProd;
  float fzProd;
  float ftProd;
  MTrack(int id, int status, int mother1, int mother2, int daughter1,
         int daughter2, int col, int acol, float px, float py, float pz,
         float e, float m, float scale, float xProd, float yProd, float zProd,
         float tProd);
  ~MTrack(){};

  TLorentzVector GetMomentum4() { return TLorentzVector(fpx, fpy, fpz, fe); };

  TLorentzVector GetVertex4() {
    return TLorentzVector(fxProd, fyProd, fzProd, ftProd);
  };
};

MTrack::MTrack(int id, int status, int mother1, int mother2, int daughter1,
               int daughter2, int col, int acol, float px, float py, float pz,
               float e, float m, float scale, float xProd, float yProd,
               float zProd, float tProd)
    : fid(id), fstatus(status), fmother1(mother1), fmother2(mother2),
      fdaughter1(daughter1), fdaughter2(daughter2), fcol(col), facol(acol),
      fpx(px), fpy(py), fpz(pz), fe(e), fm(m), fscale(scale), fxProd(xProd),
      fyProd(yProd), fzProd(zProd), ftProd(tProd){};

class MEvent {
public:
  TTree *fChain;  //! pointer to the analyzed TTree or TChain
  Int_t fCurrent; //! current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  Int_t number_tracks_event;
  Int_t process_event;
  Int_t id_track[10000];        //[number_tracks_event]
  Int_t status_track[10000];    //[number_tracks_event]
  Int_t mother1_track[10000];   //[number_tracks_event]
  Int_t mother2_track[10000];   //[number_tracks_event]
  Int_t daughter1_track[10000]; //[number_tracks_event]
  Int_t daughter2_track[10000]; //[number_tracks_event]
  Int_t col_track[10000];       //[number_tracks_event]
  Int_t acol_track[10000];      //[number_tracks_event]
  Float_t px_track[10000];      //[number_tracks_event]
  Float_t py_track[10000];      //[number_tracks_event]
  Float_t pz_track[10000];      //[number_tracks_event]
  Float_t e_track[10000];       //[number_tracks_event]
  Float_t m_track[10000];       //[number_tracks_event]
  Float_t scale_track[10000];   //[number_tracks_event]
  Float_t xProd_track[10000];   //[number_tracks_event]
  Float_t yProd_track[10000];   //[number_tracks_event]
  Float_t zProd_track[10000];   //[number_tracks_event]
  Float_t tProd_track[10000];   //[number_tracks_event]

  // List of branches
  TBranch *b_number_tracks_event; //!
  TBranch *b_process_event;       //!
  TBranch *b_id_track;            //!
  TBranch *b_status_track;        //!
  TBranch *b_mother1_track;       //!
  TBranch *b_mother2_track;       //!
  TBranch *b_daughter1_track;     //!
  TBranch *b_daughter2_track;     //!
  TBranch *b_col_track;           //!
  TBranch *b_acol_track;          //!
  TBranch *b_px_track;            //!
  TBranch *b_py_track;            //!
  TBranch *b_pz_track;            //!
  TBranch *b_e_track;             //!
  TBranch *b_m_track;             //!
  TBranch *b_scale_track;         //!
  TBranch *b_xProd_track;         //!
  TBranch *b_yProd_track;         //!
  TBranch *b_zProd_track;         //!
  TBranch *b_tProd_track;         //!

  MEvent(TTree *tree = 0);
  virtual ~MEvent();
  virtual Int_t Cut(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry);
  MTrack GetTrack(int i);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree *tree);
  virtual Bool_t Notify();
  virtual void Show(Long64_t entry = -1);
};

MEvent::MEvent(TTree *tree) : fChain(0) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject("test.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("test.root");
    }
    f->GetObject("Event", tree);
  }
  Init(tree);
}

MEvent::~MEvent() {
  if (!fChain)
    return;
  delete fChain->GetCurrentFile();
}

Int_t MEvent::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain)
    return 0;
  return fChain->GetEntry(entry);
}

MTrack MEvent::GetTrack(int i) {
  return MTrack(id_track[i], status_track[i], mother1_track[i],
                mother2_track[i], daughter1_track[i], daughter2_track[i],
                col_track[i], acol_track[i], px_track[i], py_track[i],
                pz_track[i], e_track[i], m_track[i], scale_track[i],
                xProd_track[i], yProd_track[i], zProd_track[i], tProd_track[i]);
}

Long64_t MEvent::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void MEvent::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree)
    return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("number_tracks_event", &number_tracks_event,
                           &b_number_tracks_event);
  fChain->SetBranchAddress("process_event", &process_event, &b_process_event);
  fChain->SetBranchAddress("id_track", id_track, &b_id_track);
  fChain->SetBranchAddress("status_track", status_track, &b_status_track);
  fChain->SetBranchAddress("mother1_track", mother1_track, &b_mother1_track);
  fChain->SetBranchAddress("mother2_track", mother2_track, &b_mother2_track);
  fChain->SetBranchAddress("daughter1_track", daughter1_track,
                           &b_daughter1_track);
  fChain->SetBranchAddress("daughter2_track", daughter2_track,
                           &b_daughter2_track);
  fChain->SetBranchAddress("col_track", col_track, &b_col_track);
  fChain->SetBranchAddress("acol_track", acol_track, &b_acol_track);
  fChain->SetBranchAddress("px_track", px_track, &b_px_track);
  fChain->SetBranchAddress("py_track", py_track, &b_py_track);
  fChain->SetBranchAddress("pz_track", pz_track, &b_pz_track);
  fChain->SetBranchAddress("e_track", e_track, &b_e_track);
  fChain->SetBranchAddress("m_track", m_track, &b_m_track);
  fChain->SetBranchAddress("scale_track", scale_track, &b_scale_track);
  fChain->SetBranchAddress("xProd_track", xProd_track, &b_xProd_track);
  fChain->SetBranchAddress("yProd_track", yProd_track, &b_yProd_track);
  fChain->SetBranchAddress("zProd_track", zProd_track, &b_zProd_track);
  fChain->SetBranchAddress("tProd_track", tProd_track, &b_tProd_track);
  Notify();
}

Bool_t MEvent::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void MEvent::Show(Long64_t entry) {
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain)
    return;
  fChain->Show(entry);
}
Int_t MEvent::Cut(Long64_t entry) {
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

#endif // #ifdef MEvent_cxx
