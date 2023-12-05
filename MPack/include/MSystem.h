#ifndef MSYSTEM_H
#define MSYSTEM_H

#include "TApplication.h"
#include "TBrowser.h"
#include "TCanvas.h"
#include "TClass.h"
#include "TCollection.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGButton.h"
#include "TGDNDManager.h"
#include "TGFileDialog.h"
#include "TGLabel.h"
#include "TComplex.h"
#include "TGListTree.h"
#include "TGMenu.h"
#include "TGMsgBox.h"
#include "TGPicture.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH1S.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH2Poly.h"
#include "TH3.h"
#include "TH3F.h"
#include "THashList.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLeaf.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TMVA/Reader.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMessage.h"
#include "TObjString.h"
#include "TObject.h"
#include "TPaletteAxis.h"
#include "TPaveLabel.h"
#include "TPaveStats.h"
#include "TProfile3D.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TRootEmbeddedCanvas.h"
#include "TRootHelpDialog.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTimer.h"
#include "TTree.h"
#include "TVector3.h"
#include "TProfile2D.h"
#include <TBenchmark.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TFrame.h>
#include <TH2.h>
#include <THttpServer.h>
#include <TInterpreter.h>
#include <TMemFile.h>
#include <TNtuple.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>


#include "fstream"
#include "iomanip"
#include "iosfwd"
#include "iostream"
#include "sstream"
#include "string"
#include "vector"

using namespace std;

// unique id generator
int GenerateUID() {
  static int id = 0;
  return id++;
}

TString dec2hex(int dec){
  return TString::Format("%X", dec);
}

#endif