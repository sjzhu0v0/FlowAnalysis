// Created: 2023-04-18
// Creater: zhenjun xiong
// mail: zhenjun.xiong@cern.ch
#include <TFile.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TProfile.h>
#include "string"
#include <TLatex.h>
#include <TLine.h>

#ifndef _LIBMYTOOL_H_
#define _LIBMYTOOL_H_
using namespace std;

class tool
{
    public:
        tool();
        ~tool();
        static TList* GetList(TFile* file, const string& hashlistname, const string& listname);
        static TList* GetList(TFile* file, const string& listname);
        template <class HISTO>
        static HISTO* GetHisto(TList* list, string histoname);
        static TH1* ProjectX(TH2* h2, double ylow, double yhigh, string name);
        static TH1* ProjectY(TH2* h2, double xlow, double xhigh, string name);
        static TH1* Project(TH3* h3, double ylow, double yhigh, double zlow, double zhigh);
        static TH1* Project(THn* hn, int dim, std::vector<double> low, std::vector<double> high);
        //some functions to set the histogram
        static void Standardize(TH1* h, string, string);
        static void SetMarker(TH1* h, int color, int markerstyle, int markersize);
        static void SetLine(TH1* h, int color, int linestyle, int linewidth);
        //some functions to set pad
        static TLatex* SetLatex(double x, double y, char* text, int textFont,double textsize, int textcolor);
        static TLine* SetLine(double x1, double y1, double x2, double y2, int linecolor, int linestyle, int linewidth);
        //calculate correlation coefficient between x and y
        static double CalculateCorrelationCoefficient(TH2* hist); 
    private:
};

//return the list in TFile
//due to some unknow error, the function can not be defined in HistogramManager.cpp
template <class HISTO>
HISTO* tool::GetHisto(TList* list, string histoname)
{
    // HISTO* histo = (HISTO*)list->FindObject(histoname.c_str())->Clone();
    if (!list) {
        cout << "list is nullptr!!!!" << endl;
        return nullptr;
    }
    HISTO* histo = dynamic_cast<HISTO*>(list->FindObject(histoname.c_str()));
    if (!histo) {
        return nullptr;
    }
    return histo;
};

#endif
