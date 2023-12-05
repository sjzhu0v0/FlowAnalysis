#include "libmytool.h"
#include <iostream>
#include <THashList.h>
#ifndef _HISTOGRAM_MANAGER_CPP_
#define _HISTOGRAM_MANAGER_CPP_

using namespace std;
tool ::tool()
{
    //to do
}

tool ::~tool()
{
    //to do
}
//______________________________________________________________________________
TList* tool::GetList(TFile* mfile, const string& hashlistname, const string& listname)
{
    if (!mfile) {
        // Handle the case where the file is not valid.
        return nullptr;
    }

    THashList* HashList = dynamic_cast<THashList*>(mfile->Get(hashlistname.c_str()));
    if (!HashList) {
        // Handle the case where the HashList is not found.
        return nullptr;
    }

    TList* List = dynamic_cast<TList*>(HashList->FindObject(listname.c_str()));
    if (!List) {
        // Handle the case where the List is not found.
        return nullptr;
    }
    return List;
}
//______________________________________________________________________________
TList* tool::GetList(TFile* mfile, const string& listname)
{
    if (!mfile) {
        // Handle the case where the file is not valid.
        return nullptr;
    }

    TList* List = dynamic_cast<TList*>(mfile->Get(listname.c_str()));

    if (!List) {
        // Handle the case where the List is not found.
        return nullptr;
    }
    return List;
}
//______________________________________________________________________________

//project the 2D histogram to 1D histogram
TH1* tool::ProjectX(TH2* h2, double ylow, double yhigh, string name)
{
    if (h2==nullptr) {
        cout<<"histo is nullptr, cannot project"<<endl;
        return nullptr;
    }
    // Find the bin numbers corresponding to the specified y-range.
    int binLow = h2->GetYaxis()->FindBin(ylow);
    int binHigh = h2->GetYaxis()->FindBin(yhigh);
    // Ensure binLow is less than or equal to binHigh.
    if (binLow > binHigh) {
        std::swap(binLow, binHigh);
    }
    // Project onto the X-axis within the specified y-range.
    TH1* h1 = h2->ProjectionX(name.c_str(), binLow, binHigh);
    TH1* h1Clone = dynamic_cast<TH1*>(h1->Clone());
    return h1Clone;
}
TH1* tool::ProjectY(TH2* h2, double xlow, double xhigh, string name)
{
    if (h2==nullptr) {
        cout<<"histo is nullptr, cannot project"<<endl;
        return nullptr;
    }
    // Find the bin numbers corresponding to the specified y-range.
    int binLow = h2->GetXaxis()->FindBin(xlow);
    int binHigh = h2->GetXaxis()->FindBin(xhigh);
    // Ensure binLow is less than or equal to binHigh.
    if (binLow > binHigh) {
        std::swap(binLow, binHigh);
    }
    // Project onto the X-axis within the specified y-range.
    TH1* h1 = h2->ProjectionY(name.c_str(), binLow, binHigh);
    TH1* h1Clone = dynamic_cast<TH1*>(h1->Clone());
    return h1Clone;
}

//project the 3D histogram to 1D histogram
TH1* tool::Project(TH3* h3, double ylow, double yhigh, double zlow, double zhigh)
{
    int yBinLow = h3->GetYaxis()->FindBin(ylow + 10e-6);
    int yBinHigh = h3->GetYaxis()->FindBin(yhigh - 10e-6);
    int zBinLow = h3->GetZaxis()->FindBin(zlow + 10e-6);
    int zBinHigh = h3->GetZaxis()->FindBin(zhigh - 10e-6);
    // Ensure proper bin order.
    if (yBinLow > yBinHigh) {
        std::swap(yBinLow, yBinHigh);
    }
    if (zBinLow > zBinHigh) {
        std::swap(zBinLow, zBinHigh);
    }
    // Project onto the X-axis within the specified y and z ranges.
    TH1* h1 = h3->ProjectionX("h1", yBinLow, yBinHigh, zBinLow, zBinHigh);
    // Clone the projected histogram for safe memory management.
    TH1* h1Clone = dynamic_cast<TH1*>(h1->Clone("h1Clone"));
    return h1Clone;
}

//project the nD histogram to 1D histogram
TH1* tool::Project(THn* hn, int dim, std::vector<double> low, std::vector<double> high)
{
    if (!hn || dim <= 0 || low.size() < dim || high.size() < dim) {
        return nullptr;
    }

    THn* ndimhisto = dynamic_cast<THn*>(hn->Clone());
    if (!ndimhisto) {
        return nullptr;
    }

    for (int i = 0; i < dim; ++i) {
        ndimhisto->GetAxis(i)->SetRangeUser(low[i], high[i]);
    }

    TH1* h1 = ndimhisto->Projection(0);
    TH1* h1Clone = static_cast<TH1*>(h1->Clone());

    delete ndimhisto;
    delete h1;

    return h1Clone;
}



//______________________________________________________________________________

void tool::Standardize(TH1* h1, string xaxisname, string yaxisname){
    h1->GetXaxis()->SetTitleOffset(1.1);
    h1->GetYaxis()->SetTitleOffset(1.1);
    h1->GetXaxis()->SetTitleSize(0.06);
    h1->GetYaxis()->SetTitleSize(0.06);
    h1->GetXaxis()->SetLabelSize(0.06);
    h1->GetXaxis()->SetLabelOffset(0.015);
    h1->GetYaxis()->SetLabelSize(0.06);
    h1->GetYaxis()->SetLabelOffset(0.015);
    h1->GetXaxis()->SetTitle(xaxisname.c_str());
    h1->GetYaxis()->SetTitle(yaxisname.c_str());
    h1->SetLineColor(1);
    h1->SetLineWidth(1);
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(1);
    h1->SetMarkerColor(1);
    h1->SetLineStyle(1);
    h1->SetStats(0);
    h1->SetTitle("");
}
//______________________________________________________________________________

void tool::SetMarker(TH1* h, int color, int markerstyle, int markersize){
    h->SetMarkerColor(color);
    h->SetMarkerStyle(markerstyle);
    h->SetMarkerSize(markersize);
}
//______________________________________________________________________________

void tool::SetLine(TH1* h, int color, int linestyle, int linewidth){
    h->SetLineColor(color);
    h->SetLineStyle(linestyle);
    h->SetLineWidth(linewidth);
}
//______________________________________________________________________________

TLatex* tool::SetLatex(double x, double y, char* text, int textFont,double textsize, int textcolor)
{
    TLatex* tex = new TLatex(x, y, text);
    tex->SetTextFont(textFont);
    tex->SetTextSize(textsize);
    tex->SetTextColor(textcolor);
    tex->SetNDC();
    tex->Draw();
    return tex;
};
//______________________________________________________________________________

TLine* tool::SetLine(double x1, double y1, double x2, double y2, int linecolor, int linestyle, int linewidth)
{
    TLine* line = new TLine(x1, y1, x2, y2);
    line->SetLineColor(linecolor);
    line->SetLineStyle(linestyle);
    line->SetLineWidth(linewidth);
    line->Draw();
    return line;
};

double tool::CalculateCorrelationCoefficient(TH2* hist) 
{
    // 获取直方图的数量和均值
    int nbinsX = hist->GetNbinsX();
    int nbinsY = hist->GetNbinsY();
    double meanX = hist->GetMean(1);
    double meanY = hist->GetMean(2);

    // 初始化相关系数的分子和分母
    double numerator = 0.0;
    double denominatorX = 0.0;
    double denominatorY = 0.0;

    // 遍历直方图的 bin
    for (int binX = 1; binX <= nbinsX; ++binX) {
        for (int binY = 1; binY <= nbinsY; ++binY) {
            double value = hist->GetBinContent(binX, binY);
            double x = hist->GetXaxis()->GetBinCenter(binX);
            double y = hist->GetYaxis()->GetBinCenter(binY);

            // 更新分子和分母
            numerator += (x - meanX) * (y - meanY) * value;
            denominatorX += TMath::Power(x - meanX, 2) * value;
            denominatorY += TMath::Power(y - meanY, 2) * value;
        }
    }

    // 计算相关系数
    double correlationCoefficient = numerator / (TMath::Sqrt(denominatorX) * TMath::Sqrt(denominatorY));

    return correlationCoefficient;
}

#endif