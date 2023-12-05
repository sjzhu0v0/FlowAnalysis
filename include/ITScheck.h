#include <TH1.h>
#include <TMath.h>
#include <TH2.h>
#include <iostream>

#ifndef _ALICE_RUN3_ITS_TOOL_H_
#define _ALICE_RUN3_ITS_TOOL_H_
using namespace std;

class its
{
    public:
        static string itsmatchmode[4];
        static bool ismatch(string type, int itsmap);
        static double geteff(TH1F* itsclustermap, string type);
        static double geterr(TH1F* itsclustermap, string type);
        static TH1F* getdistribution(TH2F* itsclustermap_vs_var, string type);
        static TH1F* convert_to_itseff(TH2*, double*, int, string type);
        static TH1F* convert_to_layer(TH1*);

};

string its::itsmatchmode[4] = {"ITSibAny","ITSallAny","ITSall7Layers","ITSncls2"};

bool its::ismatch(string type, int itsmap)//can be add more match mode
{
    if (type=="ITSibAny"){
        if ((itsmap & (uint8_t(1) << 0)) || (itsmap & (uint8_t(1) << 1)) || (itsmap & (uint8_t(1) << 2))){
            return true;
        }
    }
    else if (type=="ITSallAny"){
        if ((itsmap & (uint8_t(1) << 0)) || (itsmap & (uint8_t(1) << 1)) || (itsmap & (uint8_t(1) << 2)) || (itsmap & (uint8_t(1) << 3)) || (itsmap & (uint8_t(1) << 4)) || (itsmap & (uint8_t(1) << 5)) || (itsmap & (uint8_t(1) << 6))){
            return true;
        }
    }
    else if (type=="ITSall7Layers"){
        if ((itsmap & (uint8_t(1) << 0)) && (itsmap & (uint8_t(1) << 1)) && (itsmap & (uint8_t(1) << 2)) && (itsmap & (uint8_t(1) << 3)) && (itsmap & (uint8_t(1) << 4)) && (itsmap & (uint8_t(1) << 5)) && (itsmap & (uint8_t(1) << 6))){
            return true;
        }
    }
    else if (type=="ITSncls2"){
        int ncls=0;
        for (int i=0;i<7;i++){
            if (itsmap & (uint8_t(1) << i)){
                ncls++;
            }
        }
        if (ncls>=2){
            return true;
        }
    }
    return false;
};

double its::geteff(TH1F* itsclustermap, string type)
{
    double allhits=0;
    double selhits=0;
    for (int imap=0; imap<128; imap++){
        double nhits = itsclustermap->GetBinContent(imap+1);
        allhits+=nhits;        
        if (ismatch(type,imap)){
            selhits+=nhits;
        };
    }
    return selhits/allhits;
};

double its::geterr(TH1F*itsclustermap, string type)
{
    double allhits=0;
    double selhits=0;
    for (int imap=0; imap<128; imap++){
        double nhits = itsclustermap->GetBinContent(imap+1);
        allhits+=nhits;        
        if (ismatch(type,imap)){
            selhits+=nhits;
        };
    }
    return sqrt(selhits)/allhits;
}

TH1F* its::getdistribution(TH2F* itsclustermap_vs_var, string type)
{
    TH2F* itsclustermap_vs_var_clone = (TH2F*)itsclustermap_vs_var->Clone();
    TH1F* distribution = (TH1F*)itsclustermap_vs_var_clone->ProjectionX();
    distribution->Reset();
    //reset itsclustermap
    for (int imap=0; imap<128; imap++){
        if (ismatch(type,imap)){
            distribution->Add(itsclustermap_vs_var_clone->ProjectionX(Form("imap%d",imap),imap+1,imap+1));
            }
    }
    return distribution;
};

TH1F* its::convert_to_itseff(TH2* map_var,double* binning, int nbin, string type)
{
    if (map_var==nullptr){
        cout << "ITSMap is nullptr" << endl;
        return nullptr;
    }
    TH1F* eff_var = new TH1F(Form("eff_%s",type.c_str()),Form("eff_%s",type.c_str()),nbin,binning);
    time_t now = time(0);
    char* dt = ctime(&now);
    eff_var->SetName(dt);
    for (int ibin=0;ibin<nbin;ibin++){
        TH1F* itsclustermap = (TH1F*)map_var->ProjectionY(Form("itsclustermap_%d",ibin),map_var->GetXaxis()->FindBin(binning[ibin]+0.00001),map_var->GetXaxis()->FindBin(binning[ibin+1]-0.00001));
        double eff = its::geteff(itsclustermap,type);
        double err = its::geterr(itsclustermap,type);
        eff_var->SetBinContent(ibin+1,eff);
        eff_var->SetBinError(ibin+1,err);
    }
    return (TH1F*)eff_var;
};

TH1F* its::convert_to_layer(TH1* map)
{
    if (map==nullptr){
        cout << "ITSMap is nullptr" << endl;
        return nullptr;
    }
    TH1F* layer = new TH1F("layer","layer",7,0.5,7.5);
    double alltracks = map->GetEntries();
    double hits_in_layer[7]={0,0,0,0,0,0,0};

    for (int imap=0; imap<128; imap++){
        double nhits = map->GetBinContent(imap+1);
        for (int ilayer=0;ilayer<7;ilayer++){
            if (imap & (uint8_t(1) << ilayer)){
                hits_in_layer[ilayer]+=nhits;
            }
        }
    }

    for (int ilayer=0;ilayer<7;ilayer++){
        layer->SetBinContent(ilayer+1,hits_in_layer[ilayer]/alltracks);
        layer->SetBinError(ilayer+1,sqrt(hits_in_layer[ilayer])/alltracks);
    }
    return (TH1F*)layer;
};


#endif
