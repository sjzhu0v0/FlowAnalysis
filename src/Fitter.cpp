#include "Fitter.h"
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMath.h>
#include <TF1.h>
#include <TPaveText.h>
#include <iostream>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>

fitter::fitter()
{
    fGlobal = nullptr;
    fResidualBkg = nullptr;
    fSignal = nullptr;
    fCombinationBkg = nullptr;
}

fitter::~fitter()
{
    //
}
//______________________________________________________________________________

//set fit method
void fitter::set_fitmethod(string method)
{
    if (method=="unlike_like") {
        fitmethod=0;
        cout<<"------------fitting method is unlike-like sign-----------"<<endl;
    }else if (method=="fit_directly") {
        fitmethod=1;
        cout<<"------------fitting method is fit directly-----------"<<endl;
    }else if (method=="unlike_mixed") {
        fitmethod=2;
        cout<<"------------fitting method is unlike-mixed event-----------"<<endl;
    }else{
        cout<<"------------fitting method is set to default-----------"<<endl;
        cout<<"------------fitting method is unlike-like sign-----------"<<endl;
    }
};

//set fit function
void fitter::set_fitfunction(TF1* signal, TF1* residual, TF1* global, TF1* combination)
{ 
    fSignal = signal;
    fResidualBkg = residual;
    fGlobal = global;
    fCombinationBkg = combination;

    //setting
    fGlobal->SetLineColor(2);
    fGlobal->SetLineWidth(2);
    fGlobal->SetNpx(1000);

    if (fitmethod==0){
        fResidualBkg->SetLineColor(4);
        fResidualBkg->SetLineWidth(2);
        fResidualBkg->SetNpx(1000);
    }
    if (fitmethod==1){
        fSignal->SetLineColor(2);
        fSignal->SetLineWidth(2);
        fSignal->SetNpx(1000);
        fCombinationBkg->SetLineColor(4);
        fCombinationBkg->SetLineWidth(2);
        fCombinationBkg->SetNpx(1000);
    }
    if (fitmethod==2){
        fResidualBkg->SetLineColor(4);
        fResidualBkg->SetLineWidth(2);
        fResidualBkg->SetNpx(1000);
    }
};
//______________________________________________________________________________

//fit (standard mode)
void fitter::fit(TH1* unlike, TH1* like)
{
    if (fitmethod==0){//unlike sign - like sign method
        cout<<"using "<<fSignal->GetName()<<" as signal function"<<endl;
        cout<<"using "<<fResidualBkg->GetName()<<" as residual background function"<<endl;
        cout<<"---------------------------------------------------------"<<endl;
        hUnlikeSign=(TH1*)unlike->Clone();
        hLikeSign=(TH1*)like->Clone();
        Rawhisto=(TH1*)unlike->Clone();
        Rawhisto->Add(like,-1);
        set_histogram();
        set_rawhistogram();

        //prefit
        Rawhisto->Fit(fSignal,"R","",massrange[0],massrange[1]);
        for (int ipar_signal=0;ipar_signal<fSignal->GetNpar();ipar_signal++){
            fGlobal->SetParameter(ipar_signal,fSignal->GetParameter(ipar_signal));
        }
        TH1F* sideband = (TH1F*)Rawhisto->Clone();
        for (int ibin=1;ibin<=sideband->GetNbinsX();ibin++){
            if (sideband->GetBinCenter(ibin)>h_sideband[0] && sideband->GetBinCenter(ibin)<h_sideband[1]){
                sideband->SetBinContent(ibin,0);
                sideband->SetBinError(ibin,0);
            }
        }
        sideband->Fit(fResidualBkg,"R","",massrange[0],massrange[1]);
        for (int ipar_bkg=0;ipar_bkg<fResidualBkg->GetNpar();ipar_bkg++){
            // fGlobal->SetParameter(ipar_bkg+fSignal->GetNpar(),fResidualBkg->GetParameter(ipar_bkg));
            fGlobal->FixParameter(ipar_bkg+fSignal->GetNpar(),fResidualBkg->GetParameter(ipar_bkg));
        }
        sideband->Delete();
        
        //fit
        TFitResultPtr r = Rawhisto->Fit(fGlobal,"RS","",massrange[0],massrange[1]);
        
        //get fit result
        TMatrixDSym covTotal = r->GetCovarianceMatrix();
        covTotal.GetSub(fSignal->GetNpar(),fGlobal->GetNpar()-1,fSignal->GetNpar(),fGlobal->GetNpar()-1,covBG);//get the covariance matrix of background
        covTotal.GetSub(1,fGlobal->GetNpar()-1,1,fGlobal->GetNpar()-1,cov);//get the covariance matrix
        for (int ipar_bkg=0;ipar_bkg<fResidualBkg->GetNpar();ipar_bkg++){
            fResidualBkg->SetParameter(ipar_bkg,fGlobal->GetParameter(ipar_bkg+fSignal->GetNpar()));
        }
        for (int ipar_signal=0;ipar_signal<fSignal->GetNpar();ipar_signal++){
            fSignal->SetParameter(ipar_signal,fGlobal->GetParameter(ipar_signal));
        }
    }
    if (fitmethod==1){//fit unlike sign directly
        cout<<"using "<<fSignal->GetName()<<" as signal function"<<endl;
        cout<<"using "<<fCombinationBkg->GetName()<<" as combination background function"<<endl;
        hUnlikeSign=(TH1*)unlike->Clone("hUnlikeSign");
        set_histogram();

        //prefit
        hUnlikeSign->Fit(fSignal,"R","",massrange[0],massrange[1]);
        for (int ipar_signal=0;ipar_signal<fSignal->GetNpar();ipar_signal++){
            fGlobal->SetParameter(ipar_signal,fSignal->GetParameter(ipar_signal));
        }
        TH1F* sideband = (TH1F*)hUnlikeSign->Clone();
        for (int ibin=1;ibin<=sideband->GetNbinsX();ibin++){
            if (sideband->GetBinCenter(ibin)>h_sideband[0] && sideband->GetBinCenter(ibin)<h_sideband[1]){
                sideband->SetBinContent(ibin,0);
                sideband->SetBinError(ibin,0);
            }
        }
        sideband->Fit(fResidualBkg,"R","",massrange[0],massrange[1]);
        for (int ipar_bkg=0;ipar_bkg<fResidualBkg->GetNpar();ipar_bkg++){
            fGlobal->SetParameter(ipar_bkg+fSignal->GetNpar(),fResidualBkg->GetParameter(ipar_bkg));
        }
        sideband->Delete();

        //fit
        TFitResultPtr r = hUnlikeSign->Fit(fGlobal,"RSI","",massrange[0],massrange[1]);
        TMatrixDSym covTotal = r->GetCovarianceMatrix();
        covTotal.GetSub(fSignal->GetNpar(),fGlobal->GetNpar()-1,fSignal->GetNpar(),fGlobal->GetNpar()-1,covBG);//get the covariance matrix of background
        covTotal.GetSub(1,fGlobal->GetNpar()-1,1,fGlobal->GetNpar()-1,cov);//get the covariance matrix
        for (int ipar_signal=0;ipar_signal<fSignal->GetNpar();ipar_signal++){
            fSignal->SetParameter(ipar_signal,fGlobal->GetParameter(ipar_signal));
        }
        for (int ipar_bkg=0;ipar_bkg<fCombinationBkg->GetNpar();ipar_bkg++){
            fCombinationBkg->SetParameter(ipar_bkg,fGlobal->GetParameter(ipar_bkg+fSignal->GetNpar()));
        }
        //set rawhisto
        Rawhisto=(TH1*)unlike->Clone("Rawhisto");
        set_rawhistogram();
        for (int ibin=1;ibin<=Rawhisto->GetNbinsX();ibin++){
            Rawhisto->SetBinContent(ibin,Rawhisto->GetBinContent(ibin)-fCombinationBkg->Eval(Rawhisto->GetBinCenter(ibin)));
        }
        if (Rawhisto->GetMinimum()>0){
            Rawhisto->GetYaxis()->SetRangeUser(0,Rawhisto->GetMaximum());
        }else {
            Rawhisto->GetYaxis()->SetRangeUser(Rawhisto->GetMinimum()*1.25,Rawhisto->GetMaximum()*1.25);
        }
    }
    if (fitmethod==2){
        cout<<"using "<<fSignal->GetName()<<" as signal function"<<endl;
        cout<<"using "<<fResidualBkg->GetName()<<" as residual background function"<<endl;
        cout<<"---------------------------------------------------------"<<endl;
        hUnlikeSign=(TH1*)unlike->Clone();
        hLikeSign=(TH1*)like->Clone();
        hLikeSign->Scale(scale_factor);
        Rawhisto=(TH1*)unlike->Clone();
        Rawhisto->Add(hLikeSign,-1);
        set_histogram();
        set_rawhistogram();

        //prefit
        Rawhisto->Fit(fSignal,"NR","",massrange[0],massrange[1]);
        for (int ipar_signal=0;ipar_signal<fSignal->GetNpar();ipar_signal++){
            fGlobal->SetParameter(ipar_signal,fSignal->GetParameter(ipar_signal));
        }
        TH1F* sideband = (TH1F*)Rawhisto->Clone();
        for (int ibin=1;ibin<=sideband->GetNbinsX();ibin++){
            if (sideband->GetBinCenter(ibin)>h_sideband[0] && sideband->GetBinCenter(ibin)<h_sideband[1]){
                sideband->SetBinContent(ibin,0);
                sideband->SetBinError(ibin,0);
            }
        }
        sideband->Fit(fResidualBkg,"NR","",massrange[0],massrange[1]);
        for (int ipar_bkg=0;ipar_bkg<fResidualBkg->GetNpar();ipar_bkg++){
            fGlobal->SetParameter(ipar_bkg+fSignal->GetNpar(),fResidualBkg->GetParameter(ipar_bkg));
        }
        sideband->Delete();

        //fit
        TFitResultPtr r = Rawhisto->Fit(fGlobal,"RSI","",massrange[0],massrange[1]);
        
        //get fit result
        TMatrixDSym covTotal = r->GetCovarianceMatrix();
        covTotal.GetSub(fSignal->GetNpar(),fGlobal->GetNpar()-1,fSignal->GetNpar(),fGlobal->GetNpar()-1,covBG);//get the covariance matrix of background
        covTotal.GetSub(1,fGlobal->GetNpar()-1,1,fGlobal->GetNpar()-1,cov);//get the covariance matrix
        for (int ipar_bkg=0;ipar_bkg<fResidualBkg->GetNpar();ipar_bkg++){
            fResidualBkg->SetParameter(ipar_bkg,fGlobal->GetParameter(ipar_bkg+fSignal->GetNpar()));
        }
        for (int ipar_signal=0;ipar_signal<fSignal->GetNpar();ipar_signal++){
            fSignal->SetParameter(ipar_signal,fGlobal->GetParameter(ipar_signal));
        }
    }  
};
//______________________________________________________________________________
//calculate fit result
mysignal fitter::calculate(mysignal sig){
    mysignal signal;
    signal=sig;
    double rawcounts = 0;
    double rawcounts_err = 0;
    double combination_bkg = 0;
    double combination_bkg_err = 0;
    double residual_bkg = 0;
    double residual_bkg_err = 0;
    double rawsignal = 1;
    double rawsignal_err = 0;
    double totalbkg = 0;
    double totalbkg_err = 0;
    double S_over_B = 0;
    double S_over_B_err = 0;
    double significance = 0;

    //calculate the fit result
    cout<<"-----------------calculate the fit result-----------------"<<endl;
    //raw counts
    rawcounts=Rawhisto->Integral(Rawhisto->FindBin(signal.kMasswindow[0]),Rawhisto->FindBin(signal.kMasswindow[1]));
    Rawhisto->IntegralAndError(Rawhisto->FindBin(signal.kMasswindow[0]),Rawhisto->FindBin(signal.kMasswindow[1]),rawcounts_err);
    cout<<"rawcounts: "<<rawcounts<<"+-"<<rawcounts_err<<endl;
    //combination bkg
    if (fitmethod==0||fitmethod==2){
        combination_bkg=hLikeSign->Integral(hLikeSign->FindBin(signal.kMasswindow[0]),hLikeSign->FindBin(signal.kMasswindow[1]));
        hLikeSign->IntegralAndError(hLikeSign->FindBin(signal.kMasswindow[0]),hLikeSign->FindBin(signal.kMasswindow[1]),combination_bkg_err);
        cout<<"combination_bkg: "<<combination_bkg<<"+-"<<combination_bkg_err<<endl;
    }
    if (fitmethod==1){
        combination_bkg=fCombinationBkg->Integral(signal.kMasswindow[0],signal.kMasswindow[1]);
        combination_bkg_err=fCombinationBkg->IntegralError(signal.kMasswindow[0],signal.kMasswindow[1],fCombinationBkg->GetParameters(),covBG.GetMatrixArray());
        combination_bkg/=Rawhisto->GetBinWidth(1);
        combination_bkg_err/=Rawhisto->GetBinWidth(1);
        cout<<"combination_bkg: "<<combination_bkg<<"+-"<<combination_bkg_err<<endl;
    }
    //residual bkg
    if (fitmethod==0||fitmethod==2){
        residual_bkg=fResidualBkg->Integral(signal.kMasswindow[0],signal.kMasswindow[1]);
        residual_bkg_err=fResidualBkg->IntegralError(signal.kMasswindow[0],signal.kMasswindow[1],fResidualBkg->GetParameters(),covBG.GetMatrixArray());
        residual_bkg/=Rawhisto->GetBinWidth(1);
        residual_bkg_err/=Rawhisto->GetBinWidth(1);
        cout<<"residual_bkg: "<<residual_bkg<<"+-"<<residual_bkg_err<<endl;
    }
    if (fitmethod==1){
        residual_bkg=0;
        residual_bkg_err=0;
        cout<<"there is no residual bkg"<<endl;
    }
    //raw signal
    if (fitmethod==0||fitmethod==2){
        rawsignal=rawcounts-residual_bkg;
        rawsignal_err=sqrt(rawcounts_err*rawcounts_err+residual_bkg_err*residual_bkg_err);//need to be checked
        cout<<"rawsignal: "<<rawsignal<<"+-"<<rawsignal_err<<endl;
    }
    if (fitmethod==1){
        rawsignal=rawcounts;
        rawsignal_err=rawcounts_err;
        cout<<"rawsignal: "<<rawsignal<<"+-"<<rawsignal_err<<endl;
    }
    //total bkg
    totalbkg=combination_bkg+residual_bkg;
    totalbkg_err=sqrt(combination_bkg_err*combination_bkg_err+residual_bkg_err*residual_bkg_err);
    cout<<"totalbkg: "<<totalbkg<<"+-"<<totalbkg_err<<endl;
    //S/B
    S_over_B=rawsignal/totalbkg;
    S_over_B_err=rawsignal_err/totalbkg;
    cout<<"S/B: "<<S_over_B<<"+-"<<S_over_B_err<<endl;
    //significance
    if (fitmethod==0){
        significance=rawsignal/sqrt(rawsignal+2*totalbkg);
    }
    if (fitmethod==1){
        significance=rawsignal/sqrt(rawsignal+totalbkg);
    }
    if (fitmethod==2){
        significance=rawsignal/sqrt(rawsignal+totalbkg);
    }
    cout<<"significance: "<<significance<<endl;

    cout<<"-----------------------------------------------------------"<<endl;
    //save the fit result
    signal.kRawcounts=rawcounts;
    signal.kRawcounts_err=rawcounts_err;
    signal.kCombination_bkg=combination_bkg;
    signal.kCombination_bkg_err=combination_bkg_err;
    signal.kResidual_bkg=residual_bkg;
    signal.kResidual_bkg_err=residual_bkg_err;
    signal.kRawsignal=rawsignal;
    signal.kRawsignal_err=rawsignal_err;
    signal.kTotalbkg=totalbkg;
    signal.kTotalbkg_err=totalbkg_err;
    signal.kS_over_B=S_over_B;
    signal.kS_over_B_err=S_over_B_err;
    signal.kSignificance=significance;
    return signal;
}

//set histogram
void fitter::set_histogram()
{
    //rebin 
    double oldbinwidth = hUnlikeSign->GetBinWidth(1);
    int rebin = (int)(h_binwidth/oldbinwidth);
    //set unlike
    hUnlikeSign->SetName("hUnlikeSign");
    hUnlikeSign->Rebin(rebin);
    fitter::set_histogramstyle(hUnlikeSign,"#it{m}_{ee} (GeV/#it{c}^{2})",Form("Entries per %0.f MeV/#it{c}^{2}",h_binwidth*1000),2,20,1,2,1,2);
    hUnlikeSign->GetXaxis()->SetRangeUser(massrange[0],massrange[1]);
    hUnlikeSign->GetYaxis()->SetRangeUser(hUnlikeSign->GetMinimum()*0,hUnlikeSign->GetMaximum()*1.35);

    if (fitmethod==0||fitmethod==2){
        hLikeSign->SetName("hLikeSign");
        hLikeSign->Rebin(rebin);
        fitter::set_histogramstyle(hLikeSign,"#it{m}_{ee} (GeV/#it{c}^{2})",Form("Entries per %0.f MeV/#it{c}^{2}",h_binwidth*1000),4,20,1,4,1,2);
        hLikeSign->GetXaxis()->SetRangeUser(massrange[0],massrange[1]);
        hLikeSign->GetYaxis()->SetRangeUser(0,hLikeSign->GetMaximum()*1.65);
    }
};

//set raw histogram
void fitter::set_rawhistogram(){
    //set rawhisto
    Rawhisto->SetName("Rawhisto");
    double oldbinwidth = Rawhisto->GetBinWidth(1);
    int rebin = (int)(h_binwidth/oldbinwidth);
    Rawhisto->Rebin(rebin);
    fitter::set_histogramstyle(Rawhisto,"#it{m}_{ee} (GeV/#it{c}^{2})","",1,20,1,1,1,2);
    Rawhisto->GetXaxis()->SetRangeUser(massrange[0],massrange[1]);
}

//save canvas
TCanvas* fitter::save_fitcanvas(){
    cout<<"------------------standard canva-----------------"<<endl;
    TCanvas *c = new TCanvas("c","c",600,800);
    c->Divide(1,2,0,0);
    c->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.);
    gPad->SetTopMargin(0.06);
    gPad->SetRightMargin(0.);
    gPad->SetTicky();
    gPad->SetTickx();
    hUnlikeSign->Draw("E");
    if (fitmethod==0||fitmethod==2){
        hLikeSign->Draw("E same");
    }
    if (fitmethod==1){
        fGlobal->Draw("same");
        fCombinationBkg->Draw("same");
    }
    TLegend* leg = new TLegend(0.67,0.57,0.9,0.82);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->AddEntry(hUnlikeSign,"Unlike sign","p");
    if (fitmethod==0){
        leg->AddEntry(hLikeSign,"Like sign","p");
    }
    if (fitmethod==1){
        leg->AddEntry(fGlobal,"Total fit","l");
        leg->AddEntry(fCombinationBkg,"Combination bkg","l");
    }
    if (fitmethod==2){
        leg->AddEntry(hLikeSign,"mixed event","l");
    }
    leg->Draw();
    c->cd(2);
    gPad->SetTopMargin(0.);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.);
    gPad->SetTicky();
    Rawhisto->Draw("E");
    if (fitmethod==0||fitmethod==2){
        fGlobal->Draw("same");
        fResidualBkg->Draw("same");
    }
    if (fitmethod==1){
        fSignal->Draw("E same");
    }
    TLegend* leg2 = new TLegend(0.65,0.78,0.9,0.94);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.05);
    
    if (fitmethod==0||fitmethod==2)
    {
        leg2->AddEntry(Rawhisto,"Unlike - Like sign","p");
        leg2->AddEntry(fGlobal,"Total fit","l");
        leg2->AddEntry(fResidualBkg,"Residual","l");
    }
    if (fitmethod==1)
    {
        leg2->AddEntry(Rawhisto,"Opposite sign - Combinational bkg","p");
        leg2->AddEntry(fSignal,"Signal","l");
    }
    leg2->Draw();
    return c;
}
//set mass range for fitting
void fitter::set_massrange(double low, double high)
{
    massrange[0] = low;
    massrange[1] = high;
}
//______________________________________________________________________________

//set sideband range for fitting
void fitter::set_sideband(double low, double high)
{
    h_sideband[0] = low;
    h_sideband[1] = high;
}

//set bin width
void fitter::set_binwidth(double width)
{
    h_binwidth = width;
}

//______________________________________________________________________________

//set scale factor
void fitter::set_scalefactor(TH1* same_event, TH1* mixed_event)
{
    double same_event_counts, mixed_event_counts;
    double same_event_counts_err, mixed_event_counts_err;
    same_event_counts = same_event->IntegralAndError(same_event->GetXaxis()->FindBin(massrange[0]+1e-6),same_event->GetXaxis()->FindBin(massrange[1]-1e-6),same_event_counts_err);
    cout<<"same event counts = "<<same_event_counts<<endl;
    mixed_event_counts = mixed_event->IntegralAndError(mixed_event->GetXaxis()->FindBin(massrange[0]+1e-6),mixed_event->GetXaxis()->FindBin(massrange[1]-1e-6),mixed_event_counts_err);
    cout<<"mixed event counts = "<<mixed_event_counts<<endl;
    double factor = same_event_counts/mixed_event_counts;
    double factor_err = factor*TMath::Sqrt(TMath::Power(same_event_counts_err/same_event_counts,2)+TMath::Power(mixed_event_counts_err/mixed_event_counts,2));
    cout << "scale factor = " << factor << endl;
    cout << "scale factor error = " << factor_err << endl;
    scale_factor = factor;
}
//______________________________________________________________________________

//set style of histogram
void fitter::set_histogramstyle(TH1* h, string xaxisname, string yaxisname, int markercolor, int markerstyle, int markersize, int linecolor, int linestyle, int linewidth)
{
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetLabelOffset(0.015);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetYaxis()->SetLabelOffset(0.015);
    h->GetXaxis()->SetTitle(xaxisname.c_str());
    h->GetYaxis()->SetTitle(yaxisname.c_str());
    h->SetLineColor(1);
    h->SetLineWidth(1);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1);
    h->SetMarkerColor(1);
    h->SetLineStyle(1);
    h->SetStats(0);
    h->SetTitle("");
    h->SetMarkerColor(markercolor);
    h->SetMarkerStyle(markerstyle);
    h->SetMarkerSize(markersize);
    h->SetLineColor(linecolor);
    h->SetLineStyle(linestyle);
    h->SetLineWidth(linewidth);
}

//function library
double fitter::CB1(double *x, double *par)
{
    double Constant=par[0];
    double mean=par[1];
    double sigma=par[2];
    double alpha=par[3];
    double n=par[4];

    double t = (x[0]-mean)/sigma;
    if (alpha < 0) t = -t;
    double absAlpha = TMath::Abs(alpha);
    if (t >= -absAlpha) {
        return Constant*exp(-0.5*t*t);
    }
    else {
        double a = TMath::Power(n*1.0/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
        double b = n/absAlpha - absAlpha;
        return Constant*a*TMath::Power(b - t, -n);
    }
};

double fitter::CB2(double *x, double *par)
{
    double Constant=par[0];
    double mean=par[1];
    double sigma=par[2];
    double alpha=par[3];
    double n=par[4];
    double ratio=par[5];
    double mean2=par[6];

    double t = (x[0]-mean)/sigma;
    if (alpha < 0) t = -t;
    double absAlpha = TMath::Abs(alpha);
    double t2 = (x[0]-mean2)/sigma;
    if (alpha < 0) t2 = -t2;

    if (t >= -absAlpha) {
        if (t2 >= -absAlpha) {
            return Constant*exp(-0.5*t*t)+ratio*Constant*exp(-0.5*t2*t2);
        }
        else {
            double a = TMath::Power(n*1.0/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
            double b = n/absAlpha - absAlpha;
            return Constant*exp(-0.5*t*t)+ratio*Constant*a*TMath::Power(b - t2, -n);
        }
    }
    else {
        double a = TMath::Power(n*1.0/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
        double b = n/absAlpha - absAlpha;
        return Constant*a*TMath::Power(b - t, -n)+ratio*Constant*a*TMath::Power(b - t2, -n);
    }
};

double fitter::CB3(double *x, double *par)
{
    double Constant=par[0];
    double mean=par[1];
    double sigma=par[2];
    double alpha=par[3];
    double n=par[4];
    double ratio=par[5];
    double mean2=par[6];
    double ratio2=par[7];
    double mean3=par[8];

    double t = (x[0]-mean)/sigma;
    if (alpha < 0) t = -t;
    double absAlpha = TMath::Abs(alpha);
    double t2 = (x[0]-mean2)/sigma;
    if (alpha < 0) t2 = -t2;
    double t3 = (x[0]-mean3)/sigma;
    if (alpha < 0) t3 = -t3;

    if (t >= -absAlpha) {
        if (t2 >= -absAlpha) {
            if (t3 >= -absAlpha) {
                return Constant*exp(-0.5*t*t)+ratio*Constant*exp(-0.5*t2*t2)+ratio2*Constant*exp(-0.5*t3*t3);
            }
            else {
                double a = TMath::Power(n*1.0/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
                double b = n/absAlpha - absAlpha;
                return Constant*exp(-0.5*t*t)+ratio*Constant*exp(-0.5*t2*t2)+ratio2*Constant*a*TMath::Power(b - t3, -n);
            }
        }
        else {
            double a = TMath::Power(n*1.0/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
            double b = n/absAlpha - absAlpha;
            return Constant*exp(-0.5*t*t)+ratio*Constant*a*TMath::Power(b - t2, -n)+ratio2*Constant*a*TMath::Power(b - t3, -n);
        }
    }
    else {
        double a = TMath::Power(n*1.0/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
        double b = n/absAlpha - absAlpha;
        return Constant*a*TMath::Power(b - t, -n)+ratio*Constant*a*TMath::Power(b - t2, -n)+ratio2*Constant*a*TMath::Power(b - t3, -n);
    }
};

//get the function which is defined in this class
TF1* fitter::get_func(string funcname){
    TF1* func;
    if (funcname=="CB1"){
        func = new TF1("CB1",&fitter::CB1,0,15,5);
        func->SetParNames("Constant","mean","sigma","alpha","n");
    }
    if (funcname=="CB2"){
        func = new TF1("CB2",&fitter::CB2,0,15,7);
        func->SetParNames("Constant","mean","sigma","alpha","n","ratio","mean2");
    }
    if (funcname=="CB3"){
        func = new TF1("CB3",&fitter::CB3,0,15,9);
        func->SetParNames("Constant","mean","sigma","alpha","n","ratio","mean2","ratio2","mean3");
    }
    if (funcname=="pol2"){
        func = new TF1("pol2","[0]+[1]*x+[2]*x*x",0,15);
        func->SetParNames("a0","a1","a2");
    }
    if (funcname=="pol3"){
        func = new TF1("pol3","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0,15);
        func->SetParNames("a0","a1","a2","a3");
    }
    if (funcname=="expo2"){
        func = new TF1("expo2","exp([0]*x+[1])",0,15);
        func->SetParNames("a","b");
    }
    if (funcname=="expo3"){
        func = new TF1("expo3","exp([0]*x+[1])+[2]",0,15);
        func->SetParNames("a","b","c");
    }

    return func;
}