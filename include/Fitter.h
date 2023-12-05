// Created: 2023-06-29
// Creater: zhenjun xiong
// mail: zhenjun.xiong@cern.ch
// what is this: this is a class for fitting

#include <TF1.h>
#include <TH1.h>
#include <TMatrixDSym.h>
#include <TPad.h>

#ifndef _FITTER_H_
#define _FITTER_H_
using namespace std;
typedef struct {
    // string kSignalname;
    double kMasswindow[2];
    int kWindowcolor;
    double kRawcounts;
    double kRawcounts_err;
    double kCombination_bkg;
    double kCombination_bkg_err;
    double kResidual_bkg;
    double kResidual_bkg_err;
    double kRawsignal;
    double kRawsignal_err;
    double kTotalbkg;
    double kTotalbkg_err;
    double kS_over_B;
    double kS_over_B_err;
    double kSignificance;
} mysignal;

class fitter
{
    public:
        fitter();
        ~fitter();
        //setting
        TF1* fGlobal;//signal+residual
        TF1* fResidualBkg;
        TF1* fSignal;
        TF1* fCombinationBkg;
        //used histogram
        TH1* hMCsignal;
        TH1* hUnlikeSign;
        TH1* hLikeSign;
        TH1* hCombinationBkg;
        TH1* Rawhisto;
        int fitmethod=0;//0:unlike-like sign, 1:fit unlike sign directly, 2:unlike-mixed event
        double scale_factor=1.0;
        void set_fitfunction(TF1*,TF1*,TF1*,TF1*);
        void set_fitmethod(string);
        void set_massrange(double,double);
        void set_sideband(double,double);
        void set_binwidth(double);
        void set_histogram();
        void set_rawhistogram();
        void set_scalefactor(TH1*,TH1*);
        void set_histogramstyle(TH1*,string,string,int,int,int,int,int,int);
        //calculate
        // void calculate_fitresult();
        mysignal calculate(mysignal sig);
        //fit
        void fit(TH1*,TH1*);//fit the signal and background
        //save result
        TCanvas* save_fitcanvas();

        //functuin library
        TF1* get_func(string);
        static double CB1(double *x, double *par);
        static double CB2(double *x, double *par);
        static double CB3(double *x, double *par);
    private:
        //fit result and calculation
        TMatrixDSym cov;
        TMatrixDSym covBG;
        //some information
        // double masswindow[2]={2.92,3.16};
        double massrange[2]={2,4};
        double h_sideband[2]={2.7,3.2};
        double h_binwidth=0.04;//Gev/c^2

};


#endif