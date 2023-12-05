#include "libmytool.h"
#include "Fitter.h"

struct SignalExtInfo
{
    TCanvas* c;
    mysignal Jpsi;
};

SignalExtInfo SignalExt(fitter* myfunc, mysignal Jpsi,TH1* h1SE_unlike, TH1* h1SE_like, TH1* h1ME_unlike, TH1* h1ME_like){
    if (myfunc->fitmethod==0){//unlike-like sign
        myfunc->fit(h1SE_unlike,h1SE_like);
    }
    else if (myfunc->fitmethod==1){//fit unlike sign directly
        myfunc->fit(h1SE_unlike,h1SE_like);
    }
    else if (myfunc->fitmethod==2){//unlike-mixed event
        myfunc->set_scalefactor(h1SE_like,h1ME_like);
        myfunc->fit(h1SE_unlike,h1ME_like);
    }
    else{
        cout<<"fitmethod is not defined"<<endl;
    }
    return SignalExtInfo{myfunc->save_fitcanvas(), myfunc->calculate(Jpsi)};
};
