#ifndef MPHYSICS_H
#define MPHYSICS_H

#include "MSystem.h"

const double mass_proton = 0.938272;
const double mass_kaon = 0.493677;
const double mass_pion = 0.139570;
const double mass_muon = 0.105658;
const double mass_electron = 0.000511;
const double mass_d0 = 1.86483;
const double mass_lambda_c = 2.2846;

std::vector<double> CalculateCumulants(TH1 *hist, int n_required = 4) {
  std::vector<double> cumulants;

  // Calculate the mean
  double mean = hist->GetMean();

  // Calculate the first cumulant (mean)
  cumulants.push_back(mean);

  // Calculate the variance
  double variance = hist->GetStdDev() * hist->GetStdDev();

  // Calculate the second cumulant (variance)
  cumulants.push_back(variance);

  // Calculate the higher order cumulants
  for (int i = 3; i <= n_required; i++) {
    double sum = 0;
    for (int j = 1; j <= i; j++) {
      double term =
          TMath::Binomial(i, j) * TMath::Power(-1, j + 1) * cumulants[i - j];
      sum += term;
    }
    cumulants.push_back(sum / TMath::Power(variance, i / 2.0));
  }

  return cumulants;
}

double *ProcessPowerFunction(double x, double n_power, double err_x) {
  double *result = new double[2]; // [0]: result, [1]: error of result, [2]:
                                  // power, [3]: error of power
  result[0] = TMath::Power(x, n_power);
  result[1] = TMath::Abs(n_power * TMath::Power(x, n_power - 1) * err_x);
  return result;
}

double *ProcessDivision(double numerator, double denominator, double err_num,
                        double err_den) {
  double *result = new double[2]; // [0]: result, [1]: error of result
  result[0] = numerator / denominator;
  result[1] = TMath::Sqrt(TMath::Power(err_num / denominator, 2) +
                          TMath::Power(numerator * err_den, 2) /
                              TMath::Power(denominator, 4));
  return result;
}

double *ProcessSignal2Total(double signal, double total, double err_sig,
                            double err_tot) {
  double *result = new double[2]; // [0]: result, [1]: error of result
  result[0] = signal / total;
  result[1] = sqrt((total - 2 * signal) / TMath::Power(total, 3) *
                       TMath::Power(err_sig, 2) +
                   TMath::Power(signal, 2) / TMath::Power(total, 4) *
                       TMath::Power(err_tot, 2));
  return result;
}

double *ProcessMultiplication(double factor1, double factor2, double err_fac1,
                              double err_fac2) {
  double *result = new double[2]; // [0]: result, [1]: error of result
  result[0] = factor1 * factor2;
  result[1] = TMath::Sqrt(TMath::Power(err_fac1 * factor2, 2) +
                          TMath::Power(err_fac2 * factor1, 2));
  return result;
}

double *ProcessComplex(double im, double re, double err_im, double err_re) {
  double *result = new double[4]; // [0]: magnitude, [1]: error of magnitude,
                                  // [2]: phase, [3]: error of phase
  result[0] = TMath::Sqrt(im * im + re * re);
  result[1] =
      TMath::Sqrt(TMath::Power(im * err_im, 2) + TMath::Power(re * err_re, 2)) /
      result[0];
  result[2] = TMath::ATan2(im, re);
  result[3] =
      TMath::Sqrt(TMath::Power(err_im * re, 2) + TMath::Power(err_re * im, 2)) /
      (im * im + re * re);
  return result;
}

TComplex *ProcessComplexProduct(TComplex c1, TComplex c2, TComplex c1_err,
                                TComplex c2_err) {
  TComplex *result =
      new TComplex[2]; // [0]: magnitude, [1]: error of magnitude,
                       // [2]: phase, [3]: error of phase
  TComplex c = c1 * c2;
  TComplex c_err = {
      sqrt(pow(c1.Re() * c2_err.Re(), 2) + pow(c2.Re() * c1_err.Re(), 2) +
           pow(c1.Im() * c2_err.Im(), 2) + pow(c2.Im() * c1_err.Im(), 2)),
      sqrt(pow(c1.Re() * c2_err.Im(), 2) + pow(c1.Im() * c2_err.Re(), 2) +
           pow(c2.Re() * c1_err.Im(), 2) + pow(c2.Im() * c1_err.Im(), 2))};
  result[0] = c;
  result[1] = c_err;
  return result;
}

TComplex *ProcessComplexAddition(TComplex c1, TComplex c2, TComplex c1_err,
                                 TComplex c2_err) {
  TComplex *result =
      new TComplex[2]; // [0]: magnitude, [1]: error of magnitude,
                       // [2]: phase, [3]: error of phase
  TComplex c = c1 + c2;
  TComplex c_err = {sqrt(pow(c1_err.Re(), 2) + pow(c2_err.Re(), 2) +
                         pow(c1_err.Im(), 2) + pow(c2_err.Im(), 2)),
                    sqrt(pow(c1_err.Re(), 2) + pow(c2_err.Re(), 2) +
                         pow(c1_err.Im(), 2) + pow(c2_err.Im(), 2))};
  result[0] = c;
  result[1] = c_err;
  return result;
}

class MComplex {
public:
  MComplex() {
    re = 0;
    im = 0;
    re_err = 0;
    im_err = 0;
  }
  MComplex(double r, double i, double r_err, double i_err) {
    re = r;
    im = i;
    re_err = r_err;
    im_err = i_err;
  }
  MComplex(TComplex c, TComplex c_err) {
    re = c.Re();
    im = c.Im();
    re_err = c_err.Re();
    im_err = c_err.Im();
  }

  double Re() { return re; }
  double Im() { return im; }
  double ReErr() { return re_err; }
  double ImErr() { return im_err; }

  void SetRe(double r) { re = r; }
  void SetIm(double i) { im = i; }
  void SetReErr(double r_err) { re_err = r_err; }
  void SetImErr(double i_err) { im_err = i_err; }

  double Mag() { return sqrt(re * re + im * im); }
  double MagErr() {
    return sqrt(pow(re * re_err, 2) + pow(im * im_err, 2)) / Mag();
  }
  double Mag2() { return re * re + im * im; }
  double Mag2Err() { return sqrt(pow(re * re_err, 2) + pow(im * im_err, 2)); }
  double Phase() { return atan2(im, re); }
  double PhaseErr() {
    return sqrt(pow(re * im_err, 2) + pow(im * re_err, 2)) /
           (re * re + im * im);
  }

  MComplex &operator!() {
    im = -im;
    im_err = abs(im_err);
    return *this;
  }
  double operator() (int i) {
    if (i == 0) {
      return re;
    } else if (i == 1) {
      return im;
    } else if (i == 2) {
      return re_err;
    } else if (i == 3) {
      return im_err;
    } else {
      cout << "Error: MComplex::operator() (int i) - i must be 0, 1, 2, or 3."
           << endl;
      return 0;
    }
  }
  MComplex operator*(MComplex c) {
    return MComplex(re * c.Re() - im * c.Im(), re * c.Im() + im * c.Re(),
                    sqrt(pow(c.Re() * re_err, 2) + pow(c.Im() * im_err, 2) +
                         pow(re * c.ReErr(), 2) + pow(im * c.ImErr(), 2)),
                    sqrt(pow(c.Re() * im_err, 2) + pow(c.Im() * re_err, 2) +
                         pow(im * c.ReErr(), 2) + pow(re * c.ImErr(), 2)));
  }
  MComplex operator*(double c) {
    return MComplex(re * c, im * c, abs(re_err * c), abs(im_err * c));
  }
  MComplex operator/(MComplex c) {
    double err_re1 = re_err * c.Re() / c.Mag2();
    double err_re2 =
        c.ReErr() * re * (pow(c.Im(), 2) - pow(c.Re(), 2)) / c.Mag2();
    double err_re3 = im_err * c.Im() / c.Mag2();
    double err_re4 =
        c.ImErr() * im * (pow(c.Re(), 2) - pow(c.Im(), 2)) / c.Mag2();
    double err_re = sqrt(pow(err_re1, 2) + pow(err_re2, 2) + pow(err_re3, 2) +
                         pow(err_re4, 2));
    double err_im1 = im_err * c.Re() / c.Mag2();
    double err_im2 =
        c.ReErr() * im * (pow(c.Im(), 2) - pow(c.Re(), 2)) / c.Mag2();
    double err_im3 = re_err * c.Im() / c.Mag2();
    double err_im4 =
        c.ImErr() * re * (pow(c.Re(), 2) - pow(c.Im(), 2)) / c.Mag2();
    double err_im = sqrt(pow(err_im1, 2) + pow(err_im2, 2) + pow(err_im3, 2) +
                         pow(err_im4, 2));

    return MComplex(
        (re * c.Re() + im * c.Im()) / (c.Re() * c.Re() + c.Im() * c.Im()),
        (im * c.Re() - re * c.Im()) / (c.Re() * c.Re() + c.Im() * c.Im()),
        err_re, err_im);
  }
  MComplex operator/(double c) {
    return MComplex(re / c, im / c, abs(re_err / c), abs(im_err / c));
  }
  MComplex operator+(MComplex c) {
    return MComplex(re + c.Re(), im + c.Im(),
                    sqrt(re_err * re_err + c.ReErr() * c.ReErr()),
                    sqrt(im_err * im_err + c.ImErr() * c.ImErr()));
  }
  MComplex operator-(MComplex c) {
    return MComplex(re - c.Re(), im - c.Im(),
                    sqrt(re_err * re_err + c.ReErr() * c.ReErr()),
                    sqrt(im_err * im_err + c.ImErr() * c.ImErr()));
  }

  MComplex &operator*=(MComplex c) {
    double re_temp = re;
    re = re * c.Re() - im * c.Im();
    im = re_temp * c.Im() + im * c.Re();
    re_err = sqrt(pow(c.Re() * re_err, 2) + pow(c.Im() * im_err, 2) +
                  pow(re_temp * c.ReErr(), 2) + pow(im * c.ImErr(), 2));
    im_err = sqrt(pow(c.Re() * im_err, 2) + pow(c.Im() * re_err, 2) +
                  pow(im * c.ReErr(), 2) + pow(re_temp * c.ImErr(), 2));
    return *this;
  }
  MComplex &operator*=(double c) {
    re *= c;
    im *= c;
    re_err = abs(re_err * c);
    im_err = abs(im_err * c);
    return *this;
  }
  MComplex &operator/=(MComplex c) {
    double err_re1 = re_err * c.Re() / c.Mag2();
    double err_re2 =
        c.ReErr() * re * (pow(c.Im(), 2) - pow(c.Re(), 2)) / c.Mag2();
    double err_re3 = im_err * c.Im() / c.Mag2();
    double err_re4 =
        c.ImErr() * im * (pow(c.Re(), 2) - pow(c.Im(), 2)) / c.Mag2();
    re_err = sqrt(pow(err_re1, 2) + pow(err_re2, 2) + pow(err_re3, 2) +
                         pow(err_re4, 2));
    double err_im1 = im_err * c.Re() / c.Mag2();
    double err_im2 =
        c.ReErr() * im * (pow(c.Im(), 2) - pow(c.Re(), 2)) / c.Mag2();
    double err_im3 = re_err * c.Im() / c.Mag2();
    double err_im4 =
        c.ImErr() * re * (pow(c.Re(), 2) - pow(c.Im(), 2)) / c.Mag2();
    im_err = sqrt(pow(err_im1, 2) + pow(err_im2, 2) + pow(err_im3, 2) +
                         pow(err_im4, 2));

    double re_temp = re;
    re = (re * c.Re() + im * c.Im()) / (c.Re() * c.Re() + c.Im() * c.Im());
    im = (im * c.Re() - re_temp * c.Im()) / (c.Re() * c.Re() + c.Im() * c.Im());
    return *this;
  }
  MComplex &operator/=(double c) {
    re /= c;
    im /= c;
    re_err = abs(re_err / c);
    im_err = abs(im_err / c);
    return *this;
  }
  MComplex &operator+=(MComplex c) {
    re += c.Re();
    im += c.Im();
    re_err = sqrt(re_err * re_err + c.ReErr() * c.ReErr());
    im_err = sqrt(im_err * im_err + c.ImErr() * c.ImErr());
    return *this;
  }
  MComplex &operator-=(MComplex c) {
    re -= c.Re();
    im -= c.Im();
    re_err = sqrt(re_err * re_err + c.ReErr() * c.ReErr());
    im_err = sqrt(im_err * im_err + c.ImErr() * c.ImErr());
    return *this;
  }

private:
  double re;
  double im;
  double re_err;
  double im_err;
};

#endif