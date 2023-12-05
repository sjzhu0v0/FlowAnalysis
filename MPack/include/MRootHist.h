#ifndef MROOTHIST_H
#define MROOTHIST_H

#include "MSystem.h"

// get the positions of 1-D histogram peaks
vector<double> GetHistPeak(TH1D *h = nullptr, int n_rebin = 1,
                           double limit_low = 0, int n_peak_max = 3) {
  vector<double> vec_peak;
  auto hist = (TH1D *)h->Clone();
  hist->Rebin(n_rebin);
  int n_bin = hist->GetNbinsX();
  for (int i = 1; i <= n_bin; i++) {
    if (hist->GetBinContent(i) > hist->GetBinContent(i - 1) &&
        hist->GetBinContent(i) > hist->GetBinContent(i + 1) &&
        hist->GetBinCenter(i) > limit_low) {
      vec_peak.push_back(hist->GetBinCenter(i));
    }
  }
  if (vec_peak.size() <= n_peak_max)
    return vec_peak;
  else {
    sort(vec_peak.begin(), vec_peak.end());
    vector<double> vec_peak_temp;
    for (int i = 0; i < n_peak_max; i++)
      vec_peak_temp.push_back(vec_peak[i]);
    return vec_peak_temp;
  }
  return vec_peak;
}

TH1 *RegisterHistogram(
    string string_registeration_hist =
        "{h2_aa}{h2_aa;pt{312}{};y{qwe}{112}{}}{{0,0.5,1}}{10,-1,1}") {
  int n_left_bracket = 0;
  int n_right_bracket = 0;
  for (auto c : string_registeration_hist) {
    if (c == '{') {
      n_left_bracket++;
    } else if (c == '}') {
      n_right_bracket++;
    }
  }

  if (n_left_bracket != n_right_bracket) {
    cerr << "Error: the number of { and } is not equal." << endl;
    return nullptr;
  }

  if (n_left_bracket < 3) {
    cerr << "Error: the number of { and } is less than 3." << endl;
    return nullptr;
  }

  vector<string> vec_string_control;

  int n_net_bracket = 0;
  string string_temp = "";
  for (auto c : string_registeration_hist) {
    if (n_net_bracket == 0 && c == '{') {
      n_net_bracket++;
      continue;
    }
    if (c == '{') {
      n_net_bracket++;
    } else if (c == '}') {
      n_net_bracket--;
    }
    if (n_net_bracket > 0) {
      string_temp += c;
    } else if (n_net_bracket < 0) {
      cerr << "Error: mismatching { and }" << endl;
    } else {
      vec_string_control.push_back(string_temp);
      string_temp = "";
    }
  }

  if (vec_string_control.size() < 3) {
    cerr << "Error: the number of control strings is less than 3." << endl;
    return nullptr;
  }

  string string_hist_name = vec_string_control[0];
  string string_hist_title = vec_string_control[1];
  int n_dimension = vec_string_control.size() - 2;

  vector<vector<double>> vec_bin;
  for (int i = 0; i < n_dimension; i++) {
    vector<double> vec_bin_temp;
    string string_bin = vec_string_control[i + 2];
    // print string_bin
    if (string_bin == "") {
      cerr << "Error: the bin string is empty." << endl;
      return nullptr;
    }
    // two type of bin string: {x1,x2,x3...} or n,xmin,xmax
    // if the type is {n,xmin,xmax}, translate it to the first type
    if (string_bin[0] == '{') {
      // judge whether the { and } is in pair and the most outside
      if (string_bin[0] != '{' || string_bin[string_bin.size() - 1] != '}') {
        cerr << "Error: the bin string is not in pair." << endl;
        return nullptr;
      }
      // remove the most outside { and }
      string_bin = string_bin.substr(1, string_bin.size() - 2);
      // find if there is any { or } in the string
      if (string_bin.find('{') != string::npos ||
          string_bin.find('}') != string::npos) {
        cerr << "Error: the bin string contains { or }." << endl;
        return nullptr;
      }
      // split the string by ,
      for (int j = 0; j < string_bin.size(); j++) {
        string string_temp = "";
        while (string_bin[j] != ',' && j < string_bin.size()) {
          string_temp += string_bin[j];
          j++;
        }
        // cout << string_temp << endl;
        vec_bin_temp.push_back(stod(string_temp));
      }
      vec_bin.push_back(vec_bin_temp);
    } else {
      // count the number of ,
      int n_comma = 0;
      for (auto c : string_bin) {
        if (c == ',') {
          n_comma++;
        }
      }
      if (n_comma != 2) {
        cerr << "Error: the bin string is not in the form of n,xmin,xmax."
             << endl;
        return nullptr;
      }
      // get the n,xmin,xmax
      string string_n = "";
      string string_xmin = "";
      string string_xmax = "";
      int n_comma_temp = 0;
      for (auto c : string_bin) {
        if (c == ',') {
          n_comma_temp++;
          continue;
        }
        if (n_comma_temp == 0) {
          string_n += c;
        } else if (n_comma_temp == 1) {
          string_xmin += c;
        } else if (n_comma_temp == 2) {
          string_xmax += c;
        }
      }
      // calculate the bin
      int n = stoi(string_n);
      double xmin = stod(string_xmin);
      double xmax = stod(string_xmax);
      double dx = (xmax - xmin) / n;
      for (int j = 0; j <= n; j++)
        vec_bin_temp.push_back(xmin + j * dx);
      vec_bin.push_back(vec_bin_temp);
    }
  }

  if (vec_bin.size() == 1) {
    TH1D *hist = new TH1D(string_hist_name.c_str(), string_hist_title.c_str(),
                          vec_bin[0].size() - 1, &vec_bin[0][0]);
    return hist;
  }
  if (vec_bin.size() == 2) {
    TH2D *hist = new TH2D(string_hist_name.c_str(), string_hist_title.c_str(),
                          vec_bin[0].size() - 1, &vec_bin[0][0],
                          vec_bin[1].size() - 1, &vec_bin[1][0]);
    return hist;
  }
  if (vec_bin.size() == 3) {
    TH3D *hist =
        new TH3D(string_hist_name.c_str(), string_hist_title.c_str(),
                 vec_bin[0].size() - 1, &vec_bin[0][0], vec_bin[1].size() - 1,
                 &vec_bin[1][0], vec_bin[2].size() - 1, &vec_bin[2][0]);
    return hist;
  }
  if (vec_bin.size() > 3) {
    cerr << "Error: the dimension of histogram is larger than 3." << endl;
    return nullptr;
  }
  return nullptr;
}

#endif