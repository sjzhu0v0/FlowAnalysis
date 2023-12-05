#ifndef MROOTGRAPHIC_H
#define MROOTGRAPHIC_H

#include "MSystem.h"

namespace MColorSpace {
const int RGB_4_0[4][3] = {
    {223, 122, 94}, {60, 64, 91}, {130, 178, 154}, {242, 204, 142}}; // 2
const int RGB_5_0[5][3] = {{38, 70, 83},
                           {42, 157, 142},
                           {233, 196, 107},
                           {243, 162, 97},
                           {230, 111, 81}}; // 3
const int RGB_5_1[5][3] = {
    {1, 7, 19}, {1, 36, 76}, {1, 53, 101}, {255, 195, 0}, {251, 132, 2}}; // 14
const int RGB_6_0[6][3] = {{0, 0, 0},     {58, 9, 100},   {122, 27, 109},
                           {189, 55, 82}, {237, 104, 37}, {251, 180, 26}}; // 18
const int RGB_7_0[7][3] = {{144, 201, 231}, {33, 158, 188}, {19, 103, 131},
                           {2, 48, 74},     {254, 183, 5},  {255, 158, 2},
                           {250, 134, 0}}; // 5
const int RGB_7_1[7][3] = {{38, 70, 83},    {40, 114, 113},  {42, 157, 140},
                           {138, 176, 125}, {233, 196, 107}, {243, 162, 97},
                           {230, 111, 81}}; // 16
} // namespace MColorSpace

int GetColor(const int *rgb) {
  return TColor::GetColor(*rgb, *(rgb + 1), *(rgb + 2));
}

TLatex *ReturnLatex(double x, double y, const char *text, int color = 1,
                    double size = 0.05) {
  TLatex *latex = new TLatex(x, y, text);
  latex->SetNDC();
  latex->SetTextSize(size);
  latex->SetTextColor(color);
  latex->Draw();
  return latex;
}

void SetPaletteWidth(TH1 *h, double width) {
  TPaletteAxis *palette =
      (TPaletteAxis *)h->GetListOfFunctions()->FindObject("palette");
  double x_ndc_1 = palette->GetX1NDC();
  double x_ndc_2 = palette->GetX2NDC();
  double y_ndc_1 = palette->GetY1NDC();
  double y_ndc_2 = palette->GetY2NDC();
  palette->SetX1NDC(x_ndc_1 + width);
  palette->SetX2NDC(x_ndc_2 + width);
  palette->SetY1NDC(y_ndc_1);
  palette->SetY2NDC(y_ndc_2);
}

void SetPolttingStyle(TStyle *style = gStyle) {
  style->SetOptStat(0);
  style->SetPadTickX(1);
  style->SetPadTickY(1);

  style->SetPadLeftMargin(0.20);
  style->SetPadRightMargin(0.01);
  style->SetPadTopMargin(0.06);
  style->SetPadBottomMargin(0.20);

  style->SetPadBorderMode(0);
  style->SetCanvasBorderMode(0);
  style->SetFrameBorderMode(0);
  style->SetLegendBorderSize(0);
  style->SetLegendFillColor(0);
  style->SetLegendFont(42);
  style->SetLegendBorderSize(0);
  style->SetLegendFillColor(0);
  style->SetLegendTextSize(0.05);
  style->SetTitleSize(0.07, "XYZ");
  style->SetLabelSize(0.05, "XYZ");
  style->SetGridColor(16);
}

void SetH1LineStyle(TH1 *h, Color_t color = kRed, int width = 2,
                    double n_pad_x = 1) {
  h->SetLineColor(color);
  h->SetLineWidth(width);
  h->SetTitleSize(0.05 * n_pad_x, "XYZ");
  h->GetYaxis()->SetTitleOffset(0.9);
  h->GetXaxis()->SetTitleOffset(0.9);
  h->GetYaxis()->SetMaxDigits(3);
  h->GetXaxis()->SetNdivisions(504);
  h->GetYaxis()->SetNdivisions(504);
  h->SetLabelSize(0.05 * n_pad_x, "XYZ");
}

// class MCanvas {
// private:
//   TCanvas *fCanvas;
//   vector<TPad *> fVecPad;
//   vector<TPad *> fVecPadSub;

//   int fNPadX;
//   int fNPadY;

// public:
//   MCanvas(int width_x = 800, int width_y = 600, int n_pad_x = 1,
//           int n_pad_y = 1, bool is_continous_x = true,
//           bool is_continous_y = true, double margin_left = 0.2,
//           double margin_right = 0.05, double margin_top = 0.01,
//           double margin_bottom = 0.19);
//   TPad *GetPad(int i_pad);
//   TPad *GetPad(int i_pad_x, int i_pad_y);
//   ~MCanvas();
// };

// MCanvas::MCanvas(int width_x, int width_y, int n_pad_x, int n_pad_y,
//                  bool is_continous_x, bool is_continous_y, double
//                  margin_left, double margin_right, double margin_top, double
//                  margin_bottom)
//     : fNPadX(n_pad_x), fNPadY(n_pad_y) {
//   fCanvas = new TCanvas(Form("canvas_%d", GenerateUID()), "", width_x,
//   width_y); vector<double[4]> fVecPosNDCPad; // xlow,ylow,xup,yup
//   vector<double[4]> fMarginNDCPad; // left,right,top,bottom
//   fVecPosNDCPad.resize(n_pad_x * n_pad_y);
//   fMarginNDCPad.resize(n_pad_x * n_pad_y);

//   if (is_continous_x) {
//     double width =
//         (1 - margin_left / (double)n_pad_x - margin_right / (double)n_pad_x)
//         / (double)n_pad_x;
//     for (int i_pad_x = 0; i_pad_x < n_pad_x; i_pad_x++) {
//       double xlow = margin_left / (double)n_pad_x + i_pad_x * width;
//       double xup = xlow + width;
//       double left = 0;
//       double right = 0;
//       if (i_pad_x == 0) {
//         xlow = 0;
//         left = (margin_left / n_pad_x) / (margin_left / n_pad_x + width);
//       } else if (i_pad_x == n_pad_x - 1) {
//         xup = 1;
//         right = (margin_right / n_pad_x) / (margin_right / n_pad_x + width);
//       }
//       for (int i_pad_y = 0; i_pad_y < n_pad_y; i_pad_y++) {
//         fVecPosNDCPad[i_pad_x * n_pad_y + i_pad_y][0] = xlow;
//         fVecPosNDCPad[i_pad_x * n_pad_y + i_pad_y][2] = xup;
//         fMarginNDCPad[i_pad_x * n_pad_y + i_pad_y][0] = left;
//         fMarginNDCPad[i_pad_x * n_pad_y + i_pad_y][1] = right;
//       }
//     }
//   } else {
//     double width = 1. / (double)n_pad_x;
//     for (int i_pad_x = 0; i_pad_x < n_pad_x; i_pad_x++) {
//       double xlow = i_pad_x * width;
//       double xup = xlow + width;
//       double left = margin_left;
//       double right = margin_right;
//       for (int i_pad_y = 0; i_pad_y < n_pad_y; i_pad_y++) {
//         fVecPosNDCPad[i_pad_x * n_pad_y + i_pad_y][0] = xlow;
//         fVecPosNDCPad[i_pad_x * n_pad_y + i_pad_y][2] = xup;
//         fMarginNDCPad[i_pad_x * n_pad_y + i_pad_y][0] = left;
//         fMarginNDCPad[i_pad_x * n_pad_y + i_pad_y][1] = right;
//       }
//     }
//   }

//   if (is_continous_y) {
//     double width =
//         (1 - margin_bottom / (double)n_pad_y - margin_top / (double)n_pad_y)
//         / (double)n_pad_y;
//     for (int i_pad_y = 0; i_pad_y < n_pad_y; i_pad_y++) {
//       double ylow = margin_bottom / (double)n_pad_y + i_pad_y * width;
//       double yup = ylow + width;
//       double bottom = 0;
//       double top = 0;
//       if (i_pad_y == 0) {
//         ylow = 0;
//         bottom = (margin_bottom / n_pad_y) / (margin_bottom / n_pad_y +
//         width);
//       } else if (i_pad_y == n_pad_y - 1) {
//         yup = 1;
//         top = (margin_top / n_pad_y) / (margin_top / n_pad_y + width);
//       }
//       for (int i_pad_x = 0; i_pad_x < n_pad_x; i_pad_x++) {
//         fVecPosNDCPad[i_pad_x * n_pad_y + i_pad_y][1] = ylow;
//         fVecPosNDCPad[i_pad_x * n_pad_y + i_pad_y][3] = yup;
//         fMarginNDCPad[i_pad_x * n_pad_y + i_pad_y][2] = top;
//         fMarginNDCPad[i_pad_x * n_pad_y + i_pad_y][3] = bottom;
//       }
//     }
//   } else {
//     double width = 1. / (double)n_pad_y;
//     for (int i_pad_y = 0; i_pad_y < n_pad_y; i_pad_y++) {
//       double ylow = i_pad_y * width;
//       double yup = ylow + width;
//       double bottom = margin_bottom;
//       double top = margin_top;
//       for (int i_pad_x = 0; i_pad_x < n_pad_x; i_pad_x++) {
//         fVecPosNDCPad[i_pad_x * n_pad_y + i_pad_y][1] = ylow;
//         fVecPosNDCPad[i_pad_x * n_pad_y + i_pad_y][3] = yup;
//         fMarginNDCPad[i_pad_x * n_pad_y + i_pad_y][2] = top;
//         fMarginNDCPad[i_pad_x * n_pad_y + i_pad_y][3] = bottom;
//       }
//     }
//   }

//   fVecPad.resize(n_pad_x * n_pad_y);
//   fVecPadSub.resize(n_pad_x * n_pad_y);

//   for (int i_pad = 0; i_pad < n_pad_x * n_pad_y; i_pad++) {
//     fVecPad[i_pad] = new TPad(Form("pad_%d_%d", i_pad, GenerateUID()), "",
//                               fVecPosNDCPad[i_pad][0],
//                               fVecPosNDCPad[i_pad][1],
//                               fVecPosNDCPad[i_pad][2],
//                               fVecPosNDCPad[i_pad][3]);
//     fVecPad[i_pad]->SetLeftMargin(fMarginNDCPad[i_pad][0]);
//     fVecPad[i_pad]->SetRightMargin(fMarginNDCPad[i_pad][1]);
//     fVecPad[i_pad]->SetTopMargin(fMarginNDCPad[i_pad][2]);
//     fVecPad[i_pad]->SetBottomMargin(fMarginNDCPad[i_pad][3]);
//     fVecPad[i_pad]->Draw();
//   }
// }

// MCanvas::~MCanvas() {}

// TPad *MCanvas::(int i_pad) { return fVecPad[i_pad]; }

// TPad *MCanvas::(int i_pad_x, int i_pad_y) {
//   return fVecPad[i_pad_x * fNPadX + i_pad_y];
// }

#endif