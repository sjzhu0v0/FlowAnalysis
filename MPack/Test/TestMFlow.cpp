#include "MFlow.h"

void TestMFlow() {
  R__LOAD_LIBRARY(/data/Software/MPack/lib/libMFlow.so);
  MFGF::ArrayFlowInit();

  MFGF::fArrayFlow[MFGV::kHadron].Fill(0.1,1.);
  MFGF::fArrayFlow[MFGV::kHadron].Print();
  MFGF::ArrayFlowReset(MFGV::kHadron);
  MFGF::fArrayFlow[MFGV::kHadron].Print();

}