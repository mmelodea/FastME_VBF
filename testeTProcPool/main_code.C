#include <string>
#include <vector>
#include "FastME_ProcPool.C"

void main_code(void){
  string Data_Path("VBF_Samples/Vectors/VectorVBF_QED4_QCD0_BKG_FastME_1.root");  
  vector<string> MCs;
  MCs.push_back("VBF_Samples/Vectors/VectorVBF_QED4_QCD0_SIG_FastME_2.root");
  MCs.push_back("VBF_Samples/Vectors/VectorVBF_QED4_QCD0_BKG_FastME_2.root");
  MCs.push_back("VBF_Samples/Vectors/VectorVBF_QED2_QCD2_BKG_FastME_2.root");
  
  FastME_ProcPool(Data_Path,MCs,10000,3,"fme_analysis_results_Bkg");
  
}
