#include <string>
#include <vector>
#include "FastME_ProcPool.C"

void main_code(void){
  string Data_Path("VBF_Samples/Vectors/VectorVBF_QED4_QCD0_plus_VBF_QED2_QCD2_BKGs_2.root");  
  vector<string> MCs;
  MCs.push_back("VBF_Samples/Vectors/VectorVBF_QED4_QCD0_SIG_FastME_2.root");
  //MCs.push_back("VBF_Samples/Vectors/VectorVBF_QED4_QCD0_plus_VBF_QED2_QCD2_BKGs_2.root");
  MCs.push_back("VBF_Samples/Vectors/VectorVBF_QED4_QCD0_BKG_FastME_2.root");
  MCs.push_back("VBF_Samples/Vectors/VectorVBF_QED2_QCD2_BKG_FastME_2.root");
  
  ///Parameters (from left to right):
  ///1. Address to Data sample;
  ///2. Vector with address of MC samples;
  ///3. Number of Data events;
  ///4. Number of MC types (if all have different classification this is vector size);
  ///5. Number of final state particles;
  ///6. Name of the output file to store FastME analysis results.
  FastME_ProcPool(Data_Path,MCs,10000,3,6,"fme_results_bkg");
  
}
