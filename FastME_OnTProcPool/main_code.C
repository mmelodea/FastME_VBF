#include <string>
#include <vector>
#include "FastME.C"

void main_code(void){
  //string Data_Path("../VBF_Samples/Vectors/VectorVBF_QED4_QCD0_plus_VBF_QED2_QCD2_BKGs_1.root");
  string Data_Path("../VBF_Samples/Vectors/VectorVBF_QED4_QCD0_SIG_FastME_1.root"); 
 
  vector<string> MCs, MC_Names;
  MCs.push_back("../VBF_Samples/Vectors/VectorVBF_QED4_QCD0_SIG_FastME_2.root");
  //MCs.push_back("../VBF_Samples/Vectors/VectorVBF_QED4_QCD0_plus_VBF_QED2_QCD2_BKGs_2.root");
  MCs.push_back("../VBF_Samples/Vectors/VectorVBF_QED4_QCD0_BKG_FastME_2.root");
  MCs.push_back("../VBF_Samples/Vectors/VectorVBF_QED2_QCD2_BKG_FastME_2.root");
  
  MC_Names.push_back("Signal");
  MC_Names.push_back("EW Bkg");
  MC_Names.push_back("QCD Bkg");

  ///Parameters (from left to right):
  ///1. Address to Data sample;
  ///2. Vector with address of MC samples;
  ///3.	Name of tree containing	the events;
  ///4. Vector with names to each MC type;
  ///5. Number of final state particles;
  ///6. Number of cores to be used;
  ///7. Name of the output file to store FastME analysis results.
  FastME(Data_Path, MCs, "VBF", MC_Names, 6, 2, "fme_results_sig");
  
}
