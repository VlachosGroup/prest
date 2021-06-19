#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <sstream>
#include <map>
#include <vector>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include "Classheader.h"
#include "AdditionalFunctions.h"
#include "StringRegistry.h"
bool checkCCRt_O2Adsorbm0(const Molecule& m0) {
  return true;
}
bool checkCCRt_O2Adsorbm1(const Molecule& m1) {
  return true;
}
bool checkCCRt_O2Adsorbmboth(const Molecule& m0, const Molecule& m1) {
  return true;
}
bool Rt_O2Adsorbprodcons(const Molecule& m0) {
  return true;
}
bool checkCCRt_CCScissionm0(const Molecule& m0) {
  return true;
}
bool checkCCRt_CCScissionm1(const Molecule& m1) {
  return true;
}
bool checkCCRt_CCScissionmboth(const Molecule& m0, const Molecule& m1) {
  return true;
}
bool Rt_CCScissionprodcons(const Molecule& m0) {
  return true;
}
bool checkCCRt_CHScissionm0(const Molecule& m0) {
  return true;
}
bool checkCCRt_CHScissionm1(const Molecule& m1) {
  return true;
}
bool checkCCRt_CHScissionmboth(const Molecule& m0, const Molecule& m1) {
  return true;
}
bool Rt_CHScissionprodcons(const Molecule& m0) {
  return true;
}
bool checkCCRt_OHScissionm0(const Molecule& m0) {
  return true;
}
bool checkCCRt_OHScissionm1(const Molecule& m1) {
  return true;
}
bool checkCCRt_OHScissionmboth(const Molecule& m0, const Molecule& m1) {
  return true;
}
bool Rt_OHScissionprodcons(const Molecule& m0) {
  return true;
}
bool checkCCRt_HTransferCtoOm0(const Molecule& m0) {
  return true;
}
bool checkCCRt_HTransferCtoOm1(const Molecule& m1) {
  return true;
}
bool checkCCRt_HTransferCtoOmboth(const Molecule& m0, const Molecule& m1) {
  return true;
}
bool Rt_HTransferCtoOprodcons(const Molecule& m0) {
  return true;
}
bool checkCCRt_HTransferOHtoOm0(const Molecule& m0) {
  return true;
}
bool checkCCRt_HTransferOHtoOm1(const Molecule& m1) {
  return true;
}
bool checkCCRt_HTransferOHtoOmboth(const Molecule& m0, const Molecule& m1) {
  return true;
}
bool Rt_HTransferOHtoOprodcons(const Molecule& m0) {
  return true;
}
bool checkCCRt_COFormationfromOHm0(const Molecule& m0) {
  return Patternmatch(m0, Substructure("O1(-C2)", patternsize("O1(-C2)")), 1).GetDistinctMatches()==0 && true;
}
bool checkCCRt_COFormationfromOHm1(const Molecule& m1) {
  return true;
}
bool checkCCRt_COFormationfromOHmboth(const Molecule& m0, const Molecule& m1) {
  return true;
}
bool Rt_COFormationfromOHprodcons(const Molecule& m0) {
  return true;
}
bool checkCCRt_COFormationfromOm0(const Molecule& m0) {
  return Patternmatch(m0, Substructure("O1(-C2)", patternsize("O1(-C2)")), 1).GetDistinctMatches()==0 && true;
}
bool checkCCRt_COFormationfromOm1(const Molecule& m1) {
  return true;
}
bool checkCCRt_COFormationfromOmboth(const Molecule& m0, const Molecule& m1) {
  return true;
}
bool Rt_COFormationfromOprodcons(const Molecule& m0) {
  return true;
}
bool checkCCRt_H2OwCFormationm0(const Molecule& m0) {
  return true;
}
bool checkCCRt_H2OwCFormationm1(const Molecule& m1) {
  return true;
}
bool checkCCRt_H2OwCFormationmboth(const Molecule& m0, const Molecule& m1) {
  return true;
}
bool Rt_H2OwCFormationprodcons(const Molecule& m0) {
  return true;
}
bool checkCCRt_H2OwOFormationm0(const Molecule& m0) {
  return true;
}
bool checkCCRt_H2OwOFormationm1(const Molecule& m1) {
  return true;
}
bool checkCCRt_H2OwOFormationmboth(const Molecule& m0, const Molecule& m1) {
  return true;
}
bool Rt_H2OwOFormationprodcons(const Molecule& m0) {
  return true;
}
bool correctionCheck0(const Molecule& m0) {
  return true;
}
bool correctionCheck1(const Molecule& m0) {
  return true;
}
bool correctionCheck2(const Molecule& m0) {
  return true;
}
bool correctionCheck3(const Molecule& m0) {
  return true;
}
bool correctionCheck4(const Molecule& m0) {
  return true;
}
bool correctionCheck5(const Molecule& m0) {
  return true;
}
bool correctionCheck6(const Molecule& m0) {
  return true;
}
bool correctionCheck7(const Molecule& m0) {
  return true;
}
bool correctionCheck8(const Molecule& m0) {
  return true;
}
bool correctionCheck9(const Molecule& m0) {
  return true;
}
bool correctionCheck10(const Molecule& m0) {
  return true;
}
bool correctionCheck11(const Molecule& m0) {
  return true;
}
bool correctionCheck12(const Molecule& m0) {
  return true;
}
bool correctionCheck13(const Molecule& m0) {
  return true;
}
bool correctionCheck14(const Molecule& m0) {
  return true;
}
bool correctionCheck15(const Molecule& m0) {
  return true;
}
bool correctionCheck16(const Molecule& m0) {
  return true;
}
bool correctionCheck17(const Molecule& m0) {
  return true;
}
bool correctionCheck18(const Molecule& m0) {
  return true;
}
bool correctionCheck19(const Molecule& m0) {
  return true;
}
bool correctionCheck20(const Molecule& m0) {
  return true;
}
bool correctionCheck21(const Molecule& m0) {
  return true;
}
bool correctionCheck22(const Molecule& m0) {
  return true;
}
bool correctionCheck23(const Molecule& m0) {
  return true;
}
bool correctionCheck24(const Molecule& m0) {
  return true;
}
bool correctionCheck25(const Molecule& m0) {
  return true;
}
bool correctionCheck26(const Molecule& m0) {
  return true;
}
bool correctionCheck27(const Molecule& m0) {
  return true;
}
bool correctionCheck28(const Molecule& m0) {
  return true;
}
bool correctionCheck29(const Molecule& m0) {
  return true;
}
bool correctionCheck30(const Molecule& m0) {
  return true;
}
bool correctionCheck31(const Molecule& m0) {
  return true;
}
bool correctionCheck32(const Molecule& m0) {
  return true;
}
bool correctionCheck33(const Molecule& m0) {
  return true;
}
bool correctionCheck34(const Molecule& m0) {
  return true;
}
bool correctionCheck35(const Molecule& m0) {
  return true;
}
bool correctionCheck36(const Molecule& m0) {
  return true;
}
bool correctionCheck37(const Molecule& m0) {
  return true;
}
bool correctionCheck38(const Molecule& m0) {
  return true;
}
bool correctionCheck39(const Molecule& m0) {
  return true;
}
bool correctionCheck40(const Molecule& m0) {
  return true;
}
bool correctionCheck41(const Molecule& m0) {
  return true;
}
bool correctionCheck42(const Molecule& m0) {
  return true;
}
bool correctionCheck43(const Molecule& m0) {
  return true;
}
bool correctionCheck44(const Molecule& m0) {
  return true;
}
bool correctionCheck45(const Molecule& m0) {
  return true;
}
bool correctionCheck46(const Molecule& m0) {
  return true;
}
bool correctionCheck47(const Molecule& m0) {
  return true;
}
bool correctionCheck48(const Molecule& m0) {
  return true;
}
bool correctionCheck49(const Molecule& m0) {
  return true;
}
bool correctionCheck50(const Molecule& m0) {
  return true;
}
bool correctionCheck51(const Molecule& m0) {
  return true;
}
bool correctionCheck52(const Molecule& m0) {
  return true;
}
bool checkCglobalConstraints(const Molecule& m0) {
  return true;
}
bool DummyFunctionForThermoCorrelation(const Molecule& m0, double& EnthalpyA, double& EnthalpyB, double& EntropyA, double& EntropyB) { return true;}
ThermoCorrelationPtr ThermoGA::CorrelationPtr = &DummyFunctionForThermoCorrelation;
 bool ThermoGA::HasCorrelationsForThermo = false;

int main()
{
   char *programName = "combustion.txt";
   list<string> reactantlist;
  vector<string> CompositeAtomsList;
  vector<pair<string,SiteType> > CompositeSiteList;
 ifstream reactantfile("reactants.txt");
   string line;
   if (reactantfile.is_open())
   {
      while (! reactantfile.eof() )
      {
         getline (reactantfile,line);
         if (line!="")
            reactantlist.push_back(line);
      }
      reactantfile.close();
   }
   else cout << "Unable to open reactants.txt, so skipping reactant input from file"<<endl;
vector<Reactiontype> Rtypelist;
vector<KineticParamPtr> kineticFns(10);
   reactantlist.push_back("CCC");
   reactantlist.push_back("O=O");
   reactantlist.push_back("{Pt}");
CompositeAtomsList.push_back("{Pt}");

  Reactiontype Rt_O2Adsorb;
  Rt_O2Adsorb.add_reactant_pattern(Substructure("O1(=O2)", patternsize("O1(=O2)")));
  Rt_O2Adsorb.add_reactant_pattern(Substructure("{Pt}3[!|~$|>0]", patternsize("{Pt}3[!|~$|>0]")));
   map<int, int> Rt_O2Adsorbdupmap2;
  Rt_O2Adsorbdupmap2.insert(pair<int, int>(3, 4));
   Rt_O2Adsorb.AddReactantCopy(1, Rt_O2Adsorbdupmap2 );
   map<int, int> Rt_O2Adsorbdupmap3;
  Rt_O2Adsorbdupmap3.insert(pair<int, int>(3, 5));
   Rt_O2Adsorb.AddReactantCopy(1, Rt_O2Adsorbdupmap3 );
   map<int, int> Rt_O2Adsorbdupmap4;
  Rt_O2Adsorbdupmap4.insert(pair<int, int>(3, 6));
   Rt_O2Adsorb.AddReactantCopy(1, Rt_O2Adsorbdupmap4 );
   Rt_O2Adsorb.add_reactantconstraint(&checkCCRt_O2Adsorbm0);
   Rt_O2Adsorb.add_reactantconstraint(&checkCCRt_O2Adsorbm1);
   Rt_O2Adsorb.add_combined_constraint(&checkCCRt_O2Adsorbmboth);
   Rt_O2Adsorb.disconnect_bond(1,2);
   Rt_O2Adsorb.connect_bond(1,3);
   Rt_O2Adsorb.connect_bond(1,4);
   Rt_O2Adsorb.connect_bond(2,5);
   Rt_O2Adsorb.connect_bond(2,6);
  Rt_O2Adsorb.add_productconstraint(&Rt_O2Adsorbprodcons);
  Rt_O2Adsorb.setRuleName("O2Adsorb");
  Rtypelist.push_back(Rt_O2Adsorb); 

  Reactiontype Rt_CCScission;
  Rt_CCScission.add_reactant_pattern(Substructure("C1(-C2)", patternsize("C1(-C2)")));
  Rt_CCScission.add_reactant_pattern(Substructure("{Pt}3[!|~$|>0]", patternsize("{Pt}3[!|~$|>0]")));
   map<int, int> Rt_CCScissiondupmap2;
  Rt_CCScissiondupmap2.insert(pair<int, int>(3, 4));
   Rt_CCScission.AddReactantCopy(1, Rt_CCScissiondupmap2 );
   Rt_CCScission.add_reactantconstraint(&checkCCRt_CCScissionm0);
   Rt_CCScission.add_reactantconstraint(&checkCCRt_CCScissionm1);
   Rt_CCScission.add_combined_constraint(&checkCCRt_CCScissionmboth);
   Rt_CCScission.disconnect_bond(1,2);
   Rt_CCScission.connect_bond(1,3);
   Rt_CCScission.connect_bond(2,4);
  Rt_CCScission.add_productconstraint(&Rt_CCScissionprodcons);
  Rt_CCScission.setRuleName("CCScission");
  Rtypelist.push_back(Rt_CCScission); 

  Reactiontype Rt_CHScission;
  Rt_CHScission.add_reactant_pattern(Substructure("C1(-H2)", patternsize("C1(-H2)")));
  Rt_CHScission.add_reactant_pattern(Substructure("{Pt}3[!|~$|>0]", patternsize("{Pt}3[!|~$|>0]")));
   map<int, int> Rt_CHScissiondupmap2;
  Rt_CHScissiondupmap2.insert(pair<int, int>(3, 4));
   Rt_CHScission.AddReactantCopy(1, Rt_CHScissiondupmap2 );
   Rt_CHScission.add_reactantconstraint(&checkCCRt_CHScissionm0);
   Rt_CHScission.add_reactantconstraint(&checkCCRt_CHScissionm1);
   Rt_CHScission.add_combined_constraint(&checkCCRt_CHScissionmboth);
   Rt_CHScission.disconnect_bond(1,2);
   Rt_CHScission.connect_bond(1,3);
   Rt_CHScission.connect_bond(2,4);
  Rt_CHScission.add_productconstraint(&Rt_CHScissionprodcons);
  Rt_CHScission.setRuleName("CHScission");
  Rtypelist.push_back(Rt_CHScission); 

  Reactiontype Rt_OHScission;
  Rt_OHScission.add_reactant_pattern(Substructure("O1(-H2)", patternsize("O1(-H2)")));
  Rt_OHScission.add_reactant_pattern(Substructure("{Pt}3[!|~$|>0]", patternsize("{Pt}3[!|~$|>0]")));
   map<int, int> Rt_OHScissiondupmap2;
  Rt_OHScissiondupmap2.insert(pair<int, int>(3, 4));
   Rt_OHScission.AddReactantCopy(1, Rt_OHScissiondupmap2 );
   Rt_OHScission.add_reactantconstraint(&checkCCRt_OHScissionm0);
   Rt_OHScission.add_reactantconstraint(&checkCCRt_OHScissionm1);
   Rt_OHScission.add_combined_constraint(&checkCCRt_OHScissionmboth);
   Rt_OHScission.disconnect_bond(1,2);
   Rt_OHScission.connect_bond(1,3);
   Rt_OHScission.connect_bond(2,4);
  Rt_OHScission.add_productconstraint(&Rt_OHScissionprodcons);
  Rt_OHScission.setRuleName("OHScission");
  Rtypelist.push_back(Rt_OHScission); 

  Reactiontype Rt_HTransferCtoO;
  Rt_HTransferCtoO.add_reactant_pattern(Substructure("C1(-H2)", patternsize("C1(-H2)")));
  Rt_HTransferCtoO.add_reactant_pattern(Substructure("{Pt}4(-O3(-{Pt}5))", patternsize("{Pt}4(-O3(-{Pt}5))")));
   Rt_HTransferCtoO.add_reactantconstraint(&checkCCRt_HTransferCtoOm0);
   Rt_HTransferCtoO.add_reactantconstraint(&checkCCRt_HTransferCtoOm1);
   Rt_HTransferCtoO.add_combined_constraint(&checkCCRt_HTransferCtoOmboth);
   Rt_HTransferCtoO.disconnect_bond(1,2);
   Rt_HTransferCtoO.disconnect_bond(3,4);
   Rt_HTransferCtoO.connect_bond(1,4);
   Rt_HTransferCtoO.connect_bond(3,2);
  Rt_HTransferCtoO.add_productconstraint(&Rt_HTransferCtoOprodcons);
  Rt_HTransferCtoO.setRuleName("HTransferCtoO");
  Rtypelist.push_back(Rt_HTransferCtoO); 

  Reactiontype Rt_HTransferOHtoO;
  Rt_HTransferOHtoO.add_reactant_pattern(Substructure("O1(-H2)", patternsize("O1(-H2)")));
  Rt_HTransferOHtoO.add_reactant_pattern(Substructure("{Pt}4(-O3(-{Pt}5))", patternsize("{Pt}4(-O3(-{Pt}5))")));
   Rt_HTransferOHtoO.add_reactantconstraint(&checkCCRt_HTransferOHtoOm0);
   Rt_HTransferOHtoO.add_reactantconstraint(&checkCCRt_HTransferOHtoOm1);
   Rt_HTransferOHtoO.add_combined_constraint(&checkCCRt_HTransferOHtoOmboth);
   Rt_HTransferOHtoO.disconnect_bond(1,2);
   Rt_HTransferOHtoO.disconnect_bond(3,4);
   Rt_HTransferOHtoO.connect_bond(1,4);
   Rt_HTransferOHtoO.connect_bond(3,2);
  Rt_HTransferOHtoO.add_productconstraint(&Rt_HTransferOHtoOprodcons);
  Rt_HTransferOHtoO.setRuleName("HTransferOHtoO");
  Rtypelist.push_back(Rt_HTransferOHtoO); 

  Reactiontype Rt_COFormationfromOH;
  Rt_COFormationfromOH.add_reactant_pattern(Substructure("{Pt}2(-C1[!|~C|>1])", patternsize("{Pt}2(-C1[!|~C|>1])")));
  Rt_COFormationfromOH.add_reactant_pattern(Substructure("{Pt}5(-O3(-H4))", patternsize("{Pt}5(-O3(-H4))")));
   Rt_COFormationfromOH.add_reactantconstraint(&checkCCRt_COFormationfromOHm0);
   Rt_COFormationfromOH.add_reactantconstraint(&checkCCRt_COFormationfromOHm1);
   Rt_COFormationfromOH.add_combined_constraint(&checkCCRt_COFormationfromOHmboth);
   Rt_COFormationfromOH.disconnect_bond(3,5);
   Rt_COFormationfromOH.disconnect_bond(1,2);
   Rt_COFormationfromOH.connect_bond(1,3);
  Rt_COFormationfromOH.add_productconstraint(&Rt_COFormationfromOHprodcons);
  Rt_COFormationfromOH.setRuleName("COFormationfromOH");
  Rtypelist.push_back(Rt_COFormationfromOH); 

  Reactiontype Rt_COFormationfromO;
  Rt_COFormationfromO.add_reactant_pattern(Substructure("{Pt}2(-C1[!|~C|>1])", patternsize("{Pt}2(-C1[!|~C|>1])")));
  Rt_COFormationfromO.add_reactant_pattern(Substructure("{Pt}4(-O3(-{Pt}5))", patternsize("{Pt}4(-O3(-{Pt}5))")));
   Rt_COFormationfromO.add_reactantconstraint(&checkCCRt_COFormationfromOm0);
   Rt_COFormationfromO.add_reactantconstraint(&checkCCRt_COFormationfromOm1);
   Rt_COFormationfromO.add_combined_constraint(&checkCCRt_COFormationfromOmboth);
   Rt_COFormationfromO.disconnect_bond(3,4);
   Rt_COFormationfromO.disconnect_bond(1,2);
   Rt_COFormationfromO.connect_bond(1,3);
  Rt_COFormationfromO.add_productconstraint(&Rt_COFormationfromOprodcons);
  Rt_COFormationfromO.setRuleName("COFormationfromO");
  Rtypelist.push_back(Rt_COFormationfromO); 

  Reactiontype Rt_H2OwCFormation;
  Rt_H2OwCFormation.add_reactant_pattern(Substructure("C1(-H2)", patternsize("C1(-H2)")));
  Rt_H2OwCFormation.add_reactant_pattern(Substructure("{Pt}5(-O3(-H4))", patternsize("{Pt}5(-O3(-H4))")));
   Rt_H2OwCFormation.add_reactantconstraint(&checkCCRt_H2OwCFormationm0);
   Rt_H2OwCFormation.add_reactantconstraint(&checkCCRt_H2OwCFormationm1);
   Rt_H2OwCFormation.add_combined_constraint(&checkCCRt_H2OwCFormationmboth);
   Rt_H2OwCFormation.disconnect_bond(5,3);
   Rt_H2OwCFormation.disconnect_bond(2,1);
   Rt_H2OwCFormation.connect_bond(3,2);
   Rt_H2OwCFormation.connect_bond(1,5);
  Rt_H2OwCFormation.add_productconstraint(&Rt_H2OwCFormationprodcons);
  Rt_H2OwCFormation.setRuleName("H2OwCFormation");
  Rtypelist.push_back(Rt_H2OwCFormation); 

  Reactiontype Rt_H2OwOFormation;
  Rt_H2OwOFormation.add_reactant_pattern(Substructure("O1(-H2)", patternsize("O1(-H2)")));
  Rt_H2OwOFormation.add_reactant_pattern(Substructure("{Pt}5(-O3(-H4))", patternsize("{Pt}5(-O3(-H4))")));
   Rt_H2OwOFormation.add_reactantconstraint(&checkCCRt_H2OwOFormationm0);
   Rt_H2OwOFormation.add_reactantconstraint(&checkCCRt_H2OwOFormationm1);
   Rt_H2OwOFormation.add_combined_constraint(&checkCCRt_H2OwOFormationmboth);
   Rt_H2OwOFormation.disconnect_bond(5,3);
   Rt_H2OwOFormation.disconnect_bond(2,1);
   Rt_H2OwOFormation.connect_bond(3,2);
   Rt_H2OwOFormation.connect_bond(1,5);
  Rt_H2OwOFormation.add_productconstraint(&Rt_H2OwOFormationprodcons);
  Rt_H2OwOFormation.setRuleName("H2OwOFormation");
  Rtypelist.push_back(Rt_H2OwOFormation); 
   LumpingStrategy L(false);
  { GAdata gad; gad.setEnthalpy(-5.68* 4.18); gad.setEntropy(1.09* 4.18); gad.setCp(100,0.97* 4.18); gad.setCp(200,3.12* 4.18); gad.setCp(298,4.63* 4.18); gad.setCp(400,5.93* 4.18); gad.setCp(500,6.96* 4.18); gad.setCp(600,7.77* 4.18); gad.setCp(700,8.42* 4.18); gad.setCp(800,8.95* 4.18); gad.setCp(900,9.41* 4.18); gad.setCp(1000,9.8* 4.18); gad.setCp(1100,10.13* 4.18); gad.setCp(1200,10.43* 4.18); gad.setCp(1300,10.68* 4.18); gad.setCp(1400,10.9* 4.18); gad.setCp(1500,11.1* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-O6)(-C3(=O4))(-C2)(-H5)", gad); }
  { GAdata gad; gad.setEnthalpy(-3.94* 4.18); gad.setEntropy(1.3* 4.18); gad.setCp(100,0.27* 4.18); gad.setCp(200,1.75* 4.18); gad.setCp(298,3.52* 4.18); gad.setCp(400,5.07* 4.18); gad.setCp(500,6.26* 4.18); gad.setCp(600,7.15* 4.18); gad.setCp(700,7.84* 4.18); gad.setCp(800,8.39* 4.18); gad.setCp(900,8.84* 4.18); gad.setCp(1000,9.22* 4.18); gad.setCp(1100,9.54* 4.18); gad.setCp(1200,9.81* 4.18); gad.setCp(1300,10.04* 4.18); gad.setCp(1400,10.24* 4.18); gad.setCp(1500,10.41* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}6)(-C3(=O4))(-C2)(-H5)", gad); }
  { GAdata gad; gad.setEnthalpy(-7.46* 4.18); gad.setEntropy(6.01* 4.18); gad.setCp(100,2.55* 4.18); gad.setCp(200,3.78* 4.18); gad.setCp(298,5.14* 4.18); gad.setCp(400,6.71* 4.18); gad.setCp(500,8.15* 4.18); gad.setCp(600,9.36* 4.18); gad.setCp(700,10.39* 4.18); gad.setCp(800,11.27* 4.18); gad.setCp(900,12.03* 4.18); gad.setCp(1000,12.69* 4.18); gad.setCp(1100,13.26* 4.18); gad.setCp(1200,13.76* 4.18); gad.setCp(1300,14.19* 4.18); gad.setCp(1400,14.56* 4.18); gad.setCp(1500,14.89* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-C3(=O4))(-C2)(-H6)(-H5)", gad); }
  { GAdata gad; gad.setEnthalpy(-13.74* 4.18); gad.setEntropy(-5.44* 4.18); gad.setCp(100,-1.78* 4.18); gad.setCp(200,0.78* 4.18); gad.setCp(298,2.86* 4.18); gad.setCp(400,4.2* 4.18); gad.setCp(500,5.01* 4.18); gad.setCp(600,5.51* 4.18); gad.setCp(700,5.83* 4.18); gad.setCp(800,6.04* 4.18); gad.setCp(900,6.2* 4.18); gad.setCp(1000,6.32* 4.18); gad.setCp(1100,6.42* 4.18); gad.setCp(1200,6.49* 4.18); gad.setCp(1300,6.55* 4.18); gad.setCp(1400,6.6* 4.18); gad.setCp(1500,6.65* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}6)(-O5)(-C3(=O4))(-C2)", gad); }
  { GAdata gad; gad.setEnthalpy(-0.77* 4.18); gad.setEntropy(-0.19* 4.18); gad.setCp(100,-0.86* 4.18); gad.setCp(200,0.95* 4.18); gad.setCp(298,2.43* 4.18); gad.setCp(400,3.41* 4.18); gad.setCp(500,4.05* 4.18); gad.setCp(600,4.46* 4.18); gad.setCp(700,4.73* 4.18); gad.setCp(800,4.92* 4.18); gad.setCp(900,5.05* 4.18); gad.setCp(1000,5.15* 4.18); gad.setCp(1100,5.24* 4.18); gad.setCp(1200,5.31* 4.18); gad.setCp(1300,5.37* 4.18); gad.setCp(1400,5.42* 4.18); gad.setCp(1500,5.47* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}6)(-{Pt}5)(-C3(=O4))(-C2)", gad); }
  { GAdata gad; gad.setEnthalpy(-19.77* 4.18); gad.setEntropy(-0.28* 4.18); gad.setCp(100,-0.92* 4.18); gad.setCp(200,1.47* 4.18); gad.setCp(298,3.64* 4.18); gad.setCp(400,5.3* 4.18); gad.setCp(500,6.47* 4.18); gad.setCp(600,7.31* 4.18); gad.setCp(700,7.94* 4.18); gad.setCp(800,8.44* 4.18); gad.setCp(900,8.85* 4.18); gad.setCp(1000,9.2* 4.18); gad.setCp(1100,9.5* 4.18); gad.setCp(1200,9.76* 4.18); gad.setCp(1300,9.98* 4.18); gad.setCp(1400,10.18* 4.18); gad.setCp(1500,10.35* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-O4)(-C2)(-H3)", gad); }
  { GAdata gad; gad.setEnthalpy(-3.44* 4.18); gad.setEntropy(2.44* 4.18); gad.setCp(100,-0.68* 4.18); gad.setCp(200,1.48* 4.18); gad.setCp(298,3.54* 4.18); gad.setCp(400,5.05* 4.18); gad.setCp(500,6.11* 4.18); gad.setCp(600,6.85* 4.18); gad.setCp(700,7.41* 4.18); gad.setCp(800,7.85* 4.18); gad.setCp(900,8.21* 4.18); gad.setCp(1000,8.52* 4.18); gad.setCp(1100,8.78* 4.18); gad.setCp(1200,9* 4.18); gad.setCp(1300,9.2* 4.18); gad.setCp(1400,9.37* 4.18); gad.setCp(1500,9.52* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-{Pt}4)(-C2)(-H3)", gad); }
  { GAdata gad; gad.setEnthalpy(-15.18* 4.18); gad.setEntropy(7.42* 4.18); gad.setCp(100,2.47* 4.18); gad.setCp(200,4.12* 4.18); gad.setCp(298,5.58* 4.18); gad.setCp(400,7.14* 4.18); gad.setCp(500,8.53* 4.18); gad.setCp(600,9.69* 4.18); gad.setCp(700,10.66* 4.18); gad.setCp(800,11.49* 4.18); gad.setCp(900,12.21* 4.18); gad.setCp(1000,12.83* 4.18); gad.setCp(1100,13.37* 4.18); gad.setCp(1200,13.84* 4.18); gad.setCp(1300,14.25* 4.18); gad.setCp(1400,14.61* 4.18); gad.setCp(1500,14.92* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-O5)(-C2)(-H4)(-H3)", gad); }
  { GAdata gad; gad.setEnthalpy(-10.52* 4.18); gad.setEntropy(5.81* 4.18); gad.setCp(100,1.04* 4.18); gad.setCp(200,2.89* 4.18); gad.setCp(298,4.94* 4.18); gad.setCp(400,6.78* 4.18); gad.setCp(500,8.23* 4.18); gad.setCp(600,9.37* 4.18); gad.setCp(700,10.28* 4.18); gad.setCp(800,11.05* 4.18); gad.setCp(900,11.7* 4.18); gad.setCp(1000,12.27* 4.18); gad.setCp(1100,12.76* 4.18); gad.setCp(1200,13.19* 4.18); gad.setCp(1300,13.56* 4.18); gad.setCp(1400,13.89* 4.18); gad.setCp(1500,14.17* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-C2)(-H4)(-H3)", gad); }
  { GAdata gad; gad.setEnthalpy(-18.42* 4.18); gad.setEntropy(12.1* 4.18); gad.setCp(100,3.9* 4.18); gad.setCp(200,5.71* 4.18); gad.setCp(298,7.14* 4.18); gad.setCp(400,8.79* 4.18); gad.setCp(500,10.37* 4.18); gad.setCp(600,11.77* 4.18); gad.setCp(700,12.99* 4.18); gad.setCp(800,14.06* 4.18); gad.setCp(900,15.01* 4.18); gad.setCp(1000,15.84* 4.18); gad.setCp(1100,16.57* 4.18); gad.setCp(1200,17.22* 4.18); gad.setCp(1300,17.78* 4.18); gad.setCp(1400,18.28* 4.18); gad.setCp(1500,18.71* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-C2)(-H5)(-H4)(-H3)", gad); }
  { GAdata gad; gad.setEnthalpy(-18.48* 4.18); gad.setEntropy(-1.34* 4.18); gad.setCp(100,-1.82* 4.18); gad.setCp(200,0.59* 4.18); gad.setCp(298,2.49* 4.18); gad.setCp(400,3.65* 4.18); gad.setCp(500,4.34* 4.18); gad.setCp(600,4.75* 4.18); gad.setCp(700,5.01* 4.18); gad.setCp(800,5.18* 4.18); gad.setCp(900,5.31* 4.18); gad.setCp(1000,5.4* 4.18); gad.setCp(1100,5.48* 4.18); gad.setCp(1200,5.54* 4.18); gad.setCp(1300,5.59* 4.18); gad.setCp(1400,5.63* 4.18); gad.setCp(1500,5.67* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-{Pt}4)(-O3)(-C2)", gad); }
  { GAdata gad; gad.setEnthalpy(2.85* 4.18); gad.setEntropy(0.54* 4.18); gad.setCp(100,-1.15* 4.18); gad.setCp(200,1.38* 4.18); gad.setCp(298,2.63* 4.18); gad.setCp(400,3.33* 4.18); gad.setCp(500,3.77* 4.18); gad.setCp(600,4.05* 4.18); gad.setCp(700,4.23* 4.18); gad.setCp(800,4.36* 4.18); gad.setCp(900,4.45* 4.18); gad.setCp(1000,4.52* 4.18); gad.setCp(1100,4.58* 4.18); gad.setCp(1200,4.62* 4.18); gad.setCp(1300,4.66* 4.18); gad.setCp(1400,4.69* 4.18); gad.setCp(1500,4.72* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-{Pt}4)(-{Pt}3)(-C2)", gad); }
  { GAdata gad; gad.setEnthalpy(-11.83* 4.18); gad.setEntropy(3.26* 4.18); gad.setCp(100,1.2* 4.18); gad.setCp(200,2.9* 4.18); gad.setCp(298,4.45* 4.18); gad.setCp(400,5.81* 4.18); gad.setCp(500,6.87* 4.18); gad.setCp(600,7.69* 4.18); gad.setCp(700,8.33* 4.18); gad.setCp(800,8.85* 4.18); gad.setCp(900,9.3* 4.18); gad.setCp(1000,9.68* 4.18); gad.setCp(1100,10.02* 4.18); gad.setCp(1200,10.31* 4.18); gad.setCp(1300,10.57* 4.18); gad.setCp(1400,10.79* 4.18); gad.setCp(1500,10.99* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-O5)(-C3)(-C2)(-H4)", gad); }
  { GAdata gad; gad.setEnthalpy(-8.39* 4.18); gad.setEntropy(1.6* 4.18); gad.setCp(100,0.19* 4.18); gad.setCp(200,1.6* 4.18); gad.setCp(298,3.42* 4.18); gad.setCp(400,5.01* 4.18); gad.setCp(500,6.21* 4.18); gad.setCp(600,7.11* 4.18); gad.setCp(700,7.8* 4.18); gad.setCp(800,8.35* 4.18); gad.setCp(900,8.8* 4.18); gad.setCp(1000,9.17* 4.18); gad.setCp(1100,9.49* 4.18); gad.setCp(1200,9.76* 4.18); gad.setCp(1300,9.99* 4.18); gad.setCp(1400,10.19* 4.18); gad.setCp(1500,10.36* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-C3)(-C2)(-H4)", gad); }
  { GAdata gad; gad.setEnthalpy(-9.79* 4.18); gad.setEntropy(8.01* 4.18); gad.setCp(100,3* 4.18); gad.setCp(200,4.29* 4.18); gad.setCp(298,5.49* 4.18); gad.setCp(400,6.92* 4.18); gad.setCp(500,8.29* 4.18); gad.setCp(600,9.47* 4.18); gad.setCp(700,10.47* 4.18); gad.setCp(800,11.34* 4.18); gad.setCp(900,12.08* 4.18); gad.setCp(1000,12.73* 4.18); gad.setCp(1100,13.29* 4.18); gad.setCp(1200,13.78* 4.18); gad.setCp(1300,14.2* 4.18); gad.setCp(1400,14.57* 4.18); gad.setCp(1500,14.89* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-C3)(-C2)(-H5)(-H4)", gad); }
  { GAdata gad; gad.setEnthalpy(-14.1* 4.18); gad.setEntropy(-3.97* 4.18); gad.setCp(100,-1.63* 4.18); gad.setCp(200,0.87* 4.18); gad.setCp(298,2.95* 4.18); gad.setCp(400,4.27* 4.18); gad.setCp(500,5.05* 4.18); gad.setCp(600,5.51* 4.18); gad.setCp(700,5.79* 4.18); gad.setCp(800,5.98* 4.18); gad.setCp(900,6.12* 4.18); gad.setCp(1000,6.22* 4.18); gad.setCp(1100,6.3* 4.18); gad.setCp(1200,6.37* 4.18); gad.setCp(1300,6.43* 4.18); gad.setCp(1400,6.48* 4.18); gad.setCp(1500,6.53* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-O4)(-C3)(-C2)", gad); }
  { GAdata gad; gad.setEnthalpy(-2.11* 4.18); gad.setEntropy(0.89* 4.18); gad.setCp(100,-0.62* 4.18); gad.setCp(200,1.2* 4.18); gad.setCp(298,2.67* 4.18); gad.setCp(400,3.63* 4.18); gad.setCp(500,4.24* 4.18); gad.setCp(600,4.63* 4.18); gad.setCp(700,4.87* 4.18); gad.setCp(800,5.04* 4.18); gad.setCp(900,5.16* 4.18); gad.setCp(1000,5.25* 4.18); gad.setCp(1100,5.32* 4.18); gad.setCp(1200,5.38* 4.18); gad.setCp(1300,5.43* 4.18); gad.setCp(1400,5.47* 4.18); gad.setCp(1500,5.51* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-{Pt}4)(-C3)(-C2)", gad); }
  { GAdata gad; gad.setEnthalpy(-18.48* 4.18); gad.setEntropy(-0.2* 4.18); gad.setCp(100,-0.82* 4.18); gad.setCp(200,1.49* 4.18); gad.setCp(298,3.59* 4.18); gad.setCp(400,5.24* 4.18); gad.setCp(500,6.43* 4.18); gad.setCp(600,7.28* 4.18); gad.setCp(700,7.92* 4.18); gad.setCp(800,8.43* 4.18); gad.setCp(900,8.84* 4.18); gad.setCp(1000,9.19* 4.18); gad.setCp(1100,9.49* 4.18); gad.setCp(1200,9.75* 4.18); gad.setCp(1300,9.97* 4.18); gad.setCp(1400,10.17* 4.18); gad.setCp(1500,10.34* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}6)(-O5)(-C2(=O3))(-H4)", gad); }
  { GAdata gad; gad.setEnthalpy(-2.84* 4.18); gad.setEntropy(3.09* 4.18); gad.setCp(100,-0.63* 4.18); gad.setCp(200,1.22* 4.18); gad.setCp(298,3.2* 4.18); gad.setCp(400,4.72* 4.18); gad.setCp(500,5.8* 4.18); gad.setCp(600,6.58* 4.18); gad.setCp(700,7.17* 4.18); gad.setCp(800,7.64* 4.18); gad.setCp(900,8.02* 4.18); gad.setCp(1000,8.35* 4.18); gad.setCp(1100,8.63* 4.18); gad.setCp(1200,8.87* 4.18); gad.setCp(1300,9.08* 4.18); gad.setCp(1400,9.26* 4.18); gad.setCp(1500,9.42* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}6)(-{Pt}5)(-C2(=O3))(-H4)", gad); }
  { GAdata gad; gad.setEnthalpy(-10.57* 4.18); gad.setEntropy(4.05* 4.18); gad.setCp(100,0.59* 4.18); gad.setCp(200,2.43* 4.18); gad.setCp(298,4.62* 4.18); gad.setCp(400,6.57* 4.18); gad.setCp(500,8.11* 4.18); gad.setCp(600,9.31* 4.18); gad.setCp(700,10.28* 4.18); gad.setCp(800,11.09* 4.18); gad.setCp(900,11.77* 4.18); gad.setCp(1000,12.35* 4.18); gad.setCp(1100,12.85* 4.18); gad.setCp(1200,13.29* 4.18); gad.setCp(1300,13.66* 4.18); gad.setCp(1400,13.99* 4.18); gad.setCp(1500,14.27* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}6)(-C2(=O3))(-H5)(-H4)", gad); }
  { GAdata gad; gad.setEnthalpy(-14.56* 4.18); gad.setEntropy(0.21* 4.18); gad.setCp(100,-1.28* 4.18); gad.setCp(200,0.81* 4.18); gad.setCp(298,2.5* 4.18); gad.setCp(400,3.56* 4.18); gad.setCp(500,4.22* 4.18); gad.setCp(600,4.64* 4.18); gad.setCp(700,4.92* 4.18); gad.setCp(800,5.12* 4.18); gad.setCp(900,5.27* 4.18); gad.setCp(1000,5.38* 4.18); gad.setCp(1100,5.46* 4.18); gad.setCp(1200,5.53* 4.18); gad.setCp(1300,5.59* 4.18); gad.setCp(1400,5.64* 4.18); gad.setCp(1500,5.68* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}6)(-{Pt}5)(-O4)(-C2(=O3))", gad); }
  { GAdata gad; gad.setEnthalpy(3.55* 4.18); gad.setEntropy(2.09* 4.18); gad.setCp(100,-0.59* 4.18); gad.setCp(200,1.42* 4.18); gad.setCp(298,2.47* 4.18); gad.setCp(400,3.1* 4.18); gad.setCp(500,3.51* 4.18); gad.setCp(600,3.8* 4.18); gad.setCp(700,4* 4.18); gad.setCp(800,4.15* 4.18); gad.setCp(900,4.27* 4.18); gad.setCp(1000,4.36* 4.18); gad.setCp(1100,4.43* 4.18); gad.setCp(1200,4.5* 4.18); gad.setCp(1300,4.55* 4.18); gad.setCp(1400,4.59* 4.18); gad.setCp(1500,4.63* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}6)(-{Pt}5)(-{Pt}4)(-C2(=O3))", gad); }
  { GAdata gad; gad.setEnthalpy(-3.17* 4.18); gad.setEntropy(0.13* 4.18); gad.setCp(100,0.89* 4.18); gad.setCp(200,2.57* 4.18); gad.setCp(298,4.07* 4.18); gad.setCp(400,5.45* 4.18); gad.setCp(500,6.58* 4.18); gad.setCp(600,7.49* 4.18); gad.setCp(700,8.21* 4.18); gad.setCp(800,8.81* 4.18); gad.setCp(900,9.31* 4.18); gad.setCp(1000,9.74* 4.18); gad.setCp(1100,10.1* 4.18); gad.setCp(1200,10.41* 4.18); gad.setCp(1300,10.68* 4.18); gad.setCp(1400,10.92* 4.18); gad.setCp(1500,11.12* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-O7)(-C4(=O5))(-C2(=O3))(-H6)", gad); }
  { GAdata gad; gad.setEnthalpy(-4.19* 4.18); gad.setEntropy(-0.07* 4.18); gad.setCp(100,0.12* 4.18); gad.setCp(200,1.62* 4.18); gad.setCp(298,3.51* 4.18); gad.setCp(400,5.12* 4.18); gad.setCp(500,6.34* 4.18); gad.setCp(600,7.26* 4.18); gad.setCp(700,7.97* 4.18); gad.setCp(800,8.53* 4.18); gad.setCp(900,8.99* 4.18); gad.setCp(1000,9.37* 4.18); gad.setCp(1100,9.68* 4.18); gad.setCp(1200,9.95* 4.18); gad.setCp(1300,10.17* 4.18); gad.setCp(1400,10.37* 4.18); gad.setCp(1500,10.53* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}7)(-C4(=O5))(-C2(=O3))(-H6)", gad); }
  { GAdata gad; gad.setEnthalpy(-2.66* 4.18); gad.setEntropy(7.62* 4.18); gad.setCp(100,2.78* 4.18); gad.setCp(200,3.87* 4.18); gad.setCp(298,5.17* 4.18); gad.setCp(400,6.71* 4.18); gad.setCp(500,8.13* 4.18); gad.setCp(600,9.35* 4.18); gad.setCp(700,10.39* 4.18); gad.setCp(800,11.27* 4.18); gad.setCp(900,12.04* 4.18); gad.setCp(1000,12.7* 4.18); gad.setCp(1100,13.28* 4.18); gad.setCp(1200,13.77* 4.18); gad.setCp(1300,14.2* 4.18); gad.setCp(1400,14.58* 4.18); gad.setCp(1500,14.91* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-C4(=O5))(-C2(=O3))(-H7)(-H6)", gad); }
  { GAdata gad; gad.setEnthalpy(-12.48* 4.18); gad.setEntropy(-5.56* 4.18); gad.setCp(100,-1.52* 4.18); gad.setCp(200,0.65* 4.18); gad.setCp(298,2.58* 4.18); gad.setCp(400,3.87* 4.18); gad.setCp(500,4.69* 4.18); gad.setCp(600,5.24* 4.18); gad.setCp(700,5.61* 4.18); gad.setCp(800,5.88* 4.18); gad.setCp(900,6.09* 4.18); gad.setCp(1000,6.24* 4.18); gad.setCp(1100,6.36* 4.18); gad.setCp(1200,6.46* 4.18); gad.setCp(1300,6.54* 4.18); gad.setCp(1400,6.6* 4.18); gad.setCp(1500,6.65* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}7)(-O6)(-C4(=O5))(-C2(=O3))", gad); }
  { GAdata gad; gad.setEnthalpy(-13.53* 4.18); gad.setEntropy(-1.4* 4.18); gad.setCp(100,-0.99* 4.18); gad.setCp(200,0.67* 4.18); gad.setCp(298,2.1* 4.18); gad.setCp(400,3.07* 4.18); gad.setCp(500,3.74* 4.18); gad.setCp(600,4.21* 4.18); gad.setCp(700,4.55* 4.18); gad.setCp(800,4.8* 4.18); gad.setCp(900,4.99* 4.18); gad.setCp(1000,5.14* 4.18); gad.setCp(1100,5.26* 4.18); gad.setCp(1200,5.35* 4.18); gad.setCp(1300,5.43* 4.18); gad.setCp(1400,5.49* 4.18); gad.setCp(1500,5.54* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}7)(-{Pt}6)(-C4(=O5))(-C2(=O3))", gad); }
  { GAdata gad; gad.setEnthalpy(-17.06* 4.18); gad.setEntropy(0.49* 4.18); gad.setCp(100,-1.49* 4.18); gad.setCp(200,1.41* 4.18); gad.setCp(298,3.53* 4.18); gad.setCp(400,5* 4.18); gad.setCp(500,6.01* 4.18); gad.setCp(600,6.74* 4.18); gad.setCp(700,7.3* 4.18); gad.setCp(800,7.76* 4.18); gad.setCp(900,8.14* 4.18); gad.setCp(1000,8.47* 4.18); gad.setCp(1100,8.75* 4.18); gad.setCp(1200,9* 4.18); gad.setCp(1300,9.21* 4.18); gad.setCp(1400,9.39* 4.18); gad.setCp(1500,9.54* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-{Pt}4)(-O3)(-H2)", gad); }
  { GAdata gad; gad.setEnthalpy(-8.11* 4.18); gad.setEntropy(2.31* 4.18); gad.setCp(100,-1.76* 4.18); gad.setCp(200,0.58* 4.18); gad.setCp(298,3.03* 4.18); gad.setCp(400,4.64* 4.18); gad.setCp(500,5.65* 4.18); gad.setCp(600,6.32* 4.18); gad.setCp(700,6.8* 4.18); gad.setCp(800,7.19* 4.18); gad.setCp(900,7.51* 4.18); gad.setCp(1000,7.78* 4.18); gad.setCp(1100,8.01* 4.18); gad.setCp(1200,8.21* 4.18); gad.setCp(1300,8.38* 4.18); gad.setCp(1400,8.54* 4.18); gad.setCp(1500,8.67* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-{Pt}4)(-{Pt}3)(-H2)", gad); }
  { GAdata gad; gad.setEnthalpy(-18.64* 4.18); gad.setEntropy(5.42* 4.18); gad.setCp(100,0.67* 4.18); gad.setCp(200,2.96* 4.18); gad.setCp(298,5.05* 4.18); gad.setCp(400,6.83* 4.18); gad.setCp(500,8.22* 4.18); gad.setCp(600,9.32* 4.18); gad.setCp(700,10.21* 4.18); gad.setCp(800,10.96* 4.18); gad.setCp(900,11.6* 4.18); gad.setCp(1000,12.17* 4.18); gad.setCp(1100,12.66* 4.18); gad.setCp(1200,13.09* 4.18); gad.setCp(1300,13.47* 4.18); gad.setCp(1400,13.8* 4.18); gad.setCp(1500,14.1* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-O4)(-H3)(-H2)", gad); }
  { GAdata gad; gad.setEnthalpy(-2.4* 4.18); gad.setEntropy(4.7* 4.18); gad.setCp(100,-0.5* 4.18); gad.setCp(200,2.24* 4.18); gad.setCp(298,4.87* 4.18); gad.setCp(400,6.79* 4.18); gad.setCp(500,8.14* 4.18); gad.setCp(600,9.15* 4.18); gad.setCp(700,9.94* 4.18); gad.setCp(800,10.61* 4.18); gad.setCp(900,11.19* 4.18); gad.setCp(1000,11.69* 4.18); gad.setCp(1100,12.13* 4.18); gad.setCp(1200,12.52* 4.18); gad.setCp(1300,12.85* 4.18); gad.setCp(1400,13.15* 4.18); gad.setCp(1500,13.41* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-{Pt}4)(-H3)(-H2)", gad); }
  { GAdata gad; gad.setEnthalpy(-12.53* 4.18); gad.setEntropy(11.26* 4.18); gad.setCp(100,3* 4.18); gad.setCp(200,4.92* 4.18); gad.setCp(298,6.93* 4.18); gad.setCp(400,8.79* 4.18); gad.setCp(500,10.37* 4.18); gad.setCp(600,11.67* 4.18); gad.setCp(700,12.79* 4.18); gad.setCp(800,13.76* 4.18); gad.setCp(900,14.61* 4.18); gad.setCp(1000,15.37* 4.18); gad.setCp(1100,16.04* 4.18); gad.setCp(1200,16.63* 4.18); gad.setCp(1300,17.14* 4.18); gad.setCp(1400,17.6* 4.18); gad.setCp(1500,18* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-H4)(-H3)(-H2)", gad); }
  { GAdata gad; gad.setEnthalpy(-22.06* 4.18); gad.setEntropy(2.22* 4.18); gad.setCp(100,-0.68* 4.18); gad.setCp(200,1.68* 4.18); gad.setCp(298,3.1* 4.18); gad.setCp(400,3.83* 4.18); gad.setCp(500,4.21* 4.18); gad.setCp(600,4.4* 4.18); gad.setCp(700,4.5* 4.18); gad.setCp(800,4.56* 4.18); gad.setCp(900,4.61* 4.18); gad.setCp(1000,4.64* 4.18); gad.setCp(1100,4.66* 4.18); gad.setCp(1200,4.69* 4.18); gad.setCp(1300,4.71* 4.18); gad.setCp(1400,4.73* 4.18); gad.setCp(1500,4.74* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-{Pt}4)(-{Pt}3)(-O2)", gad); }
  { GAdata gad; gad.setEnthalpy(13.48* 4.18); gad.setEntropy(1.62* 4.18); gad.setCp(100,-1.8* 4.18); gad.setCp(200,-0.2* 4.18); gad.setCp(298,1.32* 4.18); gad.setCp(400,2.23* 4.18); gad.setCp(500,2.77* 4.18); gad.setCp(600,3.1* 4.18); gad.setCp(700,3.31* 4.18); gad.setCp(800,3.46* 4.18); gad.setCp(900,3.56* 4.18); gad.setCp(1000,3.64* 4.18); gad.setCp(1100,3.69* 4.18); gad.setCp(1200,3.74* 4.18); gad.setCp(1300,3.77* 4.18); gad.setCp(1400,3.8* 4.18); gad.setCp(1500,3.82* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(-{Pt}5)(-{Pt}4)(-{Pt}3)(-{Pt}2)", gad); }
  { GAdata gad; gad.setEnthalpy(-44.16* 4.18); gad.setEntropy(-9.73* 4.18); gad.setCp(100,7.96* 4.18); gad.setCp(200,10.98* 4.18); gad.setCp(298,12.82* 4.18); gad.setCp(400,14.07* 4.18); gad.setCp(500,14.96* 4.18); gad.setCp(600,15.59* 4.18); gad.setCp(700,16.06* 4.18); gad.setCp(800,16.4* 4.18); gad.setCp(900,16.66* 4.18); gad.setCp(1000,16.86* 4.18); gad.setCp(1100,17.01* 4.18); gad.setCp(1200,17.13* 4.18); gad.setCp(1300,17.23* 4.18); gad.setCp(1400,17.31* 4.18); gad.setCp(1500,17.38* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(=O2)(-C4(=O5))(-C3)", gad); }
  { GAdata gad; gad.setEnthalpy(-42.63* 4.18); gad.setEntropy(14.23* 4.18); gad.setCp(100,4.54* 4.18); gad.setCp(200,6.73* 4.18); gad.setCp(298,8.16* 4.18); gad.setCp(400,9.43* 4.18); gad.setCp(500,10.51* 4.18); gad.setCp(600,11.41* 4.18); gad.setCp(700,12.15* 4.18); gad.setCp(800,12.77* 4.18); gad.setCp(900,13.28* 4.18); gad.setCp(1000,13.72* 4.18); gad.setCp(1100,14.08* 4.18); gad.setCp(1200,14.4* 4.18); gad.setCp(1300,14.67* 4.18); gad.setCp(1400,14.9* 4.18); gad.setCp(1500,15.1* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(=O2)(-C3)(-H4)", gad); }
  { GAdata gad; gad.setEnthalpy(-52.47* 4.18); gad.setEntropy(8.73* 4.18); gad.setCp(100,2.2* 4.18); gad.setCp(200,4.62* 4.18); gad.setCp(298,6.31* 4.18); gad.setCp(400,7.4* 4.18); gad.setCp(500,8.15* 4.18); gad.setCp(600,8.68* 4.18); gad.setCp(700,9.08* 4.18); gad.setCp(800,9.38* 4.18); gad.setCp(900,9.61* 4.18); gad.setCp(1000,9.8* 4.18); gad.setCp(1100,9.95* 4.18); gad.setCp(1200,10.08* 4.18); gad.setCp(1300,10.18* 4.18); gad.setCp(1400,10.26* 4.18); gad.setCp(1500,10.34* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(=O2)(-{Pt}4)(-C3)", gad); }
  { GAdata gad; gad.setEnthalpy(-40.19* 4.18); gad.setEntropy(9.5* 4.18); gad.setCp(100,3.3* 4.18); gad.setCp(200,5.85* 4.18); gad.setCp(298,7.35* 4.18); gad.setCp(400,8.34* 4.18); gad.setCp(500,9.03* 4.18); gad.setCp(600,9.53* 4.18); gad.setCp(700,9.9* 4.18); gad.setCp(800,10.18* 4.18); gad.setCp(900,10.41* 4.18); gad.setCp(1000,10.6* 4.18); gad.setCp(1100,10.75* 4.18); gad.setCp(1200,10.88* 4.18); gad.setCp(1300,10.99* 4.18); gad.setCp(1400,11.09* 4.18); gad.setCp(1500,11.17* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(=O2)(-C4)(-C3)", gad); }
  { GAdata gad; gad.setEnthalpy(-27.89* 4.18); gad.setEntropy(37.22* 4.18); gad.setCp(100,0.68* 4.18); gad.setCp(200,1.86* 4.18); gad.setCp(298,2.73* 4.18); gad.setCp(400,3.64* 4.18); gad.setCp(500,4.52* 4.18); gad.setCp(600,5.31* 4.18); gad.setCp(700,6* 4.18); gad.setCp(800,6.59* 4.18); gad.setCp(900,7.11* 4.18); gad.setCp(1000,7.55* 4.18); gad.setCp(1100,7.93* 4.18); gad.setCp(1200,8.26* 4.18); gad.setCp(1300,8.55* 4.18); gad.setCp(1400,8.79* 4.18); gad.setCp(1500,9.01* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(=O2)(-C3(=O4))(-H5)", gad); }
  { GAdata gad; gad.setEnthalpy(-38.24* 4.18); gad.setEntropy(32.56* 4.18); gad.setCp(100,-1.18* 4.18); gad.setCp(200,0.21* 4.18); gad.setCp(298,1.15* 4.18); gad.setCp(400,1.76* 4.18); gad.setCp(500,2.2* 4.18); gad.setCp(600,2.55* 4.18); gad.setCp(700,2.83* 4.18); gad.setCp(800,3.07* 4.18); gad.setCp(900,3.28* 4.18); gad.setCp(1000,3.46* 4.18); gad.setCp(1100,3.62* 4.18); gad.setCp(1200,3.75* 4.18); gad.setCp(1300,3.87* 4.18); gad.setCp(1400,3.98* 4.18); gad.setCp(1500,4.07* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(=O2)(-{Pt}5)(-C3(=O4))", gad); }
  { GAdata gad; gad.setEnthalpy(-26.98* 4.18); gad.setEntropy(6.42* 4.18); gad.setCp(100,3.01* 4.18); gad.setCp(200,5.29* 4.18); gad.setCp(298,6.84* 4.18); gad.setCp(400,7.95* 4.18); gad.setCp(500,8.76* 4.18); gad.setCp(600,9.37* 4.18); gad.setCp(700,9.83* 4.18); gad.setCp(800,10.18* 4.18); gad.setCp(900,10.45* 4.18); gad.setCp(1000,10.67* 4.18); gad.setCp(1100,10.85* 4.18); gad.setCp(1200,10.99* 4.18); gad.setCp(1300,11.11* 4.18); gad.setCp(1400,11.2* 4.18); gad.setCp(1500,11.29* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(=O2)(-C5(=O6))(-C3(=O4))", gad); }
  { GAdata gad; gad.setEnthalpy(-48.56* 4.18); gad.setEntropy(16.05* 4.18); gad.setCp(100,4.08* 4.18); gad.setCp(200,6.1* 4.18); gad.setCp(298,7.68* 4.18); gad.setCp(400,8.95* 4.18); gad.setCp(500,9.98* 4.18); gad.setCp(600,10.81* 4.18); gad.setCp(700,11.5* 4.18); gad.setCp(800,12.07* 4.18); gad.setCp(900,12.56* 4.18); gad.setCp(1000,12.97* 4.18); gad.setCp(1100,13.32* 4.18); gad.setCp(1200,13.61* 4.18); gad.setCp(1300,13.86* 4.18); gad.setCp(1400,14.08* 4.18); gad.setCp(1500,14.27* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(=O2)(-{Pt}4)(-H3)", gad); }
  { GAdata gad; gad.setEnthalpy(-59.1* 4.18); gad.setEntropy(12.22* 4.18); gad.setCp(100,2.35* 4.18); gad.setCp(200,4.85* 4.18); gad.setCp(298,6.25* 4.18); gad.setCp(400,6.97* 4.18); gad.setCp(500,7.44* 4.18); gad.setCp(600,7.81* 4.18); gad.setCp(700,8.11* 4.18); gad.setCp(800,8.37* 4.18); gad.setCp(900,8.59* 4.18); gad.setCp(1000,8.77* 4.18); gad.setCp(1100,8.92* 4.18); gad.setCp(1200,9.05* 4.18); gad.setCp(1300,9.15* 4.18); gad.setCp(1400,9.24* 4.18); gad.setCp(1500,9.32* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("C1(=O2)(-{Pt}4)(-{Pt}3)", gad); }
  { GAdata gad; gad.setEnthalpy(-38* 4.18); gad.setEntropy(9.04* 4.18); gad.setCp(100,3.24* 4.18); gad.setCp(200,4.37* 4.18); gad.setCp(298,5.11* 4.18); gad.setCp(400,5.84* 4.18); gad.setCp(500,6.5* 4.18); gad.setCp(600,7.06* 4.18); gad.setCp(700,7.52* 4.18); gad.setCp(800,7.91* 4.18); gad.setCp(900,8.23* 4.18); gad.setCp(1000,8.51* 4.18); gad.setCp(1100,8.76* 4.18); gad.setCp(1200,8.97* 4.18); gad.setCp(1300,9.15* 4.18); gad.setCp(1400,9.32* 4.18); gad.setCp(1500,9.46* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("O1(-C2)(-H3)", gad); }
  { GAdata gad; gad.setEnthalpy(-20.24* 4.18); gad.setEntropy(2.66* 4.18); gad.setCp(100,0.63* 4.18); gad.setCp(200,1.89* 4.18); gad.setCp(298,2.77* 4.18); gad.setCp(400,3.4* 4.18); gad.setCp(500,3.86* 4.18); gad.setCp(600,4.18* 4.18); gad.setCp(700,4.42* 4.18); gad.setCp(800,4.58* 4.18); gad.setCp(900,4.7* 4.18); gad.setCp(1000,4.78* 4.18); gad.setCp(1100,4.85* 4.18); gad.setCp(1200,4.89* 4.18); gad.setCp(1300,4.92* 4.18); gad.setCp(1400,4.94* 4.18); gad.setCp(1500,4.96* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("O1(-{Pt}3)(-C2)", gad); }
  { GAdata gad; gad.setEnthalpy(0.0); gad.setEntropy(0.0); gad.setCp(298,0.0); gad.setCp(400,0.0); gad.setCp(500,0.0); gad.setCp(600,0.0); gad.setCp(800,0.0); gad.setCp(1000,0.0); gad.setCp(1500,0.0);   gad.PrepareGA(); ThermoGA::AddGA("{Pt}1", gad); }
  { GAdata gad; gad.setEnthalpy(0.0); gad.setEntropy(0.0); gad.setCp(298,0.0); gad.setCp(400,0.0); gad.setCp(500,0.0); gad.setCp(600,0.0); gad.setCp(800,0.0); gad.setCp(1000,0.0); gad.setCp(1500,0.0);   gad.PrepareGA(); ThermoGA::AddGA("{Pt}1(-O2)", gad); }
  { GAdata gad; gad.setEnthalpy(0.0); gad.setEntropy(0.0); gad.setCp(298,0.0); gad.setCp(400,0.0); gad.setCp(500,0.0); gad.setCp(600,0.0); gad.setCp(800,0.0); gad.setCp(1000,0.0); gad.setCp(1500,0.0);   gad.PrepareGA(); ThermoGA::AddGA("{Pt}1(-C2)", gad); }
  { GAdata gad; gad.setEnthalpy(0.0); gad.setEntropy(0.0); gad.setCp(298,0.0); gad.setCp(400,0.0); gad.setCp(500,0.0); gad.setCp(600,0.0); gad.setCp(800,0.0); gad.setCp(1000,0.0); gad.setCp(1500,0.0);   gad.PrepareGA(); ThermoGA::AddGA("{Pt}1(_O2)", gad); }
  { GAdata gad; gad.setEnthalpy(0.0); gad.setEntropy(0.0); gad.setCp(298,0.0); gad.setCp(400,0.0); gad.setCp(500,0.0); gad.setCp(600,0.0); gad.setCp(800,0.0); gad.setCp(1000,0.0); gad.setCp(1500,0.0);   gad.PrepareGA(); ThermoGA::AddGA("C1(-O2)", gad); }
  { GAdata gad; gad.setEnthalpy(-24.7* 4.18); gad.setEntropy(6.7* 4.18); gad.setCp(298,4.3* 4.18); gad.setCp(400,4.6* 4.18); gad.setCp(500,4.9* 4.18); gad.setCp(600,5.3* 4.18); gad.setCp(800,5.7* 4.18); gad.setCp(1000,5.8* 4.18); gad.setCp(1500,5.8* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("O1(-{Pt}3)(-C2)", gad); }
  { GAdata gad; gad.setEnthalpy(-55.19* 4.18); gad.setEntropy(-0.42* 4.18); gad.setCp(100,-0.58* 4.18); gad.setCp(200,1.11* 4.18); gad.setCp(298,2.09* 4.18); gad.setCp(400,2.8* 4.18); gad.setCp(500,3.33* 4.18); gad.setCp(600,3.72* 4.18); gad.setCp(700,4* 4.18); gad.setCp(800,4.21* 4.18); gad.setCp(900,4.36* 4.18); gad.setCp(1000,4.47* 4.18); gad.setCp(1100,4.56* 4.18); gad.setCp(1200,4.63* 4.18); gad.setCp(1300,4.68* 4.18); gad.setCp(1400,4.72* 4.18); gad.setCp(1500,4.76* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("O1(-{Pt}4)(-C2(=O3))", gad); }
  { GAdata gad; gad.setEnthalpy(-38.34* 4.18); gad.setEntropy(11.3* 4.18); gad.setCp(100,3.13* 4.18); gad.setCp(200,4.55* 4.18); gad.setCp(298,5.57* 4.18); gad.setCp(400,6.27* 4.18); gad.setCp(500,6.74* 4.18); gad.setCp(600,7.07* 4.18); gad.setCp(700,7.32* 4.18); gad.setCp(800,7.52* 4.18); gad.setCp(900,7.71* 4.18); gad.setCp(1000,7.88* 4.18); gad.setCp(1100,8.04* 4.18); gad.setCp(1200,8.19* 4.18); gad.setCp(1300,8.32* 4.18); gad.setCp(1400,8.45* 4.18); gad.setCp(1500,8.57* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("O1(-{Pt}3)(-H2)", gad); }
  { GAdata gad; gad.setEnthalpy(-61.88* 4.18); gad.setEntropy(14.76* 4.18); gad.setCp(100,4.55* 4.18); gad.setCp(200,6.95* 4.18); gad.setCp(298,8.29* 4.18); gad.setCp(400,9.1* 4.18); gad.setCp(500,9.7* 4.18); gad.setCp(600,10.17* 4.18); gad.setCp(700,10.58* 4.18); gad.setCp(800,10.96* 4.18); gad.setCp(900,11.32* 4.18); gad.setCp(1000,11.66* 4.18); gad.setCp(1100,11.98* 4.18); gad.setCp(1200,12.29* 4.18); gad.setCp(1300,12.57* 4.18); gad.setCp(1400,12.83* 4.18); gad.setCp(1500,13.07* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("O1(-H3)(-H2)", gad); }
  { GAdata gad; gad.setEnthalpy(19* 4.18); gad.setEntropy(0* 4.18); gad.setCp(298,3.7* 4.18); gad.setCp(400,3.7* 4.18); gad.setCp(500,3.7* 4.18); gad.setCp(600,3.7* 4.18); gad.setCp(800,4.2* 4.18); gad.setCp(1000,4.2* 4.18); gad.setCp(1500,4.8* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("O1(=O2)", gad); }
  { GAdata gad; gad.setEnthalpy(-31.53* 4.18); gad.setEntropy(3.04* 4.18); gad.setCp(100,-1.26* 4.18); gad.setCp(200,1.25* 4.18); gad.setCp(298,2.52* 4.18); gad.setCp(400,3.09* 4.18); gad.setCp(500,3.39* 4.18); gad.setCp(600,3.56* 4.18); gad.setCp(700,3.67* 4.18); gad.setCp(800,3.74* 4.18); gad.setCp(900,3.79* 4.18); gad.setCp(1000,3.82* 4.18); gad.setCp(1100,3.85* 4.18); gad.setCp(1200,3.87* 4.18); gad.setCp(1300,3.88* 4.18); gad.setCp(1400,3.9* 4.18); gad.setCp(1500,3.91* 4.18);   gad.PrepareGA(); ThermoGA::AddGA("O1(-{Pt}3)(-{Pt}2)", gad); }
  { GAdata gad; gad.setEnthalpy(-19.73* 4.18); gad.setEntropy(28.77* 4.18); gad.setCp(100,8.88* 4.18); gad.setCp(200,9.72* 4.18); gad.setCp(298,10.46* 4.18); gad.setCp(400,11.73* 4.18); gad.setCp(500,13.18* 4.18); gad.setCp(600,14.59* 4.18); gad.setCp(700,15.91* 4.18); gad.setCp(800,17.12* 4.18); gad.setCp(900,18.22* 4.18); gad.setCp(1000,19.2* 4.18); gad.setCp(1100,20.08* 4.18); gad.setCp(1200,20.86* 4.18); gad.setCp(1300,21.54* 4.18); gad.setCp(1400,22.15* 4.18); gad.setCp(1500,22.68* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|_{Pt}|==1](-H5)(-H4)(-H3)(-H2)", &correctionCheck0, gad); }
  { GAdata gad; gad.setEnthalpy(2.49* 4.18); gad.setEntropy(-0.6* 4.18); gad.setCp(100,-0.24* 4.18); gad.setCp(200,-0.5* 4.18); gad.setCp(298,-0.34* 4.18); gad.setCp(400,-0.18* 4.18); gad.setCp(500,-0.06* 4.18); gad.setCp(600,0.04* 4.18); gad.setCp(700,0.11* 4.18); gad.setCp(800,0.15* 4.18); gad.setCp(900,0.18* 4.18); gad.setCp(1000,0.19* 4.18); gad.setCp(1100,0.2* 4.18); gad.setCp(1200,0.19* 4.18); gad.setCp(1300,0.19* 4.18); gad.setCp(1400,0.18* 4.18); gad.setCp(1500,0.17* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==0](-C2[|-{Pt}|==0])", &correctionCheck1, gad); }
  { GAdata gad; gad.setEnthalpy(-0.84* 4.18); gad.setEntropy(0.95* 4.18); gad.setCp(100,-0.13* 4.18); gad.setCp(200,-0.85* 4.18); gad.setCp(298,-0.7* 4.18); gad.setCp(400,-0.52* 4.18); gad.setCp(500,-0.39* 4.18); gad.setCp(600,-0.32* 4.18); gad.setCp(700,-0.27* 4.18); gad.setCp(800,-0.24* 4.18); gad.setCp(900,-0.22* 4.18); gad.setCp(1000,-0.2* 4.18); gad.setCp(1100,-0.19* 4.18); gad.setCp(1200,-0.18* 4.18); gad.setCp(1300,-0.17* 4.18); gad.setCp(1400,-0.16* 4.18); gad.setCp(1500,-0.15* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|-{Pt}|==2](-O3(-H4)))", &correctionCheck2, gad); }
  { GAdata gad; gad.setEnthalpy(3.21* 4.18); gad.setEntropy(-1.04* 4.18); gad.setCp(100,-0.34* 4.18); gad.setCp(200,-0.24* 4.18); gad.setCp(298,-0.17* 4.18); gad.setCp(400,-0.1* 4.18); gad.setCp(500,-0.04* 4.18); gad.setCp(600,0.01* 4.18); gad.setCp(700,0.05* 4.18); gad.setCp(800,0.08* 4.18); gad.setCp(900,0.1* 4.18); gad.setCp(1000,0.11* 4.18); gad.setCp(1100,0.11* 4.18); gad.setCp(1200,0.11* 4.18); gad.setCp(1300,0.11* 4.18); gad.setCp(1400,0.1* 4.18); gad.setCp(1500,0.1* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==0](-C2[|-{Pt}|==1])", &correctionCheck3, gad); }
  { GAdata gad; gad.setEnthalpy(6.56* 4.18); gad.setEntropy(-1.85* 4.18); gad.setCp(100,-0.08* 4.18); gad.setCp(200,-0.03* 4.18); gad.setCp(298,-0.01* 4.18); gad.setCp(400,0.04* 4.18); gad.setCp(500,0.09* 4.18); gad.setCp(600,0.14* 4.18); gad.setCp(700,0.19* 4.18); gad.setCp(800,0.22* 4.18); gad.setCp(900,0.25* 4.18); gad.setCp(1000,0.26* 4.18); gad.setCp(1100,0.26* 4.18); gad.setCp(1200,0.26* 4.18); gad.setCp(1300,0.25* 4.18); gad.setCp(1400,0.24* 4.18); gad.setCp(1500,0.23* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==0](-C2[|-{Pt}|==2])", &correctionCheck4, gad); }
  { GAdata gad; gad.setEnthalpy(-4.69* 4.18); gad.setEntropy(-0.71* 4.18); gad.setCp(100,-0.16* 4.18); gad.setCp(200,-0.68* 4.18); gad.setCp(298,-0.38* 4.18); gad.setCp(400,-0.11* 4.18); gad.setCp(500,0.04* 4.18); gad.setCp(600,0.13* 4.18); gad.setCp(700,0.17* 4.18); gad.setCp(800,0.2* 4.18); gad.setCp(900,0.21* 4.18); gad.setCp(1000,0.21* 4.18); gad.setCp(1100,0.2* 4.18); gad.setCp(1200,0.19* 4.18); gad.setCp(1300,0.18* 4.18); gad.setCp(1400,0.17* 4.18); gad.setCp(1500,0.16* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==0](-C2[|-{Pt}|==3])", &correctionCheck5, gad); }
  { GAdata gad; gad.setEnthalpy(-106.29* 4.18); gad.setEntropy(14.89* 4.18); gad.setCp(100,3.94* 4.18); gad.setCp(200,7.95* 4.18); gad.setCp(298,10.78* 4.18); gad.setCp(400,12.87* 4.18); gad.setCp(500,14.44* 4.18); gad.setCp(600,15.62* 4.18); gad.setCp(700,16.52* 4.18); gad.setCp(800,17.23* 4.18); gad.setCp(900,17.8* 4.18); gad.setCp(1000,18.27* 4.18); gad.setCp(1100,18.67* 4.18); gad.setCp(1200,19* 4.18); gad.setCp(1300,19.29* 4.18); gad.setCp(1400,19.55* 4.18); gad.setCp(1500,19.76* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1(=O2)(-{Pt}3)(-O4(-H5))", &correctionCheck6, gad); }
  { GAdata gad; gad.setEnthalpy(-98.93* 4.18); gad.setEntropy(34.05* 4.18); gad.setCp(100,7.85* 4.18); gad.setCp(200,8.88* 4.18); gad.setCp(298,10.09* 4.18); gad.setCp(400,11.05* 4.18); gad.setCp(500,11.8* 4.18); gad.setCp(600,12.42* 4.18); gad.setCp(700,12.92* 4.18); gad.setCp(800,13.34* 4.18); gad.setCp(900,13.7* 4.18); gad.setCp(1000,13.99* 4.18); gad.setCp(1100,14.23* 4.18); gad.setCp(1200,14.43* 4.18); gad.setCp(1300,14.61* 4.18); gad.setCp(1400,14.75* 4.18); gad.setCp(1500,14.88* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|_{Pt}|==1](=O3)(=O2)", &correctionCheck7, gad); }
  { GAdata gad; gad.setEnthalpy(-4.87* 4.18); gad.setEntropy(0.77* 4.18); gad.setCp(100,0.45* 4.18); gad.setCp(200,0.53* 4.18); gad.setCp(298,0.25* 4.18); gad.setCp(400,0.05* 4.18); gad.setCp(500,-0.04* 4.18); gad.setCp(600,-0.07* 4.18); gad.setCp(700,-0.09* 4.18); gad.setCp(800,-0.09* 4.18); gad.setCp(900,-0.09* 4.18); gad.setCp(1000,-0.08* 4.18); gad.setCp(1100,-0.08* 4.18); gad.setCp(1200,-0.07* 4.18); gad.setCp(1300,-0.07* 4.18); gad.setCp(1400,-0.06* 4.18); gad.setCp(1500,-0.06* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==1](-C2[|-{Pt}|==1])", &correctionCheck8, gad); }
  { GAdata gad; gad.setEnthalpy(-9.56* 4.18); gad.setEntropy(2.9* 4.18); gad.setCp(100,1.26* 4.18); gad.setCp(200,0.95* 4.18); gad.setCp(298,0.57* 4.18); gad.setCp(400,0.33* 4.18); gad.setCp(500,0.18* 4.18); gad.setCp(600,0.1* 4.18); gad.setCp(700,0.05* 4.18); gad.setCp(800,0.02* 4.18); gad.setCp(900,0* 4.18); gad.setCp(1000,-0.01* 4.18); gad.setCp(1100,-0.02* 4.18); gad.setCp(1200,-0.03* 4.18); gad.setCp(1300,-0.03* 4.18); gad.setCp(1400,-0.03* 4.18); gad.setCp(1500,-0.03* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==1](-O2[|-{Pt}|==1])", &correctionCheck9, gad); }
  { GAdata gad; gad.setEnthalpy(-5.04* 4.18); gad.setEntropy(-1.2* 4.18); gad.setCp(100,0.1* 4.18); gad.setCp(200,0.39* 4.18); gad.setCp(298,0.25* 4.18); gad.setCp(400,0.15* 4.18); gad.setCp(500,0.09* 4.18); gad.setCp(600,0.07* 4.18); gad.setCp(700,0.06* 4.18); gad.setCp(800,0.06* 4.18); gad.setCp(900,0.06* 4.18); gad.setCp(1000,0.06* 4.18); gad.setCp(1100,0.06* 4.18); gad.setCp(1200,0.06* 4.18); gad.setCp(1300,0.06* 4.18); gad.setCp(1400,0.06* 4.18); gad.setCp(1500,0.05* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==2](-C2[|-{Pt}|==1])", &correctionCheck10, gad); }
  { GAdata gad; gad.setEnthalpy(-10.77* 4.18); gad.setEntropy(1.23* 4.18); gad.setCp(100,0.9* 4.18); gad.setCp(200,0.56* 4.18); gad.setCp(298,0.09* 4.18); gad.setCp(400,-0.16* 4.18); gad.setCp(500,-0.25* 4.18); gad.setCp(600,-0.27* 4.18); gad.setCp(700,-0.27* 4.18); gad.setCp(800,-0.25* 4.18); gad.setCp(900,-0.22* 4.18); gad.setCp(1000,-0.2* 4.18); gad.setCp(1100,-0.18* 4.18); gad.setCp(1200,-0.15* 4.18); gad.setCp(1300,-0.14* 4.18); gad.setCp(1400,-0.12* 4.18); gad.setCp(1500,-0.11* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==2](-C2[|-{Pt}|==2])", &correctionCheck11, gad); }
  { GAdata gad; gad.setEnthalpy(-22.05* 4.18); gad.setEntropy(1.74* 4.18); gad.setCp(100,1.04* 4.18); gad.setCp(200,1.09* 4.18); gad.setCp(298,0.55* 4.18); gad.setCp(400,0.17* 4.18); gad.setCp(500,-0.05* 4.18); gad.setCp(600,-0.18* 4.18); gad.setCp(700,-0.26* 4.18); gad.setCp(800,-0.31* 4.18); gad.setCp(900,-0.35* 4.18); gad.setCp(1000,-0.36* 4.18); gad.setCp(1100,-0.37* 4.18); gad.setCp(1200,-0.36* 4.18); gad.setCp(1300,-0.35* 4.18); gad.setCp(1400,-0.34* 4.18); gad.setCp(1500,-0.32* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==2](-O2[|-{Pt}|==1])", &correctionCheck12, gad); }
  { GAdata gad; gad.setEnthalpy(-20.68* 4.18); gad.setEntropy(-0.69* 4.18); gad.setCp(100,-0.09* 4.18); gad.setCp(200,-0.29* 4.18); gad.setCp(298,-0.24* 4.18); gad.setCp(400,-0.2* 4.18); gad.setCp(500,-0.17* 4.18); gad.setCp(600,-0.15* 4.18); gad.setCp(700,-0.13* 4.18); gad.setCp(800,-0.11* 4.18); gad.setCp(900,-0.09* 4.18); gad.setCp(1000,-0.08* 4.18); gad.setCp(1100,-0.07* 4.18); gad.setCp(1200,-0.06* 4.18); gad.setCp(1300,-0.05* 4.18); gad.setCp(1400,-0.04* 4.18); gad.setCp(1500,-0.04* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|-{Pt}|==1])", &correctionCheck13, gad); }
  { GAdata gad; gad.setEnthalpy(-6.51* 4.18); gad.setEntropy(0.9* 4.18); gad.setCp(100,0.94* 4.18); gad.setCp(200,0.88* 4.18); gad.setCp(298,0.46* 4.18); gad.setCp(400,0.16* 4.18); gad.setCp(500,0* 4.18); gad.setCp(600,-0.07* 4.18); gad.setCp(700,-0.1* 4.18); gad.setCp(800,-0.11* 4.18); gad.setCp(900,-0.11* 4.18); gad.setCp(1000,-0.1* 4.18); gad.setCp(1100,-0.09* 4.18); gad.setCp(1200,-0.08* 4.18); gad.setCp(1300,-0.07* 4.18); gad.setCp(1400,-0.06* 4.18); gad.setCp(1500,-0.05* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|-{Pt}|==2])", &correctionCheck14, gad); }
  { GAdata gad; gad.setEnthalpy(-11.51* 4.18); gad.setEntropy(0.8* 4.18); gad.setCp(100,-1.96* 4.18); gad.setCp(200,-1.12* 4.18); gad.setCp(298,0.19* 4.18); gad.setCp(400,1.24* 4.18); gad.setCp(500,1.97* 4.18); gad.setCp(600,2.46* 4.18); gad.setCp(700,2.81* 4.18); gad.setCp(800,3.05* 4.18); gad.setCp(900,3.22* 4.18); gad.setCp(1000,3.36* 4.18); gad.setCp(1100,3.46* 4.18); gad.setCp(1200,3.53* 4.18); gad.setCp(1300,3.6* 4.18); gad.setCp(1400,3.65* 4.18); gad.setCp(1500,3.69* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("H1[|H|==0]", &correctionCheck15, gad); }
  { GAdata gad; gad.setEnthalpy(45.56* 4.18); gad.setEntropy(4.58* 4.18); gad.setCp(100,1.96* 4.18); gad.setCp(200,0.96* 4.18); gad.setCp(298,0.56* 4.18); gad.setCp(400,0.32* 4.18); gad.setCp(500,0.2* 4.18); gad.setCp(600,0.15* 4.18); gad.setCp(700,0.15* 4.18); gad.setCp(800,0.16* 4.18); gad.setCp(900,0.17* 4.18); gad.setCp(1000,0.17* 4.18); gad.setCp(1100,0.17* 4.18); gad.setCp(1200,0.17* 4.18); gad.setCp(1300,0.17* 4.18); gad.setCp(1400,0.16* 4.18); gad.setCp(1500,0.15* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|-{Pt}|==3])", &correctionCheck16, gad); }
  { GAdata gad; gad.setEnthalpy(31.70976* 4.18); gad.setEntropy(3.18768* 4.18); gad.setCp(100,1.36416* 4.18); gad.setCp(200,0.66816* 4.18); gad.setCp(298,0.38976* 4.18); gad.setCp(400,0.22272* 4.18); gad.setCp(500,0.1392* 4.18); gad.setCp(600,0.1044* 4.18); gad.setCp(700,0.1044* 4.18); gad.setCp(800,0.11136* 4.18); gad.setCp(900,0.11832* 4.18); gad.setCp(1000,0.11832* 4.18); gad.setCp(1100,0.11832* 4.18); gad.setCp(1200,0.11832* 4.18); gad.setCp(1300,0.11832* 4.18); gad.setCp(1400,0.11136* 4.18); gad.setCp(1500,0.1044* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|-{Pt}|==2])", &correctionCheck17, gad); }
  { GAdata gad; gad.setEnthalpy(27.70048* 4.18); gad.setEntropy(2.78464* 4.18); gad.setCp(100,1.19168* 4.18); gad.setCp(200,0.58368* 4.18); gad.setCp(298,0.34048* 4.18); gad.setCp(400,0.19456* 4.18); gad.setCp(500,0.1216* 4.18); gad.setCp(600,0.0912* 4.18); gad.setCp(700,0.0912* 4.18); gad.setCp(800,0.09728* 4.18); gad.setCp(900,0.10336* 4.18); gad.setCp(1000,0.10336* 4.18); gad.setCp(1100,0.10336* 4.18); gad.setCp(1200,0.10336* 4.18); gad.setCp(1300,0.10336* 4.18); gad.setCp(1400,0.09728* 4.18); gad.setCp(1500,0.0912* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|-{Pt}|==1][!|%$|>=1])", &correctionCheck18, gad); }
  { GAdata gad; gad.setEnthalpy(27.70048* 4.18); gad.setEntropy(2.78464* 4.18); gad.setCp(100,1.19168* 4.18); gad.setCp(200,0.58368* 4.18); gad.setCp(298,0.34048* 4.18); gad.setCp(400,0.19456* 4.18); gad.setCp(500,0.1216* 4.18); gad.setCp(600,0.0912* 4.18); gad.setCp(700,0.0912* 4.18); gad.setCp(800,0.09728* 4.18); gad.setCp(900,0.10336* 4.18); gad.setCp(1000,0.10336* 4.18); gad.setCp(1100,0.10336* 4.18); gad.setCp(1200,0.10336* 4.18); gad.setCp(1300,0.10336* 4.18); gad.setCp(1400,0.09728* 4.18); gad.setCp(1500,0.0912* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-O2[|-{Pt}|==1])", &correctionCheck19, gad); }
  { GAdata gad; gad.setEnthalpy(30.38852* 4.18); gad.setEntropy(3.05486* 4.18); gad.setCp(100,1.30732* 4.18); gad.setCp(200,0.64032* 4.18); gad.setCp(298,0.37352* 4.18); gad.setCp(400,0.21344* 4.18); gad.setCp(500,0.1334* 4.18); gad.setCp(600,0.10005* 4.18); gad.setCp(700,0.10005* 4.18); gad.setCp(800,0.10672* 4.18); gad.setCp(900,0.11339* 4.18); gad.setCp(1000,0.11339* 4.18); gad.setCp(1100,0.11339* 4.18); gad.setCp(1200,0.11339* 4.18); gad.setCp(1300,0.11339* 4.18); gad.setCp(1400,0.10672* 4.18); gad.setCp(1500,0.10005* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|-{Pt}|==1][|%$|==1])", &correctionCheck20, gad); }
  { GAdata gad; gad.setEnthalpy(13.85024* 4.18); gad.setEntropy(1.39232* 4.18); gad.setCp(100,0.59584* 4.18); gad.setCp(200,0.29184* 4.18); gad.setCp(298,0.17024* 4.18); gad.setCp(400,0.09728* 4.18); gad.setCp(500,0.0608* 4.18); gad.setCp(600,0.0456* 4.18); gad.setCp(700,0.0456* 4.18); gad.setCp(800,0.04864* 4.18); gad.setCp(900,0.05168* 4.18); gad.setCp(1000,0.05168* 4.18); gad.setCp(1100,0.05168* 4.18); gad.setCp(1200,0.05168* 4.18); gad.setCp(1300,0.05168* 4.18); gad.setCp(1400,0.04864* 4.18); gad.setCp(1500,0.0456* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==2](-C2[|-{Pt}|==1][!|%$|>=1])", &correctionCheck21, gad); }
  { GAdata gad; gad.setEnthalpy(13.85024* 4.18); gad.setEntropy(1.39232* 4.18); gad.setCp(100,0.59584* 4.18); gad.setCp(200,0.29184* 4.18); gad.setCp(298,0.17024* 4.18); gad.setCp(400,0.09728* 4.18); gad.setCp(500,0.0608* 4.18); gad.setCp(600,0.0456* 4.18); gad.setCp(700,0.0456* 4.18); gad.setCp(800,0.04864* 4.18); gad.setCp(900,0.05168* 4.18); gad.setCp(1000,0.05168* 4.18); gad.setCp(1100,0.05168* 4.18); gad.setCp(1200,0.05168* 4.18); gad.setCp(1300,0.05168* 4.18); gad.setCp(1400,0.04864* 4.18); gad.setCp(1500,0.0456* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==2](-O2[|-{Pt}|==1])", &correctionCheck22, gad); }
  { GAdata gad; gad.setEnthalpy(16.53828* 4.18); gad.setEntropy(1.66254* 4.18); gad.setCp(100,0.71148* 4.18); gad.setCp(200,0.34848* 4.18); gad.setCp(298,0.20328* 4.18); gad.setCp(400,0.11616* 4.18); gad.setCp(500,0.0726* 4.18); gad.setCp(600,0.05445* 4.18); gad.setCp(700,0.05445* 4.18); gad.setCp(800,0.05808* 4.18); gad.setCp(900,0.06171* 4.18); gad.setCp(1000,0.06171* 4.18); gad.setCp(1100,0.06171* 4.18); gad.setCp(1200,0.06171* 4.18); gad.setCp(1300,0.06171* 4.18); gad.setCp(1400,0.05808* 4.18); gad.setCp(1500,0.05445* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==2](-C2[|-{Pt}|==1][|%$|==1])", &correctionCheck23, gad); }
  { GAdata gad; gad.setEnthalpy(17.85952* 4.18); gad.setEntropy(1.79536* 4.18); gad.setCp(100,0.76832* 4.18); gad.setCp(200,0.37632* 4.18); gad.setCp(298,0.21952* 4.18); gad.setCp(400,0.12544* 4.18); gad.setCp(500,0.0784* 4.18); gad.setCp(600,0.0588* 4.18); gad.setCp(700,0.0588* 4.18); gad.setCp(800,0.06272* 4.18); gad.setCp(900,0.06664* 4.18); gad.setCp(1000,0.06664* 4.18); gad.setCp(1100,0.06664* 4.18); gad.setCp(1200,0.06664* 4.18); gad.setCp(1300,0.06664* 4.18); gad.setCp(1400,0.06272* 4.18); gad.setCp(1500,0.0588* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==2](-C2[|-{Pt}|==2])", &correctionCheck24, gad); }
  { GAdata gad; gad.setEnthalpy(9.88652* 4.18); gad.setEntropy(0.99386* 4.18); gad.setCp(100,0.42532* 4.18); gad.setCp(200,0.20832* 4.18); gad.setCp(298,0.12152* 4.18); gad.setCp(400,0.06944* 4.18); gad.setCp(500,0.0434* 4.18); gad.setCp(600,0.03255* 4.18); gad.setCp(700,0.03255* 4.18); gad.setCp(800,0.03472* 4.18); gad.setCp(900,0.03689* 4.18); gad.setCp(1000,0.03689* 4.18); gad.setCp(1100,0.03689* 4.18); gad.setCp(1200,0.03689* 4.18); gad.setCp(1300,0.03689* 4.18); gad.setCp(1400,0.03472* 4.18); gad.setCp(1500,0.03255* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==1][!|%$|>=1](-C2[|-{Pt}|==1][!|%$|>=1])", &correctionCheck25, gad); }
  { GAdata gad; gad.setEnthalpy(9.88652* 4.18); gad.setEntropy(0.99386* 4.18); gad.setCp(100,0.42532* 4.18); gad.setCp(200,0.20832* 4.18); gad.setCp(298,0.12152* 4.18); gad.setCp(400,0.06944* 4.18); gad.setCp(500,0.0434* 4.18); gad.setCp(600,0.03255* 4.18); gad.setCp(700,0.03255* 4.18); gad.setCp(800,0.03472* 4.18); gad.setCp(900,0.03689* 4.18); gad.setCp(1000,0.03689* 4.18); gad.setCp(1100,0.03689* 4.18); gad.setCp(1200,0.03689* 4.18); gad.setCp(1300,0.03689* 4.18); gad.setCp(1400,0.03472* 4.18); gad.setCp(1500,0.03255* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==1][!|%$|>=1](-O2[|-{Pt}|==1])", &correctionCheck26, gad); }
  { GAdata gad; gad.setEnthalpy(12.529* 4.18); gad.setEntropy(1.2595* 4.18); gad.setCp(100,0.539* 4.18); gad.setCp(200,0.264* 4.18); gad.setCp(298,0.154* 4.18); gad.setCp(400,0.088* 4.18); gad.setCp(500,0.055* 4.18); gad.setCp(600,0.04125* 4.18); gad.setCp(700,0.04125* 4.18); gad.setCp(800,0.044* 4.18); gad.setCp(900,0.04675* 4.18); gad.setCp(1000,0.04675* 4.18); gad.setCp(1100,0.04675* 4.18); gad.setCp(1200,0.04675* 4.18); gad.setCp(1300,0.04675* 4.18); gad.setCp(1400,0.044* 4.18); gad.setCp(1500,0.04125* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==1][!|%$|>=1](-C2[|-{Pt}|==1][|%$|==1])", &correctionCheck27, gad); }
  { GAdata gad; gad.setEnthalpy(9.88652* 4.18); gad.setEntropy(0.99386* 4.18); gad.setCp(100,0.42532* 4.18); gad.setCp(200,0.20832* 4.18); gad.setCp(298,0.12152* 4.18); gad.setCp(400,0.06944* 4.18); gad.setCp(500,0.0434* 4.18); gad.setCp(600,0.03255* 4.18); gad.setCp(700,0.03255* 4.18); gad.setCp(800,0.03472* 4.18); gad.setCp(900,0.03689* 4.18); gad.setCp(1000,0.03689* 4.18); gad.setCp(1100,0.03689* 4.18); gad.setCp(1200,0.03689* 4.18); gad.setCp(1300,0.03689* 4.18); gad.setCp(1400,0.03472* 4.18); gad.setCp(1500,0.03255* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("O1[|-{Pt}|==1](-O2[|-{Pt}|==1])", &correctionCheck28, gad); }
  { GAdata gad; gad.setEnthalpy(12.529* 4.18); gad.setEntropy(1.2595* 4.18); gad.setCp(100,0.539* 4.18); gad.setCp(200,0.264* 4.18); gad.setCp(298,0.154* 4.18); gad.setCp(400,0.088* 4.18); gad.setCp(500,0.055* 4.18); gad.setCp(600,0.04125* 4.18); gad.setCp(700,0.04125* 4.18); gad.setCp(800,0.044* 4.18); gad.setCp(900,0.04675* 4.18); gad.setCp(1000,0.04675* 4.18); gad.setCp(1100,0.04675* 4.18); gad.setCp(1200,0.04675* 4.18); gad.setCp(1300,0.04675* 4.18); gad.setCp(1400,0.044* 4.18); gad.setCp(1500,0.04125* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("O1[|-{Pt}|==1](-C2[|-{Pt}|==1][|%$|==1])", &correctionCheck29, gad); }
  { GAdata gad; gad.setEnthalpy(15.17148* 4.18); gad.setEntropy(1.52514* 4.18); gad.setCp(100,0.65268* 4.18); gad.setCp(200,0.31968* 4.18); gad.setCp(298,0.18648* 4.18); gad.setCp(400,0.10656* 4.18); gad.setCp(500,0.0666* 4.18); gad.setCp(600,0.04995* 4.18); gad.setCp(700,0.04995* 4.18); gad.setCp(800,0.05328* 4.18); gad.setCp(900,0.05661* 4.18); gad.setCp(1000,0.05661* 4.18); gad.setCp(1100,0.05661* 4.18); gad.setCp(1200,0.05661* 4.18); gad.setCp(1300,0.05661* 4.18); gad.setCp(1400,0.05328* 4.18); gad.setCp(1500,0.04995* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==1][|%$|==1](-C2[|-{Pt}|==1][|%$|==1])", &correctionCheck30, gad); }
  { GAdata gad; gad.setEnthalpy(2.64248* 4.18); gad.setEntropy(0.26564* 4.18); gad.setCp(100,0.11368* 4.18); gad.setCp(200,0.05568* 4.18); gad.setCp(298,0.03248* 4.18); gad.setCp(400,0.01856* 4.18); gad.setCp(500,0.0116* 4.18); gad.setCp(600,0.0087* 4.18); gad.setCp(700,0.0087* 4.18); gad.setCp(800,0.00928* 4.18); gad.setCp(900,0.00986* 4.18); gad.setCp(1000,0.00986* 4.18); gad.setCp(1100,0.00986* 4.18); gad.setCp(1200,0.00986* 4.18); gad.setCp(1300,0.00986* 4.18); gad.setCp(1400,0.00928* 4.18); gad.setCp(1500,0.0087* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==2](-C2[|%$|==1](-C3[|-{Pt}|==2]))", &correctionCheck31, gad); }
  { GAdata gad; gad.setEnthalpy(1.32124* 4.18); gad.setEntropy(0.13282* 4.18); gad.setCp(100,0.05684* 4.18); gad.setCp(200,0.02784* 4.18); gad.setCp(298,0.01624* 4.18); gad.setCp(400,0.00928* 4.18); gad.setCp(500,0.0058* 4.18); gad.setCp(600,0.00435* 4.18); gad.setCp(700,0.00435* 4.18); gad.setCp(800,0.00464* 4.18); gad.setCp(900,0.00493* 4.18); gad.setCp(1000,0.00493* 4.18); gad.setCp(1100,0.00493* 4.18); gad.setCp(1200,0.00493* 4.18); gad.setCp(1300,0.00493* 4.18); gad.setCp(1400,0.00464* 4.18); gad.setCp(1500,0.00435* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==2](-C2[|%$|==1](-C3[|-{Pt}|==1][|%$|==1]))", &correctionCheck32, gad); }
  { GAdata gad; gad.setEnthalpy(27.70048* 4.18); gad.setEntropy(2.78464* 4.18); gad.setCp(100,1.19168* 4.18); gad.setCp(200,0.58368* 4.18); gad.setCp(298,0.34048* 4.18); gad.setCp(400,0.19456* 4.18); gad.setCp(500,0.1216* 4.18); gad.setCp(600,0.0912* 4.18); gad.setCp(700,0.0912* 4.18); gad.setCp(800,0.09728* 4.18); gad.setCp(900,0.10336* 4.18); gad.setCp(1000,0.10336* 4.18); gad.setCp(1100,0.10336* 4.18); gad.setCp(1200,0.10336* 4.18); gad.setCp(1300,0.10336* 4.18); gad.setCp(1400,0.09728* 4.18); gad.setCp(1500,0.0912* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[!|%$|>=1](-C3[|-{Pt}|==3]))", &correctionCheck33, gad); }
  { GAdata gad; gad.setEnthalpy(27.70048* 4.18); gad.setEntropy(2.78464* 4.18); gad.setCp(100,1.19168* 4.18); gad.setCp(200,0.58368* 4.18); gad.setCp(298,0.34048* 4.18); gad.setCp(400,0.19456* 4.18); gad.setCp(500,0.1216* 4.18); gad.setCp(600,0.0912* 4.18); gad.setCp(700,0.0912* 4.18); gad.setCp(800,0.09728* 4.18); gad.setCp(900,0.10336* 4.18); gad.setCp(1000,0.10336* 4.18); gad.setCp(1100,0.10336* 4.18); gad.setCp(1200,0.10336* 4.18); gad.setCp(1300,0.10336* 4.18); gad.setCp(1400,0.09728* 4.18); gad.setCp(1500,0.0912* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-O2[!|~{Pt}|>=1](-C3[|-{Pt}|==3]))", &correctionCheck34, gad); }
  { GAdata gad; gad.setEnthalpy(30.38852* 4.18); gad.setEntropy(3.05486* 4.18); gad.setCp(100,1.30732* 4.18); gad.setCp(200,0.64032* 4.18); gad.setCp(298,0.37352* 4.18); gad.setCp(400,0.21344* 4.18); gad.setCp(500,0.1334* 4.18); gad.setCp(600,0.10005* 4.18); gad.setCp(700,0.10005* 4.18); gad.setCp(800,0.10672* 4.18); gad.setCp(900,0.11339* 4.18); gad.setCp(1000,0.11339* 4.18); gad.setCp(1100,0.11339* 4.18); gad.setCp(1200,0.11339* 4.18); gad.setCp(1300,0.11339* 4.18); gad.setCp(1400,0.10672* 4.18); gad.setCp(1500,0.10005* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|%$|==1](-C3[|-{Pt}|==3]))", &correctionCheck35, gad); }
  { GAdata gad; gad.setEnthalpy(13.85024* 4.18); gad.setEntropy(1.39232* 4.18); gad.setCp(100,0.59584* 4.18); gad.setCp(200,0.29184* 4.18); gad.setCp(298,0.17024* 4.18); gad.setCp(400,0.09728* 4.18); gad.setCp(500,0.0608* 4.18); gad.setCp(600,0.0456* 4.18); gad.setCp(700,0.0456* 4.18); gad.setCp(800,0.04864* 4.18); gad.setCp(900,0.05168* 4.18); gad.setCp(1000,0.05168* 4.18); gad.setCp(1100,0.05168* 4.18); gad.setCp(1200,0.05168* 4.18); gad.setCp(1300,0.05168* 4.18); gad.setCp(1400,0.04864* 4.18); gad.setCp(1500,0.0456* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[!|%$|>=1](-C3[|-{Pt}|==2]))", &correctionCheck36, gad); }
  { GAdata gad; gad.setEnthalpy(16.53828* 4.18); gad.setEntropy(1.66254* 4.18); gad.setCp(100,0.71148* 4.18); gad.setCp(200,0.34848* 4.18); gad.setCp(298,0.20328* 4.18); gad.setCp(400,0.11616* 4.18); gad.setCp(500,0.0726* 4.18); gad.setCp(600,0.05445* 4.18); gad.setCp(700,0.05445* 4.18); gad.setCp(800,0.05808* 4.18); gad.setCp(900,0.06171* 4.18); gad.setCp(1000,0.06171* 4.18); gad.setCp(1100,0.06171* 4.18); gad.setCp(1200,0.06171* 4.18); gad.setCp(1300,0.06171* 4.18); gad.setCp(1400,0.05808* 4.18); gad.setCp(1500,0.05445* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|%$|==1](-C3[|-{Pt}|==2]))", &correctionCheck37, gad); }
  { GAdata gad; gad.setEnthalpy(9.88652* 4.18); gad.setEntropy(0.99386* 4.18); gad.setCp(100,0.42532* 4.18); gad.setCp(200,0.20832* 4.18); gad.setCp(298,0.12152* 4.18); gad.setCp(400,0.06944* 4.18); gad.setCp(500,0.0434* 4.18); gad.setCp(600,0.03255* 4.18); gad.setCp(700,0.03255* 4.18); gad.setCp(800,0.03472* 4.18); gad.setCp(900,0.03689* 4.18); gad.setCp(1000,0.03689* 4.18); gad.setCp(1100,0.03689* 4.18); gad.setCp(1200,0.03689* 4.18); gad.setCp(1300,0.03689* 4.18); gad.setCp(1400,0.03472* 4.18); gad.setCp(1500,0.03255* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[!|%$|>=1](-C3[|-{Pt}|==1][!|%$|>=1]))", &correctionCheck38, gad); }
  { GAdata gad; gad.setEnthalpy(9.88652* 4.18); gad.setEntropy(0.99386* 4.18); gad.setCp(100,0.42532* 4.18); gad.setCp(200,0.20832* 4.18); gad.setCp(298,0.12152* 4.18); gad.setCp(400,0.06944* 4.18); gad.setCp(500,0.0434* 4.18); gad.setCp(600,0.03255* 4.18); gad.setCp(700,0.03255* 4.18); gad.setCp(800,0.03472* 4.18); gad.setCp(900,0.03689* 4.18); gad.setCp(1000,0.03689* 4.18); gad.setCp(1100,0.03689* 4.18); gad.setCp(1200,0.03689* 4.18); gad.setCp(1300,0.03689* 4.18); gad.setCp(1400,0.03472* 4.18); gad.setCp(1500,0.03255* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-O2(-C3[|-{Pt}|==1][!|%$|>=1]))", &correctionCheck39, gad); }
  { GAdata gad; gad.setEnthalpy(12.529* 4.18); gad.setEntropy(1.2595* 4.18); gad.setCp(100,0.539* 4.18); gad.setCp(200,0.264* 4.18); gad.setCp(298,0.154* 4.18); gad.setCp(400,0.088* 4.18); gad.setCp(500,0.055* 4.18); gad.setCp(600,0.04125* 4.18); gad.setCp(700,0.04125* 4.18); gad.setCp(800,0.044* 4.18); gad.setCp(900,0.04675* 4.18); gad.setCp(1000,0.04675* 4.18); gad.setCp(1100,0.04675* 4.18); gad.setCp(1200,0.04675* 4.18); gad.setCp(1300,0.04675* 4.18); gad.setCp(1400,0.044* 4.18); gad.setCp(1500,0.04125* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|%$|==1](-C3[|-{Pt}|==1][!|%$|>=1]))", &correctionCheck40, gad); }
  { GAdata gad; gad.setEnthalpy(9.88652* 4.18); gad.setEntropy(0.99386* 4.18); gad.setCp(100,0.42532* 4.18); gad.setCp(200,0.20832* 4.18); gad.setCp(298,0.12152* 4.18); gad.setCp(400,0.06944* 4.18); gad.setCp(500,0.0434* 4.18); gad.setCp(600,0.03255* 4.18); gad.setCp(700,0.03255* 4.18); gad.setCp(800,0.03472* 4.18); gad.setCp(900,0.03689* 4.18); gad.setCp(1000,0.03689* 4.18); gad.setCp(1100,0.03689* 4.18); gad.setCp(1200,0.03689* 4.18); gad.setCp(1300,0.03689* 4.18); gad.setCp(1400,0.03472* 4.18); gad.setCp(1500,0.03255* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[!|%$|>=1](-O3[|-{Pt}|==1]))", &correctionCheck41, gad); }
  { GAdata gad; gad.setEnthalpy(9.88652* 4.18); gad.setEntropy(0.99386* 4.18); gad.setCp(100,0.42532* 4.18); gad.setCp(200,0.20832* 4.18); gad.setCp(298,0.12152* 4.18); gad.setCp(400,0.06944* 4.18); gad.setCp(500,0.0434* 4.18); gad.setCp(600,0.03255* 4.18); gad.setCp(700,0.03255* 4.18); gad.setCp(800,0.03472* 4.18); gad.setCp(900,0.03689* 4.18); gad.setCp(1000,0.03689* 4.18); gad.setCp(1100,0.03689* 4.18); gad.setCp(1200,0.03689* 4.18); gad.setCp(1300,0.03689* 4.18); gad.setCp(1400,0.03472* 4.18); gad.setCp(1500,0.03255* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-O2[!|~{Pt}|>=1](-O3[|-{Pt}|==1]))", &correctionCheck42, gad); }
  { GAdata gad; gad.setEnthalpy(12.529* 4.18); gad.setEntropy(1.2595* 4.18); gad.setCp(100,0.539* 4.18); gad.setCp(200,0.264* 4.18); gad.setCp(298,0.154* 4.18); gad.setCp(400,0.088* 4.18); gad.setCp(500,0.055* 4.18); gad.setCp(600,0.04125* 4.18); gad.setCp(700,0.04125* 4.18); gad.setCp(800,0.044* 4.18); gad.setCp(900,0.04675* 4.18); gad.setCp(1000,0.04675* 4.18); gad.setCp(1100,0.04675* 4.18); gad.setCp(1200,0.04675* 4.18); gad.setCp(1300,0.04675* 4.18); gad.setCp(1400,0.044* 4.18); gad.setCp(1500,0.04125* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|%$|==1](-O3[|-{Pt}|==1]))", &correctionCheck43, gad); }
  { GAdata gad; gad.setEnthalpy(12.529* 4.18); gad.setEntropy(1.2595* 4.18); gad.setCp(100,0.539* 4.18); gad.setCp(200,0.264* 4.18); gad.setCp(298,0.154* 4.18); gad.setCp(400,0.088* 4.18); gad.setCp(500,0.055* 4.18); gad.setCp(600,0.04125* 4.18); gad.setCp(700,0.04125* 4.18); gad.setCp(800,0.044* 4.18); gad.setCp(900,0.04675* 4.18); gad.setCp(1000,0.04675* 4.18); gad.setCp(1100,0.04675* 4.18); gad.setCp(1200,0.04675* 4.18); gad.setCp(1300,0.04675* 4.18); gad.setCp(1400,0.044* 4.18); gad.setCp(1500,0.04125* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-O2(-C3[|-{Pt}|==1][|%$|==1]))", &correctionCheck44, gad); }
  { GAdata gad; gad.setEnthalpy(15.17148* 4.18); gad.setEntropy(1.52514* 4.18); gad.setCp(100,0.65268* 4.18); gad.setCp(200,0.31968* 4.18); gad.setCp(298,0.18648* 4.18); gad.setCp(400,0.10656* 4.18); gad.setCp(500,0.0666* 4.18); gad.setCp(600,0.04995* 4.18); gad.setCp(700,0.04995* 4.18); gad.setCp(800,0.05328* 4.18); gad.setCp(900,0.05661* 4.18); gad.setCp(1000,0.05661* 4.18); gad.setCp(1100,0.05661* 4.18); gad.setCp(1200,0.05661* 4.18); gad.setCp(1300,0.05661* 4.18); gad.setCp(1400,0.05328* 4.18); gad.setCp(1500,0.04995* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|%$|==1](-C3[|-{Pt}|==1][|%$|==1]))", &correctionCheck45, gad); }
  { GAdata gad; gad.setEnthalpy(9.88652* 4.18); gad.setEntropy(0.99386* 4.18); gad.setCp(100,0.42532* 4.18); gad.setCp(200,0.20832* 4.18); gad.setCp(298,0.12152* 4.18); gad.setCp(400,0.06944* 4.18); gad.setCp(500,0.0434* 4.18); gad.setCp(600,0.03255* 4.18); gad.setCp(700,0.03255* 4.18); gad.setCp(800,0.03472* 4.18); gad.setCp(900,0.03689* 4.18); gad.setCp(1000,0.03689* 4.18); gad.setCp(1100,0.03689* 4.18); gad.setCp(1200,0.03689* 4.18); gad.setCp(1300,0.03689* 4.18); gad.setCp(1400,0.03472* 4.18); gad.setCp(1500,0.03255* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[!|%$|>=1](-C3[!|%$|>=1](-C4[|-{Pt}|==3])))", &correctionCheck46, gad); }
  { GAdata gad; gad.setEnthalpy(9.88652* 4.18); gad.setEntropy(0.99386* 4.18); gad.setCp(100,0.42532* 4.18); gad.setCp(200,0.20832* 4.18); gad.setCp(298,0.12152* 4.18); gad.setCp(400,0.06944* 4.18); gad.setCp(500,0.0434* 4.18); gad.setCp(600,0.03255* 4.18); gad.setCp(700,0.03255* 4.18); gad.setCp(800,0.03472* 4.18); gad.setCp(900,0.03689* 4.18); gad.setCp(1000,0.03689* 4.18); gad.setCp(1100,0.03689* 4.18); gad.setCp(1200,0.03689* 4.18); gad.setCp(1300,0.03689* 4.18); gad.setCp(1400,0.03472* 4.18); gad.setCp(1500,0.03255* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[!|%$|>=1](-O3[!|~{Pt}|>=1](-C4[|-{Pt}|==3])))", &correctionCheck47, gad); }
  { GAdata gad; gad.setEnthalpy(12.529* 4.18); gad.setEntropy(1.2595* 4.18); gad.setCp(100,0.539* 4.18); gad.setCp(200,0.264* 4.18); gad.setCp(298,0.154* 4.18); gad.setCp(400,0.088* 4.18); gad.setCp(500,0.055* 4.18); gad.setCp(600,0.04125* 4.18); gad.setCp(700,0.04125* 4.18); gad.setCp(800,0.044* 4.18); gad.setCp(900,0.04675* 4.18); gad.setCp(1000,0.04675* 4.18); gad.setCp(1100,0.04675* 4.18); gad.setCp(1200,0.04675* 4.18); gad.setCp(1300,0.04675* 4.18); gad.setCp(1400,0.044* 4.18); gad.setCp(1500,0.04125* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[!|%$|>=1](-C3[|%$|==1](-C4[|-{Pt}|==3])))", &correctionCheck48, gad); }
  { GAdata gad; gad.setEnthalpy(9.88652* 4.18); gad.setEntropy(0.99386* 4.18); gad.setCp(100,0.42532* 4.18); gad.setCp(200,0.20832* 4.18); gad.setCp(298,0.12152* 4.18); gad.setCp(400,0.06944* 4.18); gad.setCp(500,0.0434* 4.18); gad.setCp(600,0.03255* 4.18); gad.setCp(700,0.03255* 4.18); gad.setCp(800,0.03472* 4.18); gad.setCp(900,0.03689* 4.18); gad.setCp(1000,0.03689* 4.18); gad.setCp(1100,0.03689* 4.18); gad.setCp(1200,0.03689* 4.18); gad.setCp(1300,0.03689* 4.18); gad.setCp(1400,0.03472* 4.18); gad.setCp(1500,0.03255* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-O2[!|~{Pt}|>=1](-O3[!|~{Pt}|>=1](-C4[|-{Pt}|==3])))", &correctionCheck49, gad); }
  { GAdata gad; gad.setEnthalpy(12.529* 4.18); gad.setEntropy(1.2595* 4.18); gad.setCp(100,0.539* 4.18); gad.setCp(200,0.264* 4.18); gad.setCp(298,0.154* 4.18); gad.setCp(400,0.088* 4.18); gad.setCp(500,0.055* 4.18); gad.setCp(600,0.04125* 4.18); gad.setCp(700,0.04125* 4.18); gad.setCp(800,0.044* 4.18); gad.setCp(900,0.04675* 4.18); gad.setCp(1000,0.04675* 4.18); gad.setCp(1100,0.04675* 4.18); gad.setCp(1200,0.04675* 4.18); gad.setCp(1300,0.04675* 4.18); gad.setCp(1400,0.044* 4.18); gad.setCp(1500,0.04125* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-O2[!|~{Pt}|>=1](-C3[|%$|==1](-C4[|-{Pt}|==3])))", &correctionCheck50, gad); }
  { GAdata gad; gad.setEnthalpy(15.17148* 4.18); gad.setEntropy(1.52514* 4.18); gad.setCp(100,0.65268* 4.18); gad.setCp(200,0.31968* 4.18); gad.setCp(298,0.18648* 4.18); gad.setCp(400,0.10656* 4.18); gad.setCp(500,0.0666* 4.18); gad.setCp(600,0.04995* 4.18); gad.setCp(700,0.04995* 4.18); gad.setCp(800,0.05328* 4.18); gad.setCp(900,0.05661* 4.18); gad.setCp(1000,0.05661* 4.18); gad.setCp(1100,0.05661* 4.18); gad.setCp(1200,0.05661* 4.18); gad.setCp(1300,0.05661* 4.18); gad.setCp(1400,0.05328* 4.18); gad.setCp(1500,0.04995* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|%$|==1](-C3[|%$|==1](-C4[|-{Pt}|==3])))", &correctionCheck51, gad); }
  { GAdata gad; gad.setEnthalpy(1.32124* 4.18); gad.setEntropy(0.13282* 4.18); gad.setCp(100,0.05684* 4.18); gad.setCp(200,0.02784* 4.18); gad.setCp(298,0.01624* 4.18); gad.setCp(400,0.00928* 4.18); gad.setCp(500,0.0058* 4.18); gad.setCp(600,0.00435* 4.18); gad.setCp(700,0.00435* 4.18); gad.setCp(800,0.00464* 4.18); gad.setCp(900,0.00493* 4.18); gad.setCp(1000,0.00493* 4.18); gad.setCp(1100,0.00493* 4.18); gad.setCp(1200,0.00493* 4.18); gad.setCp(1300,0.00493* 4.18); gad.setCp(1400,0.00464* 4.18); gad.setCp(1500,0.00435* 4.18);   gad.PrepareGA(); ThermoGA::AddCorrections("C1[|-{Pt}|==3](-C2[|%$|==1](-C3[|%$|==1](-C4[|-{Pt}|==2])))", &correctionCheck52, gad); }
  
   if (reactantlist.size() == 0)
   {
      cout<<"No reactants specified!"<<endl;
      return -1;
   }

   try
   {
      list<string>::iterator it;
      for (it=reactantlist.begin();it!=reactantlist.end();it++)
      {
          throwsmileserror(*it);
      }
      rxn_net_gen test_network;
      test_network.AddInitialReactants(reactantlist);
      test_network.AddReactionRules(Rtypelist);
      test_network.AddGlobalConstraints(&checkCglobalConstraints);
      test_network.AddLumpingStrategy(L);
      test_network.AddCompositeSites(CompositeSiteList);
      test_network.SetCalcThermo(true);
      test_network.GenerateNetwork();
      if(L.shoudLump()) test_network.SetAllLumpedRxns();
      test_network.print_rxnlist();
      int d = 0;
   }
   catch (pair<string, int> error)
   {
      catchsmileserror (error);
   }
   catch (pair<string,string> c)
   {
      cout<<"SMILES "<<c.first<<"has an unidentified character: "<<c.second<<endl;
   }
   catch (int er)
   {
      if (er==1) cout<<"exiting execution!"<<endl;
   }
   return 0;
}
