import pyring
sf = pyring.SingleFunc()
df = pyring.DoubleFunc()
def checkCCRt_O2Adsorbm0(m0):
	return True
def checkCCRt_O2Adsorbm1(m1):
	return True
def checkCCRt_O2Adsorbmboth(m0, m1):
	return True
def Rt_O2Adsorbprodcons(m0):
	return True
def checkCCRt_CCScissionm0(m0):
	return True
def checkCCRt_CCScissionm1(m1):
	return True
def checkCCRt_CCScissionmboth(m0, m1):
	return True
def Rt_CCScissionprodcons(m0):
	return True
def checkCCRt_CHScissionm0(m0):
	return True
def checkCCRt_CHScissionm1(m1):
	return True
def checkCCRt_CHScissionmboth(m0, m1):
	return True
def Rt_CHScissionprodcons(m0):
	return True
def checkCCRt_OHScissionm0(m0):
	return True
def checkCCRt_OHScissionm1(m1):
	return True
def checkCCRt_OHScissionmboth(m0, m1):
	return True
def Rt_OHScissionprodcons(m0):
	return True
def checkCCRt_HTransferCtoOm0(m0):
	return True
def checkCCRt_HTransferCtoOm1(m1):
	return True
def checkCCRt_HTransferCtoOmboth(m0, m1):
	return True
def Rt_HTransferCtoOprodcons(m0):
	return True
def checkCCRt_HTransferOHtoOm0(m0):
	return True
def checkCCRt_HTransferOHtoOm1(m1):
	return True
def checkCCRt_HTransferOHtoOmboth(m0, m1):
	return True
def Rt_HTransferOHtoOprodcons(m0):
	return True
def checkCCRt_COFormationfromOHm0(m0):
	return pyring.Patternmatch(m0, pyring.Substructure("O1(-C2)", pyring.patternsize("O1(-C2)")), 1).GetDistinctMatches()==0 and True
def checkCCRt_COFormationfromOHm1(m1):
	return True
def checkCCRt_COFormationfromOHmboth(m0, m1):
	return True
def Rt_COFormationfromOHprodcons(m0):
	return True
def checkCCRt_COFormationfromOm0(m0):
	return pyring.Patternmatch(m0, pyring.Substructure("O1(-C2)", pyring.patternsize("O1(-C2)")), 1).GetDistinctMatches()==0 and True
def checkCCRt_COFormationfromOm1(m1):
	return True
def checkCCRt_COFormationfromOmboth(m0, m1):
	return True
def Rt_COFormationfromOprodcons(m0):
	return True
def checkCCRt_H2OwCFormationm0(m0):
	return True
def checkCCRt_H2OwCFormationm1(m1):
	return True
def checkCCRt_H2OwCFormationmboth(m0, m1):
	return True
def Rt_H2OwCFormationprodcons(m0):
	return True
def checkCCRt_H2OwOFormationm0(m0):
	return True
def checkCCRt_H2OwOFormationm1(m1):
	return True
def checkCCRt_H2OwOFormationmboth(m0, m1):
	return True
def Rt_H2OwOFormationprodcons(m0):
	return True
def checkCglobalConstraints(m0):
	return True

reactantlist = list()
CompositeAtomsList=[]
CompositeSiteList=[]
Rtypelist = []
reactantlist.append("CCC")
reactantlist.append("O=O")
reactantlist.append("{Pt}")
CompositeAtomsList.append("{Pt}")

Rt_O2Adsorb = pyring.Reactiontype()
Rt_O2Adsorb.add_reactant_pattern(pyring.Substructure("O1(=O2)", pyring.patternsize("O1(=O2)")))
Rt_O2Adsorb.add_reactant_pattern(pyring.Substructure("{Pt}3[!|~$|>0]", pyring.patternsize("{Pt}3[!|~$|>0]")))
Rt_O2Adsorbdupmap2= {}
Rt_O2Adsorbdupmap2[3] = 4
Rt_O2Adsorb.AddReactantCopy(1, Rt_O2Adsorbdupmap2)
Rt_O2Adsorbdupmap3= {}
Rt_O2Adsorbdupmap3[3] = 5
Rt_O2Adsorb.AddReactantCopy(1, Rt_O2Adsorbdupmap3)
Rt_O2Adsorbdupmap4= {}
Rt_O2Adsorbdupmap4[3] = 6
Rt_O2Adsorb.AddReactantCopy(1, Rt_O2Adsorbdupmap4)
Rt_O2Adsorb.add_reactantconstraint(sf.singlef(checkCCRt_O2Adsorbm0))
Rt_O2Adsorb.add_reactantconstraint(sf.singlef(checkCCRt_O2Adsorbm1))
Rt_O2Adsorb.add_combined_constraint(df.doublef(checkCCRt_O2Adsorbmboth))
Rt_O2Adsorb.disconnect_bond(1,2)
Rt_O2Adsorb.connect_bond(1,3)
Rt_O2Adsorb.connect_bond(1,4)
Rt_O2Adsorb.connect_bond(2,5)
Rt_O2Adsorb.connect_bond(2,6)
Rt_O2Adsorb.add_productconstraint(sf.singlef(Rt_O2Adsorbprodcons))
Rt_O2Adsorb.setRuleName("O2Adsorb")
Rtypelist.append(Rt_O2Adsorb) 

Rt_CCScission = pyring.Reactiontype()
Rt_CCScission.add_reactant_pattern(pyring.Substructure("C1(-C2)", pyring.patternsize("C1(-C2)")))
Rt_CCScission.add_reactant_pattern(pyring.Substructure("{Pt}3[!|~$|>0]", pyring.patternsize("{Pt}3[!|~$|>0]")))
Rt_CCScissiondupmap2= {}
Rt_CCScissiondupmap2[3] = 4
Rt_CCScission.AddReactantCopy(1, Rt_CCScissiondupmap2)
Rt_CCScission.add_reactantconstraint(sf.singlef(checkCCRt_CCScissionm0))
Rt_CCScission.add_reactantconstraint(sf.singlef(checkCCRt_CCScissionm1))
Rt_CCScission.add_combined_constraint(df.doublef(checkCCRt_CCScissionmboth))
Rt_CCScission.disconnect_bond(1,2)
Rt_CCScission.connect_bond(1,3)
Rt_CCScission.connect_bond(2,4)
Rt_CCScission.add_productconstraint(sf.singlef(Rt_CCScissionprodcons))
Rt_CCScission.setRuleName("CCScission")
Rtypelist.append(Rt_CCScission) 

Rt_CHScission = pyring.Reactiontype()
Rt_CHScission.add_reactant_pattern(pyring.Substructure("C1(-H2)", pyring.patternsize("C1(-H2)")))
Rt_CHScission.add_reactant_pattern(pyring.Substructure("{Pt}3[!|~$|>0]", pyring.patternsize("{Pt}3[!|~$|>0]")))
Rt_CHScissiondupmap2= {}
Rt_CHScissiondupmap2[3] = 4
Rt_CHScission.AddReactantCopy(1, Rt_CHScissiondupmap2)
Rt_CHScission.add_reactantconstraint(sf.singlef(checkCCRt_CHScissionm0))
Rt_CHScission.add_reactantconstraint(sf.singlef(checkCCRt_CHScissionm1))
Rt_CHScission.add_combined_constraint(df.doublef(checkCCRt_CHScissionmboth))
Rt_CHScission.disconnect_bond(1,2)
Rt_CHScission.connect_bond(1,3)
Rt_CHScission.connect_bond(2,4)
Rt_CHScission.add_productconstraint(sf.singlef(Rt_CHScissionprodcons))
Rt_CHScission.setRuleName("CHScission")
Rtypelist.append(Rt_CHScission) 

Rt_OHScission = pyring.Reactiontype()
Rt_OHScission.add_reactant_pattern(pyring.Substructure("O1(-H2)", pyring.patternsize("O1(-H2)")))
Rt_OHScission.add_reactant_pattern(pyring.Substructure("{Pt}3[!|~$|>0]", pyring.patternsize("{Pt}3[!|~$|>0]")))
Rt_OHScissiondupmap2= {}
Rt_OHScissiondupmap2[3] = 4
Rt_OHScission.AddReactantCopy(1, Rt_OHScissiondupmap2)
Rt_OHScission.add_reactantconstraint(sf.singlef(checkCCRt_OHScissionm0))
Rt_OHScission.add_reactantconstraint(sf.singlef(checkCCRt_OHScissionm1))
Rt_OHScission.add_combined_constraint(df.doublef(checkCCRt_OHScissionmboth))
Rt_OHScission.disconnect_bond(1,2)
Rt_OHScission.connect_bond(1,3)
Rt_OHScission.connect_bond(2,4)
Rt_OHScission.add_productconstraint(sf.singlef(Rt_OHScissionprodcons))
Rt_OHScission.setRuleName("OHScission")
Rtypelist.append(Rt_OHScission) 

Rt_HTransferCtoO = pyring.Reactiontype()
Rt_HTransferCtoO.add_reactant_pattern(pyring.Substructure("C1(-H2)", pyring.patternsize("C1(-H2)")))
Rt_HTransferCtoO.add_reactant_pattern(pyring.Substructure("{Pt}4(-O3(-{Pt}5))", pyring.patternsize("{Pt}4(-O3(-{Pt}5))")))
Rt_HTransferCtoO.add_reactantconstraint(sf.singlef(checkCCRt_HTransferCtoOm0))
Rt_HTransferCtoO.add_reactantconstraint(sf.singlef(checkCCRt_HTransferCtoOm1))
Rt_HTransferCtoO.add_combined_constraint(df.doublef(checkCCRt_HTransferCtoOmboth))
Rt_HTransferCtoO.disconnect_bond(1,2)
Rt_HTransferCtoO.disconnect_bond(3,4)
Rt_HTransferCtoO.connect_bond(1,4)
Rt_HTransferCtoO.connect_bond(3,2)
Rt_HTransferCtoO.add_productconstraint(sf.singlef(Rt_HTransferCtoOprodcons))
Rt_HTransferCtoO.setRuleName("HTransferCtoO")
Rtypelist.append(Rt_HTransferCtoO) 

Rt_HTransferOHtoO = pyring.Reactiontype()
Rt_HTransferOHtoO.add_reactant_pattern(pyring.Substructure("O1(-H2)", pyring.patternsize("O1(-H2)")))
Rt_HTransferOHtoO.add_reactant_pattern(pyring.Substructure("{Pt}4(-O3(-{Pt}5))", pyring.patternsize("{Pt}4(-O3(-{Pt}5))")))
Rt_HTransferOHtoO.add_reactantconstraint(sf.singlef(checkCCRt_HTransferOHtoOm0))
Rt_HTransferOHtoO.add_reactantconstraint(sf.singlef(checkCCRt_HTransferOHtoOm1))
Rt_HTransferOHtoO.add_combined_constraint(df.doublef(checkCCRt_HTransferOHtoOmboth))
Rt_HTransferOHtoO.disconnect_bond(1,2)
Rt_HTransferOHtoO.disconnect_bond(3,4)
Rt_HTransferOHtoO.connect_bond(1,4)
Rt_HTransferOHtoO.connect_bond(3,2)
Rt_HTransferOHtoO.add_productconstraint(sf.singlef(Rt_HTransferOHtoOprodcons))
Rt_HTransferOHtoO.setRuleName("HTransferOHtoO")
Rtypelist.append(Rt_HTransferOHtoO) 

Rt_COFormationfromOH = pyring.Reactiontype()
Rt_COFormationfromOH.add_reactant_pattern(pyring.Substructure("{Pt}2(-C1[!|~C|>1])", pyring.patternsize("{Pt}2(-C1[!|~C|>1])")))
Rt_COFormationfromOH.add_reactant_pattern(pyring.Substructure("{Pt}5(-O3(-H4))", pyring.patternsize("{Pt}5(-O3(-H4))")))
Rt_COFormationfromOH.add_reactantconstraint(sf.singlef(checkCCRt_COFormationfromOHm0))
Rt_COFormationfromOH.add_reactantconstraint(sf.singlef(checkCCRt_COFormationfromOHm1))
Rt_COFormationfromOH.add_combined_constraint(df.doublef(checkCCRt_COFormationfromOHmboth))
Rt_COFormationfromOH.disconnect_bond(3,5)
Rt_COFormationfromOH.disconnect_bond(1,2)
Rt_COFormationfromOH.connect_bond(1,3)
Rt_COFormationfromOH.add_productconstraint(sf.singlef(Rt_COFormationfromOHprodcons))
Rt_COFormationfromOH.setRuleName("COFormationfromOH")
Rtypelist.append(Rt_COFormationfromOH) 

Rt_COFormationfromO = pyring.Reactiontype()
Rt_COFormationfromO.add_reactant_pattern(pyring.Substructure("{Pt}2(-C1[!|~C|>1])", pyring.patternsize("{Pt}2(-C1[!|~C|>1])")))
Rt_COFormationfromO.add_reactant_pattern(pyring.Substructure("{Pt}4(-O3(-{Pt}5))", pyring.patternsize("{Pt}4(-O3(-{Pt}5))")))
Rt_COFormationfromO.add_reactantconstraint(sf.singlef(checkCCRt_COFormationfromOm0))
Rt_COFormationfromO.add_reactantconstraint(sf.singlef(checkCCRt_COFormationfromOm1))
Rt_COFormationfromO.add_combined_constraint(df.doublef(checkCCRt_COFormationfromOmboth))
Rt_COFormationfromO.disconnect_bond(3,4)
Rt_COFormationfromO.disconnect_bond(1,2)
Rt_COFormationfromO.connect_bond(1,3)
Rt_COFormationfromO.add_productconstraint(sf.singlef(Rt_COFormationfromOprodcons))
Rt_COFormationfromO.setRuleName("COFormationfromO")
Rtypelist.append(Rt_COFormationfromO) 

Rt_H2OwCFormation = pyring.Reactiontype()
Rt_H2OwCFormation.add_reactant_pattern(pyring.Substructure("C1(-H2)", pyring.patternsize("C1(-H2)")))
Rt_H2OwCFormation.add_reactant_pattern(pyring.Substructure("{Pt}5(-O3(-H4))", pyring.patternsize("{Pt}5(-O3(-H4))")))
Rt_H2OwCFormation.add_reactantconstraint(sf.singlef(checkCCRt_H2OwCFormationm0))
Rt_H2OwCFormation.add_reactantconstraint(sf.singlef(checkCCRt_H2OwCFormationm1))
Rt_H2OwCFormation.add_combined_constraint(df.doublef(checkCCRt_H2OwCFormationmboth))
Rt_H2OwCFormation.disconnect_bond(5,3)
Rt_H2OwCFormation.disconnect_bond(2,1)
Rt_H2OwCFormation.connect_bond(3,2)
Rt_H2OwCFormation.connect_bond(1,5)
Rt_H2OwCFormation.add_productconstraint(sf.singlef(Rt_H2OwCFormationprodcons))
Rt_H2OwCFormation.setRuleName("H2OwCFormation")
Rtypelist.append(Rt_H2OwCFormation) 

Rt_H2OwOFormation = pyring.Reactiontype()
Rt_H2OwOFormation.add_reactant_pattern(pyring.Substructure("O1(-H2)", pyring.patternsize("O1(-H2)")))
Rt_H2OwOFormation.add_reactant_pattern(pyring.Substructure("{Pt}5(-O3(-H4))", pyring.patternsize("{Pt}5(-O3(-H4))")))
Rt_H2OwOFormation.add_reactantconstraint(sf.singlef(checkCCRt_H2OwOFormationm0))
Rt_H2OwOFormation.add_reactantconstraint(sf.singlef(checkCCRt_H2OwOFormationm1))
Rt_H2OwOFormation.add_combined_constraint(df.doublef(checkCCRt_H2OwOFormationmboth))
Rt_H2OwOFormation.disconnect_bond(5,3)
Rt_H2OwOFormation.disconnect_bond(2,1)
Rt_H2OwOFormation.connect_bond(3,2)
Rt_H2OwOFormation.connect_bond(1,5)
Rt_H2OwOFormation.add_productconstraint(sf.singlef(Rt_H2OwOFormationprodcons))
Rt_H2OwOFormation.setRuleName("H2OwOFormation")
Rtypelist.append(Rt_H2OwOFormation) 
L = pyring.LumpingStrategy(False)
test_network = pyring.rxn_net_gen()
test_network.AddInitialReactants(reactantlist)
test_network.AddReactionRules(Rtypelist)
test_network.AddGlobalConstraints(sf.singlef(checkCglobalConstraints))
test_network.AddLumpingStrategy(L)
test_network.AddCompositeSites(CompositeSiteList)
test_network.SetCalcThermo(False)
test_network.GenerateNetwork()
if(L.shoudLump()):
	test_network.SetAllLumpedRxns()
test_network.print_rxnlist()
