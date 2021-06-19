<%
cfg['sources'] = ['../src/common.cpp', '../src/additionalfunc.cpp', '../src/element.cpp', '../src/atom.cpp','../src/stringreg.cpp', '../src/atomcontainer.cpp','../src/molecule.cpp','../src/substructure.cpp','../src/generated_rxn.cpp','../src/Patternmatch.cpp','../src/groupadditivity.cpp','../src/automorphs.cpp','../src/reaction.cpp','../src/lumping.cpp','../src/rng.cpp']
setup_pybind11(cfg)
%>

#include <string>
#include "../src/common.h"
#include "../src/additionalfunc.h"
#include "../src/element.h"
#include "../src/clonable.h"
#include "../src/atom.h"
#include "../src/stringreg.h"
#include "../src/atomcontainer.h"
#include "../src/molecule.h"
#include "../src/substructure.h"
#include "../src/generated_rxn.h"
#include "../src/Patternmatch.h"
#include "../src/groupadditivity.h"
#include "../src/automorphs.h"
#include "../src/reaction.h"
#include "../src/lumping.h"
#include "../src/rng.h"
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/chrono.h>

//m.def("ConstrPtr", &ConstrPtr);
//m.def("CombinedConstrPtr", &CombinedConstrPtr);

struct SingleFunc {
	std::function<bool(Molecule &)> singlef(std::function<bool(Molecule &)>func){
		return func;
	}
};

struct DoubleFunc {
	std::function<bool(Molecule &, Molecule &)> doublef(std::function<bool(Molecule &, Molecule &)>func){
		return func;
	}
};


//struct SingleMolecule {
//	bool singleconstr(std::function<bool(Molecule &)>func, Molecule &a) {
//		return func(a);
//	}
//};
//
//struct DoubleMolecule {
//	bool doubleconstr(std::function<bool(Molecule &, Molecule &)>func, Molecule &a, Molecule &b) {
//		return func(a,b);
//	}
//};

namespace py = pybind11;

PYBIND11_MODULE(pyring, m) {
	
	//common.cpp classes
	py::class_<Path>(m, "Path")
		.def(py::init<>())
		.def("add", &Path::add)
		.def("get", &Path::get)
		.def("reverse", &Path::reverse)
		.def("print", &Path::print)
		.def("size", &Path::size)
		.def("getfirst", &Path::getfirst)
		.def("getlast", &Path::getlast)
		.def("remove", &Path::remove)
		.def("join", &Path::join)
		.def("intersection", &Path::intersection)
		.def("clear", &Path::clear)
		.def("contains", &Path::contains);
	py::class_<Ringset>(m, "Ringset")
		.def("add", &Ringset::add)
		.def("get", &Ringset::get)
		.def("remove", &Ringset::remove)
		.def("print", &Ringset::print)
		.def("print_all", &Ringset::print_all)
		.def("size", &Ringset::size)
		.def("sizesort", &Ringset::sizesort)
		.def("is_present", &Ringset::is_present)
		.def("CountAtom", &Ringset::CountAtom)
		.def("SizeContainingAtom", &Ringset::SizeContainingAtom);
	py::class_<env_tuple>(m, "env_tuple")
		.def("set_flag", &env_tuple::set_flag)
		.def("get_flag", &env_tuple::get_flag)
		.def("set_type", &env_tuple::set_type)
		.def("get_type", &env_tuple::get_type)
		.def("get_freq", &env_tuple::get_freq)
		.def("set_freq", &env_tuple::set_freq)
		.def("get_comp_opr", &env_tuple::get_comp_opr)
		.def("set_comp_opr", &env_tuple::set_comp_opr);
	py::class_<env_set>(m, "env_set")
		.def("add", &env_set::add)
		.def("get", &env_set::get)
		.def("remove", &env_set::remove)
		.def("set_is_empty", &env_set::set_is_empty)
		.def("size", &env_set::size);
	
	//additionalfunc.cpp functions
	m.def("prime", &prime);
	m.def("factorial", &factorial);
	m.def("moleculesize", &moleculesize);
	m.def("patternsize", &patternsize);
	m.def("throwsmileserror", &throwsmileserror);
	m.def("catchsmileserror", &catchsmileserror);
	m.def("TokenizeIntoWords", &TokenizeIntoWords);
	m.def("StringToDouble", &StringToDouble);
	m.def("StringToInt", &StringToInt);
	m.def("DoubleToString", &DoubleToString);
	m.def("IntToStr", &IntToStr);
	m.def("WriteInfoFile", &WriteInfoFile);
	
	//element.cpp classes
	py::class_<Element>(m, "Element")
		.def(py::init<char>())
		.def("set_element_properties", &Element::set_element_properties);
	py::class_<Atomtype>(m, "Atomtype")
		.def(py::init<std::string, char>());
	
	//atom.cpp classes
	py::class_<Atom>(m, "Atom")
		.def(py::init<int>())
		.def(py::init<>())
		.def("get_nature", &Atom::get_nature)
		.def("set_nature", &Atom::set_nature)
		.def("MakeSymbolAliphatic", &Atom::MakeSymbolAliphatic)
		.def("MakeSymbolAromatic", &Atom::MakeSymbolAromatic)
		.def("IsSymbolAliphatic", &Atom::IsSymbolAliphatic);
	py::class_<SingleAtom>(m, "SingleAtom")
		.def(py::init<std::string&, std::string&, int, char>())
		.def("set_atom_symbol", &SingleAtom::set_atom_symbol)
		.def("get_atom_symbol", &SingleAtom::get_atom_symbol)
		.def("get_element_name", &SingleAtom::get_element_name)
		.def("set_isotope_number", &SingleAtom::set_isotope_number)
		.def("get_isotope_number", &SingleAtom::get_isotope_number)
		.def("set_valency_value", &SingleAtom::set_valency_value)
		.def("readjustatomproperties", &SingleAtom::readjustatomproperties)
		.def("DropRingIdentifier", &SingleAtom::DropRingIdentifier)
		.def("reset_properties", &SingleAtom::reset_properties)
		.def("set_valency", &SingleAtom::set_valency)
		.def("get_valency", &SingleAtom::get_valency)
		.def("get_charge", &SingleAtom::get_charge)
		.def("set_initial_properties", &SingleAtom::set_initial_properties)
		.def("set_atomtype_name", &SingleAtom::set_atomtype_name)
		.def("get_atomtype_name", &SingleAtom::get_atomtype_name)
		.def("get_up", &SingleAtom::get_up)
		.def("get_lp", &SingleAtom::get_lp)
		.def("get_atomic_number", &SingleAtom::get_atomic_number)
		.def("get_n_mass_number", &SingleAtom::get_n_mass_number)
		.def("get_isotope_mass_number", &SingleAtom::get_isotope_mass_number)
		.def("set_charge", &SingleAtom::set_charge);
	py::class_<CompositeAtom>(m, "CompositeAtom")
		.def(py::init<std::string, std::string>())
		.def("get_element_name", &CompositeAtom::get_element_name)
		.def("get_atom_symbol", &CompositeAtom::get_atom_symbol)
		.def("set_atom_symbol", &CompositeAtom::set_atom_symbol)
		.def("set_charge", &CompositeAtom::set_charge)
		.def("get_charge", &CompositeAtom::get_charge)
		.def("get_valency", &CompositeAtom::get_valency)
		.def("set_valency", &CompositeAtom::set_valency)
		.def("set_atomtype_name", &CompositeAtom::set_atomtype_name)
		.def("get_atomtype_name", &CompositeAtom::get_atomtype_name)
		.def("reset_properties", &CompositeAtom::reset_properties)
		.def("set_initial_properties", &CompositeAtom::set_initial_properties)
		.def("set_valency_value", &CompositeAtom::set_valency_value)
		.def("readjustatomproperties", &CompositeAtom::readjustatomproperties)
		.def("get_lp", &CompositeAtom::get_lp)
		.def("get_up", &CompositeAtom::get_up)
		.def("DropRingIdentifier", &CompositeAtom::DropRingIdentifier);
	
	//StringRegistry.cpp classes	
	py::class_<StringRegistry>(m, "StringRegistry")
		.def("InRegistry", &StringRegistry::InRegistry)
		.def("InsertIntoRegistry", &StringRegistry::InsertIntoRegistry)
		.def("getStringPointer", &StringRegistry::getStringPointer)
		.def("RemoveFromRegistry", &StringRegistry::RemoveFromRegistry);
		
	py::class_<CompositeAtomsRegistry>(m, "CompositeAtomsRegistry")
		.def("InsertIntoList", &CompositeAtomsRegistry::InsertIntoList)
		.def("getIndexOfAtom", &CompositeAtomsRegistry::getIndexOfAtom);
	
	//Atomcontainer.cpp classes
	py::class_<Atomcontainer>(m, "Atomcontainer")
		.def(py::init<int>())
		.def(py::init<>())
		.def(py::init<Atomcontainer&>())
		.def("operator=", &Atomcontainer::operator=)
		.def("setrank", &Atomcontainer::setrank)
		.def("getBplusNBatomcount", &Atomcontainer::getBplusNBatomcount)
		.def("getrank", &Atomcontainer::getrank)
		.def("getclassRank", &Atomcontainer::getclassRank)
		.def("isaromatic", &Atomcontainer::isaromatic)
		.def("get_adjacency", &Atomcontainer::get_adjacency)
		.def("get_BO", &Atomcontainer::get_BO)
		.def("getsize", &Atomcontainer::getsize)
		.def("getNN", &Atomcontainer::getNN)
		.def("setBO", &Atomcontainer::setBO)
		.def("dbcount", &Atomcontainer::dbcount)
		.def("tpcount", &Atomcontainer::tpcount)
		.def("find_BO", &Atomcontainer::find_BO)
		.def("getatomtype", &Atomcontainer::getatomtype)
		.def("print_adjacency_list", &Atomcontainer::print_adjacency_list)
		.def("getatom", &Atomcontainer::getatom)
		.def("getvalue", &Atomcontainer::getvalue)
		.def("getHydrogens", &Atomcontainer::getHydrogens)
		.def("merge", &Atomcontainer::merge)
		.def("resetlabel", &Atomcontainer::resetlabel)
		.def("findatomwithlabel", &Atomcontainer::findatomwithlabel)
		.def("getlabel", &Atomcontainer::getlabel)
		.def("setlabel", &Atomcontainer::setlabel)
		.def("formbond", &Atomcontainer::formbond)
		.def("breakbond", &Atomcontainer::breakbond)
		.def("changeBO", &Atomcontainer::changeBO)
		.def("connectedcomponents", &Atomcontainer::connectedcomponents)
		.def("setatomtypename", &Atomcontainer::setatomtypename)
		.def("setatomsymbol", &Atomcontainer::setatomsymbol)
		.def("setInitialAtomproperties", &Atomcontainer::setInitialAtomproperties)
		.def("setatomvalency", &Atomcontainer::setatomvalency)
		.def("AromaticBondCount", &Atomcontainer::AromaticBondCount)
		.def("IsAdjacentAtomWithBO", &Atomcontainer::IsAdjacentAtomWithBO)
		.def("getBatomcount", &Atomcontainer::getBatomcount)
		.def("HasNBinteractions", &Atomcontainer::HasNBinteractions)
		.def("getAromaticAtoms", &Atomcontainer::getAromaticAtoms)
		.def("getHlabel", &Atomcontainer::getHlabel)
		.def("setHlabel", &Atomcontainer::setHlabel)
		.def("findatomwithHlabel", &Atomcontainer::findatomwithHlabel)
		.def("changeHlabel", &Atomcontainer::changeHlabel)
		.def("removeHlabel", &Atomcontainer::removeHlabel)
		.def("addHlabel", &Atomcontainer::addHlabel)
		.def("setHydrogens", &Atomcontainer::setHydrogens)
		.def("getElectronicHashValue", &Atomcontainer::getElectronicHashValue)
		.def("getTotalElectronicValue", &Atomcontainer::getTotalElectronicValue)
		.def("addAtom", &Atomcontainer::addAtom)
		.def("getGroupHash", &Atomcontainer::getGroupHash)
		.def("AtomValueForGA", &Atomcontainer::AtomValueForGA)
		.def("getNNElecHash", &Atomcontainer::getNNElecHash)
		.def("getNNDoubleBondCount", &Atomcontainer::getNNDoubleBondCount)
		.def("getNNTripleBondCount", &Atomcontainer::getNNTripleBondCount)
		.def("RingCountOfAtom", &Atomcontainer::RingCountOfAtom)
		.def("getClassesFreqMapNeighboringAtom", &Atomcontainer::getClassesFreqMapNeighboringAtom)
		.def("AtomCenteredGroupForGA", &Atomcontainer::AtomCenteredGroupForGA)
		.def("getNNDoubleTripleBondFactors", &Atomcontainer::getNNDoubleTripleBondFactors)
		.def("GetElements", &Atomcontainer::GetElements);
		
	//molecule.cpp classes
	py::class_<Molecule>(m, "Molecule")
		.def(py::init<>())
		.def(py::init<std::string, int>())
		.def(py::init<Atomcontainer&>());
	
	//substructure.cpp classes
	py::class_<Substructure>(m, "Substructure")
		.def(py::init<std::string, int>())
		.def("getstring", &Substructure::getstring)
		.def("printstring", &Substructure::printstring)
		.def("isringbondcheck", &Substructure::isringbondcheck);	
	
	//generated_rxn.cpp classes
	py::class_<generated_rxn>(m, "generated_rxn")
		.def(py::init<>())
		.def("IdenticalProductsFactor", &generated_rxn::IdenticalProductsFactor)
		.def("setIntramolecularity", &generated_rxn::setIntramolecularity)
		.def("setAtomsFromOthers", &generated_rxn::setAtomsFromOthers)
		.def("set_reactionMF", &generated_rxn::set_reactionMF)
		.def("get_reactionMF", &generated_rxn::get_reactionMF)
		.def("add_reactants", &generated_rxn::add_reactants)
		.def("add_products", &generated_rxn::add_products)
		.def("get_frequency", &generated_rxn::get_frequency)
		.def("set_frequency", &generated_rxn::set_frequency);
		
	//Patternmatch.cpp classes
	py::class_<Patternmatch>(m, "Patternmatch")
		.def(py::init<Molecule&, Substructure&, int>())
		.def(py::init<>())
		.def("list_matches", &Patternmatch::list_matches)
		.def("unique_matches", &Patternmatch::unique_matches)
		.def("GetDistinctMatches", &Patternmatch::GetDistinctMatches)
		.def("print_M", &Patternmatch::print_M)
		.def("H_factor", &Patternmatch::H_factor)
		.def("H_factor_unique_match", &Patternmatch::H_factor_unique_match)
		.def("print_Ma", &Patternmatch::print_Ma)
		.def("number_of_matches", &Patternmatch::number_of_matches)
		.def("number_of_unique_matches", &Patternmatch::number_of_unique_matches)
		.def("get_unique_matches", &Patternmatch::get_unique_matches)
		.def("get_match", &Patternmatch::get_match)
		.def("getMatchFrequency", &Patternmatch::getMatchFrequency)
		.def("AtomsCoveredBySubstr", &Patternmatch::AtomsCoveredBySubstr)
		.def("IsMatchWithinRing", &Patternmatch::IsMatchWithinRing)
		.def("GetDistinctAllAndRingMatches", &Patternmatch::GetDistinctAllAndRingMatches);
		
	//groupadditivity.cpp classes
	py::class_<ThermoGA>(m, "ThermoGA")
		.def("calculateDeltaH", &ThermoGA::calculateDeltaH)
		.def("calculateDeltaS", &ThermoGA::calculateDeltaS)
		.def("calculateDeltaG", &ThermoGA::calculateDeltaG)
		.def("calculateCp", &ThermoGA::calculateCp);
		
	py::class_<ThermoValues>(m, "ThermoValues")
		.def("getThermo", &ThermoValues::getThermo)
		.def("setThermo", &ThermoValues::setThermo)
		.def("isTempAvailable", &ThermoValues::isTempAvailable);
		
	//automorphs.cpp classes
	py::class_<Automorphs>(m, "Automorphs")
		.def(py::init<Molecule&>())
		.def("NumberOfAutomorphs", &Automorphs::NumberOfAutomorphs)
		.def("SymmetryNumber", &Automorphs::SymmetryNumber);
	
	//reaction.cpp classes
	py::class_<Reactiontype>(m, "Reactiontype")
		.def(py::init<>())
		.def("add_reactant_pattern", &Reactiontype::add_reactant_pattern)
		.def("add_reactantconstraint", &Reactiontype::add_reactantconstraint)
		.def("add_productconstraint", &Reactiontype::add_productconstraint)
		.def("add_combined_constraint", &Reactiontype::add_combined_constraint)
		.def("add_mod_atomtype", &Reactiontype::add_mod_atomtype)
		.def("get_molecularity", &Reactiontype::get_molecularity)
		.def("AllowIntraMolecularRxnOnly", &Reactiontype::AllowIntraMolecularRxnOnly)
		.def("AllowIntraMolecularRxnAlso", &Reactiontype::AllowIntraMolecularRxnAlso)
		.def("isIntraMolecularAlso", &Reactiontype::isIntraMolecularAlso)
		.def("isIntraMolecularOnly", &Reactiontype::isIntraMolecularOnly)
		.def("AddReactantCopy", &Reactiontype::AddReactantCopy)
		.def("disconnect_bond", &Reactiontype::disconnect_bond)
		.def("connect_bond", &Reactiontype::connect_bond)
		.def("increaseBO", &Reactiontype::increaseBO)
		.def("decreaseBO", &Reactiontype::decreaseBO)
		.def("setCost", &Reactiontype::setCost)
		.def("getCost", &Reactiontype::getCost)
		.def("setRuleName", &Reactiontype::setRuleName)
		.def("getRuleName", &Reactiontype::getRuleName)
		.def("getFragmentCopyIndex", &Reactiontype::getFragmentCopyIndex)
		.def("BreaksAromaticity", &Reactiontype::BreaksAromaticity)
		.def("setSpeciesRank", &Reactiontype::setSpeciesRank)
		.def("getSpeciesRank", &Reactiontype::getSpeciesRank)
		.def("AllowSelfRxnOnly", &Reactiontype::AllowSelfRxnOnly)
		.def("setMinRateConst", &Reactiontype::setMinRateConst)
		.def("getMinRateConst", &Reactiontype::getMinRateConst)
		.def("shouldCheckRateConstant", &Reactiontype::shouldCheckRateConstant);
	
	py::class_<Reaction>(m, "Reaction")
		.def(py::init<std::vector<Molecule>&, Reactiontype&, std::vector<Patternmatch>&>())
		.def("get_generated_rxns", &Reaction::get_generated_rxns)
		.def("number_rxns_generated", &Reaction::number_rxns_generated);
		
	//lumping.cpp classes
	py::class_<LumpInfo>(m, "LumpInfo")
		.def(py::init<int, int, bool, int, int, std::string*>())
		.def("setMoleculeString", &LumpInfo::setMoleculeString)
		.def("getMoleculeString", &LumpInfo::getMoleculeString)
		.def("getMolStringPtr", &LumpInfo::getMolStringPtr)
		.def("getSize", &LumpInfo::getSize)
		.def("getHydrogens", &LumpInfo::getHydrogens)
		.def("getBranchValue", &LumpInfo::getBranchValue);
		
	py::class_<LumpingStrategy>(m, "LumpingStrategy")
		.def(py::init<bool>())
		.def(py::init<>())
		.def("shoudLump", &LumpingStrategy::shoudLump);
		
	
	//rng.cpp classes
	py::class_<rxn_net_gen>(m, "rxn_net_gen")
		.def(py::init<>())
		.def(py::init<rxn_net_gen*, LumpingStrategy&>())
		.def("GenerateNetwork", &rxn_net_gen::GenerateNetwork)
		.def("AddInitialReactants", &rxn_net_gen::AddInitialReactants)
		.def("AddReactionRules", &rxn_net_gen::AddReactionRules)
		.def("AddGlobalConstraints", &rxn_net_gen::AddGlobalConstraints)
		.def("AddLumpingStrategy", &rxn_net_gen::AddLumpingStrategy)
		.def("AddCompositeAtoms", &rxn_net_gen::AddCompositeAtoms)
		.def("AddCompositeSites", &rxn_net_gen::AddCompositeSites)
		.def("SetCalcThermo", &rxn_net_gen::SetCalcThermo)
		.def("print_rxnlist", &rxn_net_gen::print_rxnlist);
	
	py::class_<SingleFunc>(m, "SingleFunc")
		.def(py::init<>())
		.def("singlef", &SingleFunc::singlef);
		
    py::class_<DoubleFunc>(m, "DoubleFunc")
		.def(py::init<>())
		.def("doublef", &DoubleFunc::doublef);
	
//	py::class_<SingleMolecule>(m, "SingleMolecule")
//        .def(py::init<>())
//        .def("singleconstr", &SingleMolecule::singleconstr);
//		
//	py::class_<DoubleMolecule>(m, "DoubleMolecule")
//        .def(py::init<>())
//        .def("doubleconstr", &DoubleMolecule::doubleconstr);
	
} 