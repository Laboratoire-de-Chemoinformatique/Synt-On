import xml.etree.ElementTree as ET
from rdkit import Chem, DataStructs
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem import AddHs, AllChem
from rdkit.Chem.rdmolops import *
from rdkit.Chem.rdMolDescriptors import CalcNumRings
from rdkit.Chem.Descriptors import *
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem.Crippen import MolLogP
import datetime, os, time, random, re, resource, sys
from multiprocessing import Process, Queue
from collections import Counter
srcPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(1, srcPath)
from UsefulFunctions import *
from concurrent.futures import ProcessPoolExecutor
from functools import partial

class synthon:
    def __init__(self, smiles, cutLevel=1, directParent=None, directChildren=None, syntheticPathway=None, SynthLibProvided=False):
        self.smiles = smiles
        pat = re.compile("\[\w*:\w*\]")
        self.functionalityCount = len([self.smiles[m.start():m.start() + 2] + self.smiles[m.end() - 4:m.end()] for m in re.finditer(pat, self.smiles)])
        self.marks =  sorted(
            [self.smiles[m.start():m.start() + 2] + self.smiles[m.end() - 4:m.end()] for m in re.finditer(pat, self.smiles)])
        self.cutLevel = cutLevel
        self.directParents = []
        if directParent:
            self.directParents.extend(directParent)
        self.directChildren = []
        if directChildren:
            self.directChildren.extend(directChildren)
        self.syntheticPathway = []
        if syntheticPathway:
            self.syntheticPathway.extend(syntheticPathway)
        self.correspondingBB = None
        self.AtomNumbers = 0
        self.bbAnalogues = {}
        self.bivalentN = False
        self.SynthLibProvided = SynthLibProvided
        self.rIdsToGetIt = []

    def searchForSynthonAnalogues(self, synthLib: dict, simTh=-1):
        posAnaloguesScreeningAtomsAllowed = ["C", "F", "N", "O"]
        refMol = Chem.MolFromSmiles(self.smiles)
        refMol.UpdatePropertyCache()
        Chem.GetSymmSSSR(refMol)
        refMol.GetRingInfo().AtomRings()
        queryFP = AllChem.GetMorganFingerprintAsBitVect(refMol, radius=2, nBits=2048)
        refMarksVallences = "+".join(sorted([atom.GetSymbol() + ":" + str(atom.GetTotalDegree()) for atom in refMol.GetAtoms() if atom.GetAtomMapNum() != 0]))
        for synth in synthLib:
            if self.marks == synthLib[synth]["marks"] and synth != self.smiles and refMarksVallences == synthLib[synth]["marksVallences"]:
                if simTh != -1 and DataStructs.TanimotoSimilarity(queryFP,synthLib[synth]["fp_b"]) >= simTh :
                    self.bbAnalogues[synth] = synthLib[synth]["BBs"]
                else:
                    if self.AtomNumbers == 0:
                        self.AtomNumbers = refMol.GetNumAtoms()
                    if synthLib[synth]["n_atoms"] <= self.AtomNumbers + 1 and synthLib[synth]["n_atoms"] >= self.AtomNumbers - 1:
                        refList = [i for i in self.smiles if i.isalpha() and i!="H"]
                        qList = [i for i in synth if i.isalpha() and i!="H"]
                        qList_refList = list((Counter(qList) - Counter(refList)).elements())
                        refList_qList = list((Counter(refList) - Counter(qList)).elements())
                        if synthLib[synth]["n_atoms"] == self.AtomNumbers and sorted(refList_qList + qList_refList) == ['c', 'n']:
                            refMolRings = CalcNumRings(refMol)
                            if synthLib[synth]["n_rings"] == refMolRings:
                                self.bbAnalogues[synth] = synthLib[synth]["BBs"]
                        elif not qList_refList and not refList_qList:
                            refMolRings = CalcNumRings(refMol)
                            if synthLib[synth]["n_rings"] == refMolRings:
                                self.bbAnalogues[synth] = synthLib[synth]["BBs"]
                        elif qList_refList and len(qList_refList) == 1 and qList_refList[0] in posAnaloguesScreeningAtomsAllowed:
                            refMolRings = CalcNumRings(refMol)
                            if synthLib[synth]["n_rings"] == refMolRings:
                                analogMol = Chem.MolFromSmiles(synth)
                                analogMol.UpdatePropertyCache()
                                Chem.GetSymmSSSR(analogMol)
                                analogMol.GetRingInfo().AtomRings()
                                if analogMol.HasSubstructMatch(refMol):
                                    self.bbAnalogues[synth] = synthLib[synth]["BBs"]
                        elif refList_qList and len(refList_qList) == 0 and refList_qList[0] in posAnaloguesScreeningAtomsAllowed:
                            refMolRings = CalcNumRings(refMol)
                            if synthLib[synth]["n_rings"] == refMolRings:
                                analogMol = Chem.MolFromSmiles(synth)
                                analogMol.UpdatePropertyCache()
                                Chem.GetSymmSSSR(analogMol)
                                analogMol.GetRingInfo().AtomRings()
                                if refMol.HasSubstructMatch(analogMol):
                                    self.bbAnalogues[synth] = synthLib[synth]["BBs"]

    def printSynthonInfo(self):
        print("\n__________________________________________________________")
        print("                  Synthon Information" )
        print("__________________________________________________________\n")
        print("Synthon: " + self.smiles)
        if self.SynthLibProvided:
            if self.correspondingBB:
                print("Available. Corresponding BBs: " + self.correspondingBB)
            elif self.bbAnalogues:
                print("Synthon was not found in provided library of building blocks. " + str(len(self.bbAnalogues)) + " analog(s) has/have been found")
                print("BB analogues:\n" + "\n".join([x + " " + self.bbAnalogues[x] for x in self.bbAnalogues]))
            else:
                print("Synthon was not found in provided library of building blocks")
        if self.directParents:
            print("Parent synthons: " + ".".join(list(set([x.smiles for x in self.directParents]))))
        else:
            print("Parent synthons: - ")
        if self.directChildren:
            print("Children synthons: " + ".".join(list(set([x.smiles for x in self.directChildren]))))
        else:
            print("Children synthons: - ")

    def Ro2Filtration(self):
        mol = Chem.MolFromSmiles(self.smiles)
        mol = AddHs(mol)
        MolW = ExactMolWt(mol)
        LogP = MolLogP(mol)
        HDC = CalcNumHBD(mol)
        HAC = CalcNumHBA(mol)
        if "[NH:20]" in self.smiles or "[OH:20]" in self.smiles or "[SH:20]" in self.smiles:
            marksForCount = ["[NH:20]", "[OH:20]", "[SH:20]"]
            count = 0
            for m in marksForCount:
                count += self.smiles.count(m)
            HDC -= count
        if MolW > 200 or LogP > 2 or HDC > 2 or HAC > 4:
            return False
        else:
            return True

class syntheticPathway:
    def __init__(self, name , synthPathwayReactions,  reagentsNumber, cutLevel, directParentsSynthPathways=None, SynthLibProvided=False):
        self.name = name
        self.participatingSynthon  = []
        self.directChildrenSynthPathways = []
        self.directParentsSynthPathways = []
        if directParentsSynthPathways:
            self.directParentsSynthPathways.extend(directParentsSynthPathways)
        self.synthPathwayReactions = []
        if synthPathwayReactions:
            self.synthPathwayReactions.extend(synthPathwayReactions)
        self.reagentsNumber = reagentsNumber
        self.cutLevel = cutLevel
        self.availabilityRate = 0
        self.SynthLibProvided = SynthLibProvided

    def checkAvailability(self, SynthLib: dict, simTh=-1, FindAnaloguesOfMissingSynthons=True):
        availableAtomCounts = 0
        allAtoms = 0
        for synth in self.participatingSynthon:
            if synth.AtomNumbers==0:
                synth.AtomNumbers = Chem.MolFromSmiles(synth.smiles).GetNumAtoms()
            allAtoms += synth.AtomNumbers
            if synth.smiles in SynthLib:
                availableAtomCounts += synth.AtomNumbers
                synth.correspondingBB = SynthLib[synth.smiles]["BBs"]
            elif not synth.bbAnalogues and FindAnaloguesOfMissingSynthons:
                synth.searchForSynthonAnalogues(SynthLib, simTh)
        self.availabilityRate = round(availableAtomCounts / allAtoms, 2)
        if self.directChildrenSynthPathways:
            for comb in self.directChildrenSynthPathways:
                if comb.availabilityRate==0:
                    comb.checkAvailability(SynthLib, simTh, FindAnaloguesOfMissingSynthons)

    def printShortReagentSetInfo(self):
        if self.name=="zeroCombin":
            return
        if self.SynthLibProvided:
            print(self.name + " " + ".".join([x.smiles for x in self.participatingSynthon]) + " " + "Availability rate (% of atoms of fragmented molecule coming from available synthons): " + str(self.availabilityRate))
        else:
            print(self.name + " " + ".".join([x.smiles for x in self.participatingSynthon]))

    def printDetailedReagentsSetInfo(self):
        if self.name!="zeroCombin":
            print("\n**********************************************************")
            print("Reagent set Information " + self.name )
            print("**********************************************************")
            print("Reactions: " + "->".join(self.synthPathwayReactions))
            print("Required Synthons: " + ".".join([x.smiles for x in self.participatingSynthon]))
            print("Number of reagents: " + str(self.reagentsNumber))
            print("Number of stages: " + str(self.cutLevel))
            if self.SynthLibProvided:
                print("Availability rate (% of atoms of fragmented molecule coming from available synthons): " + str(self.availabilityRate))
            if self.directChildrenSynthPathways:
                print("Children reagent sets: ")
                for comb in self.directChildrenSynthPathways:
                    print("Reactions: " + comb.name + " " +  "->".join(comb.synthPathwayReactions) + " ||| " + "Participating synthons: " + ".".join([x.smiles for x in comb.participatingSynthon]))
            if self.directParentsSynthPathways:
                print("Parent reagent sets: ")
                for comb in self.directParentsSynthPathways:
                    print("Reactions: " + comb.name + " " +  "->".join(comb.synthPathwayReactions) + " ||| " + "Participating synthons: " + ".".join([x.smiles for x in comb.participatingSynthon]))
            for synth in self.participatingSynthon:
                synth.printSynthonInfo()

    def getSynthonsForAnaloguesGeneration(self, SynthLib, simTh, strictAvailabilityMode = True):
        synthonsDict = {}
        SynthonsForAnaloguesSynthesis = set()
        if not strictAvailabilityMode:
            for ind, synth in enumerate(self.participatingSynthon):

                if "Reagent_" + str(ind + 1) not in synthonsDict:
                    synthonsDict["Reagent_" + str(ind + 1)] = { "synthons": [], "bivalentN": synth.bivalentN }
                synthonsDict["Reagent_" + str(ind + 1)]["synthons"].append(synth.smiles)
                if not synth.correspondingBB:
                    SynthonsForAnaloguesSynthesis.add(synth.smiles + " missingBB originalBB\n")
                else:
                    SynthonsForAnaloguesSynthesis.add(synth.smiles + " " + synth.correspondingBB + " originalBB\n")
                if not synth.bbAnalogues and synth.correspondingBB:
                    synth.searchForSynthonAnalogues(SynthLib, simTh)
                if synth.bbAnalogues:
                    for analog in synth.bbAnalogues:
                        SynthonsForAnaloguesSynthesis.add(analog + " " + synth.bbAnalogues[analog] + " " + synth.smiles + " analog\n")
                        synthonsDict["Reagent_" + str(ind + 1)]["synthons"].append(analog)
            return synthonsDict, SynthonsForAnaloguesSynthesis
        else:
            for ind, synth in enumerate(self.participatingSynthon):
                if "Reagent_" + str(ind + 1) not in synthonsDict:
                    synthonsDict["Reagent_" + str(ind + 1)] = { "synthons": [], "bivalentN": synth.bivalentN }
                if not synth.bbAnalogues and synth.correspondingBB:
                    synth.searchForSynthonAnalogues(SynthLib, simTh)
                if synth.correspondingBB:
                    SynthonsForAnaloguesSynthesis.add(synth.smiles + " " + synth.correspondingBB + " originalBB\n")
                    synthonsDict["Reagent_" + str(ind + 1)]["synthons"].append(synth.smiles)
                if synth.bbAnalogues:
                    for analog in synth.bbAnalogues:
                        SynthonsForAnaloguesSynthesis.add(analog + " " + synth.bbAnalogues[analog] + " " + synth.smiles + " analog\n")
                        synthonsDict["Reagent_" + str(ind + 1)]["synthons"].append(analog)
                if len(synthonsDict["Reagent_" + str(ind + 1)]["synthons"]) == 0:
                    return None, None
            return synthonsDict, SynthonsForAnaloguesSynthesis

class enumeration:
    def __init__(self, outDir, Synthons = None, reactionSMARTS = None, maxNumberOfReactedSynthons=6, MWupperTh=None, MWlowerTh=None,
                  desiredNumberOfNewMols = 1000, nCores=1, analoguesEnumeration=False):
        if analoguesEnumeration and Synthons!=None:
            self.__Synthons = None
            self.__monoFuncBB = None
            self.__poliFuncBB = None
            self.__biFuncBB = None
            self.__maxNumberOfReactedSynthons = len(Synthons)
        else:
            self.__Synthons = []
            self.__monoFuncBB = []
            self.__poliFuncBB = []
            self.__biFuncBB = []
            self.__maxNumberOfReactedSynthons = maxNumberOfReactedSynthons
        self.__MWfiltration=False
        if MWupperTh:
            self.__MWupperTh = MWupperTh
            self.__MWfiltration = True
        if  MWlowerTh:
            self.__MWlowerTh = MWlowerTh
            self.__MWfiltration = True
        self.__reactions = []
        self.__perpSynthonsAndReactions(Synthons, reactionSMARTS, analoguesEnumeration)
        self.__desiredNumberOfNewMols = desiredNumberOfNewMols
        self.__nCores = nCores
        self.__outDir = outDir
        self.__genNonUniqMols = 0
        self.results = set()
        self.allReconstructedMols = []
        self.__marksCombinations = {'C:10': ['N:20', 'O:20', 'C:20', 'c:20', 'n:20', 'S:20'],
                                    'c:10': ['N:20', 'O:20', 'C:20', 'c:20', 'n:20', 'S:20'],
                                    'c:20': ['N:11', 'C:10', 'c:10'], 'C:20': ['C:10', 'c:10'],
                                    'c:21': ['N:20', 'O:20', 'n:20'], 'C:21': ['N:20', 'n:20'],
                                    'N:20': ['C:10', 'c:10', 'C:21', 'c:21', 'S:10'], 'N:11': ['c:20'],
                                    'n:20': ['C:10', 'c:10', 'C:21', 'c:21'], 'O:20': ['C:10', 'c:10', 'c:21'],
                                    'S:20': ['C:10', 'c:10'], 'S:10': ['N:20'], 'C:30': ['C:40', 'N:40'],
                                    'C:40': ['C:30'], 'C:50': ['C:50'], 'C:70': ['C:60', 'c:60'],
                                    'c:60':['C:70'], 'C:60': ['C:70'], 'N:40': ['C:30'] }


    def getReconstructedMols(self, allowedToRunSubprocesses=False, randomSeed=None, seed = (0,0), mainRun = True):
        pat = re.compile("\[\w*:\w*\]")
        if randomSeed==None:
            seed, randomSeed = self.__getRandomSeed(seed)
        if randomSeed:
            allowedMarks = []
            for key in [Chem.MolToSmiles(randomSeed, canonical=True)[m.start() + 1] + ":" + Chem.MolToSmiles(randomSeed, canonical=True)[m.end() - 3:m.end() - 1]
                        for m in re.finditer(pat, Chem.MolToSmiles(randomSeed, canonical=True))]:
                allowedMarks.extend(self.__marksCombinations[key])
            allowedMarks = set(allowedMarks)
            Pool = []
            for partner in self.__poliFuncBB + self.__biFuncBB + self.__monoFuncBB:
                partnerMarks = set([Chem.MolToSmiles(partner, canonical=True)[m.start() + 1] + ":" + Chem.MolToSmiles(partner, canonical=True)[m.end() - 3:m.end() - 1]
                        for m in re.finditer(pat, Chem.MolToSmiles(partner, canonical=True))])
                if allowedMarks.intersection(partnerMarks):
                    if allowedToRunSubprocesses:
                        Pool, nAlive = self.__countAndMergeActiveThreads(Pool)
                        while nAlive >= self.__nCores:
                            time.sleep(30)
                            Pool, nAlive = self.__countAndMergeActiveThreads(Pool)
                        queue = Queue()
                        process = Process(target=self.__molReconsrtuction, args=(randomSeed, partner, 1, None, queue,))
                        process.start()
                        Pool.append([process, queue, False]) # [the process, queue, was it already joined]
                        if self.__genNonUniqMols >= self.__desiredNumberOfNewMols:
                            break
                    else:
                        self.results.update(self.__molReconsrtuction(randomSeed, partner, 1, queue=None))
                        print("Number of so far reconstructed unique molecules = " + str(len(self.results)))
                        if len(self.results) >= self.__desiredNumberOfNewMols:
                            break
            if allowedToRunSubprocesses:
                Pool, nAlive = self.__countAndMergeActiveThreads(Pool)
                while nAlive > 0:
                    time.sleep(5)
                    Pool, nAlive = self.__countAndMergeActiveThreads(Pool)

        whileCount = 0
        while mainRun:
            if not allowedToRunSubprocesses and len(self.results) >= self.__desiredNumberOfNewMols :
                break
            elif allowedToRunSubprocesses and self.__genNonUniqMols >= self.__desiredNumberOfNewMols :
                break
            whileCount += 1
            seed = list(seed)
            seed[0] +=1
            seed = tuple(seed)
            seed, randomSeed = self.__getRandomSeed(seed)
            if randomSeed is None:
                break
            self.results.update(self.getReconstructedMols(allowedToRunSubprocesses, randomSeed, seed, mainRun = False))
        d_names, f_names, main_dir = listDir(self.__outDir)
        with open(os.path.join(self.__outDir, "FinalOut_allEnumeratedCompounds_DuplicatesCanBePresent.smi"), "ab") as out:
            for file in f_names:
                if "temp_" in file:
                    with open(os.path.join(self.__outDir, file), "rb") as f:
                        out.write(f.read())
                        out.flush()
                    os.remove(os.path.join(self.__outDir, file))
        return self.results

    def AnaloguesGeneration(self):
        if type(self.__Synthons) == dict:
            for rid in range(len(self.__reactions)):
                reaction = self.__reactions[rid]
                for ind1 in range(1, len(self.__Synthons) + 1):
                    for ind2 in range(ind1 + 1, len(self.__Synthons) + 1):
                        SynthonsetPartner = 'Reagent_' + str(ind2)
                        Synthonset = 'Reagent_' + str(ind1)
                        if self.__checkBB_reactionCombination(self.__Synthons[Synthonset][0], self.__Synthons[SynthonsetPartner][0], reaction):
                            for firstBB in self.__Synthons[Synthonset]:
                                for secondBB in self.__Synthons[SynthonsetPartner]:
                                    SynthonsetsUsed = set()
                                    SynthonsetsUsed.add(Synthonset)
                                    SynthonsetsUsed.add(SynthonsetPartner)
                                    self.results.update(self.__molAnaloguesLibEnumeration(reagent=firstBB, partner=secondBB,
                                      numberOfBBalreadyReacted=1, SynthonsetsUsed=SynthonsetsUsed, reactionToUse = rid))
                                    if len(self.results) >= self.__desiredNumberOfNewMols:
                                        return self.results
            return self.results
        else:
            print("Separate reconstructur should be evocken for Molecule analogues generation "
                  "from available BBs (set analoguesEnumeration=True) and for "
                  "simple molecules enumeration from the list of BBs (set analoguesEnumeration=False).")
            exit()

    def __molAnaloguesLibEnumeration(self, reagent, partner, numberOfBBalreadyReacted,
                                                    SynthonsetsUsed, reactionToUse, usedReactions = None, firstLaunch=True):
        if not usedReactions:
            usedReactions = []
        allProducts = set()
        products = self.__reactions[reactionToUse].RunReactants((reagent, partner))
        if not products:
            products = self.__reactions[reactionToUse].RunReactants((partner, reagent))
        if products:
            for prodSet in products:
                for prod in prodSet:
                    prodSMILES = Chem.MolToSmiles(prod, canonical=True)
                    functionality = prodSMILES.count(":")
                    prod.UpdatePropertyCache()
                    Chem.GetSymmSSSR(prod)
                    prod.GetRingInfo().NumRings()
                    newUsedReactions = usedReactions.copy()
                    newUsedReactions.append(reactionToUse)
                    if functionality:
                        for rid in range(len(self.__reactions)):
                            if rid not in newUsedReactions:
                                reaction = self.__reactions[rid]
                                for SynthonsetPartner in self.__Synthons:
                                    if SynthonsetPartner not in SynthonsetsUsed and self.__checkBB_reactionCombination(prod,
                                                                              self.__Synthons[SynthonsetPartner][0], reaction):
                                        for secondBB in self.__Synthons[SynthonsetPartner]:
                                            newSynthonsetsUsed = set()
                                            for i in SynthonsetsUsed:
                                                newSynthonsetsUsed.add(i)
                                            newSynthonsetsUsed.add(SynthonsetPartner)
                                            subResults = self.__molAnaloguesLibEnumeration(reagent=prod,
                                                  partner=secondBB, numberOfBBalreadyReacted=numberOfBBalreadyReacted +1,
                                                  SynthonsetsUsed=newSynthonsetsUsed, reactionToUse=rid, usedReactions = newUsedReactions,
                                                                                           firstLaunch=False)
                                            if subResults:
                                                allProducts.update(subResults)
                                            if len(allProducts)>=self.__desiredNumberOfNewMols and firstLaunch:
                                                return list(allProducts)
                    elif not functionality and numberOfBBalreadyReacted + 1 >= len(self.__Synthons):
                        allProducts.add(prodSMILES)
                        if len(allProducts) >= self.__desiredNumberOfNewMols and firstLaunch:
                            return list(allProducts)
        return list(allProducts)

    def __checkBB_reactionCombination(self, reagent, partner, reaction):
        products = reaction.RunReactants((reagent, partner))
        if products:
            return True
        else:
            products = reaction.RunReactants((partner, reagent))
            if products:
                return True
            else:
                return False

    def __molReconsrtuction(self, reagent, partner, numberOfBBalreadyReacted,  alreadyUsedReactionsInds=None, queue=None):
        if not alreadyUsedReactionsInds:
            alreadyUsedReactionsInds = set()
        pat = re.compile("\[\w*:\w*\]")
        allProducts = []
        for id,reaction in enumerate(self.__reactions):
            if id in alreadyUsedReactionsInds:
                continue
            products = reaction.RunReactants((reagent, partner))
            if not products:
                products = reaction.RunReactants((partner, reagent))
            if products:
                for prodSet in products:
                    for prod in prodSet:
                        prodSMILES = Chem.MolToSmiles(prod, canonical=True)
                        functionality = prodSMILES.count(":")
                        prod.UpdatePropertyCache()
                        Chem.GetSymmSSSR(prod)
                        prod.GetRingInfo().NumRings()
                        if functionality:
                            allowedMarks = []
                            for key in [prodSMILES[m.start() + 1] + ":" + prodSMILES[m.end() - 3:m.end() - 1] for m in
                                        re.finditer(pat, prodSMILES)]:
                                allowedMarks.extend(self.__marksCombinations[key])
                            allowedMarks = set(allowedMarks)
                            if numberOfBBalreadyReacted + 1 <= self.__maxNumberOfReactedSynthons - functionality:
                                if numberOfBBalreadyReacted + 1 == self.__maxNumberOfReactedSynthons - functionality:
                                    for newPartner in self.__monoFuncBB:
                                        partnerMarks = set([
                                            Chem.MolToSmiles(newPartner, canonical=True)[m.start() + 1] + ":" + Chem.MolToSmiles(newPartner, canonical=True)[
                                            m.end() - 3:m.end() - 1] for m in re.finditer(pat, Chem.MolToSmiles(newPartner, canonical=True))])
                                        if allowedMarks.intersection(partnerMarks):
                                            alreadyUsedReactionsInds.add(id)
                                            subResults = self.__molReconsrtuction(prod, newPartner, numberOfBBalreadyReacted+1,
                                                                                  alreadyUsedReactionsInds=alreadyUsedReactionsInds)
                                            if subResults:
                                                allProducts.extend(subResults)
                                elif functionality>3:
                                    for newPartner in self.__monoFuncBB + self.__biFuncBB:
                                        partnerMarks = set([
                                            Chem.MolToSmiles(newPartner, canonical=True)[m.start() + 1] + ":" + Chem.MolToSmiles(
                                                newPartner, canonical=True)[m.end() - 3:m.end() - 1]
                                            for m in re.finditer(pat, Chem.MolToSmiles(newPartner, canonical=True))])
                                        if allowedMarks.intersection(partnerMarks):
                                            alreadyUsedReactionsInds.add(id)
                                            subResults = self.__molReconsrtuction(prod, newPartner, numberOfBBalreadyReacted+1,
                                                                                  alreadyUsedReactionsInds=alreadyUsedReactionsInds)
                                            if subResults:
                                                allProducts.extend(subResults)
                                                allProducts = list(set(allProducts))
                                else:
                                    for newPartner in self.__monoFuncBB + self.__biFuncBB + self.__poliFuncBB:
                                        partnerMarks = set([
                                            Chem.MolToSmiles(newPartner, canonical=True)[m.start() + 1] + ":" + Chem.MolToSmiles(
                                                newPartner, canonical=True)[
                                                                                                m.end() - 3:m.end() - 1]
                                            for m in re.finditer(pat, Chem.MolToSmiles(newPartner, canonical=True))])
                                        if allowedMarks.intersection(partnerMarks):
                                            alreadyUsedReactionsInds.add(id)
                                            subResults = self.__molReconsrtuction(prod, newPartner, numberOfBBalreadyReacted+1,
                                                                                  alreadyUsedReactionsInds=alreadyUsedReactionsInds)
                                            if subResults:
                                                allProducts.extend(subResults)
                                                allProducts = list(set(allProducts))
                        else:
                            if prodSMILES not in allProducts:
                                if self.__MWfiltration:
                                    molW = ExactMolWt(prod)
                                    if self.__MWlowerTh and self.__MWupperTh:
                                        if molW <= self.__MWupperTh and molW >= self.__MWlowerTh:
                                            allProducts.append(prodSMILES)
                                    elif self.__MWlowerTh:
                                        if molW >= self.__MWlowerTh:
                                            allProducts.append(prodSMILES)
                                    elif self.__MWupperTh:
                                        if molW <= self.__MWupperTh:
                                            allProducts.append(prodSMILES)
                                else:
                                    allProducts.append(prodSMILES)

        if queue is None:
            return list(set(allProducts))
        else:
            queue.put(list(set(allProducts)))

    def __perpSynthonsAndReactions(self, Synthons, reactionSMARTS, analoguesEnumeration):
        if Synthons == None:
            raise ValueError("No building blocks are provided for enumeration")
        else:
            if analoguesEnumeration:
                self.__Synthons = dict()
                for subset in Synthons:
                    bbsList = [Chem.MolFromSmiles(x) for x in Synthons[subset]["synthons"]]
                    self.__Synthons[subset] = self.__PrepMolForReconstruction(bbsList, Synthons[subset]["bivalentN"])
            else:
                monoFuncBB = []
                biFuncBB = []
                poliFuncBB = []
                for bb in Synthons:
                    if bb.count(":") == 1:
                        monoFuncBB.append(bb)
                    elif bb.count(":") == 2:
                        biFuncBB.append(bb)
                    else:
                        poliFuncBB.append(bb)
                if monoFuncBB == None:
                    raise ValueError("No monofunctional building blocks are provided for enumeration")
                else:
                    self.__monoFuncBB = self.__PrepMolForReconstruction([Chem.MolFromSmiles(x) for x in monoFuncBB])
                    random.shuffle(self.__monoFuncBB)
                if biFuncBB:
                    self.__biFuncBB = self.__PrepMolForReconstruction([Chem.MolFromSmiles(x) for x in biFuncBB])
                    random.shuffle(self.__biFuncBB)
                if poliFuncBB:
                    self.__poliFuncBB = self.__PrepMolForReconstruction([Chem.MolFromSmiles(x) for x in poliFuncBB])
                    random.shuffle(self.__poliFuncBB)
                self.__Synthons = (self.__poliFuncBB, self.__biFuncBB, self.__monoFuncBB)
            if reactionSMARTS == None:
                raise ValueError("No reactions are provided for enumeration")
            else:
                for react in reactionSMARTS:
                    if react:
                        reaction = Reactions.ReactionFromSmarts(react)
                        self.__reactions.append(reaction)

    def __PrepMolForReconstruction(sefl, molList, bivalentN=False): # list of Mol objects
        newList = []
        labels = [10, 20, 30, 40, 50, 60, 70, 21, 11]
        atomsForMarking = [23, 74, 72, 104, 105, 106, 107, 108, 109]
        atomsForMarkingForDoubleBonds = [72, 104, 105]
        for mol in molList:
            mol = AddHs(mol)
            for atom in mol.GetAtoms():
                if atom.GetAtomMapNum() != 0:
                    repl = atomsForMarking[labels.index(atom.GetAtomMapNum())]
                    replCount = 0
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 1:
                            mol.GetAtomWithIdx(neighbor.GetIdx()).SetAtomicNum(repl)
                            replCount += 1
                            if repl not in atomsForMarkingForDoubleBonds and replCount == 1:
                                break
                            elif replCount == 2:
                                break
            mol = RemoveHs(mol)
            newList.append(mol)
        if bivalentN:
            for mol in molList:
                mol = AddHs(mol)
                for atom in mol.GetAtoms():
                    if atom.GetAtomMapNum() != 0:
                        repl = atomsForMarking[labels.index(atom.GetAtomMapNum())]
                        replCount = 0
                        if atom.GetSymbol() == "N" and atom.GetAtomMapNum() == 20:
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetAtomicNum() == 1:
                                    mol.GetAtomWithIdx(neighbor.GetIdx()).SetAtomicNum(repl)
                                    replCount += 1
                                    if replCount == 2:
                                        break
                        else:
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetAtomicNum() == 1:
                                    mol.GetAtomWithIdx(neighbor.GetIdx()).SetAtomicNum(repl)
                                    replCount += 1
                                    if repl not in atomsForMarkingForDoubleBonds and replCount == 1:
                                        break
                                    elif replCount == 2:
                                        break
                mol = RemoveHs(mol)
                newList.append(mol)
        return newList # list of Mol objects

    def __getRandomSeed(self, seed):
        randomSeed = None
        while not randomSeed:
            if seed[1] > 2:
                break
            if seed[0] < len(self.__Synthons[seed[1]]):
                randomSeed = self.__Synthons[seed[1]][seed[0]]
            else:
                seed = list(seed)
                seed[1] += 1
                seed[0] = 0
                seed = tuple(seed)
                continue
        return seed, randomSeed

    def __countAndMergeActiveThreads(self, Pool):
        nAlive = 0
        for p in Pool:
            if p[1].empty():
                nAlive += 1
            else:
                reconstructedMols = p[1].get()
                # self.results.extend(p[1].get())
                # self.results = list(set(self.results))
                self.__genNonUniqMols += len(reconstructedMols)
                print("Number of so far reconstructed molecules (may contain duplicates) = " + str(self.__genNonUniqMols))
                d_names, f_names, main_dir = listDir(self.__outDir)
                n = 0
                for file in f_names:
                    if "temp_" in file:
                        if int(file.split("_")[1]) > n:
                            n = int(file.split("_")[1])
                with open(os.path.join(self.__outDir, "temp_" + str(n + 1)), "w") as out:
                    out.writelines("%s\n" % mols for mols in reconstructedMols)
                p[0].join()
                p[-1] = True
                p[1].close()
                p[1].join_thread()
        for i in range(len(Pool) - 1, -1, -1):
            if Pool[i][-1]:
                Pool.pop(i)
        return Pool, nAlive

class fragmentation:

    def __init__(self, fragmentationMode="use_all", reactionsToWorkWith = "R1-R13", maxNumberOfReactionCentersPerFragment = 3,
                 MaxNumberOfStages = 5, FragmentsToIgnore = None,
                 setupFile = os.path.join(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0], "config" , "Setup.xml"),
                 macroCycleSetupFile = os.path.join(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0], "config" , "SetupForMacrocycles.xml"),
                 FindAnaloguesOfMissingSynthons = False, parsedSynthLib = False, SynthLibrary=None, Ro2SynthonsFiltration = False):
        max_rec = 0x100000
        # May segfault without this line. 0x100 is a guess at the size of each stack frame.
        resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
        sys.setrecursionlimit(max_rec)
        sys.getrecursionlimit()
        self.__fragmentationMode = fragmentationMode
        if FragmentsToIgnore:
            self.__SmilesToIgnore = FragmentsToIgnore
        else:
            self.__SmilesToIgnore = []
        self.__maxNumberOfReactionCentersPerFragment = maxNumberOfReactionCentersPerFragment
        self.__MaxNumberOfStages = MaxNumberOfStages
        self.__setupFile = setupFile
        self.__reactionSetup, self.__reactionsToWorkWith = self.__getSetup(reactionsToWorkWith)
        self.__macroCycleSetupFile = macroCycleSetupFile
        self.__reactionMacroCycleSetup, self.__macroCyclicReactionsToWorkWith = self.__getMacroCycleSetup(reactionsToWorkWith.replace("R", "MR"))
        forbiddenMarks = [{'[N:11]', '[c:10]'}, {'[N:11]', '[c:20]'}, {'[N:11]', '[O:20]'},{'[N:11]', '[C:30]'},
                          {'[N:11]', '[C:10]'} , {'[N:11]', '[N:11]'}, {'[N:11]', '[c:70]'},
                          {'[N:11]', '[C:40]'}, {'[N:11]', '[C:50]'}, {'[N:11]', '[S:10]'},
                                {'[c:11]', '[C:20]'}, {'[c:70]', '[c:20]'},
            {'[c:21]', '[C:60]'},{'[c:21]', '[N:11]'},  {'[c:20]', '[c:21]'}, {'[c:70]', '[C:21]'}, {'[c:20]', '[C:21]'}]
        self.forbiddenMarks = set([frozenset(x) for x in forbiddenMarks])
        if SynthLibrary:
            self.FindAnaloguesOfMissingSynthons = FindAnaloguesOfMissingSynthons
            if parsedSynthLib:
                self.SynthLib = SynthLibrary
            else:
                print("Processing BB library. It may take a few minutes, depending on the library size")
                self.SynthLib = self.__readSynthLib(SynthLibrary, Ro2SynthonsFiltration)
        else:
            self.SynthLib = None

    def cutWithHierarchyStorred(self, mol):
        max_rec = 0x100000
        # May segfault without this line. 0x100 is a guess at the size of each stack frame.
        resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
        sys.setrecursionlimit(max_rec)
        sys.getrecursionlimit()
        if self.__fragmentationMode == "one_by_one":
            allSynthons, allSyntheticPathways, ind = self.__firstMolCut(mol)
        else:
            allSynthons, allSyntheticPathways = self.__firstMolCut(mol)
        if allSyntheticPathways:
            startingCombintaions = allSyntheticPathways.copy()
            for comb in startingCombintaions:
                if allSyntheticPathways[comb].cutLevel == 1:
                    for synth in allSyntheticPathways[comb].participatingSynthon:
                        if synth != allSynthons["InitMol"]:
                            self.__cutOneSynthonHierarchically(allSynthons[synth.smiles], allSyntheticPathways[comb],
                                                                   allSynthons, allSyntheticPathways, cutLevel=2)
                            if self.__fragmentationMode == "one_by_one" and allSyntheticPathways!=startingCombintaions:
                                break

        return allSyntheticPathways,allSynthons

    def getReactionForReconstruction(self, reactionList=None):
        reactionForReconstruction = []
        if reactionList:
            for key in reactionList:
                if "R14" not in key:
                    reactionForReconstruction.append(self.__reactionSetup[key]["ReconstructionSMARTS"])
        else:
            for key in self.__reactionSetup:
                if "R14" not in key:
                    reactionForReconstruction.append(self.__reactionSetup[key]["ReconstructionSMARTS"])
        return reactionForReconstruction

    def __readSynthLib(self, SynthLib, Ro2Filtration):
        fragBegTime = datetime.datetime.now()
        availableSynthons = {}
        pat = re.compile("\[\w*:\w*\]")
        for line in open(SynthLib):
            sline = line.strip()
            if sline:
                mol = Chem.MolFromSmiles(sline.split()[0])
                if Ro2Filtration:
                    mol = AddHs(mol)
                    MolW = ExactMolWt(mol)
                    LogP = MolLogP(mol)
                    HDC = CalcNumHBD(mol)
                    HAC = CalcNumHBA(mol)
                    if MolW>200 or LogP > 2 or HDC > 2 or HAC > 4:
                        continue
                availableSynthons[sline.split()[0]] = {}
                availableSynthons[sline.split()[0]]["BBs"] = sline.split()[1]
                if self.FindAnaloguesOfMissingSynthons:
                    availableSynthons[sline.split()[0]]["fp_b"] = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
                availableSynthons[sline.split()[0]]["n_atoms"] = mol.GetNumAtoms()
                availableSynthons[sline.split()[0]]["n_rings"] = CalcNumRings(mol)
                availableSynthons[sline.split()[0]]["marks"] = sorted(
                    [sline.split()[0][m.start():m.start() + 2] + sline.split()[0][m.end() - 4:m.end()] for m in
                     re.finditer(pat, sline.split()[0])])
                availableSynthons[sline.split()[0]]["marksVallences"] = "+".join(sorted([atom.GetSymbol() + ":" + 
                            str(atom.GetTotalDegree()) for atom in mol.GetAtoms() if atom.GetAtomMapNum() != 0]))
        print("Lib BB reading time:")
        print(datetime.datetime.now() - fragBegTime)
        return availableSynthons

    def __getSetup(self, reactionsToWorkWith):
        if self.__setupFile:
            tree = ET.parse(self.__setupFile)
        else:
            tree = ET.parse("Setup.xml")
        N_SyntOn_setup = tree.getroot()
        allReactions, reactionSetup = self.__getReactionSMARTS(N_SyntOn_setup)
        if self.__fragmentationMode == "use_all":
            reaction_list = allReactions[:-1]
        elif self.__fragmentationMode == "include_only" or self.__fragmentationMode == "one_by_one":
            reaction_list = self.__getReactionList(reactionsToWorkWith, allReactions)
        elif self.__fragmentationMode == "block by block":
            reaction_list = []
            for block in reactionsToWorkWith.split(";"):
                reaction_list.append(self.__getReactionList(block, allReactions))
        elif self.__fragmentationMode == "exclude_some":
            reaction_list = [reaction for reaction in allReactions if
                             reaction not in self.__getReactionList(reactionsToWorkWith, allReactions)]
        #reaction_list -> Ids of reactions that were specified by user to be used for the fragmentation
        return reactionSetup, reaction_list

    def __getMacroCycleSetup(self, reactionsToWorkWith):
        if self.__macroCycleSetupFile:
            tree = ET.parse(self.__macroCycleSetupFile)
        else:
            tree = ET.parse("Setup.xml")
        N_SyntOn_setup = tree.getroot()
        allReactions, reactionSetup = self.__getReactionSMARTS(N_SyntOn_setup)
        if self.__fragmentationMode == "use_all":
            reaction_list = allReactions[:-1]
        elif self.__fragmentationMode == "include_only" or self.__fragmentationMode == "one_by_one":
            reaction_list = self.__getReactionList(reactionsToWorkWith, allReactions, MacroCycles=True)
        elif self.__fragmentationMode == "block by block":
            reaction_list = []
            for block in self.__macroCyclicReactionsToWorkWith.split(";"):
                reaction_list.append(self.__getReactionList(block, allReactions, MacroCycles=True))
        elif self.__fragmentationMode == "exclude_some":
            reaction_list = [reaction for reaction in allReactions if
                             reaction not in self.__getReactionList(reactionsToWorkWith, allReactions, MacroCycles=True)]
        #reaction_list -> Ids of reactions that were specified by user to be used for the fragmentation
        return reactionSetup, reaction_list

    def __getReactionSMARTS(self, N_SyntOn_setup: ET.Element):
        reactionSetup = {}
        for child in N_SyntOn_setup:
            if child.tag == "AvailableReactions":
                for ch in child:
                    for subCh in ch:
                        if subCh.get('SMARTS'):
                            reactionSetup[subCh.tag] = {}
                            reactionSetup[subCh.tag]["Name"] = subCh.get('name')
                            reactionSetup[subCh.tag]["SMARTS"] = subCh.get('SMARTS')
                            reactionSetup[subCh.tag]["Labels"] = subCh.get('Labels')
                            reactionSetup[subCh.tag]["ReconstructionSMARTS"] = subCh.get('ReconstructionReaction')
        allReactions = list(reactionSetup.keys()) #allReactions -> list of all available reactions
        return allReactions, reactionSetup

    def __getReactionList(self, reactList: str, allReactions: list, MacroCycles=False) -> list:
        reaction_list = []
        for reactCode in reactList.strip().split(","):
            reactCode = reactCode.strip()
            if "-" in reactCode:
                l = reactCode.split("-")[0]
                u = reactCode.split("-")[1]
                if "." not in l:
                    li = allReactions.index(l + ".1")
                else:
                    li = allReactions.index(l)
                if "." not in u and not MacroCycles:
                    ui = allReactions.index("R" + str(int(u[1:]) + 1) + ".1")
                elif "." not in u:
                    ui = allReactions.index("MR" + str(int(u[2:]) + 1) + ".1")
                else:
                    ui = allReactions.index(u) + 1
                reaction_list.extend(allReactions[li:ui])
            else:
                if "." in reactCode:
                    reaction_list.append(reactCode)
                else:
                    li = allReactions.index(reactCode + ".1")
                    if reactCode == allReactions[-1].split(".")[0]:
                        reaction_list.extend(allReactions[li:])
                    else:
                        if not MacroCycles:
                            ui = allReactions.index("R" + str(int(reactCode[1:]) + 1) + ".1")
                        else:
                            ui = allReactions.index("MR" + str(int(reactCode[2:]) + 1) + ".1")
                        reaction_list.extend(allReactions[li:ui])
        return reaction_list

    def __firstMolCut(self, mol ):
        mol.UpdatePropertyCache()
        Chem.GetSymmSSSR(mol)
        rings = mol.GetRingInfo().AtomRings()
        macroCycle = False
        for sets in rings:
            if len(sets) > 11:
                macroCycle = True
                break
        if macroCycle:
            reactionsToWorkWith = self.__macroCyclicReactionsToWorkWith
            reactionSetup = self.__reactionMacroCycleSetup
        else:
            reactionsToWorkWith = self.__reactionsToWorkWith
            reactionSetup = self.__reactionSetup
        allSynthons = {}
        allSyntheticPathways = {}
        allSyntheticPathways["zeroCombin"] = syntheticPathway(name="zeroCombin", reagentsNumber=0,
                                                            synthPathwayReactions=None, cutLevel=0)
        allSynthons["InitMol"] = synthon(Chem.MolToSmiles(mol, canonical=True),
                                         syntheticPathway={allSyntheticPathways["zeroCombin"]}, cutLevel=0)
        allSyntheticPathways["zeroCombin"].participatingSynthon.append(allSynthons["InitMol"])
        cutLevel = 1
        cutCount = 1
        for ind,rId in enumerate(reactionsToWorkWith):
            cuttingRule = Reactions.ReactionFromSmarts(reactionSetup[rId]['SMARTS'])
            products = cuttingRule.RunReactants((mol,))
            if products:
                setsCount = len(products)
                for i in range(setsCount):
                    prodSet = products[i]
                    if not prodSet:
                        continue
                    synthonsProdSet = []
                    Ignore = False
                    for product in prodSet:
                        if product.GetNumHeavyAtoms() < 3 or Chem.MolToSmiles(product, canonical=True) in self.__SmilesToIgnore:
                            Ignore = True
                            break
                    if Ignore:
                        continue
                    # print(rId)
                    for product in prodSet:

                        Labels = reactionSetup[rId]['Labels']
                        for lablesSet in Labels.split("|"):
                            if macroCycle:
                                labledSynthon =self.__getMacroCycleLabledSmiles(product, lablesSet.split(";"), cutCount=None)
                            else:
                                labledSynthon = self.__getLabledSmiles(product, lablesSet.split(";"), cutCount=None)
                        if not labledSynthon or labledSynthon.count(":") > self.__maxNumberOfReactionCentersPerFragment or (
                                CalcNumRings(Chem.MolFromSmiles(labledSynthon)) == 0 and
                                labledSynthon.count(":") > (Chem.MolFromSmiles(labledSynthon).GetNumHeavyAtoms()+1) / 3) or \
                                (CalcNumRings(Chem.MolFromSmiles(labledSynthon)) != 0 and labledSynthon.count(
                                    ":") > (Chem.MolFromSmiles(labledSynthon).GetNumHeavyAtoms()+1) / 2):
                        #if labledSynthon and labledSynthon.count(":") > self.__maxNumberOfReactionCentersPerFragment:
                            Ignore = True
                            break
                        if labledSynthon:
                            synthonsProdSet.append(labledSynthon)
                    if Ignore:
                        continue
                    if synthonsProdSet:
                        reagSetName = rId + "_" + str(i)
                        reagSetKey = ".".join(synthonsProdSet)
                        if reagSetKey not in allSyntheticPathways:
                            if self.SynthLib:
                                allSyntheticPathways[reagSetKey] = syntheticPathway(name=reagSetName,
                                    reagentsNumber=len(prodSet), synthPathwayReactions=[reactionSetup[rId]['Name']],
                                                                              cutLevel=1, SynthLibProvided = True)
                            else:
                                allSyntheticPathways[reagSetKey] = syntheticPathway(name=reagSetName,
                                        reagentsNumber=len(prodSet), synthPathwayReactions=[reactionSetup[rId]['Name']],
                                                                               cutLevel=1)
                            allSyntheticPathways[reagSetKey].directParentsSynthPathways.append(allSyntheticPathways["zeroCombin"])

                        else:
                            continue
                        for synth in synthonsProdSet:
                            if synth not in allSynthons:
                                if self.SynthLib:
                                    allSynthons[synth] = synthon(synth,
                                                                 syntheticPathway=[allSyntheticPathways[reagSetKey]],
                                                                 cutLevel=cutLevel, SynthLibProvided = True)
                                else:
                                    allSynthons[synth] = synthon(synth, syntheticPathway=[allSyntheticPathways[reagSetKey]],
                                                             cutLevel=cutLevel)
                            elif allSyntheticPathways[reagSetKey] not in set(allSynthons[synth].syntheticPathway):
                                allSynthons[synth].syntheticPathway.append(allSyntheticPathways[reagSetKey])
                            allSyntheticPathways[reagSetKey].participatingSynthon.append(allSynthons[synth])
                            if self.__fragmentationMode == "one_by_one":
                                allSynthons[synth].rIdsToGetIt.append(ind)
                        #allSyntheticPathways[reagSetKey].printDetailedReagentsSetInfo()
                if self.__fragmentationMode == "one_by_one":
                    return allSynthons, allSyntheticPathways, ind
        return allSynthons, allSyntheticPathways

    def __cutOneSynthonHierarchically(self, synthOld, oldComb, allSynthons, allSyntheticPathways, cutLevel):
        if synthOld.smiles == "CC(Oc1cc(-c2c(C[NH:20]C)nn(C)c2C#N)cnc1N)c1cc(F)ccc1[CH:10]=O":
            print("*************************************************")
        mol = Chem.MolFromSmiles(synthOld.smiles)
        mol.UpdatePropertyCache()
        Chem.GetSymmSSSR(mol)
        mol.GetRingInfo().NumRings()
        synthOld.AtomNumbers = mol.GetNumAtoms()
        successfulCut = False
        for ind,rId in enumerate(self.__reactionsToWorkWith):
            if self.__fragmentationMode == "one_by_one" and ind < synthOld.rIdsToGetIt[-1]:
                continue
            cuttingRule = Reactions.ReactionFromSmarts(self.__reactionSetup[rId]['SMARTS'])
            products = cuttingRule.RunReactants((mol,))
            if products:
                setsCount = len(products)
                for i in range(setsCount):
                    prodSet = products[i]
                    synthonsProdSet = []
                    Ignore = False
                    for product in prodSet:
                        if product.GetNumHeavyAtoms() < 3 or ("[V]" in Chem.MolToSmiles(product, canonical=True) and product.GetNumHeavyAtoms() < 4)\
                                or Chem.MolToSmiles(product, canonical=True) in self.__SmilesToIgnore:
                            Ignore = True
                            break
                    if Ignore:
                        continue
                    for product in prodSet:
                        Labels = self.__reactionSetup[rId]['Labels']
                        for lablesSet in Labels.split("|"):
                            labledSynthon = self.__getLabledSmiles(product, lablesSet.split(";"), cutCount=None)
                        if Chem.MolToSmiles(product) == "*c1c(C[NH:20]C)nn(C)c1C#N" or Chem.MolToSmiles(product) == "*c1cnc(N)c(OC(C)c2cc(F)ccc2[CH:10]=O)c1":
                            print("#################################")
                            print(labledSynthon)
                        if not labledSynthon:
                            print(Chem.MolToSmiles(mol, canonical=True))
                            print(".".join([Chem.MolToSmiles(x, canonical=True) for x in prodSet]))
                        elif labledSynthon.count(":") > self.__maxNumberOfReactionCentersPerFragment or (
                                CalcNumRings(Chem.MolFromSmiles(labledSynthon)) == 0 and
                                labledSynthon.count(":") > (Chem.MolFromSmiles(labledSynthon).GetNumHeavyAtoms() +1) / 3) or \
                                (CalcNumRings(Chem.MolFromSmiles(labledSynthon)) != 0 and labledSynthon.count(
                                    ":") > (Chem.MolFromSmiles(labledSynthon).GetNumHeavyAtoms() +1) / 2):
                        #if labledSynthon and labledSynthon.count(":") > self.__maxNumberOfReactionCentersPerFragment:
                            if Chem.MolToSmiles(product) == "*c1c(C[NH:20]C)nn(C)c1C#N" or Chem.MolToSmiles(
                                    product) == "*c1cnc(N)c(OC(C)c2cc(F)ccc2[CH:10]=O)c1":
                                print("Ignore1")
                            Ignore = True
                            break
                        pat = re.compile("\[\w*:\w*\]")
                        marks = frozenset(
                            [labledSynthon[m.start():m.start() + 2] + labledSynthon[m.end() - 4:m.end()] for m in
                             re.finditer(pat, labledSynthon)])
                        if marks in self.forbiddenMarks:
                            if Chem.MolToSmiles(product) == "*c1c(C[NH:20]C)nn(C)c1C#N" or Chem.MolToSmiles(
                                    product) == "*c1cnc(N)c(OC(C)c2cc(F)ccc2[CH:10]=O)c1":
                                print("Ignore2")
                            Ignore = True
                            break
                        if labledSynthon:
                            synthonsProdSet.append(labledSynthon)
                    if Ignore:
                        continue
                    if synthonsProdSet:
                        if synthOld.marks:
                            marksCheckPrev = []
                            allMarks = []
                            oldMarksPresentInNewSynthons = []
                            for labledSynthon in synthonsProdSet:
                                marksNew = sorted([labledSynthon[m.start():m.start() + 2] + labledSynthon[m.end() - 4:m.end()] for m
                                     in re.finditer(pat, labledSynthon)])
                                allMarks.extend([labledSynthon[m.start():m.start() + 2] + labledSynthon[m.end() - 4:m.end()] for m
                                     in re.finditer(pat, labledSynthon)])
                                for m in synthOld.marks:
                                    if m in marksNew:
                                        oldMarksPresentInNewSynthons.append(m)
                                if len(marksNew) == len(synthOld.marks) and marksNew != synthOld.marks:
                                    marksCheckPrev.append(0)
                                else:
                                    marksCheckPrev.append(1)
                            if sorted(oldMarksPresentInNewSynthons) == synthOld.marks:
                                marksCheckPrev.append(1)
                            if 1 not in set(marksCheckPrev):
                                if labledSynthon=="Cn1nc(C[NH:20]C)[cH:21]c1C#N":
                                    print(synthOld.marks)
                                    print(marksNew)
                                    print("Continue 1")
                                continue
                            else:
                                allNewMarks = set(allMarks)

                                checkTotalMarks = [1 for mark in synthOld.marks if mark not in allNewMarks]
                                if checkTotalMarks:
                                    if labledSynthon == "Cn1nc(C[NH:20]C)[cH:21]c1C#N":
                                        print("Continue 2")
                                    continue
                        successfulCut = True
                        reagSetName1 = oldComb.name.split("|")
                        reagSetName1.append(rId +  "_" + str(i))
                        reagSetName1.sort()
                        reagSetName = "|".join(reagSetName1)
                        reagSetKeyList = []
                        reagSetKeyList.extend(synthonsProdSet)
                        reagSetKeyList.extend([synth.smiles for synth in oldComb.participatingSynthon if synthOld != synth])
                        reagSetKeyList.sort()
                        reagSetKey = ".".join(reagSetKeyList)
                        if reagSetKey not in allSyntheticPathways:
                            if self.SynthLib:
                                allSyntheticPathways[reagSetKey] = syntheticPathway(name = reagSetName,
                                reagentsNumber=len(prodSet)+oldComb.reagentsNumber -1,
                                synthPathwayReactions=oldComb.synthPathwayReactions,
                                directParentsSynthPathways=[oldComb], cutLevel=cutLevel,  SynthLibProvided = True)
                            else:
                                allSyntheticPathways[reagSetKey] = syntheticPathway(name=reagSetName,
                                                                                  reagentsNumber=len(
                                                                                      prodSet) + oldComb.reagentsNumber - 1,
                                                                                  synthPathwayReactions=oldComb.synthPathwayReactions,
                                                                                  directParentsSynthPathways=[oldComb],
                                                                                  cutLevel=cutLevel)
                            if self.__reactionSetup[rId]['Name'] not in set(allSyntheticPathways[reagSetKey].synthPathwayReactions):
                                allSyntheticPathways[reagSetKey].synthPathwayReactions.append(self.__reactionSetup[rId]['Name'])
                            oldComb.directChildrenSynthPathways.append(allSyntheticPathways[reagSetKey])
                        elif oldComb.name != reagSetName and allSyntheticPathways[reagSetKey] not in oldComb.directParentsSynthPathways:
                            allSyntheticPathways[reagSetKey].directParentsSynthPathways.append(oldComb)
                            oldComb.directChildrenSynthPathways.append(allSyntheticPathways[reagSetKey])
                            continue
                        else:
                            continue
                        for synth in synthonsProdSet:
                            if synth not in allSynthons:
                                if self.SynthLib:
                                    allSynthons[synth] = synthon(synth,
                                                                 syntheticPathway=[allSyntheticPathways[reagSetKey]],
                                                                 cutLevel=cutLevel, SynthLibProvided = True)
                                else:
                                    allSynthons[synth] = synthon(synth,
                                            syntheticPathway=[allSyntheticPathways[reagSetKey]], cutLevel=cutLevel)

                            elif allSyntheticPathways[reagSetKey] not in set(allSynthons[synth].syntheticPathway):
                                allSynthons[synth].syntheticPathway.append(allSyntheticPathways[reagSetKey])

                            allSyntheticPathways[reagSetKey].participatingSynthon.append(allSynthons[synth])
                            allSynthons[synth].directParents.append(synthOld)
                            synthOld.directChildren.append(allSynthons[synth])
                            if self.__fragmentationMode == "one_by_one":
                                allSynthons[synth].rIdsToGetIt.append(ind)
                        parentMarksList = [m[0] for m in re.finditer(pat, synthOld.smiles)]
                        kidsMarksList = [m[0] for m in re.finditer(pat, "".join(synthonsProdSet))]
                        if '[NH2:20]' in kidsMarksList and parentMarksList.count('[NH:20]') == kidsMarksList.count('[NH:20]') + 1:
                            for synthToMark in synthonsProdSet:
                                if '[NH2:20]' in synthToMark:
                                    allSynthons[synthToMark].bivalentN = True
                        for oldSynth in oldComb.participatingSynthon:
                            if oldSynth != synthOld:
                                allSyntheticPathways[reagSetKey].participatingSynthon.append(oldSynth)
                        if cutLevel < self.__MaxNumberOfStages:
                            for synth in allSyntheticPathways[reagSetKey].participatingSynthon:
                                self.__cutOneSynthonHierarchically(allSynthons[synth.smiles],
                                                                   allSyntheticPathways[reagSetKey], allSynthons,
                                                                           allSyntheticPathways, cutLevel + 1)
                        #allSyntheticPathways[reagSetKey].printDetailedReagentsSetInfo()
                if self.__fragmentationMode == "one_by_one" and successfulCut:
                    break

    def __getMacroCycleLabledSmiles(self, productMolecule: Chem.rdchem.Mol, Labels:list, cutCount):
        productSmiles = Chem.MolToSmiles(productMolecule, canonical=True)
        for label in Labels:
            if productSmiles.find(label.split("->")[0]) != -1:
                labeledSmiles = checkLable(productSmiles, label)
                if labeledSmiles:
                    if "*" in labeledSmiles:
                        labeledSmiles = self.__getMacroCycleLabledSmiles(Chem.MolFromSmiles(labeledSmiles), Labels, cutCount)
                    return labeledSmiles
        print("WARNING! No lable were assigned to the smiles: " + productSmiles)
        return False

    def __getLabledSmiles(self, productMolecule: Chem.rdchem.Mol, Labels:list, cutCount):
        productSmiles = Chem.MolToSmiles(productMolecule, canonical=True)
        for label in Labels:
            for sublabel in label.split(","):
                if productSmiles.find(sublabel.split("->")[0]) != -1:
                    labeledSmiles = checkLable(productSmiles, sublabel)
                    if labeledSmiles:
                        return labeledSmiles
        print("WARNING! No lable were assigned to the smiles: " + productSmiles)
        return False

def fragmentMolecule(smiles, SyntOnfragmentor, simTh=-1):
    mol = readMol(smiles)
    if mol:
        RemoveStereochemistry(mol)
        rule = "[O;-1;D1;$([O-]C=O):1]>>[OH;+0;D1:1]"
        cuttingRule = Reactions.ReactionFromSmarts(rule)
        products = cuttingRule.RunReactants((mol,))
        if products:
            mol = products[0][0]
        allSyntheticPathways, allSynthons = SyntOnfragmentor.cutWithHierarchyStorred(mol)
        if SyntOnfragmentor.SynthLib != None:
            firstSynthons = getShortestSyntheticPathways(allSyntheticPathways)
            for comb in firstSynthons:
                comb.checkAvailability(SyntOnfragmentor.SynthLib, simTh, SyntOnfragmentor.FindAnaloguesOfMissingSynthons)
        return allSyntheticPathways, allSynthons
    else:
        return None,None

def analoguesLibraryGeneration(Smiles_molNameTuple, SyntOnfragmentor, outDir, simTh=-1, strictAvailabilityMode=False, desiredNumberOfNewMols=1000):
    with open(os.path.join(outDir, "SynthonsForAnalogsGenerationForMol" + str(Smiles_molNameTuple[1]) + ".smi"), "w") as outSynthons:
        allSyntheticPathways, allSynthons = fragmentMolecule(Smiles_molNameTuple[0], SyntOnfragmentor, simTh=simTh)
        """fsynthonsAfterOneCut = getShortestSyntheticPathways(allSyntheticPathways)
        shortestSynthesis = findShortestSynthPathWithAvailableSynthLib(fsynthonsAfterOneCut, showAll=False,
                                                                    firstLaunch=True)"""
        if allSyntheticPathways and len(allSyntheticPathways) > 1:
            CompletePath = False
            leafsComb = getLongestSyntheticPathways(allSyntheticPathways)
            reconstructedMols = set()
            for comb in leafsComb:
                if comb.availabilityRate == 1.0:
                    CompletePath = True
                if "MR" in comb.name:
                    continue
                synthonsDict, SynthonsForAnaloguesSynthesisLocal = comb.getSynthonsForAnaloguesGeneration(SyntOnfragmentor.SynthLib,
                                                             simTh, strictAvailabilityMode = strictAvailabilityMode)
                if synthonsDict and SynthonsForAnaloguesSynthesisLocal:
                    outSynthons.write("****************************************** " + comb.name + " ******************************************\n")
                    outSynthons.writelines(SynthonsForAnaloguesSynthesisLocal)
                    reactionsUsedInFragmentationReactions = [rid.split("_")[0] for rid in comb.name.split("|")]
                    reactionForReconstruction = SyntOnfragmentor.getReactionForReconstruction(reactionsUsedInFragmentationReactions)
                    enumerator = enumeration(outDir=outDir, Synthons=synthonsDict,
                                                   reactionSMARTS=reactionForReconstruction, maxNumberOfReactedSynthons=len(synthonsDict),
                                                   desiredNumberOfNewMols=desiredNumberOfNewMols, nCores=1, analoguesEnumeration=True)
                    #reconstructedMols.update(enumerator.newAnaloguesGeneration())
                    reconstructedMols.update(enumerator.AnaloguesGeneration())
                    if len(reconstructedMols)>=desiredNumberOfNewMols:
                        break
                    #print(comb.name)
            if not CompletePath:
                synthonsAfterOneCut = getShortestSyntheticPathways(allSyntheticPathways)
                shortestSynthesis = findShortestSynthPathWithAvailableBBlib(synthonsAfterOneCut, showAll=False,
                                                                            firstLaunch=True)
                for comb in shortestSynthesis:
                    if comb.availabilityRate == 1.0:
                        if "MR" in comb.name:
                            continue
                        #comb.printDetailedReagentsSetInfo()
                        synthonsDict, SynthonsForAnaloguesSynthesisLocal = comb.getSynthonsForAnaloguesGeneration(
                            SyntOnfragmentor.SynthLib, simTh, strictAvailabilityMode=strictAvailabilityMode)
                        #comb.printDetailedReagentsSetInfo()
                        if synthonsDict and SynthonsForAnaloguesSynthesisLocal:
                            outSynthons.write(
                                "****************************************** " + comb.name + " ******************************************\n")
                            outSynthons.writelines(SynthonsForAnaloguesSynthesisLocal)
                            reactionsUsedInFragmentationReactions = [rid.split("_")[0] for rid in comb.name.split("|")]
                            reactionForReconstruction = SyntOnfragmentor.getReactionForReconstruction(
                                reactionsUsedInFragmentationReactions)
                            enumerator = enumeration(outDir=outDir, Synthons=synthonsDict,
                                                           reactionSMARTS=reactionForReconstruction,
                                                           maxNumberOfReactedSynthons=len(synthonsDict),
                                                           desiredNumberOfNewMols=1000000, nCores=1,
                                                           analoguesEnumeration=True)
                            # reconstructedMols.update(enumerator.newAnaloguesGeneration())
                            reconstructedMols.update(enumerator.AnaloguesGeneration())

            with open(os.path.join(outDir, "AnalogsForMol" + str(Smiles_molNameTuple[1]) + ".smi"), "w") as out:
                for recMolsSmiles in reconstructedMols:
                    out.write(recMolsSmiles + "\n")

def findShortestSynthPathWithAvailableBBlib(firstSynthons:list, allCombSynthons = None, firstLaunch=True, showAll=True):
    max_rec = 0x100000
    # May segfault without this line. 0x100 is a guess at the size of each stack frame.
    resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
    sys.setrecursionlimit(max_rec)
    sys.getrecursionlimit()
    if firstLaunch:
        allCombSynthons = []
    for comb in firstSynthons:
        allCombSynthons.append(comb)
        if comb.directChildrenSynthPathways:
            findShortestSynthPathWithAvailableBBlib(comb.directChildrenSynthPathways, allCombSynthons, firstLaunch=False)
    if firstLaunch:
        #allCombSynthons.sort(key=lambda x: (x.availabilityRate, x.reagentsNumber*-1), reverse=True)
        allCombSynthons.sort(key=lambda x: (x.availabilityRate, x.reagentsNumber), reverse=True)
        rate = allCombSynthons[0].availabilityRate
        reagents = allCombSynthons[0].reagentsNumber
        if not showAll:
            shortestAvailablePathways = [allCombSynthons[0]]
            return shortestAvailablePathways
        shortestAvailablePathways = []
        add = True
        for comb in allCombSynthons:
            if comb.availabilityRate == rate and comb.reagentsNumber == reagents:
                for savedComb in shortestAvailablePathways:
                    if set([x.smiles for x in comb.participatingSynthon]) == set([x.smiles for x in savedComb.participatingSynthon]):
                        add = False
                        break
                if add:
                    shortestAvailablePathways.append(comb)
                #else:
                    #savedComb.synthPathwayReactions.extend(comb.synthPathwayReactions)
                    #savedComb.synthPathwayReactions = list(set(savedComb.synthPathwayReactions))
        print(str(len(shortestAvailablePathways)) + " equivalent synthetic pathway(s) have been found.")
        return shortestAvailablePathways

def getShortestSyntheticPathways(allSyntheticPathways):
    firstSynthons = []
    for comb in allSyntheticPathways:
        if allSyntheticPathways[comb].cutLevel == 1:
            firstSynthons.append(allSyntheticPathways[comb])
    return firstSynthons

def getLongestSyntheticPathways(allSyntheticPathways):
    leafs = []
    for comb in allSyntheticPathways:
        if not allSyntheticPathways[comb].directChildrenSynthPathways and comb != "zeroCombin":

            #print(".................----------------------------------........................")
            #print(str(allSyntheticPathways[comb].name) + " " + ".".join([x.smiles for x in allSyntheticPathways[comb].participatingSynthon]))
            leafs.append(allSyntheticPathways[comb])
    return leafs
