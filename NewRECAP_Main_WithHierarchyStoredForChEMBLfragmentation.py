import xml.etree.ElementTree as ET
from rdkit import Chem, DataStructs
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem import AddHs, AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdmolops import *
from rdkit.Chem.rdMolDescriptors import CalcNumRings
from rdkit.Chem.Descriptors import *
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem.Crippen import MolLogP
import datetime, os, time, random, re, resource, sys
from multiprocessing import Process, Queue
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from collections import Counter

# check if Setup file is specified OC=O.[O-]C=O
"""class SynthISynthons:
    def __init__(self, , ):"""

def main(inp, BBlib, outDir, simTh, strictAvailabilityMode, nCores=-1, analoguesLibGen=False, simBBselection=False, Ro2Filtration=False):

    SmilesToIgnore = ["*C(C)C", "*C(=O)C", "*C=O", "*[V]C=O", "*[V]C(C)C", "*[V]C(=O)C"]
    if simTh == -1:
        simBBselection = False
    SynthIfragmentor = fragmentation(mode="use all",
                    maxNumderOfReactionCentersPerFragment=3, MaxNumberOfStages = 5, FragmentsToIgnore=SmilesToIgnore,
                    BBlibrary=BBlib, FindAnaloguesOfMissingBBs=simBBselection, Ro2SynthonsFiltration=Ro2Filtration) # try applying weak refferences here for the analogues generation
    if analoguesLibGen and nCores==-1:
        for molNumb, line in enumerate(open(inp)):
            if line.strip():
                Smiles = line.strip().split()[0]
                analoguesLibraryGeneration((Smiles, molNumb+1), SynthIfragmentor, outDir, simTh = simTh,
                                           strictAvailabilityMode=strictAvailabilityMode)
    elif analoguesLibGen:
        Smiles_molNumb_List = []
        for molNumb, line in enumerate(open(inp)):
            if line.strip():
                Smiles = line.strip().split()[0]
                Smiles_molNumb_List.append((Smiles, molNumb+1))
        fixed_analogsGenerationFunction = partial(analoguesLibraryGeneration, outDir=outDir, simTh=simTh,
                             SynthIfragmentor=SynthIfragmentor, strictAvailabilityMode=strictAvailabilityMode)
        nCores = args.nCores
        finalLog = []
        with ProcessPoolExecutor(max_workers=nCores) as executor:
            for out in executor.map(fixed_analogsGenerationFunction, Smiles_molNumb_List):
                finalLog.append(out)
        print(finalLog)
    else:
        outOnePath = open(inp + "_out", "w")
        outAllSynthons = open("allSythons_" + inp + "_out", "w")
        for line in open(inp):
            if line.strip():
                Smiles = line.strip().split()[0]
                allReagentSets = fragmentMolecule(Smiles, SynthIfragmentor, simTh, storeHierarchy=True)

                if allReagentSets and len(allReagentSets)>1:
                    fsynthonsAfterOneCut = getFirstLevelSynthons(allReagentSets)
                    shortestSynthesis = findShortestSynthPathWithAvailableBBlib(fsynthonsAfterOneCut,  showAll=False,
                                                                                firstLaunch = True)
                    outOnePath.write(line.strip() + " " + ".".join([x.smiles for x in shortestSynthesis[0].participatingSynthon]) + " " +shortestSynthesis[0].name +
                              " " +  str(shortestSynthesis[0].reagentsNumber) + " " + str(shortestSynthesis[0].availabilityRate) + " ")
                    outOnePath.flush()
                    availableSynthons = {}
                    notAvaialableSynthons = []
                    for synth in shortestSynthesis[0].participatingSynthon:
                        if synth.correspondingBB != None:
                            availableSynthons[synth.smiles] = synth.correspondingBB
                        else:
                            notAvaialableSynthons.append(synth.smiles)
                    notAvaialableSynthons = list(set(notAvaialableSynthons))
                    outOnePath.write("AvailableSynthons:" + "|".join([synth + "->" + availableSynthons[synth] for synth in availableSynthons]))
                    outOnePath.flush()
                    if notAvaialableSynthons:
                        outOnePath.write(" NotAvailableSynthons:" + "|".join(notAvaialableSynthons))
                        outOnePath.flush()
                    outOnePath.write("\n")
                    outOnePath.flush()
                    """# save all combinations
                    leafsComb = getLastLevelSynthons(allReagentSets)
                    n=0
                    allWrittenCombinations = set()
                    outAllSynthons.write(line.strip())
                    for comb in leafsComb:
                        allSynthons = []
                        for synth in comb.participatingSynthon:
                            allSynthons.append(synth.smiles)
                        allSynthons = frozenset(allSynthons)
                        if allSynthons not in allWrittenCombinations:
                            n+=1
                            allWrittenCombinations.add(allSynthons)
                            allSynthons = list(allSynthons)
                            outAllSynthons.write(" " + ".".join(allSynthons))
                    outAllSynthons.write(" " + str(n) +" \n")
                    outAllSynthons.flush()"""
                    # save all synthons
                    leafsComb = getLastLevelSynthons(allReagentSets)
                    availableSynthons = []
                    notAvaialableSynthons = []
                    for comb in leafsComb:
                        for synth in comb.participatingSynthon:
                            if synth.correspondingBB != None:
                                availableSynthons.append(synth.smiles)
                            else:
                                notAvaialableSynthons.append(synth.smiles)
                    availableSynthons = list(set(availableSynthons))
                    notAvaialableSynthons = list(set(notAvaialableSynthons))
                    outAllSynthons.write(line.strip() +" AvailableSynthons:" + ".".join(availableSynthons))
                    if notAvaialableSynthons:
                        outAllSynthons.write(" notAvaialableSynthons:" + ".".join(notAvaialableSynthons))
                    outAllSynthons.write("\n")
                    outAllSynthons.flush()
                    shortestSynthesis.clear()
                    availableSynthons.clear()
                    notAvaialableSynthons.clear()
                    leafsComb.clear()
                    allReagentSets.clear()
                    fsynthonsAfterOneCut.clear()
                else:
                    outOnePath.write(line.strip() + "\n")
                    outOnePath.flush()
                    outAllSynthons.write(line)
                    outAllSynthons.flush()
        outOnePath.close()
        outAllSynthons.close()
        return "Finished"
        """print("findShortestSynthPathWithAvailableBBlib time:")
        print(datetime.datetime.now() - fragBegTime)
        for comb in shortestSynthesis:
            print(comb.printReagentsSetInfo())
        for comb in getLastLevelSynthons(allReagentSets):
            print(comb.printReagentsSetInfo())
    
        lastSynthons = getLastLevelSynthons(allReagentSets)
        synthons = []
        for comb in lastSynthons:
            synthons.extend([x.smiles for x in comb.participatingSynthon])
        print(synthons)
    
        reactionForReconstruction = SynthIfragmentor.getReactionForReconstruction()
        reconstructor = reconstruction(outDir="/data/scaffolds/SynthI", BBs = list(set(synthons)), reactionSMARTS=reactionForReconstruction, maxBlocks=6,
                                       MWth=750, minNumberOfNewMols = 1000)
        reconstructBegTime = datetime.datetime.now()
        reconstructedMols = reconstructor.getReconstructedMols(allowedToRunSubprocesses=True)
        print("Reconstructing time:")
        print(datetime.datetime.now() - reconstructBegTime)
        print(len(reconstructedMols))
        print(reconstructedMols)
        if Chem.MolToSmiles(Chem.MolFromSmiles(TestSmiles)) in reconstructedMols:
            print("Reconstructed!")
        else:
            print("Not reconstructed!")"""

class synthon:
    def __init__(self, smiles, cutLevel=1, directParent=None, directChildren=None, synthonsCombinations=None):
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
        self.synthonsCombinations = []
        if synthonsCombinations:
            self.synthonsCombinations.extend(synthonsCombinations)
        self.correspondingBB = None
        self.AtomNumbers = 0
        self.bbAnalogues = {}
        self.bivalentN = False

    def searchForSynthonAnalogues(self, bbLib, simTh):
        posAnaloguesScreeningAtomsAllowed = ["C", "F", "N", "O"]
        refMol = Chem.MolFromSmiles(self.smiles)
        refMol.UpdatePropertyCache()
        Chem.GetSymmSSSR(refMol)
        refMol.GetRingInfo().AtomRings()
        queryFP = AllChem.GetMorganFingerprintAsBitVect(refMol, radius=2, nBits=2048)
        refMarksVallences = "+".join(sorted([atom.GetSymbol() + ":" + str(atom.GetTotalDegree()) for atom in refMol.GetAtoms() if atom.GetAtomMapNum() != 0]))
        for synth in bbLib:
            if self.marks == bbLib[synth]["marks"] and synth != self.smiles and refMarksVallences == bbLib[synth]["marksVallences"]:
                if simTh != -1 and DataStructs.TanimotoSimilarity(queryFP,bbLib[synth]["fp_b"]) >= simTh :
            #if self.marks == bbLib[synth]["marks"] and DataStructs.TverskySimilarity(queryFP,bbLib[synth]["fp_b"],  0.9, 0.1) >= simTh:
                    self.bbAnalogues[synth] = bbLib[synth]["BBs"]
            # Experimental block including Positional Variational screening approach
                else:
                    """analogMolRings = CalcNumRings(analogMol)
                    refMolRings = CalcNumRings(refMol)
                    if analogMolRings == refMolRings:
                        if self.AtomNumbers == 0:
                            self.AtomNumbers = refMol.GetNumAtoms()
                        synthAtomNumbers = analogMol.GetNumAtoms()
                        if synthAtomNumbers <= self.AtomNumbers + 1 and synthAtomNumbers >= self.AtomNumbers - 1:
                            refList = [i for i in self.smiles if i.isalpha() and i!="H"]
                            qList = [i for i in synth if i.isalpha() and i!="H"]
                            qList_refList = list((Counter(qList) - Counter(refList)).elements())
                            refList_qList = list((Counter(refList) - Counter(qList)).elements())
                            if synthAtomNumbers == self.AtomNumbers and sorted(refList_qList + qList_refList) == ['c', 'n']:
                                self.bbAnalogues[synth] = bbLib[synth]["BBs"]
                            if not qList_refList and not refList_qList:
                                if analogMol.HasSubstructMatch(refMol) or refMol.HasSubstructMatch(analogMol):
                                    self.bbAnalogues[synth] = bbLib[synth]["BBs"]
                            elif qList_refList and qList_refList[0] in posAnaloguesScreeningAtomsAllowed:
                                if analogMol.HasSubstructMatch(refMol):
                                    self.bbAnalogues[synth] = bbLib[synth]["BBs"]
                            elif refList_qList and refList_qList[0] in posAnaloguesScreeningAtomsAllowed:
                                if refMol.HasSubstructMatch(analogMol):
                                    self.bbAnalogues[synth] = bbLib[synth]["BBs"]"""

                    if self.AtomNumbers == 0:
                        self.AtomNumbers = refMol.GetNumAtoms()
                    if bbLib[synth]["n_atoms"] <= self.AtomNumbers + 1 and bbLib[synth]["n_atoms"] >= self.AtomNumbers - 1:
                        refList = [i for i in self.smiles if i.isalpha() and i!="H"]
                        qList = [i for i in synth if i.isalpha() and i!="H"]
                        qList_refList = list((Counter(qList) - Counter(refList)).elements())
                        refList_qList = list((Counter(refList) - Counter(qList)).elements())
                        if bbLib[synth]["n_atoms"] == self.AtomNumbers and sorted(refList_qList + qList_refList) == ['c', 'n']:
                            refMolRings = CalcNumRings(refMol)
                            if bbLib[synth]["n_rings"] == refMolRings:
                                self.bbAnalogues[synth] = bbLib[synth]["BBs"]
                        elif not qList_refList and not refList_qList:
                            refMolRings = CalcNumRings(refMol)
                            if bbLib[synth]["n_rings"] == refMolRings:
                                self.bbAnalogues[synth] = bbLib[synth]["BBs"]
                        elif qList_refList and len(qList_refList) == 1 and qList_refList[0] in posAnaloguesScreeningAtomsAllowed:
                            refMolRings = CalcNumRings(refMol)
                            if bbLib[synth]["n_rings"] == refMolRings:
                                analogMol = Chem.MolFromSmiles(synth)
                                analogMol.UpdatePropertyCache()
                                Chem.GetSymmSSSR(analogMol)
                                analogMol.GetRingInfo().AtomRings()
                                if analogMol.HasSubstructMatch(refMol):
                                    self.bbAnalogues[synth] = bbLib[synth]["BBs"]
                        elif refList_qList and len(refList_qList) == 0 and refList_qList[0] in posAnaloguesScreeningAtomsAllowed:
                            refMolRings = CalcNumRings(refMol)
                            if bbLib[synth]["n_rings"] == refMolRings:
                                analogMol = Chem.MolFromSmiles(synth)
                                analogMol.UpdatePropertyCache()
                                Chem.GetSymmSSSR(analogMol)
                                analogMol.GetRingInfo().AtomRings()
                                if refMol.HasSubstructMatch(analogMol):
                                    self.bbAnalogues[synth] = bbLib[synth]["BBs"]

    def printSynthonInfo(self):
        print("__________________________________________________________")
        print("Synthons Information")
        print("__________________________________________________________")
        print("Synthon: " + self.smiles)
        if self.correspondingBB:
            print("Available. Corresponding BBs: " + self.correspondingBB)
        elif self.bbAnalogues:
            print("Synthon was not found in provided library of building blocks. " + str(len(self.bbAnalogues)) + " analog(s) has/have been found")
            print("\n".join([x + " " + self.bbAnalogues[x] for x in self.bbAnalogues]))
        else:
            print("Synthon was not found in provided library of building blocks")
        if self.directParents:
            print("Parent synthons: " + ".".join(list(set([x.smiles for x in self.directParents]))))
        if self.directChildren:
            print("Children synthons: " + ".".join(list(set([x.smiles for x in self.directChildren]))))
        print("BB analogues: " + str(len(self.bbAnalogues)))
        print(self.bbAnalogues)

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

class reagentsCombinations:
    def __init__(self, name , reactionsToGetIt,  reagentsNumber, cutLevel, directParentsCombinations=None, ):
        self.name = name
        self.participatingSynthon  = []
        self.directChildrenCombinations = []
        self.directParentsCombinations = []
        if directParentsCombinations:
            self.directParentsCombinations.extend(directParentsCombinations)
        self.reactionsToGetIt = []
        if reactionsToGetIt:
            self.reactionsToGetIt.extend(reactionsToGetIt)
        self.reagentsNumber = reagentsNumber
        self.cutLevel = cutLevel
        self.availabilityRate = 0

    def checkAvailability(self, BBlib, simTh, FindAnaloguesOfMissingBBs=True):
        #availableSynthons = 0
        availableAtomCounts = 0
        allAtoms = 0
        for synth in self.participatingSynthon:
            if synth.AtomNumbers==0:
                synth.AtomNumbers = Chem.MolFromSmiles(synth.smiles).GetNumAtoms()
            allAtoms += synth.AtomNumbers
            if synth.smiles in BBlib:
                #availableSynthons += 1
                availableAtomCounts += synth.AtomNumbers
                synth.correspondingBB = BBlib[synth.smiles]["BBs"]
            elif not synth.bbAnalogues and FindAnaloguesOfMissingBBs:
                synth.searchForSynthonAnalogues(BBlib, simTh)
        #self.availabilityRate = round(availableSynthons / self.reagentsNumber, 2)
        self.availabilityRate = round(availableAtomCounts / allAtoms, 2)
        if self.directChildrenCombinations:
            for comb in self.directChildrenCombinations:
                if comb.availabilityRate==0:
                    comb.checkAvailability(BBlib, simTh, FindAnaloguesOfMissingBBs)

    def printReagentsSetInfo(self):
        if self.name!="zeroCombin":
            print("**********************************************************")
            print("Reagent set Information " + self.name )
            print("**********************************************************")
            print("Reactions: " + "->".join(self.reactionsToGetIt))
            print("Required Synthons: " + ".".join([x.smiles for x in self.participatingSynthon]))
            print("Number of reagents: " + str(self.reagentsNumber))
            print("Number of stages: " + str(self.cutLevel))
            print("Availability rate (% of mol. atoms coming from available synthons): " + str(self.availabilityRate))
            if self.directChildrenCombinations:
                print("Children reagent sets: ")
                for comb in self.directChildrenCombinations:
                    print("Reactions: " + "->".join(self.reactionsToGetIt) + "|||" + "Participating synthons: " + ".".join([x.smiles for x in comb.participatingSynthon]))
            for synth in self.participatingSynthon:
                synth.printSynthonInfo()

    def getBBsForAnaloguesGeneration(self, BBlib, simTh, strictAvailabilityMode = True):
        synthonsDict = {}
        BBsForAnaloguesSynthesis = set()
        if not strictAvailabilityMode:
            for ind, synth in enumerate(self.participatingSynthon):

                if "Reagent_" + str(ind + 1) not in synthonsDict:
                    synthonsDict["Reagent_" + str(ind + 1)] = { "synthons": [], "bivalentN": synth.bivalentN }
                synthonsDict["Reagent_" + str(ind + 1)]["synthons"].append(synth.smiles)
                if not synth.correspondingBB:
                    BBsForAnaloguesSynthesis.add(synth.smiles + " missingBB originalBB\n")
                else:
                    BBsForAnaloguesSynthesis.add(synth.smiles + " " + synth.correspondingBB + " originalBB\n")
                if not synth.bbAnalogues and synth.correspondingBB:
                    synth.searchForSynthonAnalogues(BBlib, simTh)
                if synth.bbAnalogues:
                    for analog in synth.bbAnalogues:
                        BBsForAnaloguesSynthesis.add(analog + " " + synth.bbAnalogues[analog] + " " + synth.smiles + " analog\n")
                        synthonsDict["Reagent_" + str(ind + 1)]["synthons"].append(analog)
            return synthonsDict, BBsForAnaloguesSynthesis
        else:
            for ind, synth in enumerate(self.participatingSynthon):
                if "Reagent_" + str(ind + 1) not in synthonsDict:
                    synthonsDict["Reagent_" + str(ind + 1)] = { "synthons": [], "bivalentN": synth.bivalentN }
                if not synth.bbAnalogues and synth.correspondingBB:
                    synth.searchForSynthonAnalogues(BBlib, simTh)
                if synth.correspondingBB:
                    BBsForAnaloguesSynthesis.add(synth.smiles + " " + synth.correspondingBB + " originalBB\n")
                    synthonsDict["Reagent_" + str(ind + 1)]["synthons"].append(synth.smiles)
                if synth.bbAnalogues:
                    for analog in synth.bbAnalogues:
                        BBsForAnaloguesSynthesis.add(analog + " " + synth.bbAnalogues[analog] + " " + synth.smiles + " analog\n")
                        synthonsDict["Reagent_" + str(ind + 1)]["synthons"].append(analog)
                if len(synthonsDict["Reagent_" + str(ind + 1)]["synthons"]) == 0:
                    return None, None
            return synthonsDict, BBsForAnaloguesSynthesis

class reconstruction:
    def __init__(self, outDir, BBs = None, reactionSMARTS = None, maxBlocks=6, MWth=750,
                  minNumberOfNewMols = 1000, nCores=1, analoguesEnumeration=False):

        if analoguesEnumeration and BBs!=None:
            self.__BBs = None
            self.__monoFuncBB = None
            self.__poliFuncBB = None
            self.__biFuncBB = None
            self.__maxBlocks = len(BBs)
        else:
            self.__BBs = []
            self.__monoFuncBB = []
            self.__poliFuncBB = []
            self.__biFuncBB = []
            self.__maxBlocks = maxBlocks
        self.__reactions = []
        self.__perpBBsAndReactions(BBs, reactionSMARTS, analoguesEnumeration)
        for key in self.__BBs:
            print([Chem.MolToSmiles(x) for x in self.__BBs[key]])
        self.__MWth = MWth
        self.__minNumberOfNewMols = minNumberOfNewMols
        self.__nCores = nCores
        self.__outDir = outDir
        self.__genNonUniqMols = 0
        self.results = set()
        self.allReconstructedMols = []
        self.__marksCombinationsOld = {'V': {'W'}, 'W': {'Mt', 'V', 'Hs'}, 'Hs': {'W'}, 'Mt': {'W'}, 'Hf': {'Rf'}, 'Rf': {'Hf'},
                                    'Db': {'Db'},  'Bh': {'Sg'}, 'Sg': {'Bh'}}
        self.__marksCombinations = {'C:10': ['N:20', 'O:20', 'C:20', 'c:20', 'n:20', 'S:20'],
                                    'c:10': ['N:20', 'O:20', 'C:20', 'c:20', 'n:20', 'S:20'],
                                    'c:20': ['N:11', 'C:10', 'c:10'], 'C:20': ['C:10', 'c:10'],
                                    'c:21': ['N:20', 'O:20', 'n:20'], 'C:21': ['N:20', 'n:20'],
                                    'N:20': ['C:10', 'c:10', 'C:21', 'c:21', 'S:10'], 'N:11': ['c:20'],
                                    'n:20': ['C:10', 'c:10', 'C:21', 'c:21'], 'O:20': ['C:10', 'c:10', 'c:21'],
                                    'S:20': ['C:10', 'c:10'], 'S:10': ['N:20'], 'C:30': ['C:40', 'N:40'],
                                    'C:40': ['C:30'], 'C:50': ['C:50'], 'C:60': ['C:70', 'c:70'],
                                    'c:70':['C:60'], 'C:70': ['C:60'], 'N:40': ['C:30'] }


    def getReconstructedMols(self, allowedToRunSubprocesses, randomSeed=None, seed = (0,0), mainRun = True):
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
                    else:
                        self.results.update(self.__molReconsrtuction(randomSeed, partner, 1, queue=None))
                        print("Number of so far reconstructed unique molecules = " + str(len(self.results)))
            if allowedToRunSubprocesses:
                Pool, nAlive = self.__countAndMergeActiveThreads(Pool)
                while nAlive > 0:
                    time.sleep(5)
                    Pool, nAlive = self.__countAndMergeActiveThreads(Pool)

        whileCount = 0
        while mainRun:
            if not allowedToRunSubprocesses and len(self.results) >= self.__minNumberOfNewMols :
                break
            elif allowedToRunSubprocesses and self.__genNonUniqMols >= self.__minNumberOfNewMols :
                break
            whileCount += 1
            seed = list(seed)
            seed[0] +=1
            seed = tuple(seed)
            seed, randomSeed = self.__getRandomSeed(seed)
            if randomSeed is None:
                break
            self.results.update(self.getReconstructedMols(allowedToRunSubprocesses, randomSeed, seed, mainRun = False))
        d_names, f_names, main_dir = self.__listDir(self.__outDir)
        with open(os.path.join(self.__outDir, "FinalOut_withDuplicates.smi"), "ab") as out:
            for file in f_names:
                if "temp_" in file:
                    with open(os.path.join(self.__outDir, file), "rb") as f:
                        out.write(f.read())
                        out.flush()
                    os.remove(os.path.join(self.__outDir, file))
        return self.results

    def ALTERNATIVEnewAnaloguesGeneration(self):
        if type(self.__BBs) == dict:
            for rid in range(len(self.__reactions)):
                reaction = self.__reactions[rid]
                for ind1 in range(1, len(self.__BBs) + 1):
                    for ind2 in range(ind1 + 1, len(self.__BBs) + 1):
                        BBsetPartner = 'Reagent_' + str(ind2)
                        BBset = 'Reagent_' + str(ind1)
                        if self.__checkBB_reactionCombination(self.__BBs[BBset][0], self.__BBs[BBsetPartner][0], reaction):
                            for firstBB in self.__BBs[BBset]:
                                for secondBB in self.__BBs[BBsetPartner]:
                                    BBsetsUsed = set()
                                    BBsetsUsed.add(BBset)
                                    BBsetsUsed.add(BBsetPartner)
                                    self.results.update(self.__ALTERNATIVE_NEWmolAnaloguesLibEnumeration(reagent=firstBB, partner=secondBB,
                                      numberOfBBalreadyReacted=1, BBsetsUsed=BBsetsUsed, reactionToUse = rid))
            return self.results
        else:
            print("Separate reconstructur should be evocken for Molecule analogues generation "
                  "from available BBs (set analoguesEnumeration=True) and for "
                  "simple molecules enumeration from the list of BBs (set analoguesEnumeration=False).")
            exit()

    def newAnaloguesGeneration(self,  launchCount=1): # instead of getMolsAnalogues
        if type(self.__BBs) == dict:
            rid = 0
            reaction = self.__reactions[rid]
            for ind1 in range(1, len(self.__BBs) + 1):
                for ind2 in range(ind1 + 1, len(self.__BBs) + 1):
                    BBsetPartner = 'Reagent_' + str(ind2)
                    BBset = 'Reagent_' + str(ind1)
                    if self.__checkBB_reactionCombination(self.__BBs[BBset][0], self.__BBs[BBsetPartner][0], reaction):
                        for firstBB in self.__BBs[BBset]:
                            for secondBB in self.__BBs[BBsetPartner]:
                                BBsetsUsed = set()
                                BBsetsUsed.add(BBset)
                                BBsetsUsed.add(BBsetPartner)
                                self.results.update(self.__NEWmolAnaloguesLibEnumeration(reagent=firstBB, partner=secondBB,
                                  numberOfBBalreadyReacted=1, BBsetsUsed=BBsetsUsed, reactionToUse = rid))
            if not len(self.results) and launchCount <= 3:
                random.shuffle(self.__reactions)
                self.newAnaloguesGeneration(launchCount=launchCount + 1)
            return self.results
        else:
            print("Separate reconstructur should be evocken for Molecule analogues generation "
                  "from available BBs (set analoguesEnumeration=True) and for "
                  "simple molecules enumeration from the list of BBs (set analoguesEnumeration=False).")
            exit()

    def __ALTERNATIVE_NEWmolAnaloguesLibEnumeration(self, reagent, partner, numberOfBBalreadyReacted,
                                                    BBsetsUsed, reactionToUse, usedReactions = None):
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
                                for BBsetPartner in self.__BBs:
                                    if BBsetPartner not in BBsetsUsed and self.__checkBB_reactionCombination(prod,
                                                                              self.__BBs[BBsetPartner][0], reaction):
                                        for secondBB in self.__BBs[BBsetPartner]:
                                            newBBsetsUsed = set()
                                            for i in BBsetsUsed:
                                                newBBsetsUsed.add(i)
                                            newBBsetsUsed.add(BBsetPartner)
                                            subResults = self.__ALTERNATIVE_NEWmolAnaloguesLibEnumeration(reagent=prod,
                                                  partner=secondBB, numberOfBBalreadyReacted=numberOfBBalreadyReacted +1,
                                                  BBsetsUsed=newBBsetsUsed, reactionToUse=rid, usedReactions = newUsedReactions)
                                            if subResults:
                                                allProducts.update(subResults)
                    elif not functionality and numberOfBBalreadyReacted + 1 >= len(self.__BBs):
                        allProducts.add(prodSMILES)
        return list(allProducts)

    def __NEWmolAnaloguesLibEnumeration(self, reagent, partner, numberOfBBalreadyReacted, BBsetsUsed, reactionToUse):
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
                    if functionality and reactionToUse + 1< len(self.__reactions):
                        reaction = self.__reactions[reactionToUse + 1]
                        for BBsetPartner in self.__BBs:
                            if BBsetPartner not in BBsetsUsed and self.__checkBB_reactionCombination(prod,
                                                                      self.__BBs[BBsetPartner][0], reaction):
                                for secondBB in self.__BBs[BBsetPartner]:
                                    newBBsetsUsed = set()
                                    for i in BBsetsUsed:
                                        newBBsetsUsed.add(i)
                                    newBBsetsUsed.add(BBsetPartner)
                                    subResults = self.__NEWmolAnaloguesLibEnumeration(reagent=prod,
                                          partner=secondBB, numberOfBBalreadyReacted=numberOfBBalreadyReacted +1,
                                          BBsetsUsed=newBBsetsUsed, reactionToUse=reactionToUse+1)
                                    if subResults:
                                        allProducts.update(subResults)
                    elif not functionality:
                        allProducts.add(prodSMILES)
        return list(allProducts)

    def __checkBB_reactionCombination(self, reagent, partner, reaction):
        """if Chem.MolToSmiles(reagent) == "CN1CCNCCN(c2ccc(NC3CCCOCC3)c(C(=O)Nc3[nH]nc4cc[c:10]([V])cc34)c2)CC1":
            print("DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD")
            print(Chem.MolToSmiles(partner))
            print(Chem.MolToSmiles(reagent))
            print(Reactions.ReactionToSmarts(reaction))"""
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
                            if numberOfBBalreadyReacted + 1 <= self.__maxBlocks - functionality:
                                if numberOfBBalreadyReacted + 1 == self.__maxBlocks - functionality:
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

                            if prodSMILES not in allProducts: # and ExactMolWt(prod) <= self.__MWth
                                allProducts.append(prodSMILES)

        if queue is None:
            return list(set(allProducts))
        else:
            queue.put(list(set(allProducts)))

    def __perpBBsAndReactions(self, BBs, reactionSMARTS, analoguesEnumeration):
        if BBs == None:
            raise ValueError("No building blocks are provided for reconstruction")
        else:
            if analoguesEnumeration:
                self.__BBs = dict()
                for subset in BBs:
                    bbsList = [Chem.MolFromSmiles(x) for x in BBs[subset]["synthons"]]
                    self.__BBs[subset] = self.__PrepMolForReconstruction(bbsList, BBs[subset]["bivalentN"])
            else:
                monoFuncBB = []
                biFuncBB = []
                poliFuncBB = []
                for bb in BBs:
                    if bb.count(":") == 1:
                        monoFuncBB.append(bb)
                    elif bb.count(":") == 2:
                        biFuncBB.append(bb)
                    else:
                        poliFuncBB.append(bb)
                if monoFuncBB == None:
                    raise ValueError("No monofunctional building blocks are provided for reconstruction")
                else:
                    self.__monoFuncBB = self.__PrepMolForReconstruction([Chem.MolFromSmiles(x) for x in monoFuncBB])
                    random.shuffle(self.__monoFuncBB)
                if biFuncBB:
                    self.__biFuncBB = self.__PrepMolForReconstruction([Chem.MolFromSmiles(x) for x in biFuncBB])
                    random.shuffle(self.__biFuncBB)
                if poliFuncBB:
                    self.__poliFuncBB = self.__PrepMolForReconstruction([Chem.MolFromSmiles(x) for x in poliFuncBB])
                    random.shuffle(self.__poliFuncBB)
                self.__BBs = (self.__poliFuncBB, self.__biFuncBB, self.__monoFuncBB)
            if reactionSMARTS == None:
                raise ValueError("No reactions are provided for reconstruction")
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
            if seed[0] < len(self.__BBs[seed[1]]):
                randomSeed = self.__BBs[seed[1]][seed[0]]
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
                d_names, f_names, main_dir = self.__listDir(self.__outDir)
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

    def __listDir(self, path):
        d_names = []
        f_names = []
        for a, b, c in os.walk(path):
            main_dir = str(a)
            d_names = b
            f_names = c
            break
        return d_names, f_names, main_dir

class fragmentation:
    def __init__(self, mode="use all",maxNumderOfReactionCentersPerFragment = 3, MaxNumberOfStages = 5, FragmentsToIgnore = None,
                 setupFile = os.path.join(os.path.split(os.path.realpath(__file__))[0], "Setup.xml"), BBlibrary=None,
                 macroCycleSetupFile = os.path.join(os.path.split(os.path.realpath(__file__))[0], "SetupForMacrocycles.xml"),
                 FindAnaloguesOfMissingBBs = False, parsedBBlib = False, Ro2SynthonsFiltration = False):
        max_rec = 0x100000
        # May segfault without this line. 0x100 is a guess at the size of each stack frame.
        resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
        sys.setrecursionlimit(max_rec)
        sys.getrecursionlimit()
        self.__mode = mode
        if FragmentsToIgnore:
            self.__SmilesToIgnore = FragmentsToIgnore
        else:
            self.__SmilesToIgnore = []
        self.__maxNumderOfReactionCentersPerFragment = maxNumderOfReactionCentersPerFragment
        self.__MaxNumberOfStages = MaxNumberOfStages
        self.__setupFile = setupFile
        self.__reactionSetup, self.__reactionsToWorkWith = self.__getSetup()
        self.__macroCycleSetupFile = macroCycleSetupFile
        self.__reactionMacroCycleSetup, self.__macroCyclicReactionsToWorkWith = self.__getMacroCycleSetup()
        forbiddenMarks = [{'[N:11]', '[c:10]'}, {'[N:11]', '[c:20]'}, {'[N:11]', '[O:20]'},{'[N:11]', '[C:30]'},
                          {'[N:11]', '[C:10]'} , {'[N:11]', '[N:11]'}, {'[N:11]', '[c:70]'},
                          {'[N:11]', '[C:40]'}, {'[N:11]', '[C:50]'}, {'[N:11]', '[S:10]'},
                                {'[c:11]', '[C:20]'}, {'[c:70]', '[c:20]'},
            {'[c:21]', '[C:60]'},{'[c:21]', '[N:11]'},  {'[c:20]', '[c:21]'}, {'[c:70]', '[C:21]'}, {'[c:20]', '[C:21]'}]
        self.forbiddenMarks = set([frozenset(x) for x in forbiddenMarks])
        if BBlibrary:
            self.FindAnaloguesOfMissingBBs = FindAnaloguesOfMissingBBs
            if parsedBBlib:
                self.BBlib = BBlibrary
            else:
                print("Processing BB library. It may take a few minutes, depending on the library size")
                self.BBlib = self.__readBBlib(BBlibrary, Ro2SynthonsFiltration)
        else:
            self.BBlib = None

    def __readBBlib(self, BBlib, Ro2Filtration):
        fragBegTime = datetime.datetime.now()
        availableBBs = {}
        pat = re.compile("\[\w*:\w*\]")
        for line in open(BBlib):
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
                availableBBs[sline.split()[0]] = {}
                availableBBs[sline.split()[0]]["BBs"] = sline.split()[1]
                if self.FindAnaloguesOfMissingBBs:
                    availableBBs[sline.split()[0]]["fp_b"] = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
                availableBBs[sline.split()[0]]["n_atoms"] = mol.GetNumAtoms()
                availableBBs[sline.split()[0]]["n_rings"] = CalcNumRings(mol)
                availableBBs[sline.split()[0]]["marks"] = sorted(
                    [sline.split()[0][m.start():m.start() + 2] + sline.split()[0][m.end() - 4:m.end()] for m in
                     re.finditer(pat, sline.split()[0])])
                availableBBs[sline.split()[0]]["marksVallences"] = "+".join(sorted([atom.GetSymbol() + ":" + 
                            str(atom.GetTotalDegree()) for atom in mol.GetAtoms() if atom.GetAtomMapNum() != 0]))
        print("Lib BB reading time:")
        print(datetime.datetime.now() - fragBegTime)
        return availableBBs

    def cutByMultipleReactions(self, molSmiles,  cutCount=None ,
                               allSynthons = None, leaves = None, cutLevel = 1):
        if allSynthons is None:
            allSynthons = []
            leaves = []
        synthons = []
        mol = Chem.MolFromSmiles(molSmiles)
        pat = re.compile("\[\w*:\w*\]")
        marksPrevious = set([molSmiles[m.start():m.start()+2]+molSmiles[m.end()-4:m.end()] for m in re.finditer(pat, molSmiles)])
        for rId in self.__reactionsToWorkWith:
            """if cutLevel == 1 and len(leaves)>=10:
                break"""
            cuttingRule = Reactions.ReactionFromSmarts(self.__reactionSetup[rId]['SMARTS'])
            products = cuttingRule.RunReactants((mol,))
            if products:
                for prodSet in products:
                    if not prodSet:
                        continue
                    """if cutLevel == 1 and len(leaves) >= 10:
                        break"""
                    synthonsProdSet = []
                    Ignore = False
                    for product in prodSet:
                        if product.GetNumHeavyAtoms() < 3 or Chem.MolToSmiles(product, canonical=True) in self.__SmilesToIgnore:
                            Ignore = True
                            break
                    if Ignore:
                        continue
                    #print(rId)
                    for product in prodSet:
                        Labels = self.__reactionSetup[rId]['Labels']
                        for lablesSet in Labels.split("|"):
                            labledSynthon = self.__getLabledSmiles(product, lablesSet.split(";"), cutCount=None)
                        """if not labledSynthon:
                            print(Chem.MolToSmiles(mol, canonical=True))
                            print(".".join([Chem.MolFromSmiles(x) for x in prodSet]))
                        if labledSynthon.count(":") > self.__maxNumderOfReactionCentersPerFragment or (CalcNumRings(Chem.MolFromSmiles(labledSynthon)) ==0 and
                                        labledSynthon.count(":") > Chem.MolFromSmiles(labledSynthon).GetNumHeavyAtoms()/3) or \
                            (CalcNumRings(Chem.MolFromSmiles(labledSynthon)) !=0 and labledSynthon.count(":") > Chem.MolFromSmiles(labledSynthon).GetNumHeavyAtoms()/2):"""
                        if labledSynthon and labledSynthon and labledSynthon.count(":") > self.__maxNumderOfReactionCentersPerFragment:
                            Ignore = True
                            break
                        marks = set([labledSynthon[m.start():m.start() + 2] + labledSynthon[m.end() - 4:m.end()] for m in
                                     re.finditer(pat, labledSynthon)])
                        if marks in self.forbiddenMarks:
                            Ignore = True
                            break
                        if labledSynthon:
                            synthonsProdSet.append(labledSynthon)

                    if Ignore:
                        continue

                    if synthonsProdSet:
                        if marksPrevious:
                            marksCheckPrev = []
                            for labledSynthon in synthonsProdSet:
                                marksNew = set([labledSynthon[m.start():m.start() + 2] + labledSynthon[m.end() - 4:m.end()] for m in
                                     re.finditer(pat, labledSynthon)])
                                if len(marksNew) == len(marksPrevious) and marksNew != marksPrevious:
                                    marksCheckPrev.append(0)
                                else:
                                    marksCheckPrev.append(1)

                            if 1 not in set(marksCheckPrev):
                                continue
                        synthons.extend(synthonsProdSet)
                        synthons = list(set(synthons))
                        allSynthons.extend(synthons)
                        allSynthons = list(set(allSynthons))
                        for synth in synthonsProdSet:
                            if cutCount:  # For labeling
                                self.cutByMultipleReactions(synth,
                                                       cutCount + 1, allSynthons=allSynthons, leaves=leaves,
                                                       cutLevel=cutLevel + 1)
                            else:
                                self.cutByMultipleReactions(synth,
                                                       allSynthons=allSynthons, leaves=leaves, cutLevel=cutLevel + 1)
        if not synthons:
            leaves.append(molSmiles)
            leaves = list(set(leaves))
        return set(allSynthons), set(leaves)

    def cutWithHierarchyStorred(self, mol, cutLevel = 1, allSynthons=None, allReagentSets=None):
        max_rec = 0x100000
        # May segfault without this line. 0x100 is a guess at the size of each stack frame.
        resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
        sys.setrecursionlimit(max_rec)
        sys.getrecursionlimit()
        if cutLevel==1 and allSynthons==None and allSynthons==None:
            allSynthons, allReagentSets = self.__firstMolCut(mol)
        startingCombintaions = allReagentSets.copy()
        for comb in startingCombintaions:
            if allReagentSets[comb].cutLevel == 1:
                for synth in allReagentSets[comb].participatingSynthon:
                    if synth != allSynthons["InitMol"]:
                        self.__cutOneSynthonHierarchically(allSynthons[synth.smiles], allReagentSets[comb], allSynthons, allReagentSets, cutLevel+1)

        return allReagentSets

    def getReactionForReconstruction(self, reactionList=None):
        reactionForReconstruction = []
        if reactionList:
            for key in reactionList:
                reactionForReconstruction.append(self.__reactionSetup[key]["ReconstructionSMARTS"])
        else:
            for key in self.__reactionSetup:
                reactionForReconstruction.append(self.__reactionSetup[key]["ReconstructionSMARTS"])
        return reactionForReconstruction

    def __getSetup(self):
        if self.__setupFile:
            tree = ET.parse(self.__setupFile)
        else:
            tree = ET.parse("Setup.xml")
        N_SynthI_setup = tree.getroot()
        SubReactions, reactionSetup = self.__getReactionSMARTS(N_SynthI_setup)
        if self.__mode == "use all":
            reaction_list = SubReactions
        elif self.__mode == "include only" or self.__mode == "one by one":
            reaction_list = self.__getReactionList(self.__reactionsToWorkWith, SubReactions)
        elif self.__mode == "block by block":
            reaction_list = []
            for block in self.__reactionsToWorkWith.split(";"):
                reaction_list.append(self.__getReactionList(block, SubReactions))
        elif self.__mode == "exclude some":
            reaction_list = [reaction for reaction in SubReactions if
                             reaction not in self.__getReactionList(self.__reactionsToWorkWith, SubReactions)]
        #reaction_list -> Ids of reactions that were specified by user to be used for the fragmentation
        return reactionSetup, reaction_list

    def __getMacroCycleSetup(self):
        if self.__macroCycleSetupFile:
            tree = ET.parse(self.__macroCycleSetupFile)
        else:
            tree = ET.parse("Setup.xml")
        N_SynthI_setup = tree.getroot()
        SubReactions, reactionSetup = self.__getReactionSMARTS(N_SynthI_setup)
        if self.__mode == "use all":
            reaction_list = SubReactions
        elif self.__mode == "include only" or self.__mode == "one by one":
            reaction_list = self.__getReactionList(self.__macroCyclicReactionsToWorkWith, SubReactions)
        elif self.__mode == "block by block":
            reaction_list = []
            for block in self.self.__macroCyclicReactionsToWorkWith.split(";"):
                reaction_list.append(self.__getReactionList(block, SubReactions))
        elif self.__mode == "exclude some":
            reaction_list = [reaction for reaction in SubReactions if
                             reaction not in self.__getReactionList(self.__macroCyclicReactionsToWorkWith, SubReactions)]
        #reaction_list -> Ids of reactions that were specified by user to be used for the fragmentation
        return reactionSetup, reaction_list

    def __getReactionSMARTS(self, N_SynthI_setup: ET.Element):
        reactionSetup = {}
        for child in N_SynthI_setup:
            if child.tag == "AvailableReactions":
                for ch in child:
                    for subCh in ch:
                        if subCh.get('SMARTS'):
                            reactionSetup[subCh.tag] = {}
                            reactionSetup[subCh.tag]["Name"] = subCh.get('name')
                            reactionSetup[subCh.tag]["SMARTS"] = subCh.get('SMARTS')
                            reactionSetup[subCh.tag]["Labels"] = subCh.get('Labels')
                            reactionSetup[subCh.tag]["ReconstructionSMARTS"] = subCh.get('ReconstructionReaction')
        SubReactions = list(reactionSetup.keys()) #SubReactions -> list of all available reactions
        return SubReactions, reactionSetup

    def __getReactionList(self, reactList: str, SubReactions: list) -> list:
        reaction_list = []
        for reactCode in reactList.strip().split(","):
            reactCode = reactCode.strip()
            if "-" in reactCode:
                l = reactCode.split("-")[0]
                u = reactCode.split("-")[1]
                if "." not in l:
                    li = SubReactions.index(l + ".1")
                else:
                    li = SubReactions.index(l)
                if "." not in u:
                    ui = SubReactions.index("R" + str(int(u[1:]) + 1) + ".1")
                else:
                    ui = SubReactions.index(u) + 1
                reaction_list.extend(SubReactions[li:ui])
            else:
                if "." in reactCode:
                    reaction_list.append(reactCode)
                else:
                    li = SubReactions.index(reactCode + ".1")
                    if reactCode == SubReactions[-1].split(".")[0]:
                        reaction_list.extend(SubReactions[li:])
                    else:
                        ui = SubReactions.index("R" + str(int(reactCode[1:]) + 1) + ".1")
                        reaction_list.extend(SubReactions[li:ui])
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
        allReagentSets = {}
        allReagentSets["zeroCombin"] = reagentsCombinations(name="zeroCombin", reagentsNumber=0,
                                                            reactionsToGetIt=None, cutLevel=0)
        allSynthons["InitMol"] = synthon(Chem.MolToSmiles(mol, canonical=True),
                                         synthonsCombinations={allReagentSets["zeroCombin"]}, cutLevel=0)
        allReagentSets["zeroCombin"].participatingSynthon.append(allSynthons["InitMol"])
        cutLevel = 1
        cutCount = 1
        for rId in reactionsToWorkWith:
            """if cutLevel == 1 and len(leaves)>=10:
                break"""
            cuttingRule = Reactions.ReactionFromSmarts(reactionSetup[rId]['SMARTS'])
            products = cuttingRule.RunReactants((mol,))
            if products:
                setsCount = len(products)
                for i in range(setsCount):
                    prodSet = products[i]
                    if not prodSet:
                        continue
                    """if cutLevel == 1 and len(leaves) >= 10:
                        break"""
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
                        """if not labledSynthon:
                            print(Chem.MolToSmiles(mol, canonical=True))
                            print(".".join([Chem.MolFromSmiles(x) for x in prodSet]))"""
                        if labledSynthon.count(":") > self.__maxNumderOfReactionCentersPerFragment or (
                                CalcNumRings(Chem.MolFromSmiles(labledSynthon)) == 0 and
                                labledSynthon.count(":") > (Chem.MolFromSmiles(labledSynthon).GetNumHeavyAtoms()+1) / 3) or \
                                (CalcNumRings(Chem.MolFromSmiles(labledSynthon)) != 0 and labledSynthon.count(
                                    ":") > (Chem.MolFromSmiles(labledSynthon).GetNumHeavyAtoms()+1) / 2):
                        #if labledSynthon and labledSynthon.count(":") > self.__maxNumderOfReactionCentersPerFragment:
                            Ignore = True
                            break
                        if labledSynthon:
                            synthonsProdSet.append(labledSynthon)

                    if Ignore:
                        continue
                    if synthonsProdSet:
                        reagSetName = rId + "_" + str(i)
                        reagSetKey = ".".join(synthonsProdSet)
                        if reagSetKey not in allReagentSets:
                            allReagentSets[reagSetKey] = reagentsCombinations(name=reagSetName,
                                                                               reagentsNumber=len(prodSet),
                                                                               reactionsToGetIt=[
                                                                                   reactionSetup[rId]['Name']],
                                                                               cutLevel=1)
                            allReagentSets[reagSetKey].directParentsCombinations.append(allReagentSets["zeroCombin"])

                        else:
                            continue
                        for synth in synthonsProdSet:
                            if synth not in allSynthons:
                                allSynthons[synth] = synthon(synth, synthonsCombinations=[allReagentSets[reagSetKey]],
                                                             cutLevel=cutLevel)
                            elif allReagentSets[reagSetKey] not in set(allSynthons[synth].synthonsCombinations):
                                allSynthons[synth].synthonsCombinations.append(allReagentSets[reagSetKey])
                            allReagentSets[reagSetKey].participatingSynthon.append(allSynthons[synth])
                        #allReagentSets[reagSetKey].printReagentsSetInfo()
        return allSynthons, allReagentSets

    def __cutOneSynthonHierarchically(self, synthOld, oldComb, allSynthons, allReagentSets, cutLevel):
        mol = Chem.MolFromSmiles(synthOld.smiles)
        mol.UpdatePropertyCache()
        Chem.GetSymmSSSR(mol)
        mol.GetRingInfo().NumRings()
        synthOld.AtomNumbers = mol.GetNumAtoms()
        for rId in self.__reactionsToWorkWith:
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
                        if not labledSynthon:
                            print(Chem.MolToSmiles(mol, canonical=True))
                            print(".".join([Chem.MolToSmiles(x, canonical=True) for x in prodSet]))
                        elif labledSynthon.count(":") > self.__maxNumderOfReactionCentersPerFragment or (
                                CalcNumRings(Chem.MolFromSmiles(labledSynthon)) == 0 and
                                labledSynthon.count(":") > (Chem.MolFromSmiles(labledSynthon).GetNumHeavyAtoms() +1) / 3) or \
                                (CalcNumRings(Chem.MolFromSmiles(labledSynthon)) != 0 and labledSynthon.count(
                                    ":") > (Chem.MolFromSmiles(labledSynthon).GetNumHeavyAtoms() +1) / 2):
                        #if labledSynthon and labledSynthon.count(":") > self.__maxNumderOfReactionCentersPerFragment:
                            Ignore = True
                            break
                        pat = re.compile("\[\w*:\w*\]")
                        marks = frozenset(
                            [labledSynthon[m.start():m.start() + 2] + labledSynthon[m.end() - 4:m.end()] for m in
                             re.finditer(pat, labledSynthon)])
                        if marks in self.forbiddenMarks:
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
                            for labledSynthon in synthonsProdSet:
                                marksNew = set(
                                    [labledSynthon[m.start():m.start() + 2] + labledSynthon[m.end() - 4:m.end()] for m
                                     in re.finditer(pat, labledSynthon)])
                                allMarks.extend([labledSynthon[m.start():m.start() + 2] + labledSynthon[m.end() - 4:m.end()] for m
                                     in re.finditer(pat, labledSynthon)])
                                if len(marksNew) == len(synthOld.marks) and marksNew != synthOld.marks:
                                    marksCheckPrev.append(0)
                                else:

                                    marksCheckPrev.append(1)
                            if 1 not in set(marksCheckPrev):
                                continue
                            else:
                                allNewMarks = set(allMarks)

                                checkTotalMarks = [1 for mark in synthOld.marks if mark not in allNewMarks]
                                if checkTotalMarks:
                                    continue
                        reagSetName1 = oldComb.name.split("|")
                        reagSetName1.append(rId +  "_" + str(i))
                        reagSetName1.sort()
                        reagSetName = "|".join(reagSetName1)
                        reagSetKeyList = []
                        reagSetKeyList.extend(synthonsProdSet)
                        reagSetKeyList.extend([synth.smiles for synth in oldComb.participatingSynthon if synthOld != synth])
                        reagSetKeyList.sort()
                        reagSetKey = ".".join(reagSetKeyList)
                        if reagSetKey not in allReagentSets:
                            allReagentSets[reagSetKey] = reagentsCombinations(name = reagSetName,
                                reagentsNumber=len(prodSet)+oldComb.reagentsNumber -1,
                                reactionsToGetIt=oldComb.reactionsToGetIt,
                                directParentsCombinations=[oldComb], cutLevel=cutLevel)
                            if self.__reactionSetup[rId]['Name'] not in set(allReagentSets[reagSetKey].reactionsToGetIt):
                                allReagentSets[reagSetKey].reactionsToGetIt.append(self.__reactionSetup[rId]['Name'])
                            oldComb.directChildrenCombinations.append(allReagentSets[reagSetKey])
                        elif oldComb.name != reagSetName and allReagentSets[reagSetKey] not in oldComb.directParentsCombinations:
                            allReagentSets[reagSetKey].directParentsCombinations.append(oldComb)
                            oldComb.directChildrenCombinations.append(allReagentSets[reagSetKey])
                            continue
                        else:
                            continue
                        for synth in synthonsProdSet:
                            if synth not in allSynthons:
                                allSynthons[synth] = synthon(synth,
                                            synthonsCombinations=[allReagentSets[reagSetKey]], cutLevel=cutLevel)

                            elif allReagentSets[reagSetKey] not in set(allSynthons[synth].synthonsCombinations):
                                allSynthons[synth].synthonsCombinations.append(allReagentSets[reagSetKey])

                            allReagentSets[reagSetKey].participatingSynthon.append(allSynthons[synth])
                            allSynthons[synth].directParents.append(synthOld)
                            synthOld.directChildren.append(allSynthons[synth])

                        parentMarksList = [m[0] for m in re.finditer(pat, synthOld.smiles)]
                        kidsMarksList = [m[0] for m in re.finditer(pat, "".join(synthonsProdSet))]
                        if '[NH2:20]' in kidsMarksList and parentMarksList.count('[NH:20]') == kidsMarksList.count('[NH:20]') + 1:
                            for synthToMark in synthonsProdSet:
                                if '[NH2:20]' in synthToMark:
                                    allSynthons[synthToMark].bivalentN = True
                        for oldSynth in oldComb.participatingSynthon:
                            if oldSynth != synthOld:
                                allReagentSets[reagSetKey].participatingSynthon.append(oldSynth)
                        if cutLevel < self.__MaxNumberOfStages:
                            for synth in allReagentSets[reagSetKey].participatingSynthon:
                                self.__cutOneSynthonHierarchically(allSynthons[synth.smiles],
                                                                   allReagentSets[reagSetKey], allSynthons,
                                                                           allReagentSets, cutLevel + 1)
                        #allReagentSets[reagSetKey].printReagentsSetInfo()
                        if self.__mode == "one by one":
                            break

    def __checkLable(self, productSmiles:str, Label:str):
        goodValenceSmiles = None
        initialProductSmiles = productSmiles
        if Label.split("->")[0][1] == "S":
            hCount = 1
            out = productSmiles.replace(Label.split("->")[0],
                                            "[" + Label.split("->")[1].split(":")[0] + "H" + str(hCount) + ":" +
                                            Label.split("->")[1].split(":")[1] + "]")
            goodValenceSmiles = out
        else:
            for hCount in range(1, 5):
                if hCount == 1:
                    if "+" in Label and "H" not in Label:
                        out = productSmiles.replace(Label.split("->")[0],
                                                    "[" + Label.split("->")[1].split(":")[0] + "+:" +
                                                    Label.split("->")[1].split(":")[1] + "]")

                    else:
                        out = productSmiles.replace(Label.split("->")[0],
                                                    "[" + Label.split("->")[1].split(":")[0] + ":" +
                                                    Label.split("->")[1].split(":")[1] + "]")

                    check = self.__CheckMolStructure(out, Label)
                    if check:
                        goodValenceSmiles = out
                        break
                if "+" in Label and "H" not in Label:
                    out = productSmiles.replace(Label.split("->")[0],
                                                "[" + Label.split("->")[1].split(":")[0] + "H" + str(hCount) + "+:" +
                                                Label.split("->")[1].split(":")[1] + "]")
                else:
                    out = productSmiles.replace(Label.split("->")[0],
                                                "[" + Label.split("->")[1].split(":")[0] + "H" + str(hCount) + ":" +
                                                Label.split("->")[1].split(":")[1] + "]")

                newMol = Chem.MolFromSmiles(out)
                if not newMol:
                    break
                else:
                    goodValenceSmiles = out
                    check = self.__CheckMolStructure(goodValenceSmiles, Label)
                    if check:
                        break
        if not goodValenceSmiles:
            print("Problem with structure check: " + productSmiles + " " + out)
        else:
            return self.__generateMajorTautFromSynthonSmiles(goodValenceSmiles)

    def __generateMajorTautFromSynthonSmiles(self, initSmiles):
        enumerator = rdMolStandardize.TautomerEnumerator()
        initMol = Chem.MolFromSmiles(initSmiles)
        nHinit = set()
        for atom in initMol.GetAtoms():
            if atom.GetAtomMapNum() != 0:
                nHinit.add(atom.GetTotalNumHs())
        initMol.UpdatePropertyCache()
        Chem.GetSymmSSSR(initMol)
        tautMol = enumerator.Canonicalize(initMol)
        tautSmiles = Chem.MolToSmiles(tautMol, canonical=True)
        initSmiles = Chem.MolToSmiles(Chem.MolFromSmiles(initSmiles), canonical=True)
        if tautSmiles == initSmiles:
            return tautSmiles
        nHtaut = set()
        for atom in tautMol.GetAtoms():
            if atom.GetAtomMapNum() != 0:
                nHtaut.add(atom.GetTotalNumHs())
        if nHinit == nHtaut:
            return tautSmiles
        else:
            return initSmiles

    def __getMacroCycleLabledSmiles(self, productMolecule: Chem.rdchem.Mol, Labels:list, cutCount):
        productSmiles = Chem.MolToSmiles(productMolecule, canonical=True)
        for label in Labels:
            if productSmiles.find(label.split("->")[0]) != -1:
                labeledSmiles = self.__checkLable(productSmiles, label)
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
                    labeledSmiles = self.__checkLable(productSmiles, sublabel)
                    if labeledSmiles:
                        return labeledSmiles
        print("WARNING! No lable were assigned to the smiles: " + productSmiles)
        return False

    def __CheckMolStructure(self, goodValenceSmiles, label):
        vallences = {"C": 4, "N": 3, "N+": 4, "O": 2, "S:10": 6, "S:20": 2}
        try:
            mol = Chem.MolFromSmiles(goodValenceSmiles)
        except:
            return False
        else:
            if mol:
                for atom in mol.GetAtoms():
                    if atom.GetAtomMapNum() != 0:
                        symbol = atom.GetSymbol()
                        if symbol == "C" or symbol == "O" or (symbol == "N" and "+" not in label):
                            if atom.GetTotalValence() < vallences[symbol]:
                                return False
                        elif symbol == "S":
                            if atom.GetTotalValence() < vallences[symbol + ":" + str(atom.GetAtomMapNum())]:
                                return False
                        elif atom.GetTotalValence() < vallences["N+"]:
                                return False
                return True
            else:
                return False

def fragmentMolecule(smiles, SynthIfragmentor, simTh, storeHierarchy = True):
    mol = readMol(smiles)
    if mol:
        RemoveStereochemistry(mol)
        rule = "[O;-1;D1;$([O-]C=O):1]>>[OH;+0;D1:1]"
        cuttingRule = Reactions.ReactionFromSmarts(rule)
        products = cuttingRule.RunReactants((mol,))
        if products:
            mol = products[0][0]
        if storeHierarchy:
            allReagentSets = SynthIfragmentor.cutWithHierarchyStorred(mol)
        else:
            allReagentSets = SynthIfragmentor.cutByMultipleReactions(mol)
        if SynthIfragmentor.BBlib != None:
            firstSynthons = getFirstLevelSynthons(allReagentSets)
            for comb in firstSynthons:
                comb.checkAvailability(SynthIfragmentor.BBlib, simTh, SynthIfragmentor.FindAnaloguesOfMissingBBs)
        return allReagentSets
    else:
        return None

def analoguesLibraryGeneration(Smiles_molNumbTuple, SynthIfragmentor, outDir, simTh, strictAvailabilityMode):
    with open(os.path.join(outDir, "BBsForAnalogsGenerationForMol" + str(Smiles_molNumbTuple[1]) + ".smi"),
              "w") as outBBs:
        allReagentSets = fragmentMolecule(Smiles_molNumbTuple[0], SynthIfragmentor, storeHierarchy=True, simTh=simTh)
        """fsynthonsAfterOneCut = getFirstLevelSynthons(allReagentSets)
        shortestSynthesis = findShortestSynthPathWithAvailableBBlib(fsynthonsAfterOneCut, showAll=False,
                                                                    firstLaunch=True)"""
        print(Smiles_molNumbTuple)
        if allReagentSets and len(allReagentSets) > 1:
            CompletePath = False
            leafsComb = getLastLevelSynthons(allReagentSets)
            reconstructedMols = set()
            for comb in leafsComb:
                if comb.availabilityRate == 1.0:
                    CompletePath = True
                if "MR" in comb.name:
                    continue
                synthonsDict, BBsForAnaloguesSynthesisLocal = comb.getBBsForAnaloguesGeneration(SynthIfragmentor.BBlib,
                                                             simTh, strictAvailabilityMode = strictAvailabilityMode)
                comb.printReagentsSetInfo()
                if synthonsDict and BBsForAnaloguesSynthesisLocal:
                    outBBs.write("****************************************** " + comb.name + " ******************************************\n")
                    outBBs.writelines(BBsForAnaloguesSynthesisLocal)
                    reactionsUsedInFragmentationReactions = [rid.split("_")[0] for rid in comb.name.split("|")]
                    reactionForReconstruction = SynthIfragmentor.getReactionForReconstruction(reactionsUsedInFragmentationReactions)
                    print("_______________-----------------------------_______________")
                    print(synthonsDict)
                    print("_______________-----------------------------_______________")
                    reconstructor = reconstruction(outDir=outDir, BBs=synthonsDict,
                                                   reactionSMARTS=reactionForReconstruction, maxBlocks=len(synthonsDict),
                                                   MWth=750, minNumberOfNewMols=1000000, nCores=1, analoguesEnumeration=True)
                    #reconstructedMols.update(reconstructor.newAnaloguesGeneration())
                    reconstructedMols.update(reconstructor.ALTERNATIVEnewAnaloguesGeneration())
                    print(comb.name)
            if not CompletePath:
                synthonsAfterOneCut = getFirstLevelSynthons(allReagentSets)
                shortestSynthesis = findShortestSynthPathWithAvailableBBlib(synthonsAfterOneCut, showAll=False,
                                                                            firstLaunch=True)
                for comb in shortestSynthesis:
                    if comb.availabilityRate == 1.0:
                        if "MR" in comb.name:
                            continue
                        #comb.printReagentsSetInfo()
                        synthonsDict, BBsForAnaloguesSynthesisLocal = comb.getBBsForAnaloguesGeneration(
                            SynthIfragmentor.BBlib, simTh, strictAvailabilityMode=strictAvailabilityMode)
                        comb.printReagentsSetInfo()
                        if synthonsDict and BBsForAnaloguesSynthesisLocal:
                            outBBs.write(
                                "****************************************** " + comb.name + " ******************************************\n")
                            outBBs.writelines(BBsForAnaloguesSynthesisLocal)
                            reactionsUsedInFragmentationReactions = [rid.split("_")[0] for rid in comb.name.split("|")]
                            reactionForReconstruction = SynthIfragmentor.getReactionForReconstruction(
                                reactionsUsedInFragmentationReactions)
                            print("_______________-----------------------------_______________")
                            print(synthonsDict)
                            print("_______________-----------------------------_______________")
                            reconstructor = reconstruction(outDir=outDir, BBs=synthonsDict,
                                                           reactionSMARTS=reactionForReconstruction,
                                                           maxBlocks=len(synthonsDict),
                                                           MWth=750, minNumberOfNewMols=1000000, nCores=1,
                                                           analoguesEnumeration=True)
                            # reconstructedMols.update(reconstructor.newAnaloguesGeneration())
                            reconstructedMols.update(reconstructor.ALTERNATIVEnewAnaloguesGeneration())

            with open(os.path.join(outDir, "AnalogsForMol" + str(Smiles_molNumbTuple[1]) + ".smi"), "w") as out:
                for recMolsSmiles in reconstructedMols:
                    out.write(recMolsSmiles + "\n")

def findShortestSynthPathWithAvailableBBlib(firstSynthons:list, allCombSynthons = None, firstLaunch=False, showAll=True):
    max_rec = 0x100000
    # May segfault without this line. 0x100 is a guess at the size of each stack frame.
    resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
    sys.setrecursionlimit(max_rec)
    sys.getrecursionlimit()
    if firstLaunch:
        allCombSynthons = []
    for comb in firstSynthons:
        allCombSynthons.append(comb)
        if comb.directChildrenCombinations:
            findShortestSynthPathWithAvailableBBlib(comb.directChildrenCombinations, allCombSynthons)
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
                    #savedComb.reactionsToGetIt.extend(comb.reactionsToGetIt)
                    #savedComb.reactionsToGetIt = list(set(savedComb.reactionsToGetIt))
        print(str(len(shortestAvailablePathways)) + " equivalent synthetic pathways has been found.")
        return shortestAvailablePathways

def getFirstLevelSynthons(allReagentSets):
    firstSynthons = []
    for comb in allReagentSets:
        if allReagentSets[comb].cutLevel == 1:
            firstSynthons.append(allReagentSets[comb])
    return firstSynthons

def getLastLevelSynthons(allReagentSets):
    leafs = []
    for comb in allReagentSets:
        if not allReagentSets[comb].directChildrenCombinations and comb != "zeroCombin":

            #print(".................----------------------------------........................")
            #print(str(allReagentSets[comb].name) + " " + ".".join([x.smiles for x in allReagentSets[comb].participatingSynthon]))
            leafs.append(allReagentSets[comb])
    return leafs

def readMol(smiles):
    try:
        targetMol = Chem.MolFromSmiles(smiles)
    except:
        if "[nH]" in smiles:
            modifiedSmiles = smiles.replace("[nH]", "n", 1)

        elif "n" in smiles:
            modifiedSmiles = smiles.replace("n", "[nH]", 1)
        elif "[B]" in smiles:
            modifiedSmiles = smiles.replace("[B]", "[B-]")
        else:
            modifiedSmiles = smiles
        try:
            targetMol = Chem.MolFromSmiles(modifiedSmiles)
        except:
            if "[nH]" in modifiedSmiles:
                try:
                    modifiedSmiles = smiles.replace("[nH]", "n")
                    targetMol = Chem.MolFromSmiles(modifiedSmiles)
                except:
                    print(smiles + " was not processed by rdkit")
            else:
                print(smiles + " was not processed by rdkit")
    else:
        if targetMol == None:
            if "[nH]" in smiles:
                modifiedSmiles = smiles.replace("[nH]", "n", 1)
            elif "n" in smiles:
                modifiedSmiles = smiles.replace("n", "[nH]", 1)
            elif "[B]" in smiles:
                modifiedSmiles = smiles.replace("[B]", "[B-]")
            else:
                modifiedSmiles = smiles
            try:
                targetMol = Chem.MolFromSmiles(modifiedSmiles)
            except:
                if "[nH]" in modifiedSmiles:
                    try:
                        modifiedSmiles = smiles.replace("[nH]", "n")
                        targetMol = Chem.MolFromSmiles(modifiedSmiles)
                    except:
                        print(smiles + " was not processed by rdkit")
                else:
                    print(smiles + " was not processed by rdkit")
    if targetMol != None:
        return targetMol
    else:
        return None

def countLines(file):
    lines = 0
    f = open(file, "rb")
    buf_size = 1024*1024
    read_f = f.raw.read
    buf = read_f(buf_size)
    while buf:
        lines += buf.count(b'\n')
        buf = read_f(buf_size)
    f.close()
    return lines

def splitFileByLines(inpFile, outName, linesPerFile):
    n = 0
    smallFile = None
    outNamesList = []
    for ind, line in enumerate(open(inpFile)):
        if ind % linesPerFile == 0:
            if smallFile:
                smallFile.close()
            n+=1
            smallFile = open(outName + "_" + str(n), "w")
            outNamesList.append(outName + "_" + str(n))
        smallFile.write(line)
    if smallFile:
        smallFile.close()
    return outNamesList

#print(list(N_SynthI_setup))
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="N-SynthI", epilog="Yuliana Zabolotna 2020",
                                     prog="N-SynthI")
    parser.add_argument("-i", "--input", type=str, help="input file")
    parser.add_argument("-oD", "--outDir", type=str, help="output directory")
    parser.add_argument("-bbL", "--BBlib", type=str, default=None, help="library of BBs to compare to")
    parser.add_argument("--nCores", default=-1, type=int, help="Number of available cores.")
    parser.add_argument("--simTh", default=-1, type=float, help="Similarity threshold for BB analogues search. "
                    "If not specified, only positional variational approach will be used for BBs search")
    parser.add_argument("--analoguesLibGen",  action="store_true", help="Generate library of analogues from input mol")
    parser.add_argument("--strictAvailabilityMode", action="store_true", help="Only fully synthesizable analogues are generated. "
                                    "Alternatively, unavailable synthons resulted from compound fragmentation"
                                    " will still be used for its analogues generation.")
    parser.add_argument("--simBBselection", action="store_true", help=" Used always with analoguesLibGen."
        "For library generation will be used not only synthons from input molecules but also their closest analogues")
    parser.add_argument("--Ro2Filtration", action="store_true", help="Filter input synthons library by Ro2 (MW <= 200, logP <= 2, H-bond donors count <= 2 and H-bond acceptors count <= 4)")

    args = parser.parse_args()
    if args.nCores == -1 or args.analoguesLibGen:
        main(args.input, args.BBlib, args.outDir, args.simTh, args.strictAvailabilityMode, args.nCores,
             args.analoguesLibGen, args.simBBselection, args.Ro2Filtration)
    else:
        wc = countLines(args.input)
        if wc<args.nCores:
            nCores = wc
        else:
            nCores = args.nCores
        linesPerFile = wc//nCores
        outNamesList = splitFileByLines(args.input, args.input, linesPerFile)
        fixed_main = partial(main, outDir=args.outDir, simTh=args.simTh, BBlib=args.BBlib,
                             strictAvailabilityMode=args.strictAvailabilityMode, Ro2Filtration=args.Ro2Filtration)
        finalLog = []
        with ProcessPoolExecutor(max_workers=nCores) as executor:
            for out in executor.map(fixed_main, outNamesList):
                finalLog.append(out)
        with open(args.input + "_out", "w") as outOnePath, open("allSythons_" + args.input + "_out", "w") as outAllSynthons:
            for inp in outNamesList:
                for line in open(inp + "_out"):
                    if line.strip():
                        outOnePath.write(line)
                os.remove(inp + "_out")
                for line in open("allSythons_" + inp + "_out"):
                    if line.strip():
                        outAllSynthons.write(line)
                os.remove("allSythons_" + inp + "_out")
        for file in outNamesList:
            os.remove(file)






