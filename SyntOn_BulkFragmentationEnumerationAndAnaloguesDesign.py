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
srcPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(1, srcPath)
from src.UsefulFunctions import *
from src.SyntOn import *


def main(inp, SynthLibrary, outDir, simTh, strictAvailabilityMode, nCores=-1, analoguesLibGen=False,
         Ro2Filtration=False, fragmentationMode="use_all", reactionsToWorkWith = "R1-R13", MaxNumberOfStages=5,
         maxNumberOfReactionCentersPerFragment=3, desiredNumberOfNewMols = 1000, enumerationMode=False, MWupperTh=1000,
         MWlowerTh=100):

    SmilesToIgnore = ["*C(C)C", "*C(=O)C", "*C=O", "*[V]C=O", "*[V]C(C)C", "*[V]C(=O)C"]
    if simTh == -1:
        simBBselection = False
    else:
        simBBselection = True
    SyntOnfragmentor = fragmentation(fragmentationMode=fragmentationMode, reactionsToWorkWith=reactionsToWorkWith,
                    maxNumberOfReactionCentersPerFragment=maxNumberOfReactionCentersPerFragment,
                    MaxNumberOfStages = MaxNumberOfStages, FragmentsToIgnore=SmilesToIgnore,
                    SynthLibrary=SynthLibrary, FindAnaloguesOfMissingSynthons=simBBselection,
                    Ro2SynthonsFiltration=Ro2Filtration)
    if analoguesLibGen and nCores==-1:
        for molNumb, line in enumerate(open(inp)):
            if line.strip():
                Smiles = line.strip().split()[0]
                analoguesLibraryGeneration((Smiles, molNumb+1), SyntOnfragmentor, outDir, simTh = simTh,
                                           strictAvailabilityMode=strictAvailabilityMode,
                                           desiredNumberOfNewMols=desiredNumberOfNewMols)
    elif analoguesLibGen:
        Smiles_molNumb_List = []
        for molNumb, line in enumerate(open(inp)):
            if line.strip():
                Smiles = line.strip().split()[0]
                Smiles_molNumb_List.append((Smiles, molNumb+1))
        fixed_analogsGenerationFunction = partial(analoguesLibraryGeneration, outDir=outDir, simTh=simTh,
                             SyntOnfragmentor=SyntOnfragmentor, strictAvailabilityMode=strictAvailabilityMode,
                                                  desiredNumberOfNewMols=desiredNumberOfNewMols)
        nCores = args.nCores
        finalLog = []
        with ProcessPoolExecutor(max_workers=nCores) as executor:
            for out in executor.map(fixed_analogsGenerationFunction, Smiles_molNumb_List):
                finalLog.append(out)
        print(finalLog)
    elif not enumerationMode:
        outOnePath = open(inp + "_out", "w")
        outAllSynthons = open("allSythons_" + inp + "_out", "w")
        for line in open(inp):
            if line.strip():
                Smiles = line.strip().split()[0]
                allReagentSets, allSynthons = fragmentMolecule(Smiles, SyntOnfragmentor, simTh)

                if allReagentSets and len(allReagentSets)>1:
                    fsynthonsAfterOneCut = getShortestSyntheticPathways(allReagentSets)
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
                    leafsComb = getLongestSyntheticPathways(allReagentSets)
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
    else:
        reactionForReconstruction = SyntOnfragmentor.getReactionForReconstruction()
        synthons = []
        for line in open(inp):
            sline = line.strip()
            if sline and "*" not in sline:
                synthons.append(sline.split()[0])
        enumerator = enumeration(outDir=outDir, Synthons=list(set(synthons)),
                                 reactionSMARTS=reactionForReconstruction, maxNumberOfReactedSynthons=MaxNumberOfStages+1,
                                 MWupperTh=MWupperTh, MWlowerTh = MWlowerTh,
                                 desiredNumberOfNewMols = desiredNumberOfNewMols, nCores = nCores)
        reconstructedMols = enumerator.getReconstructedMols(allowedToRunSubprocesses=True)
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Compound fragmentaitiona and analogues generation.",
                                     epilog="_________________________________________________________________________________________________________________________"
                                            "\n\nCode implementation:                Yuliana Zabolotna, Alexandre Varnek\n"
                                            "                                    Laboratoire de Chémoinformatique, Université de Strasbourg.\n\n"
                                            "Knowledge base (SMARTS library):    Dmitriy M.Volochnyuk, Sergey V.Ryabukhin, Kostiantyn Gavrylenko, Olexandre Oksiuta\n"
                                            "                                    Institute of Organic Chemistry, National Academy of Sciences of Ukraine\n"
                                            "                                    Kyiv National Taras Shevchenko University\n"
                                            "2021 Strasbourg, Kiev",
                                     prog="SyntOn_BulkFragmentationEnumerationAndAnaloguesDesign", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", type=str, help="input file")
    parser.add_argument("-oD", "--outDir", type=str, help="Output directory to write analogues.")
    parser.add_argument("--SynthLibrary", type=str, default=None, help="Library of available synthons. Generated from avaialable BBs using SyntOn_BBsBulkClassificationAndSynthonization.py")
    parser.add_argument("--nCores", default=-1, type=int, help="Number of CPUs available for parallelization.")
    parser.add_argument("--simTh", default=-1, type=float, help="Similarity threshold for BB analogues search. "
                    "If not specified, only positional variational approach will be used for BBs search")
    parser.add_argument("--analoguesLibGen",  action="store_true", help="Generate library of analogues from input mol")
    parser.add_argument("--strictAvailabilityMode", action="store_true", help="Only fully synthesizable analogues are generated. "
                                    "Alternatively, unavailable synthons resulted from compound fragmentation"
                                    " will still be used for its analogues generation.")
    parser.add_argument("--Ro2Filtration", action="store_true", help="Filter input synthons library by Ro2 (MW <= 200, logP <= 2, H-bond donors count <= 2 and H-bond acceptors count <= 4)")
    parser.add_argument("--fragmentationMode", default="use_all", type=str, help="Mode of fragmentation (defines how the reaction list is specified)"
                                                                    "\nPossible options: use_all, include_only, exclude_some, one_by_one"
                                                                    "\n(default: use_all)")
    parser.add_argument("--reactionsToWorkWith", default="R1-R13", type=str, help="List of RiDs to be used."
                                                                                  "\n(default: R1-R13 (all reactions) ")
    parser.add_argument("--desiredNumberOfNewMols", default=1000, type=int,
                        help="Desired number of new compounds to be generated (in case of anaogues generation - number of analogues per compound)."
                             "\n(default: 1000)")

    parser.add_argument("--MaxNumberOfStages", default=5, type=int, help="Maximal number of stages during fragmentation or enumeration."
                                                                         "\n(default: 5)")
    parser.add_argument("--maxNumberOfReactionCentersPerFragment", default=3, type=int, help="Maximal number of reaction centers per fragment."
                                                                                             "\n(default: 3)")
    parser.add_argument("--enumerationMode",  action="store_true", help="Enumerate library using input synthons")
    parser.add_argument("--MWupperTh", default=1000, type=int,
                        help="Maximum molecular weight allowed for generated compounds."
                             "\n(default: 1000)")
    parser.add_argument("--MWlowerTh", default=100, type=int,
                        help="Minimum molecular weight allowed for generated compounds."
                             "\n(default: 100)")

    args = parser.parse_args()
    if args.nCores == -1 or args.analoguesLibGen or args.enumerationMode:
        main(args.input, args.SynthLibrary, args.outDir, args.simTh, args.strictAvailabilityMode, args.nCores,
             args.analoguesLibGen, args.Ro2Filtration, fragmentationMode=args.fragmentationMode,
             reactionsToWorkWith = args.reactionsToWorkWith,
             MaxNumberOfStages=args.MaxNumberOfStages,
             maxNumberOfReactionCentersPerFragment=args.maxNumberOfReactionCentersPerFragment,
             desiredNumberOfNewMols=args.desiredNumberOfNewMols,enumerationMode=args.enumerationMode,
             MWupperTh=args.MWupperTh, MWlowerTh=args.MWlowerTh)
    else:
        wc = countLines(args.input)
        if wc<args.nCores:
            nCores = wc
        else:
            nCores = args.nCores
        linesPerFile = wc//nCores
        outNamesList = splitFileByLines(args.input, args.input, linesPerFile)
        fixed_main = partial(main, outDir=args.outDir, simTh=args.simTh, SynthLibrary=args.SynthLibrary,
                             strictAvailabilityMode=args.strictAvailabilityMode, Ro2Filtration=args.Ro2Filtration,
                             fragmentationMode = args.fragmentationMode, reactionsToWorkWith = args.reactionsToWorkWith,
                             MaxNumberOfStages = args.MaxNumberOfStages,
                             maxNumberOfReactionCentersPerFragment = args.maxNumberOfReactionCentersPerFragment,
                             desiredNumberOfNewMols=args.desiredNumberOfNewMols)
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
