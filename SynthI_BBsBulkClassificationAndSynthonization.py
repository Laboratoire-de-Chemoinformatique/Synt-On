import os,sys
srcPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(1, srcPath)
from src.SynthI_Classifier import *
from src.SynthI_BBs import *
from src.UsefulFunctions import *
from rdkit import Chem
"""from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem.rdmolops import *"""
from concurrent.futures import ProcessPoolExecutor
"""from rdkit.Chem.Descriptors import *
from rdkit.Chem.rdMolDescriptors import *"""
from functools import partial


def main(inp,keepPG, output=None, Ro2Filtr=False):
    if not output:
        output=inp
    with open(inp) as inp, open(output + "_BBmode.smi", "w") as out, open(output + "_Synthmode.smi", "w") as outS,\
        open(output + "_NotProcessed", "w") as notProc, open(output + "_NotClassified", "w") as notClassified:
        for line in inp:
            sline = line.strip()
            if sline:
                finalSynthon = {}
                classified = False
                withoutSynthons = True
                allClasses = []
                for molSmiles in line.split()[0].split("."):
                    print(molSmiles)
                    if "[nH+]" in molSmiles:
                        molSmiles = molSmiles.replace("[nH+]", "[nH]:", 1)
                    initMol = readMol(molSmiles)
                    if initMol == None:
                        notProc.write(line)
                        notProc.flush()
                        continue
                    Classes = BBClassifier(mol=initMol)
                    allClasses.extend(Classes)
                    if Classes:
                        classified = True
                        for clas in Classes:
                            if "MedChemHighlights" not in clas and "DEL" not in clas:
                                withoutSynthons = False
                                break
                    ind2del = []
                    for ind,Class in enumerate(Classes):
                        if "MedChemHighlights" in Class or "DEL" in Class:
                            ind2del.append(ind)
                    if ind2del:
                        for ind in ind2del[::-1]:
                            Classes.pop(ind)
                    azoles,fSynt = mainSynthonsGenerator( Chem.MolToSmiles(initMol), keepPG, Classes,returnBoolAndDict=True)
                    for synth1 in fSynt:
                        if synth1 not in finalSynthon:
                            finalSynthon[synth1] = fSynt[synth1].copy()
                        else:
                            finalSynthon[synth1].update(fSynt[synth1].copy())
                if not classified:
                    notClassified.write(line)
                    notClassified.flush()
                    continue
                if not finalSynthon and withoutSynthons:
                    out.write(sline + " " + ",".join(allClasses) + " -\n")
                    out.flush()
                    continue
                if finalSynthon:
                    if azoles:
                        allClasses.append("nHAzoles_nHAzoles")
                    out.write(sline + " " + ",".join(allClasses) + " " + ".".join(list(finalSynthon)) + " " + str(len(finalSynthon)) + "\n")
                    for synth2 in finalSynthon:
                        if Ro2Filtr:
                            goodSynth, paramsValList = Ro2Filtration(synth2)
                            if not goodSynth:
                                continue
                            else:
                                outS.write(synth2 + " " + "+".join(finalSynthon[synth2]) + " " + sline + "\n")
                        else:
                            outS.write(synth2 + " " + "+".join(finalSynthon[synth2]) + " " + sline + "\n")
                elif not withoutSynthons:
                    out.write(sline + " " + ",".join(Classes) + " NoSynthonsWereGenerated\n")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="BBs classification and Synthons generation for large BBs libraries",
                                     epilog="Code implementation:                Yuliana Zabolotna, Alexandre Varnek\n"
                                            "                                    Laboratoire de Chémoinformatique, Université de Strasbourg.\n\n"
                                            "Knowledge base (SMARTS library):    Dmitriy M.Volochnyuk, Sergey V.Ryabukhin, Kostiantyn Gavrylenko, Olexandre Oksiuta\n"
                                            "                                    Institute of Organic Chemistry, National Academy of Sciences of Ukraine\n"
                                            "                                    Kyiv National Taras Shevchenko University\n"
                                            "2021 Strasbourg, Kiev",
                                     prog="SynthI_BBsBulkClassificationAndSynthonization", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", type=str, help="Input file containing 2 columns building blocks smiles and ids.")
    parser.add_argument("-o", "--output", type=str, help="Output files suffix name.")
    parser.add_argument("--keepPG", action="store_true", help="Write both protected and unprotected "
                                            "synthons to the output (concerns Boc, Bn, Fmoc, Cbz and Esters protections).")
    parser.add_argument("--Ro2Filtr", action="store_true", help="Write only synthons satisfying Ro2 (MW <= 200, logP <= 2, H-bond donors count <= 2 and H-bond acceptors count <= 4)")

    parser.add_argument("--nCores", default=-1, type=int, help="Number of available cores for parallel calculations. Memory usage is optimized, so maximal number of parallel processes can be launched.")
    args = parser.parse_args()
    if args.nCores == -1:
        main(args.input, args.keepPG, args.output, args.Ro2Filtr)
    else:
        wc = countLines(args.input)
        linesPerFile = wc // args.nCores
        outNamesList = splitFileByLines(args.input, args.input, linesPerFile)
        fixed_main = partial(main,  keepPG=args.keepPG,
                             Ro2Filtr=args.Ro2Filtr)
        nCores = args.nCores
        finalLog = []
        with ProcessPoolExecutor(max_workers=nCores) as executor:
            for out in executor.map(fixed_main, outNamesList):
                finalLog.append(out)
        with open(args.output + "_BBmode.smi", "w") as out, open(args.output + "_Synthmode.smi", "w") as outS, \
                open(args.output + "_NotProcessed", "w") as notProc, open(args.output + "_NotClassified", "w") as notClassified:
            for inp in outNamesList:
                for line in open(inp + "_BBmode.smi"):
                    if line.strip():
                        out.write(line)
                os.remove(inp + "_BBmode.smi")
                for line in open( inp + "_Synthmode.smi"):
                    if line.strip():
                        outS.write(line)
                os.remove(inp + "_Synthmode.smi")
                for line in open( inp + "_NotProcessed"):
                    if line.strip():
                        notProc.write(line)
                os.remove(inp + "_NotProcessed")
                for line in open( inp + "_NotClassified"):
                    if line.strip():
                        notClassified.write(line)
                os.remove(inp + "_NotClassified")
        for file in outNamesList:
            os.remove(file)
    with open(args.output + "_Synthmode.smi") as inpS, open("ttt", "w") as outS:
        synthons = {}
        for line in inpS:
            sline = line.strip()
            if sline and len(sline.split()) == 4:
                if sline.split()[0] not in synthons:
                    synthons[sline.split()[0]] = {"Class": set(sline.split()[1].split(",")),
                                                  "BBs": {sline.split()[2]}, "BB_Ids": {sline.split()[3]},
                                                  "Count": 1}
                else:
                    synthons[sline.split()[0]]["Class"].update(sline.split()[1].split(","))
                    synthons[sline.split()[0]]["BBs"].add(sline.split()[2])
                    synthons[sline.split()[0]]["BB_Ids"].add(sline.split()[3])
                    synthons[sline.split()[0]]["Count"] += 1
        for ind, synth in enumerate(synthons):
            outS.write(synth + " " + "+".join(synthons[synth]['BB_Ids']) + " " +
                ",".join(synthons[synth]['Class']) + " " + ";".join(synthons[synth]['BBs']) + " " + str(
                synthons[synth]['Count']) + " " + args.output + "_" + str(ind+1) + "\n")

    os.remove(args.output + "_Synthmode.smi")
    os.rename("ttt", args.output + "_Synthmode.smi")
