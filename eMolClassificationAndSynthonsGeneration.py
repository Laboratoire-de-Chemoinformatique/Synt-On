from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import *
from BBsClasification import BBClassifier
from SynthI_BBs import mainSynthonsGenerator, readMol
from rdkit.Chem.rdmolops import *
from concurrent.futures import ProcessPoolExecutor
from rdkit.Chem import AddHs
from rdkit.Chem.Descriptors import *
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem.Crippen import MolLogP
from functools import partial
import os

def main(inp, outputClasses, BB_mode, Synth_mode, SMARTSLib, Marks, keepPG, output=None, Ro2Filtration=False):
    if not output:
        output=inp
    with open(inp) as inp, open(output + "BBmode.smi", "w") as out, open(output + "Synthmode.smi", "w") as outS,\
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
                    print(molSmiles)
                    query = Chem.MolFromSmarts("[#6]-[#6]-[#8]-[#6].[#6]-[#8]-[#6](-[#6])=O.[#6]-[#8]-[#6](-[#6])=O.[#6]-[#8]-[#6](-[#6])=O")
                    if initMol == None or initMol.HasSubstructMatch(query):
                        notProc.write(line)
                        notProc.flush()
                        continue
                    Classes = BBClassifier(SMARTSLib, mol=initMol)
                    allClasses.extend(Classes)
                    if Classes:
                        classified = True
                    elif Classes and outputClasses:
                        for clas in Classes:
                            if "MedChemHighlights" not in clas and "DEL" not in clas:
                                withoutSynthons = False
                    ind2del = []
                    for ind,Class in enumerate(Classes):
                        if "MedChemHighlights" in Class or "DEL" in Class:
                            ind2del.append(ind)
                    if ind2del:
                        for ind in ind2del[::-1]:
                            Classes.pop(ind)
                    if len(line.split()[0].split("."))>1:
                        azoles,fSynt = mainSynthonsGenerator(Classes, Chem.MolToSmiles(initMol), keepPG, Marks, multiFrag=True)
                    else:
                        azoles,fSynt = mainSynthonsGenerator(Classes, Chem.MolToSmiles(initMol), keepPG, Marks)
                    print(fSynt)
                    for synth1 in fSynt:
                        if synth1 not in finalSynthon:
                            finalSynthon[synth1] = fSynt[synth1].copy()
                        else:
                            finalSynthon[synth1].update(fSynt[synth1].copy())
                    print(finalSynthon)
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
                    if outputClasses and BB_mode:
                        out.write(sline + " " + ",".join(allClasses) + " " + ".".join(list(finalSynthon)) + "\n")
                    elif BB_mode:
                        out.write(sline + " " + ".".join(list(finalSynthon)) + "\n")
                    if Synth_mode:
                        for synth2 in finalSynthon:
                            if Ro2Filtration:
                                mol = Chem.MolFromSmiles(synth2)
                                mol = AddHs(mol)
                                MolW = ExactMolWt(mol)
                                LogP = MolLogP(mol)
                                HDC = CalcNumHBD(mol)
                                HAC = CalcNumHBA(mol)
                                if "[NH:20]" in synth2 or "[OH:20]" in synth2 or "[SH:20]" in synth2:
                                    marksForCount = ["[NH:20]", "[OH:20]", "[SH:20]"]
                                    count = 0
                                    for m in marksForCount:
                                        count += synth2.count(m)
                                    HDC -= count
                                if MolW > 200 or LogP > 2 or HDC > 2 or HAC > 4:
                                    continue
                            outS.write(synth2 + " " + "+".join(finalSynthon[synth2]) + " " + sline + "\n")
                elif not withoutSynthons:
                    out.write(sline + " " + ",".join(Classes) + " NoSynthonsWereGenerated\n")

    """if Synth_mode:
        os.rename(output, output + "_temp")
        synthDict = {}
        for line in open(output + "_temp"):
            sline = line.strip()
            if sline and sline.split()[0] not in synthDict:
                synthDict[sline.split()[0]] = {}
                synthDict[sline.split()[0]]["BBs"] = sline.split()[1]
                synthDict[sline.split()[0]]["Ids"] = sline.split()[2]
            elif sline:
                synthDict[sline.split()[0]]["BBs"] += "+" + sline.split()[1]
                synthDict[sline.split()[0]]["Ids"] += "+" + sline.split()[2]
        with open(output, "w") as out:
            for synth in synthDict:
                out.write(synth + " " + synthDict[synth]["BBs"] + " " + synthDict[synth]["Ids"] + "\n")"""

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

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="BBs classification and Synthons generation for eMol",
                                     epilog="Yuliana Zabolotna 2020",
                                     prog="BBs Project")
    parser.add_argument("-i", "--input", type=str, help="input building blocks smiles")
    parser.add_argument("-o", "--output", type=str, help="output")
    parser.add_argument("-oC", "--outputClasses", action="store_true", help="add classes to the output")
    parser.add_argument("-BB", "--BB_mode", action="store_true", help="BB mode for the output (each line correspond to "
                                                                      "the BB, last col contain all corresponding synthons)")
    parser.add_argument("-synth", "--Synth_mode", action="store_true", help="Synthons mode for the output "
                                                                            "(each line correspond to the synthon, "
                                                                            "last col contain all corresponding BBs)")
    parser.add_argument("-sL", "--SMARTSLib", type=str, default="SMARTSLibNew.json", help="SMARTSLibrary")
    parser.add_argument("--Marks", type=str, default="BB_Marks.xml", help="BB marks")
    parser.add_argument("--keepPG", action="store_true", help="The script will keep both protected and unprotected "
                                            "synthons (concerns Boc, Bn, Fmoc, Cbz and Esters protections).")
    parser.add_argument("--Ro2Filtration", action="store_true", help="Write only synthons satisfying Ro2 (MW <= 200, logP <= 2, H-bond donors count <= 2 and H-bond acceptors count <= 4)")

    parser.add_argument("--nCores", default=-1, type=int, help="Number of available cores.")
    args = parser.parse_args()
    if args.nCores == -1:
        main(args.input, args.outputClasses, args.BB_mode, args.Synth_mode, args.SMARTSLib, args.Marks,
              args.keepPG, args.output, args.Ro2Filtration)
    else:
        wc = countLines(args.input)
        linesPerFile = wc // args.nCores
        outNamesList = splitFileByLines(args.input, args.input, linesPerFile)
        fixed_main = partial(main, outputClasses=args.outputClasses, BB_mode=args.BB_mode, Synth_mode=args.Synth_mode,
                             SMARTSLib=args.SMARTSLib, Marks=args.Marks, keepPG=args.keepPG,
                             Ro2Filtration=args.Ro2Filtration)
        nCores = args.nCores
        finalLog = []
        with ProcessPoolExecutor(max_workers=nCores) as executor:
            for out in executor.map(fixed_main, outNamesList):
                finalLog.append(out)
        with open(args.output + "BBmode.smi", "w") as out, open(args.output + "Synthmode.smi", "w") as outS, \
                open(args.output + "_NotProcessed", "w") as notProc, open(args.output + "_NotClassified", "w") as notClassified:
            for inp in outNamesList:
                for line in open(inp + "BBmode.smi"):
                    if line.strip():
                        out.write(line)
                os.remove(inp + "BBmode.smi")
                for line in open( inp + "Synthmode.smi"):
                    if line.strip():
                        outS.write(line)
                os.remove(inp + "Synthmode.smi")
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
    with open(args.output + "Synthmode.smi") as inpS, open("ttt", "w") as outS:
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

    os.remove(args.output + "Synthmode.smi")
    os.rename("ttt", args.output + "Synthmode.smi")
