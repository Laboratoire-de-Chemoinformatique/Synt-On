import xml.etree.ElementTree as ET
from rdkit import Chem, DataStructs
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem import AddHs, AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdmolops import *                                                                       #is needed
from rdkit.Chem.rdMolDescriptors import CalcNumRings
from rdkit.Chem.Descriptors import *
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem.Crippen import MolLogP
import datetime, os, time, random, re, resource, sys
from multiprocessing import Process, Queue
from collections import Counter

def readMol(smiles):
    try:
        targetMol = Chem.MolFromSmiles(smiles)
        RemoveStereochemistry(targetMol)
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
                    RemoveStereochemistry(targetMol)
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
                RemoveStereochemistry(targetMol)
            except:
                if "[nH]" in modifiedSmiles:
                    try:
                        modifiedSmiles = smiles.replace("[nH]", "n")
                        targetMol = Chem.MolFromSmiles(modifiedSmiles)
                        RemoveStereochemistry(targetMol)
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

def listDir(path):
    d_names = []
    f_names = []
    for a, b, c in os.walk(path):
        main_dir = str(a)
        d_names = b
        f_names = c
        break
    return d_names, f_names, main_dir

def Ro2Filtration(synthonSmiles):
    mol = Chem.MolFromSmiles(synthonSmiles)
    functionality = synthonSmiles.count(":")
    Bivalent_electrophilic = synthonSmiles.count(":30")
    Bivalent_nucleophilic = synthonSmiles.count(":40")
    Bivalent_neutral = synthonSmiles.count(":50")
    MolW = ExactMolWt(mol) - functionality - Bivalent_electrophilic - Bivalent_nucleophilic - Bivalent_neutral
    LogP = MolLogP(mol)
    HDC = CalcNumHBD(mol)
    HAC = CalcNumHBA(mol)
    if "[NH:20]" in synthonSmiles or "[OH:20]" in synthonSmiles or "[SH:20]" in synthonSmiles:
        marksForCount = ["[NH:20]", "[OH:20]", "[SH:20]"]
        count = 0
        for m in marksForCount:
            count += synthonSmiles.count(m)
        HDC -= count
    if MolW > 200 or LogP > 2 or HDC > 2 or HAC > 4:
        return False, ["MolW=" + str(MolW), "LogP=" + str(LogP), "HDC=" + str(HDC), "HAC=" + str(HAC)]
    else:
        return True, ["MolW=" + str(MolW), "LogP=" + str(LogP), "HDC=" + str(HDC), "HAC=" + str(HAC)]

def CheckMolStructure(goodValenceSmiles, label):
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

def generateMajorTautFromSynthonSmiles(initSmiles):
    enumerator = rdMolStandardize.TautomerEnumerator()
    initMol = Chem.MolFromSmiles(initSmiles)
    nHinit = []
    for atom in initMol.GetAtoms():
        if atom.GetAtomMapNum() != 0:
            nHinit.append(atom.GetTotalNumHs())
    initMol.UpdatePropertyCache()
    Chem.GetSymmSSSR(initMol)
    tautMol = enumerator.Canonicalize(initMol)
    tautSmiles = Chem.MolToSmiles(tautMol, canonical=True)
    initSmiles = Chem.MolToSmiles(Chem.MolFromSmiles(initSmiles), canonical=True)
    if tautSmiles == initSmiles:
        return tautSmiles
    nHtaut = []
    for atom in tautMol.GetAtoms():
        if atom.GetAtomMapNum() != 0:
            nHtaut.append(atom.GetTotalNumHs())
    if nHinit == nHtaut:
        return tautSmiles
    else:
        return initSmiles

def checkLable(productSmiles:str, Label:str):
        goodValenceSmiles = None
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

                    check = CheckMolStructure(out, Label)
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
                    check = CheckMolStructure(goodValenceSmiles, Label)
                    if check:
                        break
        if not goodValenceSmiles:
            print("Problem with structure check: " + productSmiles + " " + out)
        else:
            return generateMajorTautFromSynthonSmiles(goodValenceSmiles)

def readSyntonLib(synthLibFile, Ro2Filtration=False, FindAnaloguesOfMissingBBs=False):
    fragBegTime = datetime.datetime.now()
    availableBBs = {}
    pat = re.compile("\[\w*:\w*\]")
    for line in open(synthLibFile):
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
            if FindAnaloguesOfMissingBBs:
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
