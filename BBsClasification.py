from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import *
import os, json

def main(args):
    for ind, line in enumerate(open(args.input)):
        sline = line.strip()
        if sline and ind != 0:
            if ind % 1000 == 0:
                print(str(ind))
            Classes = BBClassifier(args.SMARTSLib, molSmiles=sline.split()[0])
            if Classes:
                for Class in Classes:
                    with open(Class + ".smi", "a") as out:
                        out.write(line)
            elif Classes == None:
                with open("Not_Processed.smi", "a") as out:
                    out.write(line)
            else:
                with open("Not_Classified.smi", "a") as out:
                    out.write(line)


def BBClassifier(SMARTSLib="SMARTSLibNew.json", molSmiles=None, mol=None):
    """d_names,f_names,main_dir = listDir(args.SMARTSLib)
    SmartsLib = {}
    for dir in d_names:
        SmartsLib[dir] = {}
        d_n, f_n, m_d = listDir(os.path.join(args.SMARTSLib, dir))
        for file in f_n:
            f = open(os.path.join(args.SMARTSLib, dir, file))
            ClassName = f.readline().strip()
            SmartsLib[dir][ClassName] = {"ShouldContainAtLeastOne": [], "ShouldAlsoContain": [], "shouldNotContain": []}

            SmartsLib[dir][ClassName]["ShouldContainAtLeastOne"], SmartsLib[dir][ClassName]["ShouldAlsoContain"], \
            SmartsLib[dir][ClassName]["shouldNotContain"] = parseClassSpecification(f)
            f.close()
    with open('SMARTSLib.json', 'w') as outfile:
        json.dump(SmartsLib, outfile)"""
    Classes = []
    with open(SMARTSLib) as input:
        SmartsLib = json.load(input)
    if molSmiles!=None and mol==None:
        mol = __ChangeMolToBeReadalbeByRDKIT(molSmiles)
        if mol != None:
            mol.UpdatePropertyCache()
        else:
            print(molSmiles + " was not processed by rdkit")
            return None
    elif molSmiles==None and mol==None:
        print("ERROR! Input Smiles or Mol object should be provided")
        exit()

    for bigClass in SmartsLib:
        for subClass in SmartsLib[bigClass]:
            if classChecker(SmartsLib[bigClass][subClass]["ShouldContainAtLeastOne"],
                            SmartsLib[bigClass][subClass]["ShouldAlsoContain"],
                            SmartsLib[bigClass][subClass]["shouldNotContain"], mol):
                Classes.append(bigClass + "_" + subClass)
    return Classes

def parseClassSpecification(file):
    shouldNotContain = []
    ShouldContainAtLeastOne = []
    ShouldAlsoContain = []
    sline = file.readline().strip()
    while sline:
        if sline=="@OR@":
            sline = file.readline().strip()
            while sline and sline != "$$$$":
                ShouldContainAtLeastOne.append(sline.split()[0])
                sline = file.readline().strip()
        elif sline == "@AND@":
            sline = file.readline().strip()
            while sline and sline != "$$$$":
                ShouldAlsoContain.append(sline.split()[0])
                sline = file.readline().strip()
        elif sline == "@NOT@":
            sline = file.readline().strip()
            while sline and sline != "$$$$":
                shouldNotContain.append(sline.split()[0])
                sline = file.readline().strip()
        else:
            sline = file.readline().strip()
    return ShouldContainAtLeastOne, ShouldAlsoContain, shouldNotContain

def classChecker(ShouldContainAtLeastOne, ShouldAlsoContain, shouldNotContain, mol):
    match = False
    for query1 in ShouldContainAtLeastOne:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(query1)):
            match = True
            break
    if match and ShouldAlsoContain:
        for query2 in ShouldAlsoContain:
            if not mol.HasSubstructMatch(Chem.MolFromSmarts(query2)):
                match = False
                break
    if match and shouldNotContain:
        for query3 in shouldNotContain:
            q3 =Chem.MolFromSmarts(query3)
            #q3.UpdatePropertyCache()
            ttt = mol.HasSubstructMatch(q3)
            if ttt:
                match = False
    return match

def listDir(path):
    d_names = []
    f_names = []
    for a,b,c in os.walk(path):
        main_dir = str(a)
        d_names = b
        f_names = c
        break
    return d_names,f_names,main_dir

def __ChangeMolToBeReadalbeByRDKIT(smiles):
    try:
        targetMol = Chem.MolFromSmiles(smiles)
    except:
        modifiedSmiles = smiles
        if "[nH]" in smiles:
            modifiedSmiles = smiles.replace("[nH]", "n", 1)
        elif "n" in smiles:
            modifiedSmiles = smiles.replace("n", "[nH]", 1)
        elif "[B]" in smiles:
            modifiedSmiles = smiles.replace("[B]", "[B-]")
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
        modifiedSmiles = smiles
        if targetMol==None:
            if "[nH]" in smiles:
                modifiedSmiles = smiles.replace("[nH]", "n", 1)
            elif "n" in smiles:
                modifiedSmiles = smiles.replace("n", "[nH]", 1)
            elif "[B]" in smiles:
                modifiedSmiles = smiles.replace("[B]", "[B-]")
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
    if targetMol!=None:
        return targetMol
    else:
        return None

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="ChEMBL compounds retreival by target", epilog="Yuliana Zabolotna 2020",
                                     prog="ChEMBL compounds retreival by target")
    parser.add_argument("-sL", "--SMARTSLib", type=str, default="SMARTSLibNew.json", help="SMARTSLibrary")
    parser.add_argument("-i", "--input", type=str, help="Input SMILES file")
    args = parser.parse_args()
    main(args)

