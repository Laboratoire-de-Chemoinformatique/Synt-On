from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import *
import os, json,sys
srcPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(1, srcPath)
from UsefulFunctions import *

def main(args):
    for ind, line in enumerate(open(args.input)):
        sline = line.strip()
        if sline and ind != 0:
            if ind % 1000 == 0:
                print(str(ind))
            Classes = BBClassifier(molSmiles=sline.split()[0])
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


def BBClassifier(molSmiles=None, mol=None):
    SMARTSLib = os.path.join(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0], "config" ,"SMARTSLibNew.json")
    Classes = []
    with open(SMARTSLib) as input:
        SmartsLib = json.load(input)
    if molSmiles!=None and mol==None:
        mol = readMol(molSmiles)
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
            if __classChecker(SmartsLib[bigClass][subClass]["ShouldContainAtLeastOne"],
                            SmartsLib[bigClass][subClass]["ShouldAlsoContain"],
                            SmartsLib[bigClass][subClass]["shouldNotContain"], mol):
                Classes.append(bigClass + "_" + subClass)
    return Classes

def __classChecker(ShouldContainAtLeastOne, ShouldAlsoContain, shouldNotContain, mol):
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



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Classification of building blocks. Separates provided library into several sublibraries according to the reagents classess.",
                                     epilog="Code implementation:                Yuliana Zabolotna, Alexandre Varnek\n"
                                            "                                    Laboratoire de Chémoinformatique, Université de Strasbourg.\n\n"
                                            "Knowledge base (SMARTS library):    Dmitriy M.Volochnyuk, Sergey V.Ryabukhin, Kostiantyn Gavrylenko, Olexandre Oksiuta\n"
                                            "                                    Institute of Organic Chemistry, National Academy of Sciences of Ukraine\n"
                                            "                                    Kyiv National Taras Shevchenko University\n"
                                            "2021 Strasbourg, Kiev",
                                     prog="SyntOn_Classifier", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", type=str, help="Input SMILES file")
    args = parser.parse_args()
    main(args)

