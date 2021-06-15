from rdkit import Chem
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem.Scaffolds.MurckoScaffold import *



def main(args):
    PGdict = {"NCbz": "[N:1][C;$(C(=O)O[CH2]c1[cH][cH][cH][cH][cH]1):2]>>[N:1]",
              "NFmoc":"[N:1][C;$(C(=O)O[CH2][CH]1c2[cH][cH][cH][cH]c2-c3[cH][cH][cH][cH]c13):2]>>[N:1]",
              "NBnz":"[N;+0;$(N[CH2]c1[cH][cH][cH][cH][cH]1);!$(N[C,S,P]=[O,S,N]):1][C;$([CH2]c1[cH][cH][cH][cH][cH]1):2]>>[N:1]",
              "COOBnz": "[O;$(O(C)C([#6])=O):1][C;$([CH2]c1[cH][cH][cH][cH][cH]1):2]>>[OH:1]",
              "Boronics":"[B;$(B(O@C)O@C):1][#6:2]>>[#6:2]",
              "Oxiranes":"[C:1]1[O:2][C:3]1>>[C:1]([OH:2])[C;+0:3]"}
    with open(args.output + ".smi", "w") as out:
        scaffoldsCount = {}
        for line in open(args.input):
            sline = line.strip()
            if sline:
                mol = readMol(sline.split()[0])
                if mol:
                    for pg in PGdict:
                        mol = removePG(PGdict[pg], mol)
                    scaffold = MurckoScaffoldSmiles(mol=mol)
                    out.write(sline + " " + Chem.MolToSmiles(mol))
                    if scaffold:
                        if scaffold not in scaffoldsCount:
                            scaffoldsCount[scaffold] = 0
                        scaffoldsCount[scaffold] += 1
                        out.write(" " + scaffold + "\n")
                    else:
                        out.write(" linearMolecule\n")
    scaffoldsCountSorted = {r: scaffoldsCount[r] for r in sorted(scaffoldsCount, key=scaffoldsCount.get, reverse=True)}
    scaffoldsCount.clear()
    with open(args.output + "_scaffoldsCounts.smi", "w") as outCounts:
        for scaffold in scaffoldsCountSorted:
            outCounts.write(scaffold + " " + str(scaffoldsCountSorted[scaffold]) + "\n")
    with open(args.output + "_cumulativeprecentage.smi", "w") as outCumPer:
        cumSum = 0
        TotalCompNumb = sum(scaffoldsCountSorted.values())
        TotalScaffNumb = len(scaffoldsCountSorted)
        for ind,scaff in enumerate(scaffoldsCountSorted):
            cumSum += scaffoldsCountSorted[scaff]
            outCumPer.write(str(int(round((ind + 1) / TotalScaffNumb * 100))) + " " + str(
                int(round(cumSum / TotalCompNumb * 100))) + "\n")
    scaffoldsCountSorted.clear()
    scaffoldPlot(args.output + "_cumulativeprecentage.smi", args.output)


def scaffoldPlot(cumPercentageFile, outName):
    from matplotlib import pyplot as plt
    from numpy import genfromtxt
    Data = genfromtxt(cumPercentageFile, delimiter=' ', names=['x', 'y'])
    fig, ax = plt.subplots()
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.plot(Data['x'], Data['y'], color="darkgreen")
    plt.ylim(ymin=0, ymax=100)
    plt.xlim(xmin=0, xmax=100)
    plt.ylabel("Fraction of BBs, %", fontweight='bold', fontsize=14)
    plt.xlabel("Fraction of scaffolds, %", fontweight='bold', fontsize=14)
    plt.title("Cumulative Scaffold Frequency Plot", fontweight='bold', fontsize=14)
    plt.savefig("Scaffolds_FreqPlot_" + outName + ".png")


def removePG(reactionRule, mol):
    q = Chem.MolFromSmarts(reactionRule.split(">>")[0])
    cuttingRule = Reactions.ReactionFromSmarts(reactionRule)
    while mol.HasSubstructMatch(q):
        products = cuttingRule.RunReactants((mol,))
        mol = products[0][0]
        mol.UpdatePropertyCache()
        Chem.GetSymmSSSR(mol)
        mol.GetRingInfo().NumRings()
    return mol

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
        if targetMol==None:
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
    if targetMol!=None:
        try:
            RemoveStereochemistry(targetMol)
        except:
            print("Molecule has been processed, but stereo could not been removed")
        return targetMol
    else:
        return None


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="ChemSpace Atlas Data collection", epilog="Yuliana Zabolotna, Starsbourg 2020", prog="ChemSpace ResAnalyser")
    parser.add_argument("-i", "--input", type=str, help="Input")
    parser.add_argument("-o", "--output", type=str, default="eMol_Scaffolds.smi", help="Output")
    args = parser.parse_args()
    main(args)