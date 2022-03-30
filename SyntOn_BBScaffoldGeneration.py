import sys,os
from rdkit import Chem
srcPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(1, srcPath)
from src.UsefulFunctions import *
from src.SyntOn_BBs import *

def main(args):
    with open(args.output + "_Scaffolds.smi", "w") as out:
        scaffoldsCount = {}
        for line in open(args.input):
            sline = line.strip()
            if sline:
                scaffold, mol = generateScaffoldForBB(sline.split()[0], returnObjects=True)
                if mol:
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


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="BBs Scaffold analysis. Generates meaningful BBs scaffolds after removing ring-containing leaving and protective groups. Count scaffolds occurrence in the provided collection of BBs, and construct cumulative scaffold frequency plot.",
                                     epilog="Code implementation:                Yuliana Zabolotna, Alexandre Varnek\n"
                                            "                                    Laboratoire de Chémoinformatique, Université de Strasbourg.\n\n"
                                            "Knowledge base (SMARTS library):    Dmitriy M.Volochnyuk, Sergey V.Ryabukhin, Kostiantyn Gavrylenko, Olexandre Oksiuta\n"
                                            "                                    Institute of Organic Chemistry, National Academy of Sciences of Ukraine\n"
                                            "                                    Kyiv National Taras Shevchenko University\n"
                                            "2021 Strasbourg, Kiev",
                                     prog="SyntOn_BBScaffoldGeneration", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", type=str, help="Input BBs file.")
    parser.add_argument("-o", "--output", type=str, help="Output files suffix name.")
    args = parser.parse_args()
    main(args)