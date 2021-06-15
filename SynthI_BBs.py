import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import rdChemReactions as Reactions
from BBsClasification import BBClassifier
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem.rdmolops import *
from rdkit.Chem import AddHs
from rdkit.Chem.Descriptors import *
from rdkit.Chem.Crippen import MolLogP
import re
from rdkit.Chem.MolStandardize import rdMolStandardize


def main(args):
    molSmiles = args.smiles
    if "[nH+]" in molSmiles:
        molSmiles = molSmiles.replace("[nH+]", "[nH]:", 1)
    initMol = readMol(molSmiles)
    if initMol == None:
        print("Mol " + molSmiles + " was not processed")
        return
    Classes = BBClassifier(args.SMARTSLib, mol=initMol)
    print(Classes)
    ClassesWithSynthons = [clas for clas in Classes if "MedChemHighlights" not in clas and "DEL" not in clas ]
    if len(molSmiles.split("."))>1:
        azoles,finalSynthon = mainSynthonsGenerator(ClassesWithSynthons, Chem.MolToSmiles(initMol), args.keepPG, args.bbM, multiFrag=True)
    else:
        azoles,finalSynthon = mainSynthonsGenerator(ClassesWithSynthons, Chem.MolToSmiles(initMol), args.keepPG, args.bbM)

    print(".".join([x for x in finalSynthon]))
    for synth in finalSynthon:
        print(synth + " " + "+".join(finalSynthon[synth]))

def mainSynthonsGenerator(Classes, initSmi, keepPG, BBmarks="BB_Marks.xml", multiFrag=False, SMARTSLib="SMARTSLibNew.json"):
    tree = ET.parse(BBmarks)
    BB_Marks = tree.getroot()
    MarksSetup = __getReactionSMARTS(BB_Marks)
    polyfunc = False
    polyfuncName = []
    keepSynthonsWithPG = keepPG
    molsToWorkWith = {initSmi: set()}
    ind = 0
    finalSynthon = {}
    polyfuncInd = []
    synthonsAfterMonofuncClasses = {}
    for Cls in Classes:
        if "Bifunctional" in Cls or "Trifunctional" in Cls:
            polyfunc = True
            polyfuncName.append(Cls)
            if not ("Nboc" in Cls or "Ncbz" in Cls or "Nfmoc" in Cls  or "Ester" in Cls or "TFAc" in Cls):
                keepSynthonsWithPG = True
                break
    while ind < len(Classes):
        if "Bifunctional" in Classes[ind] or "Trifunctional" in Classes[ind]:
            polyfuncInd.append(ind)
            ind += 1
            continue
        else:
            if polyfunc and "Bifunctional" not in Classes[ind]:
                ignoreThisClass = False
                for subName in [y for x in polyfuncName for y in x.split("_")]:
                    if subName in Classes[ind] and "Bifunctional_NbnDi_Amines" not in polyfuncName:
                        ignoreThisClass = True
                        ind += 1
                        break
                if ignoreThisClass:
                    continue
            for mol in molsToWorkWith:
                synthons = classesAssignement(Classes[ind], molsToWorkWith[mol], mol, MarksSetup, keepSynthonsWithPG)
                if synthons:
                    for synth in synthons:
                        if synth not in synthonsAfterMonofuncClasses:
                            synthonsAfterMonofuncClasses[synth] = synthons[synth].copy()
                        else:
                            synthonsAfterMonofuncClasses[synth].update(synthons[synth])
            for synth in synthonsAfterMonofuncClasses:
                if synth not in molsToWorkWith:
                    molsToWorkWith[synth] = synthonsAfterMonofuncClasses[synth].copy()
                else:
                    molsToWorkWith[synth].update(synthonsAfterMonofuncClasses[synth])

        ind += 1
    if keepSynthonsWithPG or not polyfunc:
        if synthonsAfterMonofuncClasses:
            for synth in synthonsAfterMonofuncClasses:
                if synth not in finalSynthon:
                    finalSynthon[synth] = synthonsAfterMonofuncClasses[synth].copy()
                else:
                    finalSynthon[synth].update(synthonsAfterMonofuncClasses[synth])
    if polyfunc:
        for i, ind in enumerate(polyfuncInd):
            if i < len(polyfuncInd) - 1 :
                extraMols = {}
            print("dddddddddddddddddddddddddddddddddd")
            print(str(ind))
            print(polyfuncInd)
            print(Classes[ind])
            print(molsToWorkWith)
            for mol in molsToWorkWith:
                synthons = classesAssignement(Classes[ind], molsToWorkWith[mol], mol, MarksSetup, keepSynthonsWithPG)
                if synthons:
                    for synth in synthons:
                        if synth not in finalSynthon:
                            finalSynthon[synth] = synthons[synth].copy()
                        else:
                            finalSynthon[synth].update(synthons[synth])
                        if i < len(polyfuncInd) - 1 and synth not in molsToWorkWith:
                            if synth not in extraMols:
                                extraMols[synth] = synthons[synth].copy()
                            else:
                                extraMols[synth].update(synthons[synth])
            if i < len(polyfuncInd) - 1 and extraMols:
                for synth in extraMols:
                    if synth not in molsToWorkWith:
                        molsToWorkWith[synth] = extraMols[synth].copy()
        """for mol in extraMols:
            for ind in polyfuncInd:
                synthons = classesAssignement(Classes[ind], extraMols[mol], mol, MarksSetup, keepSynthonsWithPG)
                if synthons:
                    for synth in synthons:
                        if synth not in finalSynthon:
                            finalSynthon[synth] = synthons[synth].copy()"""

    for clas in Classes:
        if "Trifunctional" in clas and "Ester" in clas and "Acid" in clas:
            additionalSynthons = generateBiacideSynthonForTrifunctional(finalSynthon, clas)
            if additionalSynthons:
                for synth in additionalSynthons:
                    if synth not in finalSynthon:
                        finalSynthon[synth] = additionalSynthons[synth].copy()
                    else:
                        finalSynthon[synth].update(additionalSynthons[synth])

    additionalSynthons = azolesSynthonPostGeneration(finalSynthon)
    if additionalSynthons:
        azoles = True
        for synth in additionalSynthons:
            if synth not in finalSynthon:
                finalSynthon[synth] = additionalSynthons[synth].copy()
            else:
                finalSynthon[synth].update(additionalSynthons[synth])
    else:
        azoles = False
    if not finalSynthon and not multiFrag and "Esters_Esters" in Classes:
        ReactionLIST = "[C;$(C(=O)[#6]):1][O:2]>>*[C;+0:1]|[O;!R;$(O(C(=O)[#6])[CX4,c]):1][C;$(C(=O)):2]>>*[O;+0:1]"
        LabelList = "*C->C:10,*[13C]->13C:10,*[13CH]->13C:10|*O->O:20"
        finalSynthon = NormalSynthonsGenerator(LabelList, ReactionLIST,
                                set(), "Esters_Esters", Chem.MolFromSmiles(initSmi),
                                func=1)
    if "Ketones_Ketones" in Classes:
        newSynthToAdd = {}
        for synthon in finalSynthon:
            if "Ketones_Ketones" in finalSynthon[synthon]:
                synthMol = readMol(synthon)
                newClasses = BBClassifier(SMARTSLib, mol=synthMol)
                for cls in newClasses:
                    if "Alcohols" in cls:
                        nAzoles, nFinalSynthon = mainSynthonsGenerator([cls], synthon, keepPG, BBmarks)
                        for newSynth in nFinalSynthon:
                            if newSynth not in finalSynthon and newSynth not in newSynthToAdd:
                                newSynthToAdd[newSynth] = set()
                                newSynthToAdd[newSynth].update(finalSynthon[synthon])
                                newSynthToAdd[newSynth].update(nFinalSynthon[newSynth])
        if newSynthToAdd:
            for newSynth in newSynthToAdd:
                finalSynthon[newSynth] = newSynthToAdd[newSynth]
    return azoles, finalSynthon

def classesAssignement(CurrentClass, PreviousClasses, molSmi, MarksSetup, keepSynthonsWithPG=True):
    additionalBifuncClasses = ["Aminoacids_N-AliphaticAmino_Acid", "Aminoacids_N-AromaticAmino_Acid", "Reagents_DiAmines"]

    PGBifunctional = ["Bifunctional_Acid_Ester","Bifunctional_Acid_Nitro", "Bifunctional_Aldehyde_Ester", "Bifunctional_Amine_Ester",
                      "Bifunctional_Ester_Isocyanates", "Bifunctional_Ester_SO2X", "Bifunctional_Aldehyde_Nitro",

                      "Bifunctional_NbocAmino_Acid", "Bifunctional_NcbzAmino_Acid", "Bifunctional_Isothiocyanates_Acid", "Bifunctional_NfmocAmino_Acid",
                      "Bifunctional_Aldehyde_Nboc", "Bifunctional_NTFAcAmino_Acid",

                      "Bifunctional_Boronics_Ncbz", "Bifunctional_Boronics_Nfmoc",

                      "Bifunctional_NbnDi_Amines", "Bifunctional_NbocDi_Amines", 'Bifunctional_NcbzDi_Amines', "Bifunctional_NfmocDi_Amines",
                      "Bifunctional_NTFAcDi_Amines", 'Bifunctional_Di_Amines_NotherCarbamates',

                      'Trifunctional_Acid_Aldehyde_Nitro', "Trifunctional_Acid_ArylHalide_Ester",
                      "Trifunctional_Acid_ArylHalide_Nitro", 'Trifunctional_Amines_ArylHalide_Nitro',
                      "Trifunctional_NbocAmino_Acid_AlkyneCH", "Trifunctional_NbocAmino_Acid_ArylHalide",
                      "Trifunctional_NfmocAmino_Acid_AlkyneCH", "Trifunctional_NfmocAmino_Acid_ArylHalide"]

    FirstReactionAsPreparation = ["Bifunctional_Acid_Aldehyde", "Bifunctional_Aldehyde_ArylHalide",
                             "Bifunctional_Aldehyde_SO2X", "Bifunctional_Boronics_Acid",
                             "Bifunctional_Boronics_Aldehyde", "Bifunctional_Hydroxy_Aldehyde",

                             'Trifunctional_Acid_Aldehyde_ArylHalide', "Trifunctional_Acid_Aldehyde_Acetylenes",
                                  'Trifunctional_Acid_Aldehyde_Nitro', 'Trifunctional_Amines_ArylHalide_Nitro',
                                  "Trifunctional_NbocAmino_Acid_AlkyneCH", "Trifunctional_NfmocAmino_Acid_AlkyneCH",
                                  "Trifunctional_Di_Esters_Amino"]

    PolymerReagents = ["Reagents_PoliOxiranes", "Esters_PoliEsters", "Reagents_PoliIsocyanates", "SulfonylHalides_Poli_Sulfonylhalides"]

    trifuncClassesWithTwoPGs = ['Trifunctional_Acid_Ester_Nitro', "Trifunctional_NbocAmino_Acid_Ester",
                                    "Trifunctional_NbocAmino_Acid_Nitro", "Trifunctional_Amines_Nboc_Ester",
                                "Trifunctional_Nboc_NCbz_Amino_Acid", "Trifunctional_Nboc_Nfmoc_Amino_Acid",
                                "Trifunctional_NfmocAmino_Acid_Ester", "Trifunctional_NfmocAmino_Acid_Nitro",
                                "Trifunctional_Di_Esters_Amino"]

    if CurrentClass in trifuncClassesWithTwoPGs or "Trifunctional" in CurrentClass:
        func = 3
    elif "Bifunctional" in CurrentClass or CurrentClass in additionalBifuncClasses:
        func = 2
    else:
        func = 1
    if CurrentClass in FirstReactionAsPreparation:
        firstReactionAsPrep = True
    else:
        firstReactionAsPrep=False
    if CurrentClass in trifuncClassesWithTwoPGs:
        twoPGs=True
    else:
        twoPGs=False
    labledSynthons = {}
    if CurrentClass in PolymerReagents:
        MolsToWorkWith = {molSmi: PreviousClasses}
        for i in range(len(MarksSetup[CurrentClass]['Labels'].split("|"))):
            synthons = SynthonsGeneratorsForPolymerReagents(MarksSetup[CurrentClass]['Labels'].split("|")[i],
                                    MarksSetup[CurrentClass]['SMARTS'].split("|")[i], CurrentClass, MolsToWorkWith)
            if synthons:
                for synth in synthons:
                    if synth not in labledSynthons:
                        labledSynthons[synth] = synthons[synth].copy()
                    else:
                        labledSynthons[synth].update(synthons[synth])
        return labledSynthons
    elif CurrentClass in PGBifunctional or twoPGs:
        print("*********************************************************************")
        print(molSmi)
        print("*********************************************************************")
        synthons = ProtectiveGroupRemoval(MarksSetup[CurrentClass]['Labels'], MarksSetup[CurrentClass]['SMARTS'],
                                          Chem.MolFromSmiles(molSmi), keepSynthonsWithPG,
                                          firstReactionAsPrep, func, PreviousClasses, CurrentClass, twoPGs)

    elif firstReactionAsPrep:
        synthons = FirstReactionAsPrep(MarksSetup[CurrentClass]['Labels'], MarksSetup[CurrentClass]['SMARTS'],
                                        PreviousClasses, CurrentClass,Chem.MolFromSmiles(molSmi), func)
    else:
        synthons = NormalSynthonsGenerator(MarksSetup[CurrentClass]['Labels'], MarksSetup[CurrentClass]['SMARTS'],
                                           PreviousClasses, CurrentClass, Chem.MolFromSmiles(molSmi),
                                                 func=func)
    if synthons:
        for synth in synthons:
            if synth not in labledSynthons:
                labledSynthons[synth] = synthons[synth].copy()
            else:
                labledSynthons[synth].update(synthons[synth])

    return labledSynthons

def azolesSynthonPostGeneration(labledSynthons):
    pat = re.compile("\[\w*:\w*\]")
    Class = "nHAzoles_nHAzoles"
    additionalSynthons = {}
    maxMark = 0
    for molSmiles in labledSynthons:
        marksPrevious = [molSmiles[m.start():m.start() + 2] + molSmiles[m.end() - 4:m.end()] for m in
                         re.finditer(pat, molSmiles)]
        if len(marksPrevious)>maxMark:
            maxMark=len(marksPrevious)
    for molSmiles in labledSynthons:
        query = Chem.MolFromSmarts("[nHr5;!$(nc=O)]")
        mol = Chem.MolFromSmiles(molSmiles)
        marksPrevious = [molSmiles[m.start():m.start() + 2] + molSmiles[m.end() - 4:m.end()] for m in re.finditer(pat, molSmiles)]
        if len(marksPrevious)==maxMark and mol.HasSubstructMatch(query):
            cuttingRule = Reactions.ReactionFromSmarts("[nH;r5:1]>>*[n:1]")
            label = "*n->n:20"
            products = cuttingRule.RunReactants((mol,))
            for productSet in products:
                for product in productSet:
                    labledSynthon = __getLabledSmiles(product, label)
                    if marksPrevious:
                        for synth in labledSynthon:
                            marksNew = [synth[m.start():m.start() + 2] + synth[m.end() - 4:m.end()] for m in
                                     re.finditer(pat, synth)]
                            if len(marksNew)>len(marksPrevious):
                                if synth not in additionalSynthons:
                                    additionalSynthons[synth] = labledSynthons[molSmiles].copy()
                                else:
                                    additionalSynthons[synth].update(labledSynthons[molSmiles])
                                additionalSynthons[synth].add(Class)
                    else:
                        for synth in labledSynthon:
                            if synth not in additionalSynthons:
                                additionalSynthons[synth] = labledSynthons[molSmiles].copy()
                            else:
                                additionalSynthons[synth].update(labledSynthons[molSmiles])
                            additionalSynthons[synth].add(Class)
    return additionalSynthons

def generateBiacideSynthonForTrifunctional(labledSynthons, Class):
    additionalSynthons = {}
    for molSmiles in labledSynthons:
        query = Chem.MolFromSmarts(
            "[O;$(O=C([#6])[OD1])].[O;$(O([CH3])C([#6])=O),$(O([CH2][CH3])C([#6])=O),$(O([CH2]c1[cH][cH][cH][cH][cH]1)C([#6])=O),$(O(C([CH3])([CH3])[CH3])C([#6])=O),$(O([CH2][CH]=[CH2])C([#6])=O)]")
        mol = Chem.MolFromSmiles(molSmiles)
        pat = re.compile("\[\w*:\w*\]")
        marksPrevious = [molSmiles[m.start():m.start() + 2] + molSmiles[m.end() - 4:m.end()] for m in re.finditer(pat, molSmiles)]
        if mol.HasSubstructMatch(query):
            cuttingRule = Reactions.ReactionFromSmarts(
                "[O;$(O(C)C([#6])=O):1][C;$([CH3]),$([CH2][CH3]),$([CH2]c1[cH][cH][cH][cH][cH]1),$(C([CH3])([CH3])[CH3]),$([CH2][CH]=[CH2]):2]>>[OH:1]")
            label = "No"
            products = cuttingRule.RunReactants((mol,))
            for productSet in products:
                for product in productSet:
                    labledSynthon = __getLabledSmiles(product, label)
                    if marksPrevious:
                        for synth in labledSynthon:
                            marksNew = [synth[m.start():m.start() + 2] + synth[m.end() - 4:m.end()] for m in
                                     re.finditer(pat, synth)]
                            if len(marksNew)>len(marksPrevious):
                                if synth not in additionalSynthons:
                                    additionalSynthons[synth] = labledSynthons[molSmiles].copy()
                                else:
                                    additionalSynthons[synth].update(labledSynthons[molSmiles])
                                additionalSynthons[synth].add(Class)
                    else:
                        for synth in labledSynthon:
                            if synth not in additionalSynthons:
                                additionalSynthons[synth] = labledSynthons[molSmiles].copy()
                            else:
                                additionalSynthons[synth].update(labledSynthons[molSmiles])
                            additionalSynthons[synth].add(Class)
    return additionalSynthons

def SynthonsGeneratorsForPolymerReagents(Label, rule, Class, MolsToWorkWith, finalSynthons=None, firstLaunch=True, Deprotection=False):
    if finalSynthons==None:
        finalSynthons = {}
    newMolsToWorkWith = {}
    cuttingRule = Reactions.ReactionFromSmarts(rule)
    for molSmiles in MolsToWorkWith:
        pat = re.compile("\[\w*:\w*\]")
        marksPrevious = [molSmiles[m.start():m.start() + 2] + molSmiles[m.end() - 4:m.end()] for m in re.finditer(pat, molSmiles)]
        products = cuttingRule.RunReactants((Chem.MolFromSmiles(molSmiles),))
        if not products and not firstLaunch:
            if molSmiles not in finalSynthons:
                finalSynthons[molSmiles] = MolsToWorkWith[molSmiles].copy()
            else:
                finalSynthons[molSmiles].update(MolsToWorkWith[molSmiles])
            continue
        for productSet in products:
            for product in productSet:
                labledSynthons = __getLabledSmiles(product, Label)
                if marksPrevious and not Deprotection:
                    for synth in labledSynthons:
                        marksNew = [synth[m.start():m.start() + 2] + synth[m.end() - 4:m.end()] for m in
                                        re.finditer(pat, synth)]
                        if len(marksNew) > len(marksPrevious):
                            if synth not in newMolsToWorkWith:
                                newMolsToWorkWith[synth] = MolsToWorkWith[molSmiles].copy()
                            else:
                                newMolsToWorkWith[synth].update(MolsToWorkWith[molSmiles])
                                newMolsToWorkWith[synth].add(Class)
                else:
                    for synth in labledSynthons:
                        if synth not in newMolsToWorkWith:
                            newMolsToWorkWith[synth] = MolsToWorkWith[molSmiles].copy()
                        else:
                            newMolsToWorkWith[synth].update(MolsToWorkWith[molSmiles])
                            newMolsToWorkWith[synth].add(Class)
    if newMolsToWorkWith:
        print("#############################################")
        print(newMolsToWorkWith)
        SynthonsGeneratorsForPolymerReagents(Label, rule, Class, newMolsToWorkWith, finalSynthons, firstLaunch=False,
                                             Deprotection=Deprotection)
    if firstLaunch:
        return finalSynthons

def ProtectiveGroupRemoval(LabelsLIST, ReactionLIST, mol, keepSynthonsWithPG, firstReactionAsPrep, func, PreviousClasses,
                           CurrentClass, twoPGs=False):
    LabelsLISTBeforePGRemoval = LabelsLIST.split("|No|")[0]
    ReactionLISTBeforePGRemoval = ReactionLIST.split("|")[:len(LabelsLISTBeforePGRemoval.split("|"))]
    firstStopInd = len(LabelsLISTBeforePGRemoval.split("|"))
    finalSynthons = {}
    if firstReactionAsPrep:
        synthonsBeforeFirstPGremoval = FirstReactionAsPrep(LabelsLISTBeforePGRemoval, "|".join(ReactionLISTBeforePGRemoval),
                                                           PreviousClasses, CurrentClass, mol, func)

    else:
        if "Ester" in CurrentClass and "Acid" in CurrentClass and func==3 and not twoPGs:
            synthonsBeforeFirstPGremoval = NormalSynthonsGenerator(LabelsLISTBeforePGRemoval,
                                                                   "|".join(ReactionLISTBeforePGRemoval),
                                                                PreviousClasses, CurrentClass, mol, func=3)
        elif (func==3 and twoPGs) or func==2:

            synthonsBeforeFirstPGremoval = NormalSynthonsGenerator(LabelsLISTBeforePGRemoval, "|".join(ReactionLISTBeforePGRemoval),
                                                                   PreviousClasses, CurrentClass, mol, func=1)

        elif func==3:
            synthonsBeforeFirstPGremoval = NormalSynthonsGenerator(LabelsLISTBeforePGRemoval,
                                                                   "|".join(ReactionLISTBeforePGRemoval),
                                                                   PreviousClasses, CurrentClass, mol, func=2)

    if CurrentClass == "Trifunctional_Di_Esters_Amino":

        SynthonsWithoutPG = SynthonsGeneratorsForPolymerReagents(LabelsLIST.split("|")[2],
                                            ReactionLIST.split("|")[2],  CurrentClass,
                                             synthonsBeforeFirstPGremoval, Deprotection=True)

        if SynthonsWithoutPG:
            for synth in SynthonsWithoutPG:
                if synth not in finalSynthons:
                    finalSynthons[synth] = SynthonsWithoutPG[synth].copy()
                else:
                    finalSynthons[synth].update(SynthonsWithoutPG[synth])
        lastSynthons = SynthonsGeneratorsForPolymerReagents(LabelsLIST.split("|")[3],
                                                        ReactionLIST.split("|")[3],
                                                                        CurrentClass, SynthonsWithoutPG)
        if lastSynthons:
            for synth in lastSynthons:
                if synth not in finalSynthons:
                    finalSynthons[synth] = lastSynthons[synth].copy()
                else:
                    finalSynthons[synth].update(lastSynthons[synth])
        return finalSynthons
    PGremovalRule = ReactionLIST.split("|")[firstStopInd]
    if keepSynthonsWithPG or PGremovalRule =="[N;+0,+1;$([N+](=O)([#6])[O-]),$(N(=O)([#6])=O):1](=[O:2])=,-[O;+0,-1:3]>>[NH2,+0:1]":
        for synth in synthonsBeforeFirstPGremoval:
            if synth not in finalSynthons:
                finalSynthons[synth] = synthonsBeforeFirstPGremoval[synth].copy()
            else:
                finalSynthons[synth].update(synthonsBeforeFirstPGremoval[synth])
    cuttingRule = Reactions.ReactionFromSmarts(PGremovalRule)
    PGlable = LabelsLIST.split("|")[firstStopInd]
    SynthonsWithoutPG = {}
    for smi in synthonsBeforeFirstPGremoval:
        if func==3 and not twoPGs and smi.count(":")<2 and "Ester" not in CurrentClass and "AlkyneCH" not in CurrentClass: #this force algorythm to use for the PG removal step only synthons where all unprotected functional groups has been transformed into synthons
            continue
        products = cuttingRule.RunReactants((Chem.MolFromSmiles(smi),))
        for productSet in products:
            for product in productSet:
                labledSynthon = __getLabledSmiles(product, PGlable)
                for synth in labledSynthon:
                    if synth not in SynthonsWithoutPG:
                        SynthonsWithoutPG[synth] = set()
                    SynthonsWithoutPG[synth].update(PreviousClasses)
                    SynthonsWithoutPG[synth].add(CurrentClass)
    LabelsLISTBetweenPGRemoval = LabelsLIST.split("|No|")[1]
    ReactionLISTBetweenPGRemoval = ReactionLIST.split("|")[len(LabelsLISTBeforePGRemoval.split("|")) + 1:len(
        LabelsLISTBeforePGRemoval.split("|")) + len(LabelsLISTBetweenPGRemoval.split("|")) + 1]
    synthonsBetweenPGremoval = SynthonsWithoutPG.copy()
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
    print(SynthonsWithoutPG)
    for newSynthon in SynthonsWithoutPG:
        newSynthonsBetweenPGremoval = NormalSynthonsGenerator(LabelsLISTBetweenPGRemoval,
                                                              "|".join(ReactionLISTBetweenPGRemoval),
                                                              PreviousClasses, CurrentClass, Chem.MolFromSmiles(newSynthon), func=1)
        if newSynthonsBetweenPGremoval:
            for synth in newSynthonsBetweenPGremoval:
                if synth not in synthonsBetweenPGremoval:
                    synthonsBetweenPGremoval[synth] = newSynthonsBetweenPGremoval[synth].copy()
                else:
                    synthonsBetweenPGremoval[synth].update(newSynthonsBetweenPGremoval[synth])
    if len(LabelsLIST.split("|No|"))==3:
        secondStopInd = len(LabelsLISTBeforePGRemoval.split("|")) + len(LabelsLISTBetweenPGRemoval.split("|")) + 1
        PGremovalRule = ReactionLIST.split("|")[secondStopInd]
    if len(LabelsLIST.split("|No|"))==2 or keepSynthonsWithPG or \
            PGremovalRule =="[N;+0,+1;$([N+](=O)([#6])[O-]),$(N(=O)([#6])=O):1](=[O:2])=,-[O;+0,-1:3]>>[NH2,+0:1]":
        for synth in SynthonsWithoutPG:
            if synth not in finalSynthons:
                finalSynthons[synth] = SynthonsWithoutPG[synth].copy()
            else:
                finalSynthons[synth].update(SynthonsWithoutPG[synth])
        for synth in synthonsBetweenPGremoval:
            if synth not in finalSynthons:
                finalSynthons[synth] = synthonsBetweenPGremoval[synth].copy()
            else:
                finalSynthons[synth].update(synthonsBetweenPGremoval[synth])
    if len(LabelsLIST.split("|No|"))==3:
        cuttingRule = Reactions.ReactionFromSmarts(PGremovalRule)
        PGlable = LabelsLIST.split("|")[secondStopInd]
        SynthonsWithout2PG = {}
        for mol in synthonsBetweenPGremoval:
            products = cuttingRule.RunReactants((Chem.MolFromSmiles(mol),))
            for productSet in products:
                for product in productSet:
                    labledSynthon = __getLabledSmiles(product, PGlable)
                    if labledSynthon:
                        for synth in labledSynthon:
                            if synth not in SynthonsWithout2PG:
                                SynthonsWithout2PG[synth] = set()
                            SynthonsWithout2PG[synth].update(PreviousClasses)
                            SynthonsWithout2PG[synth].add(CurrentClass)

        for synth in SynthonsWithout2PG:
            if synth not in finalSynthons:
                finalSynthons[synth] = SynthonsWithout2PG[synth].copy()
            else:
                finalSynthons[synth].update(SynthonsWithout2PG[synth])
        LabelsLast = LabelsLIST.split("|No|")[2]
        ReactionLast = ReactionLIST.split("|")[len(LabelsLISTBeforePGRemoval.split("|")) + 1 + len(LabelsLISTBetweenPGRemoval.split("|"))+1:]
        for newSynthon in SynthonsWithout2PG:
            lastSynthons = NormalSynthonsGenerator(LabelsLast,"|".join(ReactionLast),
                                                        PreviousClasses, CurrentClass, Chem.MolFromSmiles(newSynthon), func=1)
            #lastSynthons.extend(FirstReactionAsPrep(LabelsLast, "|".join(ReactionLast),Chem.MolFromSmiles(newSynthon), func=1, Class=Class)
            if lastSynthons:
                for synth in lastSynthons:
                    if synth not in finalSynthons:
                        finalSynthons[synth] = lastSynthons[synth].copy()
                    else:
                        finalSynthons[synth].update(lastSynthons[synth])
    return finalSynthons

def NormalSynthonsGenerator(LabelsLIST, ReactionLIST, PreviousClasses, CurrentClass, mol, func=1, usedInds = None):
    if usedInds == None:
        usedInds = []
    labledSynthons = {}
    pat = re.compile("\[\w*:\w*\]")
    molSmiles = Chem.MolToSmiles(mol, canonical=True)
    marksPrevious = [molSmiles[m.start():m.start() + 2] + molSmiles[m.end() - 4:m.end()] for m in re.finditer(pat, molSmiles)]
    for ind, rule in enumerate(ReactionLIST.split("|")):
        if ind not in usedInds:
            try:
                cuttingRule = Reactions.ReactionFromSmarts(rule)
            except:
                print("########################")
                print(Chem.MolToSmiles(mol, canonical=True))
                print(rule)
                exit()
            products = cuttingRule.RunReactants((mol,))
            Label = LabelsLIST.split("|")[ind]
            for productSet in products:
                for product in productSet:
                    labledSynthon = __getLabledSmiles(product, Label)
                    if labledSynthon==None:
                        print(Chem.MolToSmiles(mol, canonical=True))
                        exit()
                    if marksPrevious:
                        for synth in labledSynthon:
                            marksNew = [synth[m.start():m.start() + 2] + synth[m.end() - 4:m.end()] for m in
                                     re.finditer(pat, synth)]
                            if len(marksNew)>len(marksPrevious):
                                if synth not in labledSynthons:
                                    labledSynthons[synth] = set()
                                labledSynthons[synth].update(PreviousClasses)
                                labledSynthons[synth].add(CurrentClass)
                    else:
                        for synth in labledSynthon:
                            if synth not in labledSynthons:
                                labledSynthons[synth] = set()
                            labledSynthons[synth].update(PreviousClasses)
                            labledSynthons[synth].add(CurrentClass)
                    newSynthons = None
                    if func == 2:
                        newMol = Chem.MolFromSmiles(labledSynthon[0])
                        usedInds.append(ind)
                        newSynthons = NormalSynthonsGenerator(LabelsLIST, ReactionLIST, PreviousClasses, CurrentClass, newMol, func=1,
                                                              usedInds=usedInds)
                    elif func == 3:
                        newMol = Chem.MolFromSmiles(labledSynthon[0])
                        usedInds.append(ind)
                        newSynthons = NormalSynthonsGenerator(LabelsLIST, ReactionLIST, PreviousClasses, CurrentClass, newMol, func=2,
                                                              usedInds=usedInds)
                    if newSynthons:
                        for synth in newSynthons:
                            if synth not in labledSynthons:
                                labledSynthons[synth] = newSynthons[synth].copy()
                            else:
                                labledSynthons[synth].update(newSynthons[synth])
    return labledSynthons

def FirstReactionAsPrep(LabelsLIST, ReactionLIST, PreviousClasses, CurrentClass, mol, func):
    rule = ReactionLIST.split("|")[0]
    cuttingRule = Reactions.ReactionFromSmarts(rule)
    products = cuttingRule.RunReactants((mol,))
    if not products:
        if len(LabelsLIST.split("|"))>1:
            return FirstReactionAsPrep("|".join(LabelsLIST.split("|")[1:]), "|".join(ReactionLIST.split("|")[1:]),
                                       PreviousClasses, CurrentClass, mol, func)
        else:
            return {Chem.MolToSmiles(mol, canonical=True): PreviousClasses}
    else:
        Label = LabelsLIST.split("|")[0]
        synthonsAsInpForTheNextStep = {}
        for productSet in products:
            for product in productSet:
                labledSynthon = __getLabledSmiles(product, Label)
                if labledSynthon:
                    for synth in labledSynthon:
                        if synth not in synthonsAsInpForTheNextStep:
                            synthonsAsInpForTheNextStep[synth] = set()
                        synthonsAsInpForTheNextStep[synth].update(PreviousClasses)
                        synthonsAsInpForTheNextStep[synth].add(CurrentClass)

        if len(ReactionLIST.split("|"))==1:
            return synthonsAsInpForTheNextStep
        if not synthonsAsInpForTheNextStep and "Boronics" in CurrentClass:
            rule = ReactionLIST.split("|")[1]
            cuttingRule = Reactions.ReactionFromSmarts(rule)
            products = cuttingRule.RunReactants((mol,))
            Label = LabelsLIST.split("|")[1]
            for productSet in products:
                for product in productSet:
                    labledSynthon = __getLabledSmiles(product, Label)
                    if labledSynthon:
                        for synth in labledSynthon:
                            if synth not in synthonsAsInpForTheNextStep:
                                synthonsAsInpForTheNextStep[synth] = set()
                            synthonsAsInpForTheNextStep[synth].update(PreviousClasses)
                            synthonsAsInpForTheNextStep[synth].add(CurrentClass)
            for synth in synthonsAsInpForTheNextStep:
                labledSynthons = NormalSynthonsGenerator("|".join(LabelsLIST.split("|")[1:]), "|".join(ReactionLIST.split("|")[1:]),
                                            PreviousClasses, CurrentClass, Chem.MolFromSmiles(synth), func=func - 1)
                if labledSynthons:
                    for synth in labledSynthons:
                        if synth not in synthonsAsInpForTheNextStep:
                            synthonsAsInpForTheNextStep[synth] = labledSynthons[synth].copy()
                        else:
                            synthonsAsInpForTheNextStep[synth].update(labledSynthons[synth])
            return synthonsAsInpForTheNextStep
        lastSynthons = {}
        for synth in synthonsAsInpForTheNextStep:
            labledSynthons = NormalSynthonsGenerator("|".join(LabelsLIST.split("|")[1:]), "|".join(ReactionLIST.split("|")[1:]),
                                                          PreviousClasses, CurrentClass, Chem.MolFromSmiles(synth), func=func-1)
            if labledSynthons:
                for synth in labledSynthons:
                    if synth not in lastSynthons:
                        lastSynthons[synth] = labledSynthons[synth].copy()
                    else:
                        lastSynthons[synth].update(labledSynthons[synth])
        for synth in lastSynthons:
            if synth not in synthonsAsInpForTheNextStep:
                synthonsAsInpForTheNextStep[synth] = lastSynthons[synth]
            else:
                synthonsAsInpForTheNextStep[synth].update(lastSynthons[synth])
        return synthonsAsInpForTheNextStep

def __getLabledSmiles(productMolecule: Chem.rdchem.Mol, Label:str):
    productSmiles = Chem.MolToSmiles(productMolecule, canonical=True)
    labeledSmilesList = []
    if Label != "No":
        for sublabel in Label.split(","):
            if productSmiles.find(sublabel.split("->")[0]) != -1:
                labeledSmiles = __checkLable(productSmiles, sublabel)
                if labeledSmiles==None:
                    return None
                if "*" in labeledSmiles:
                    productSmiles = labeledSmiles
                    continue
                elif labeledSmiles:
                    labeledSmilesList.append(labeledSmiles)
        if labeledSmilesList:
            return list(set(labeledSmilesList))
    print("WARNING! No lable was assigned to the smiles: " + productSmiles)
    return [productSmiles]

def __generateMajorTautFromSynthonSmiles(initSmiles):
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

def __checkLable(productSmiles:str, Label:str):
    goodValenceSmiles = None
    if Label.split("->")[0][1] == "S":
        hCount = 1
        out = productSmiles.replace(Label.split("->")[0],
                                        "[" + Label.split("->")[1].split(":")[0] + "H" + str(hCount) + ":" +
                                        Label.split("->")[1].split(":")[1] + "]")
        goodValenceSmiles = out
    else:
        for hCount in range(1, 5):
            if "+" in Label:
                out = productSmiles.replace(Label.split("->")[0],
                                                "[" + Label.split("->")[1].split(":")[0] + "H" + str(hCount) + "+:" +
                                                Label.split("->")[1].split(":")[1] + "]")
            else:
                out = productSmiles.replace(Label.split("->")[0],
                                                "[" + Label.split("->")[1].split(":")[0] + "H" + str(hCount) + ":" +
                                                Label.split("->")[1].split(":")[1] + "]")

            newMol = Chem.MolFromSmiles(out)
            if not newMol:
                if "[nH1:" in out:
                    modifiedSmiles = out.replace("[nH1:", "[n:", 1)
                    check2 = __CheckMolStructure(modifiedSmiles, Label)
                    if check2:
                        break
                    else:
                        goodValenceSmiles = modifiedSmiles
                else:
                    break
            else:
                goodValenceSmiles = out
                check = __CheckMolStructure(goodValenceSmiles, Label)
                if check:
                    break
    if not goodValenceSmiles:
        print("Problem with structure check: " + productSmiles + " " + out)
        return None
    return __generateMajorTautFromSynthonSmiles(goodValenceSmiles)

def __CheckMolStructure(goodValenceSmiles, label):
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

def __getReactionSMARTS(BB_Marks: ET.Element):
    MarksSetup = {}
    for child in BB_Marks:
        for subCh in child:
            if subCh.get('SMARTS'):
                MarksSetup[child.tag + "_" + subCh.tag] = {}
                MarksSetup[child.tag + "_" + subCh.tag]["SMARTS"] = subCh.get('SMARTS')
                MarksSetup[child.tag + "_" + subCh.tag]["Labels"] = subCh.get('Labels')
    return MarksSetup

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

def Ro2Filtration(synthonSmiles):
    mol = Chem.MolFromSmiles(synthonSmiles)
    mol = AddHs(mol)
    MolW = ExactMolWt(mol)
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
        return False, "MolW=" + str(MolW), "LogP=" + str(LogP), "HDC=" + str(HDC), "HAC=" + str(HAC)
    else:
        return True, "MolW=" + str(MolW), "LogP=" + str(LogP), "HDC=" + str(HDC), "HAC=" + str(HAC)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Synthons Identification - BBs module - synthons generation from BBs",
                                     epilog="Yuliana Zabolotna 2020",
                                     prog="SynthI-BBs")
    parser.add_argument("-sL", "--SMARTSLib", type=str, default="SMARTSLibNew.json", help="SMARTSLibrary")
    parser.add_argument("-sm", "--smiles", type=str, help="input Molecule")
    parser.add_argument("-bbM", "--bbM", type=str, default="BB_Marks.xml", help="BB marks")
    parser.add_argument("-o", "--output", type=str, help="output file")
    parser.add_argument("--keepPG", action="store_true", help="The script will keep both protected and unprotected "
                                            "synthons (concerns Boc, Bn, Fmoc, Cbz and Esters protections).")

    args = parser.parse_args()
    main(args)
