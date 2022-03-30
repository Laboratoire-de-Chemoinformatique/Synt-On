import re,os,sys
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem.Scaffolds.MurckoScaffold import *
srcPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(1, srcPath)
from SyntOn_Classifier import BBClassifier
from UsefulFunctions import *


def mainSynthonsGenerator(initSmi, keepPG=False, Classes=None, returnDict=False, returnBoolAndDict=False):
    solventsToIgnore = ["OC(=O)C(=O)O", "CC(=O)O", "OS(=O)(=O)O", "[O-]Cl(=O)(=O)=O", "OP(=O)(O)O", "OC(=O)C(F)(F)F",
                        "OS(=O)(=O)C(F)(F)F", "OC(=O)O", "[O-]S(=O)(=O)C(F)(F)F", "OC=O", "OC(=O)/C=C\C(=O)O", "[O-]C(=O)C(F)(F)F",
                        "OC(=O)/C=C/C(=O)O"]
    canonicalSmilesOfSolventsToIgnore = set([Chem.MolToSmiles(Chem.MolFromSmiles(x), canonical=True) for x in solventsToIgnore])
    initMol = readMol(initSmi)
    query = Chem.MolFromSmarts(
        "[#6]-[#6]-[#8]-[#6].[#6]-[#8]-[#6](-[#6])=O.[#6]-[#8]-[#6](-[#6])=O.[#6]-[#8]-[#6](-[#6])=O")
    if initMol == None or initMol.HasSubstructMatch(query):
        finalSynthon = {}
        azoles = False
    elif len(initSmi.split(".")) > 1: # case of input mixtures
        finalSynthon = {}
        azoles = False
        for smi in initSmi.split("."):
            mol = readMol(smi)
            if Chem.MolToSmiles(mol, canonical=True) not in canonicalSmilesOfSolventsToIgnore:
                nAzoles, nFinalSynthon = mainSynthonsGenerator(smi, keepPG, returnBoolAndDict=True)
                if nAzoles:
                    azoles = True
                if nFinalSynthon:
                    for newSynth in nFinalSynthon:
                        if newSynth not in finalSynthon:
                            finalSynthon[newSynth] = nFinalSynthon[newSynth].copy()
    else:
        if Classes == None:
            AllClasses = BBClassifier(mol=initMol)
            Classes = [clas for clas in AllClasses if "MedChemHighlights" not in clas and "DEL" not in clas]
        BBmarks = os.path.join(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0], "config" , "BB_Marks.xml")
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
                    synthons = __synthonsAssignement(Classes[ind], molsToWorkWith[mol], mol, MarksSetup, keepSynthonsWithPG)
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
                for mol in molsToWorkWith:
                    synthons = __synthonsAssignement(Classes[ind], molsToWorkWith[mol], mol, MarksSetup, keepSynthonsWithPG)
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
                    synthons = __synthonsAssignement(Classes[ind], extraMols[mol], mol, MarksSetup, keepSynthonsWithPG)
                    if synthons:
                        for synth in synthons:
                            if synth not in finalSynthon:
                                finalSynthon[synth] = synthons[synth].copy()"""

        for clas in Classes:
            if "Trifunctional" in clas and "Ester" in clas and "Acid" in clas:
                additionalSynthons = __generateBiacideSynthonForTrifunctional(finalSynthon, clas)
                if additionalSynthons:
                    for synth in additionalSynthons:
                        if synth not in finalSynthon:
                            finalSynthon[synth] = additionalSynthons[synth].copy()
                        else:
                            finalSynthon[synth].update(additionalSynthons[synth])

        additionalSynthons = __azolesSynthonPostGeneration(finalSynthon)
        if additionalSynthons:
            azoles = True
            for synth in additionalSynthons:
                if synth not in finalSynthon:
                    finalSynthon[synth] = additionalSynthons[synth].copy()
                else:
                    finalSynthon[synth].update(additionalSynthons[synth])
        else:
            azoles = False
        if not finalSynthon and "Esters_Esters" in Classes:
            ReactionLIST = "[C;$(C(=O)[#6]):1][O:2]>>*[C;+0:1]|[O;!R;$(O(C(=O)[#6])[CX4,c]):1][C;$(C(=O)):2]>>*[O;+0:1]"
            LabelList = "*C->C:10,*[13C]->13C:10,*[13CH]->13C:10|*O->O:20"
            finalSynthon = __NormalSynthonsGenerator(LabelList, ReactionLIST,
                                    set(), "Esters_Esters", initMol,
                                    func=1)
        if "Ketones_Ketones" in Classes:
            newSynthToAdd = {}
            for synthon in finalSynthon:
                if "Ketones_Ketones" in finalSynthon[synthon]:
                    synthMol = readMol(synthon)
                    newClasses = BBClassifier(mol=synthMol)
                    for cls in newClasses:
                        if "Alcohols" in cls:
                            nAzoles, nFinalSynthon = mainSynthonsGenerator( synthon, keepPG, [cls], returnBoolAndDict=True)
                            for newSynth in nFinalSynthon:
                                if newSynth not in finalSynthon and newSynth not in newSynthToAdd:
                                    newSynthToAdd[newSynth] = set()
                                    newSynthToAdd[newSynth].update(finalSynthon[synthon])
                                    newSynthToAdd[newSynth].update(nFinalSynthon[newSynth])
            if newSynthToAdd:
                for newSynth in newSynthToAdd:
                    finalSynthon[newSynth] = newSynthToAdd[newSynth]
    if returnDict:
        return finalSynthon
    elif returnBoolAndDict:
        return azoles, finalSynthon
    else:
        print("\n\n\n___________________________________________________________________________________________________")
        print("All generated synthons (" + str(len(finalSynthon)) +"): " + ".".join([x for x in finalSynthon]))
        print("Col1-Synton Col2-RespectiveBBsClass")
        for synth in finalSynthon:
            print(synth + " " + "+".join(finalSynthon[synth]))

def __synthonsAssignement(CurrentClass, PreviousClasses, molSmi, MarksSetup, keepSynthonsWithPG=True):
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

    __FirstReactionAsPreparation = ["Bifunctional_Acid_Aldehyde", "Bifunctional_Aldehyde_ArylHalide",
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
    if CurrentClass in __FirstReactionAsPreparation:
        firstReactionAsPrep = True
    else:
        firstReactionAsPrep=False
    if CurrentClass in trifuncClassesWithTwoPGs:
        twoPGs=True
    else:
        twoPGs=False
    labledSynthons = {}
    mol = readMol(molSmi)
    if CurrentClass in PolymerReagents:
        MolsToWorkWith = {Chem.MolToSmiles(mol, canonical=True): PreviousClasses}
        for i in range(len(MarksSetup[CurrentClass]['Labels'].split("|"))):
            synthons = __SynthonsGeneratorsForPolymerReagents(MarksSetup[CurrentClass]['Labels'].split("|")[i],
                                    MarksSetup[CurrentClass]['SMARTS'].split("|")[i], CurrentClass, MolsToWorkWith)
            if synthons:
                for synth in synthons:
                    if synth not in labledSynthons:
                        labledSynthons[synth] = synthons[synth].copy()
                    else:
                        labledSynthons[synth].update(synthons[synth])
        return labledSynthons
    elif CurrentClass in PGBifunctional or twoPGs:
        synthons = __ProtectiveGroupRemoval(MarksSetup[CurrentClass]['Labels'], MarksSetup[CurrentClass]['SMARTS'],
                                          mol, keepSynthonsWithPG,
                                          firstReactionAsPrep, func, PreviousClasses, CurrentClass, twoPGs)

    elif firstReactionAsPrep:
        synthons = __FirstReactionAsPrep(MarksSetup[CurrentClass]['Labels'], MarksSetup[CurrentClass]['SMARTS'],
                                        PreviousClasses, CurrentClass, mol, func)
    else:
        synthons = __NormalSynthonsGenerator(MarksSetup[CurrentClass]['Labels'], MarksSetup[CurrentClass]['SMARTS'],
                                           PreviousClasses, CurrentClass, mol,
                                                 func=func)
    if synthons:
        for synth in synthons:
            if synth not in labledSynthons:
                labledSynthons[synth] = synthons[synth].copy()
            else:
                labledSynthons[synth].update(synthons[synth])

    return labledSynthons

def __azolesSynthonPostGeneration(labledSynthons):
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
        mol = readMol(molSmiles)
        marksPrevious = [molSmiles[m.start():m.start() + 2] + molSmiles[m.end() - 4:m.end()] for m in re.finditer(pat, molSmiles)]
        if len(marksPrevious)==maxMark and mol.HasSubstructMatch(query):
            cuttingRule = Reactions.ReactionFromSmarts("[nH;r5:1]>>*[n:1]")
            label = "*n->n:20"
            products = cuttingRule.RunReactants((mol,))
            for productSet in products:
                for product in productSet:
                    labledSynthon = __getBBLabledSmiles(product, label)
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

def __generateBiacideSynthonForTrifunctional(labledSynthons, Class):
    additionalSynthons = {}
    for molSmiles in labledSynthons:
        query = Chem.MolFromSmarts(
            "[O;$(O=C([#6])[OD1])].[O;$(O([CH3])C([#6])=O),$(O([CH2][CH3])C([#6])=O),$(O([CH2]c1[cH][cH][cH][cH][cH]1)C([#6])=O),$(O(C([CH3])([CH3])[CH3])C([#6])=O),$(O([CH2][CH]=[CH2])C([#6])=O)]")
        mol = readMol(molSmiles)
        pat = re.compile("\[\w*:\w*\]")
        marksPrevious = [molSmiles[m.start():m.start() + 2] + molSmiles[m.end() - 4:m.end()] for m in re.finditer(pat, molSmiles)]
        if mol.HasSubstructMatch(query):
            cuttingRule = Reactions.ReactionFromSmarts(
                "[O;$(O(C)C([#6])=O):1][C;$([CH3]),$([CH2][CH3]),$([CH2]c1[cH][cH][cH][cH][cH]1),$(C([CH3])([CH3])[CH3]),$([CH2][CH]=[CH2]):2]>>[OH:1]")
            label = "No"
            products = cuttingRule.RunReactants((mol,))
            for productSet in products:
                for product in productSet:
                    labledSynthon = __getBBLabledSmiles(product, label)
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

def __SynthonsGeneratorsForPolymerReagents(Label, rule, Class, MolsToWorkWith, finalSynthons=None, firstLaunch=True, Deprotection=False):
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
                labledSynthons = __getBBLabledSmiles(product, Label)
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
        __SynthonsGeneratorsForPolymerReagents(Label, rule, Class, newMolsToWorkWith, finalSynthons, firstLaunch=False,
                                             Deprotection=Deprotection)
    if firstLaunch:
        return finalSynthons

def __ProtectiveGroupRemoval(LabelsLIST, ReactionLIST, mol, keepSynthonsWithPG, firstReactionAsPrep, func, PreviousClasses,
                           CurrentClass, twoPGs=False):
    LabelsLISTBeforePGRemoval = LabelsLIST.split("|No|")[0]
    ReactionLISTBeforePGRemoval = ReactionLIST.split("|")[:len(LabelsLISTBeforePGRemoval.split("|"))]
    firstStopInd = len(LabelsLISTBeforePGRemoval.split("|"))
    finalSynthons = {}
    if firstReactionAsPrep:
        synthonsBeforeFirstPGremoval = __FirstReactionAsPrep(LabelsLISTBeforePGRemoval, "|".join(ReactionLISTBeforePGRemoval),
                                                           PreviousClasses, CurrentClass, mol, func)

    else:
        if "Ester" in CurrentClass and "Acid" in CurrentClass and func==3 and not twoPGs:
            synthonsBeforeFirstPGremoval = __NormalSynthonsGenerator(LabelsLISTBeforePGRemoval,
                                                                   "|".join(ReactionLISTBeforePGRemoval),
                                                                PreviousClasses, CurrentClass, mol, func=3)
        elif (func==3 and twoPGs) or func==2:

            synthonsBeforeFirstPGremoval = __NormalSynthonsGenerator(LabelsLISTBeforePGRemoval, "|".join(ReactionLISTBeforePGRemoval),
                                                                   PreviousClasses, CurrentClass, mol, func=1)

        elif func==3:
            synthonsBeforeFirstPGremoval = __NormalSynthonsGenerator(LabelsLISTBeforePGRemoval,
                                                                   "|".join(ReactionLISTBeforePGRemoval),
                                                                   PreviousClasses, CurrentClass, mol, func=2)

    if CurrentClass == "Trifunctional_Di_Esters_Amino":

        SynthonsWithoutPG = __SynthonsGeneratorsForPolymerReagents(LabelsLIST.split("|")[2],
                                            ReactionLIST.split("|")[2],  CurrentClass,
                                             synthonsBeforeFirstPGremoval, Deprotection=True)

        if SynthonsWithoutPG:
            for synth in SynthonsWithoutPG:
                if synth not in finalSynthons:
                    finalSynthons[synth] = SynthonsWithoutPG[synth].copy()
                else:
                    finalSynthons[synth].update(SynthonsWithoutPG[synth])
        lastSynthons = __SynthonsGeneratorsForPolymerReagents(LabelsLIST.split("|")[3],
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
                labledSynthon = __getBBLabledSmiles(product, PGlable)
                for synth in labledSynthon:
                    if synth not in SynthonsWithoutPG:
                        SynthonsWithoutPG[synth] = set()
                    SynthonsWithoutPG[synth].update(PreviousClasses)
                    SynthonsWithoutPG[synth].add(CurrentClass)
    LabelsLISTBetweenPGRemoval = LabelsLIST.split("|No|")[1]
    ReactionLISTBetweenPGRemoval = ReactionLIST.split("|")[len(LabelsLISTBeforePGRemoval.split("|")) + 1:len(
        LabelsLISTBeforePGRemoval.split("|")) + len(LabelsLISTBetweenPGRemoval.split("|")) + 1]
    synthonsBetweenPGremoval = SynthonsWithoutPG.copy()
    for newSynthon in SynthonsWithoutPG:
        newSynthonsBetweenPGremoval = __NormalSynthonsGenerator(LabelsLISTBetweenPGRemoval,
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
                    labledSynthon = __getBBLabledSmiles(product, PGlable)
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
            lastSynthons = __NormalSynthonsGenerator(LabelsLast,"|".join(ReactionLast),
                                                        PreviousClasses, CurrentClass, Chem.MolFromSmiles(newSynthon), func=1)
            #lastSynthons.extend(__FirstReactionAsPrep(LabelsLast, "|".join(ReactionLast),Chem.MolFromSmiles(newSynthon), func=1, Class=Class)
            if lastSynthons:
                for synth in lastSynthons:
                    if synth not in finalSynthons:
                        finalSynthons[synth] = lastSynthons[synth].copy()
                    else:
                        finalSynthons[synth].update(lastSynthons[synth])
    return finalSynthons

def __NormalSynthonsGenerator(LabelsLIST, ReactionLIST, PreviousClasses, CurrentClass, mol, func=1, usedInds = None):
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
                    labledSynthon = __getBBLabledSmiles(product, Label)
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
                        newSynthons = __NormalSynthonsGenerator(LabelsLIST, ReactionLIST, PreviousClasses, CurrentClass, newMol, func=1,
                                                              usedInds=usedInds)
                    elif func == 3:
                        newMol = Chem.MolFromSmiles(labledSynthon[0])
                        usedInds.append(ind)
                        newSynthons = __NormalSynthonsGenerator(LabelsLIST, ReactionLIST, PreviousClasses, CurrentClass, newMol, func=2,
                                                              usedInds=usedInds)
                    if newSynthons:
                        for synth in newSynthons:
                            if synth not in labledSynthons:
                                labledSynthons[synth] = newSynthons[synth].copy()
                            else:
                                labledSynthons[synth].update(newSynthons[synth])
    return labledSynthons

def __FirstReactionAsPrep(LabelsLIST, ReactionLIST, PreviousClasses, CurrentClass, mol, func):
    rule = ReactionLIST.split("|")[0]
    cuttingRule = Reactions.ReactionFromSmarts(rule)
    products = cuttingRule.RunReactants((mol,))
    if not products:
        if len(LabelsLIST.split("|"))>1:
            return __FirstReactionAsPrep("|".join(LabelsLIST.split("|")[1:]), "|".join(ReactionLIST.split("|")[1:]),
                                       PreviousClasses, CurrentClass, mol, func)
        else:
            return {Chem.MolToSmiles(mol, canonical=True): PreviousClasses}
    else:
        Label = LabelsLIST.split("|")[0]
        synthonsAsInpForTheNextStep = {}
        for productSet in products:
            for product in productSet:
                labledSynthon = __getBBLabledSmiles(product, Label)
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
                    labledSynthon = __getBBLabledSmiles(product, Label)
                    if labledSynthon:
                        for synth in labledSynthon:
                            if synth not in synthonsAsInpForTheNextStep:
                                synthonsAsInpForTheNextStep[synth] = set()
                            synthonsAsInpForTheNextStep[synth].update(PreviousClasses)
                            synthonsAsInpForTheNextStep[synth].add(CurrentClass)
            for synth in synthonsAsInpForTheNextStep:
                labledSynthons = __NormalSynthonsGenerator("|".join(LabelsLIST.split("|")[1:]), "|".join(ReactionLIST.split("|")[1:]),
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
            labledSynthons = __NormalSynthonsGenerator("|".join(LabelsLIST.split("|")[1:]), "|".join(ReactionLIST.split("|")[1:]),
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

def __getBBLabledSmiles(productMolecule: Chem.rdchem.Mol, Label:str):
    productSmiles = Chem.MolToSmiles(productMolecule, canonical=True)
    labeledSmilesList = []
    if Label != "No":
        for sublabel in Label.split(","):
            if productSmiles.find(sublabel.split("->")[0]) != -1:
                labeledSmiles = checkLable(productSmiles, sublabel)
                if labeledSmiles==None:
                    return None
                if "*" in labeledSmiles:
                    productSmiles = labeledSmiles
                    continue
                elif labeledSmiles:
                    labeledSmilesList.append(labeledSmiles)
        if labeledSmilesList:
            return list(set(labeledSmilesList))
    #print("WARNING! No lable was assigned to the smiles: " + productSmiles)
    return [productSmiles]

def __getReactionSMARTS(BB_Marks: ET.Element):
    MarksSetup = {}
    for child in BB_Marks:
        for subCh in child:
            if subCh.get('SMARTS'):
                MarksSetup[child.tag + "_" + subCh.tag] = {}
                MarksSetup[child.tag + "_" + subCh.tag]["SMARTS"] = subCh.get('SMARTS')
                MarksSetup[child.tag + "_" + subCh.tag]["Labels"] = subCh.get('Labels')
    return MarksSetup

def generateScaffoldForBB(smiles, returnObjects=False):
    scaffold = None
    mol = None
    PGdict = {"NCbz": "[N:1][C;$(C(=O)O[CH2]c1[cH][cH][cH][cH][cH]1):2]>>[N:1]",
              "NFmoc": "[N:1][C;$(C(=O)O[CH2][CH]1c2[cH][cH][cH][cH]c2-c3[cH][cH][cH][cH]c13):2]>>[N:1]",
              "NBnz": "[N;+0;$(N[CH2]c1[cH][cH][cH][cH][cH]1);!$(N[C,S,P]=[O,S,N]):1][C;$([CH2]c1[cH][cH][cH][cH][cH]1):2]>>[N:1]",
              "COOBnz": "[O;$(O(C)C([#6])=O):1][C;$([CH2]c1[cH][cH][cH][cH][cH]1):2]>>[OH:1]",
              "Boronics": "[B;$(B(O@C)O@C):1][#6:2]>>[#6:2]",
              "Oxiranes": "[C:1]1[O:2][C:3]1>>[C:1]([OH:2])[C;+0:3]"}
    mol = readMol(smiles)
    if mol:
        for pg in PGdict:
            mol = __removePGforScaffolds(PGdict[pg], mol)
        scaffold = MurckoScaffoldSmiles(mol=mol)
    if returnObjects:
        return scaffold,mol
    else:
        return scaffold

def __removePGforScaffolds(reactionRule, mol):
    q = Chem.MolFromSmarts(reactionRule.split(">>")[0])
    cuttingRule = Reactions.ReactionFromSmarts(reactionRule)
    while mol.HasSubstructMatch(q):
        products = cuttingRule.RunReactants((mol,))
        mol = products[0][0]
        mol.UpdatePropertyCache()
        Chem.GetSymmSSSR(mol)
        mol.GetRingInfo().NumRings()
    return mol

"""def __checkLable(productSmiles:str, Label:str):
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
                    check2 = CheckMolStructure(modifiedSmiles, Label)
                    if check2:
                        break
                    else:
                        goodValenceSmiles = modifiedSmiles
                else:
                    break
            else:
                goodValenceSmiles = out
                check = CheckMolStructure(goodValenceSmiles, Label)
                if check:
                    break
    if not goodValenceSmiles:
        print("Problem with structure check: " + productSmiles + " " + out)
        return None
    return generateMajorTautFromSynthonSmiles(goodValenceSmiles)"""

