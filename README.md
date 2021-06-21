# SynthI
Open-source tool for synthons-based library design.

# Table of Contents
* [Prerequisites](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#prerequisites)
    * [System requirements](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#system-requirements)
    * [Comopounds preprocessing](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#comopounds-preprocessing)
* [SynthI-Classifier](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#synthi-classifier)
    * [Bulk BBs classification](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#bulk-bbs-classification) 
* [SynthI-BBs](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#synthi-bbs)
    * [Scaffold analysis](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#scaffold-analysis)
    * [BBs synthonization](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#bbs-synthonization)
    * [Bulk synthons generation for the large BBs library](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#bulk-synthons-generation-for-the-large-bbs-library)
* [SynthI-Fragmentation](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#synthi-fragmentation)
    * [Fragmentation with a default setup](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#fragmentation-with-a-default-setup)
    * [Bulk compounds fragmentation](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#bulk-compounds-fragmentation)
    * [Detailed classes description](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#detailed-classes-description)

## Prerequisites

### System requirements
SynthI is a suit of scripts written in python. It should be used with the following dependencies: python 3.9.0, rdkit 2021.03.1, matplotlib 3.4.2  and numpy 1.20.2. 

Several build-in python modules are also used, but they are usually installed by default (*datetime, os, time, random, re, resource, sys, multiprocessing, collections, xml*). All other modules are custom written and provided within the package. 

The scripts were run in a linux workstation with 15 processors.

### Comopounds preprocessing

All BBs structures need to be sanitized and standartized idependently by user prior to SynthI usage. Solvents and contriones should be deleted. There is no need to generate major tautomer form as soon as SynthI will do it for each generated synthons separately.

## SynthI-Classifier
[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

This module returns the list of classes assigned to the given BB. 

```python
>>> from SynthI.src.SynthI_Classifier import BBClassifier
>>> BBClassifier(molSmiles="CCOC(=O)C1=C(N)SC=C1C2CC2")
```
```text
['Bifunctional_Amine_Ester', 'PrimaryAmines_PriAmines_Het-Anilines']
```
If the SMILES cannot be processed by RdKit, the following messages will appear:

```python
>>> BBClassifier(molSmiles="C1CCC(CC1)P([Ir][N]2=CC=CC=C2)(C3CCCCC3)C4CCCCC4.C/1C/C=C\CC/C=C1")
```
```text
[23:37:46] Explicit valence for atom # 8 N, 4, is greater than permitted
[23:37:46] Explicit valence for atom # 8 N, 4, is greater than permitted
C1CCC(CC1)P([Ir][N]2=CC=CC=C2)(C3CCCCC3)C4CCCCC4.C/1C/C=C\CC/C=C1 was not processed by rdkit
```
As soon as heterocyclization reaction will be available only in SynhtI.2.0, some of the heterocyclization reagents will not be classified in SynthI.1.0:

```python
>>> BBClassifier(molSmiles="CCOC=1C=C(CC#N)C=CC1OCC(F)(F)F")
```
```text
[]
```

### Bulk BBs classification 

[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

Also *SynthI_Classifier.py* can be launched as a comand line tool for the BBs library separation into several sublibraries according to the BBs classes:

```shell script
$ python3 SynthI/src/SynthI_Classifier.py -h
```
```text
usage: SynthI_Classifier [-h] [-i INPUT]

Classification of building blocks. Separates provided library into several sublibraries according to the reagents classess.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input SMILES file

Code implementation:                Yuliana Zabolotna, Alexandre Varnek
                                    Laboratoire de Chémoinformatique, Université de Strasbourg.

Knowledge base (SMARTS library):    Dmitriy M.Volochnyuk, Sergey V.Ryabukhin, Kostiantyn Gavrylenko, Olexandre Oksiuta 
                                    Institute of Organic Chemistry, National Academy of Sciences of Ukraine
                                    Kyiv National Taras Shevchenko University
2021 Strasbourg, Kiev 

```
As a result separate files for each BB class found in the provided library file will be created. The name contains both class and subclass names separated by underscore - *SecondaryAmines_Cyc-Aliphatic.smi* or *Acid_Aliphatic_Acid.smi*.

## SynthI-BBs 

[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

SynthI-BBs allows to perform scaffold analysis of BBs and generate exhaustively all possible synthons from a given BB. The position of the functional groups as well as type of the resulting intermediate product (cation, anion, radical etc.) is encoded in synthon’s SMILES by introducing special system of labels. 

|Synthon label | Synthon example | Nature of the reaction center | Example of corresponding reagent classes |
| :-: | :-: | :-: | :-: |
| AHn:10 | C1CC1[CH3:10] | Electrophilic | Acyl, aryl and alkyl halides, sulfonylhalides, anhydrides,  acides, aminoacids, esters, alcohols, aldehydes, ketones, Weinreb amides, acylated azides, iso(thio)cyanates, oxiranes |
| AHn:20 | C1CC1[NH2:20]  | Nucleophilic | Alcohols, thiols, amines, amides, NH-azoles, hydrazines, hydrazides, hydroxylamines, oximes, esters, element organics, metal organics, ketones, aryl and allyl sulphones, alkenes for Heck couplings |
| CHn:30 | C1CC1C[CH3:30]  | Bivalent electrophilic | Aldehydes, ketones |
| AHn:40 | C1CC1C[NH2:40]  | Bivalent nucleophilic | Ketones, primary amines, hydrazines, hydroxylamines, reagents for olefination (Jullia-Kocienski, Wittig, Horner-Wadsworth-Emmons) |
| CH3:50 | C1CC1[CH3:50]  | Bivalent neutral | Terminal alkenes (for metathesis) |
| CHn:60 | c1cc[cH:60]nc1  | Electrophilic radical | Minisci CH-partners, Michael acceptors |
| CHn:70 | C1CC1[CH3:70]  | Nucleophilic radical | BF3 and MIDA boronates, oxalate alkyl esters, NOPhtal alkyl esters, sulphinates |
| CHn:21 | C1CC1C[CH3:21]  | Boronics-derived nucleophilic | Boronic reagents |
| NH:11 | R1[NH:11]R2 | Electrophilic nitrogen | Benzoyl O-acylated hydroxilamines |

### Scaffold analysis

[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

SynthI-BBs allow to generate meaningful scaffolds from BBs by removing any ring-containing moieties that will not be kept in the reaction product and thus are irrelevant in BB analysis (e.g. protective (Bnz, Cbz, Fmoc) and leaving groups (boronics, oxiranes, etc.))

```python
>>> from SynthI.src.SynthI_BBs import generateScaffoldForBB
>>> generateScaffoldForBB("OC(=O)C=1C=CC=C(NC(=O)OCC2C=3C=CC=CC3C=4C=CC=CC24)C1")
```
```text
# here some RdKit messages can appear, they can be ignored

'c1ccccc1'
```

Also *SynthI_BBScaffoldGeneration.py* can be launched as a comand line tool for the scaffold analysis of large BBs library:

```shell script
$ python3 SynthI/SynthI_BBScaffoldGeneration.py -h 

usage: SynthI_BBScaffoldGeneration [-h] [-i INPUT] [-o OUTPUT]

BBs Scaffold analysis. Generates meaningful BBs scaffolds after removing ring-containing leaving and protective groups. Count scaffolds occurrence in the provided collection of BBs, and construct cumulative scaffold frequency plot

optional arguments: 

  -h, --help            show this help message and exit
  -i INPUT, --input INPUT 
                        Input BBs file.
  -o OUTPUT, --output OUTPUT
                        Output files suffix name. 

Code implementation:                Yuliana Zabolotna, Alexandre Varnek
                                    Laboratoire de Chémoinformatique, Université de Strasbourg.

Knowledge base (SMARTS library):    Dmitriy M.Volochnyuk, Sergey V.Ryabukhin, Kostiantyn Gavrylenko, Olexandre Oksiuta 
                                    Institute of Organic Chemistry, National Academy of Sciences of Ukraine
                                    Kyiv National Taras Shevchenko University
2021 Strasbourg, Kiev 

```

It generates four files: *outSuffixName_Scaffolds.smi*, *outSuffixName_scaffoldsCounts.smi*, *outSuffixName_cumulativeprecentage.smi*, *Scaffolds_FreqPlot_outSuffixName.png*.

### BBs synthonization

[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

In case if BB contain protective groups, protected synthons will be discarded by default

```python
>>> from SynthI.src.SynthI_BBs import mainSynthonsGenerator
>>> mainSynthonsGenerator("CCOC(=O)C1=C(N)SC=C1C2CC2")
```
```text
# here some RdKit messages can appear, they can be ignored 
___________________________________________________________________________________________________ 
All generated synthons (4): O=C(O)c1c(C2CC2)csc1[NH2:20].O=C(O)c1c(C2CC2)csc1[NH2:40].O=[CH:10]c1c(C2CC2)csc1[NH2:20].O=[CH:10]c1c(C2CC2)csc1[NH2:40]
Col1-Synton Col2-RespectiveBBsClass
O=C(O)c1c(C2CC2)csc1[NH2:20] Bifunctional_Amine_Ester
O=C(O)c1c(C2CC2)csc1[NH2:40] Bifunctional_Amine_Ester
O=[CH:10]c1c(C2CC2)csc1[NH2:20] Bifunctional_Amine_Ester 
O=[CH:10]c1c(C2CC2)csc1[NH2:40] Bifunctional_Amine_Ester 
```
They can be kept if `keepPG=True` is specified: 

```python
>>> mainSynthonsGenerator("CCOC(=O)C1=C(N)SC=C1C2CC2", keepPG=True)
```
```text
# here some RdKit messages can appear, they can be ignored
___________________________________________________________________________________________________ 
All generated synthons (6): CCOC(=O)c1c(C2CC2)csc1[NH2:20].CCOC(=O)c1c(C2CC2)csc1[NH2:40].O=C(O)c1c(C2CC2)csc1[NH2:20].O=C(O)c1c(C2CC2)csc1[NH2:40].O=[CH:10]c1c(C2CC2)csc1[NH2:20].O=[CH:10]c1c(C2CC2)csc1[NH2:40]
Col1-Synton Col2-RespectiveBBsClass 
CCOC(=O)c1c(C2CC2)csc1[NH2:20] Bifunctional_Amine_Ester
CCOC(=O)c1c(C2CC2)csc1[NH2:40] Bifunctional_Amine_Ester
O=C(O)c1c(C2CC2)csc1[NH2:20] Bifunctional_Amine_Ester
O=C(O)c1c(C2CC2)csc1[NH2:40] Bifunctional_Amine_Ester
O=[CH:10]c1c(C2CC2)csc1[NH2:20] Bifunctional_Amine_Ester
O=[CH:10]c1c(C2CC2)csc1[NH2:40] Bifunctional_Amine_Ester
```
As it was mentioned before, solvents and contriones should be removed before using SynthI. In case if two moieties will be present in input, syntons for both of them will be generated:
```python
>>> mainSynthonsGenerator("NCCN1CCOCC1.CC1=NN=C(N1)SCC(O)=O")
```
```text
# here some RdKit messages can appear, they can be ignored
___________________________________________________________________________________________________
All generated synthons (4): C1CN(CC[NH2:20])CCO1.C1CN(CC[NH2:40])CCO1.Cc1n[nH]c(SC[CH:10]=O)n1.Cc1nc(SC[CH:10]=O)[nH:20]n1 
Col1-Synton Col2-RespectiveBBsClass 
C1CN(CC[NH2:20])CCO1 Aminoacids_N-AliphaticAmino_Acid 
C1CN(CC[NH2:40])CCO1 Aminoacids_N-AliphaticAmino_Acid  
Cc1n[nH]c(SC[CH:10]=O)n1 Aminoacids_N-AliphaticAmino_Acid
Cc1nc(SC[CH:10]=O)[nH:20]n1 Aminoacids_N-AliphaticAmino_Acid+nHAzoles_nHAzoles
```
As an exception, several the most popularly occured solvents and contrions in BBs libraries are always ignored if present in a mixture (e.g. acetic, carbonic, formic, oxalic, trifluoroacetic acid etc.):
```python
>>> mainSynthonsGenerator("OC(=O)C(F)(F)F.O=C1CCCC2(CCCC2)C1")
```
```text
# here some RdKit messages can appear, they can be ignored
___________________________________________________________________________________________________
All generated synthons (6): C1CCC2(C1)CCC[CH2:10]C2.C1CCC2(C1)CCC[CH2:30]C2.O=C1CC2(CCCC2)CC[CH2:40]1.O=C1CCCC2(CCCC2)[CH2:40]1.O[CH:10]1CCCC2(CCCC2)C1.C1CCC2(C1)CCC[CH:10]([OH:20])C2  
Col1-Synton Col2-RespectiveBBsClass 
C1CCC2(C1)CCC[CH2:10]C2 Ketones_Ketones
C1CCC2(C1)CCC[CH2:30]C2 Ketones_Ketones
O=C1CC2(CCCC2)CC[CH2:40]1 Ketones_Ketones
O=C1CCCC2(CCCC2)[CH2:40]1 Ketones_Ketones 
O[CH:10]1CCCC2(CCCC2)C1 Ketones_Ketones
C1CCC2(C1)CCC[CH:10]([OH:20])C2 Alcohols_Aliphatic_alcohols+Ketones_Ketones
```

If used inside a custom forkflow user may want to get a python object that will be used latter in the script. In this case specify option `returnDict=True`:
```python
>>> synthonsDictionary = mainSynthonsGenerator("CCOC(=O)C1=C(N)SC=C1C2CC2", returnDict=True)

# here some RdKit messages can appear, they can be ignored

>>> synthonsDictionary

{'O=C(O)c1c(C2CC2)csc1[NH2:20]': {'Bifunctional_Amine_Ester'}, 
'O=C(O)c1c(C2CC2)csc1[NH2:40]': {'Bifunctional_Amine_Ester'},
'O=[CH:10]c1c(C2CC2)csc1[NH2:20]': {'Bifunctional_Amine_Ester'},
'O=[CH:10]c1c(C2CC2)csc1[NH2:40]': {'Bifunctional_Amine_Ester'}}
```

Resulted synthons can be filtered according to the Ro2:

```python
>>> from SynthI.src.UsefulFunctions import Ro2Filtration
>>> Ro2Filtration("O=C(O)c1c(C2CC2)csc1[NH2:20]")

(True, ['MolW=183.035399528', 'LogP=1.9059000000000001', 'HDC=2', 'HAC=4']) 

>>> Ro2Filtration("C1CCC2(C1)CCC[CH2:10]C2")

(False, ['MolW=138.140850576', 'LogP=3.510900000000002', 'HDC=0', 'HAC=0'])
```

### Bulk synthons generation for the large BBs library

[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

In case of large BB library, classification and synthonization of BBs can be performed using command line tool *SynthI/SynthI_BBsBulkClassificationAndSynthonization.py*:

```shell script
$ python3 SynthI/SynthI_BBsBulkClassificationAndSynthonization.py -h

usage: SynthI_BBsBulkClassificationAndSynthonization [-h] [-i INPUT] [-o OUTPUT] [--keepPG] [--Ro2Filtr] [--nCores NCORES] 

BBs classification and Synthons generation for large BBs libraries

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT        
                        Input file containing building blocks smiles and ids.
  -o OUTPUT, --output OUTPUT
                        Output files suffix name.
  --keepPG              Write both protected and unprotected synthons to the output (concerns Boc, Bn, Fmoc, Cbz and Esters protections).  
  --Ro2Filtr            Write only synthons satisfying Ro2 (MW <= 200, logP <= 2, H-bond donors count <= 2 and H-bond acceptors count <= 4)  
  --nCores NCORES       Number of available cores for parallel calculations. Memory usage is optimized, so maximal number of parallel processes can be launched.

Code implementation:                Yuliana Zabolotna, Alexandre Varnek
                                    Laboratoire de Chémoinformatique, Université de Strasbourg.

Knowledge base (SMARTS library):    Dmitriy M.Volochnyuk, Sergey V.Ryabukhin, Kostiantyn Gavrylenko, Olexandre Oksiuta 
                                    Institute of Organic Chemistry, National Academy of Sciences of Ukraine
                                    Kyiv National Taras Shevchenko University
2021 Strasbourg, Kiev 

```
As a result four files will be generated: *outSuffixName_BBmode.smi*, *outSuffixName_Synthmode.smi*, *outSuffixName_NotClassified*, *outSuffixName_NotProcessed*. 

Each line of the file *outSuffixName_BBmode.smi* contains SMILES and ID of the classified BB, assigned classes and list of generated synthons:
```text
1-BB 2-ID 3-BBClasses 4-Synthons 5-NumberOfGeneratedSynthons

CNCC=1C=CC=C(OC)C1 EN300-06971 SecondaryAmines_secBenzylic,SecondaryAmines_secAliphatic COc1cccc(C[NH:20]C)c1 1 
CC(C)(N)C#N EN300-29407 PrimaryAmines_PriAmines_Aliphatic CC(C)(C#N)[NH2:40].CC(C)(C#N)[NH2:20] 2
Cl.CCCOC=1C=CC=CC1N EN300-08901 PrimaryAmines_PriAmines_Anilines CCCOc1ccccc1[NH2:40].CCCOc1ccccc1[NH2:20] 2 
OCC1CCCCC1 EN300-21580 Alcohols_Aliphatic_alcohols C1CCC(C[OH:20])CC1.C1CCC([CH3:10])CC1 2 
``` 
File *outSuffixName_Synthmode.smi* contains the same information, but the focus is now on the unique synthons:
```text
1-synthon 2-IDsOfRespectiveBBs 3-BBClasses 4-RespectiveBBs 5-NumberOfRespectiveBBs 6-UniqSynthonsID

Clc1ccccc1C(CBr)[OH:20] EN300-43119 Alcohols_Aliphatic_alcohols O[C@H](CBr)C=1C=CC=CC1Cl 1 outSuffixName_1371
Clc1ccccc1[CH2:10]CBr EN300-43119 Alcohols_Aliphatic_alcohols O[C@H](CBr)C=1C=CC=CC1Cl 1 outSuffixName_1372
OC(c1ccccc1Cl)[CH3:10] EN300-43119 AlkylHalides_Alkyl_halides O[C@H](CBr)C=1C=CC=CC1Cl 1 outSuffixName_1373 
Clc1ccccc1C([CH3:10])[OH:20] EN300-43119 AlkylHalides_Alkyl_halides+Alcohols_Aliphatic_alcohols O[C@H](CBr)C=1C=CC=CC1Cl 1 outSuffixName_1374 
N[CH:10]=O EN300-50197 Reagents_Isocyanates ClC(Cl)(Cl)C(=O)N=C=O 1 outSuffixName_12733  
``` 
Files *outSuffixName_NotProcessed* and *outSuffixName_NotClassified* contain not processed by RdKit or processed but not classified BBs.

## SynthI-Fragmentation

[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

This module allows to fragment given molecule and generate synthons that correspond to particular BBs, needed to easily synthesize input compound.

### Fragmentation with a default setup

[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

In order to perform compoound fragmentation, first, the Fragmentor (Instant of the class fragmentation) should be initialized

```python
>>> from SynthI.src.SynthI import * 

>>> SynthLibrary = "/pathToTheSynthonsLib/outENSynthmode.smi"
>>> FragmentsToIgnore = ["*C(C)C", "*C(=O)C", "*C=O", "*[V]C=O", "*[V]C(C)C", "*[V]C(=O)C"]

>>> SynthIfragmentor = fragmentation(mode="use_all", maxNumderOfReactionCentersPerFragment=3, MaxNumberOfStages = 5, 
...                                     SynthLibrary=SynthLibrary, FragmentsToIgnore=FragmentsToIgnore, 
...                                     FindAnaloguesOfMissingSynthons=True)
```
```text
Processing BB library. It may take a few minutes, depending on the library size

Lib BB reading time:
0:01:33.865876
```
Latter `SynthIfragmentor` cna be used to fragment different molecules with a help of function `fragmentMolecule (smiles, SynthIfragmentor, simTh=-1)`

```python
>>> smi = "NC(=O)OC(CN1N=CN=N1)C1=CC=CC=C1Cl"
>>> allSyntheticPathways, allSynthons = fragmentMolecule(smi, SynthIfragmentor)
# here some RdKit messages can appear, they can be ignored
```

Both `allSyntheticPathways` and `allSynthons` are dictionaries as values containing instances of the classes *synthon* and *syntheticPathway* ( [see detailes here](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#detailed-classes-description) )
 
```python
>>> allSynthons
 
{'InitMol': <SynthI.src.SynthI.synthon at 0x7fcb4588b430>,
 'N[CH:10]=O': <SynthI.src.SynthI.synthon at 0x7fcb4588b550>,
 'Clc1ccccc1C(Cn1ncnn1)[OH:20]': <SynthI.src.SynthI.synthon at 0x7fcb4588bbb0>,
 'c1nn[nH:20]n1': <SynthI.src.SynthI.synthon at 0x7fcb4588b6a0>,
 'NC(=O)OC(c1ccccc1Cl)[CH3:10]': <SynthI.src.SynthI.synthon at 0x7fcb4588b3d0>,
 'NC(=O)OC(c1ccccc1Cl)[CH3:21]': <SynthI.src.SynthI.synthon at 0x7fcb4588beb0>,
 'Clc1ccccc1C([CH3:10])[OH:20]': <SynthI.src.SynthI.synthon at 0x7fcb3007b520>,
 'Clc1ccccc1C([OH:20])[CH3:21]': <SynthI.src.SynthI.synthon at 0x7fcb3007be80>,
 'c1nnn(C[CH2:10][OH:20])n1': <SynthI.src.SynthI.synthon at 0x7fcb3007ba00>, 
 'Clc1cccc[cH:20]1': <SynthI.src.SynthI.synthon at 0x7fcb3007ba30>,
 'Clc1ccccc1[CH2:10][OH:20]': <SynthI.src.SynthI.synthon at 0x7fcb3007bdf0>,
 'c1nnn([CH3:20])n1': <SynthI.src.SynthI.synthon at 0x7fcb3007bc70>} 
```
The detailed information about each synthon can be retreived:
```python
>>> allSynthons['c1nnn(C[CH2:10][OH:20])n1'].printSynthonInfo()    
```
```text
__________________________________________________________ 
                  Synthon Information         
__________________________________________________________

Synthon: c1nnn(C[CH2:10][OH:20])n1   
Synthon was not found in provided library of building blocks. 2 analog(s) has/have been found 
BB analogues: 
C[CH:10](c1nc[nH]n1)[OH:20] EN300-137277 
C[CH:10](c1c[nH]nn1)[OH:20] EN300-7472008 
Parent synthons: Clc1ccccc1C(Cn1ncnn1)[OH:20]
Children synthons: -
```
For the synthetic pathways there is a way to display short and detailed information
```python
>>> for key in allSyntheticPathways:
...     allSyntheticPathways[key].printShortReagentSetInfo()
```
```text
R2.2_0 N[CH:10]=O.Clc1ccccc1C(Cn1ncnn1)[OH:20] Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17
R5.1_0 c1nn[nH:20]n1.NC(=O)OC(c1ccccc1Cl)[CH3:10] Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.0 
R5.2_0 c1nn[nH:20]n1.NC(=O)OC(c1ccccc1Cl)[CH3:21] Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.0
R2.2_0|R5.1_0 c1nn[nH:20]n1.Clc1ccccc1C([CH3:10])[OH:20].N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.72 
R2.2_0|R5.2_0 c1nn[nH:20]n1.Clc1ccccc1C([OH:20])[CH3:21].N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17 
R10.1_0|R2.2_0 c1nnn(C[CH2:10][OH:20])n1.Clc1cccc[cH:20]1.N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17
R10.1_1|R2.2_0 Clc1ccccc1[CH2:10][OH:20].c1nnn([CH3:20])n1.N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17 
```
This synthetic pathways are organized in a disconnection hierarchy, that can be navigated with the help of several functions: 
* *getLongestSyntheticPathways()* - creates a list of the synthetic pathways including the largest number of stages (leafs of the hierarchy)

    ```python
    >>> LongestSyntheticPathways = getLongestSyntheticPathways(allSyntheticPathways)
    >>> for ind,reagentSet in enumerate(LongestSyntheticPathways):
    ...     print("reagentSet " + str(ind) + " :")
    ...     reagentSet.printShortReagentSetInfo()
    ```  
    ```text
    reagentSet 0 :
    R2.2_0|R5.1_0 c1nn[nH:20]n1.Clc1ccccc1C([CH3:10])[OH:20].N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.72
    reagentSet 1 :
    R2.2_0|R5.2_0 c1nn[nH:20]n1.Clc1ccccc1C([OH:20])[CH3:21].N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17
    reagentSet 2 :
    R10.1_0|R2.2_0 c1nnn(C[CH2:10][OH:20])n1.Clc1cccc[cH:20]1.N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17 
    reagentSet 3 : 
    R10.1_1|R2.2_0 Clc1ccccc1[CH2:10][OH:20].c1nnn([CH3:20])n1.N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17
    ```
  
* *getShortestSyntheticPathways()* - creates a list of the synthetic pathways including only one stage (roots of the hierarchy)
   
    ```python
    >>> firstLevelSynthonsCombinations = getShortestSyntheticPathways(allSyntheticPathways) 
    >>> for ind,reagentSet in enumerate(firstLevelSynthonsCombinations):
    ...     print("reagentSet " + str(ind) + " :")
    ...     reagentSet.printShortReagentSetInfo()
     ```
     ```text
    reagentSet 0 :
    R2.2_0 N[CH:10]=O.Clc1ccccc1C(Cn1ncnn1)[OH:20] Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17
    reagentSet 1 :   
    R5.1_0 c1nn[nH:20]n1.NC(=O)OC(c1ccccc1Cl)[CH3:10] Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.0
    reagentSet 2 : 
    R5.2_0 c1nn[nH:20]n1.NC(=O)OC(c1ccccc1Cl)[CH3:21] Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.0
     ```
* *findShortestSynthPathWithAvailableBBlib()* - creates a list of the synthetic pathways having the highest value of Availability rate (% of atoms of fragmented molecule coming from available synthons)

     ```python
    >>> shortestSynthesis = findShortestSynthPathWithAvailableBBlib(firstLevelSynthonsCombinations ,  showAll=True) 

    1 equivalent synthetic pathway(s) have been found.

    >>> for ind,reagentSet in enumerate(shortestSynthesis):
    ...     print("reagentSet " + str(ind) + " :")
    ...     reagentSet.printShortReagentSetInfo()
    ```
    ```text
    reagentSet 0 :
    R2.2_0|R5.1_0 c1nn[nH:20]n1.Clc1ccccc1C([CH3:10])[OH:20].N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.72 
    ```
Method *printDetailedReagentsSetInfo()* of the class *syntheticPathway* ( [see detailes here](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#detailed-classes-description) ) can be used for retreiving detailed information about the selected synthetic pathway: 

```python
>>> shortestSynthesis[0].printDetailedReagentsSetInfo()  
```
```text
**********************************************************
Reagent set Information R2.2_0|R5.1_0
********************************************************** 
Reactions: O-Acylation by O=C(+)-X reagents->nH-SN alkylation  
Required Synthons: c1nn[nH:20]n1.Clc1ccccc1C([CH3:10])[OH:20].N[CH:10]=O 
Number of reagents: 3  
Number of stages: 2    
Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.72
Parent reagent sets:  
Reactions: R2.2_0 O-Acylation by O=C(+)-X reagents ||| Participating synthons: N[CH:10]=O.Clc1ccccc1C(Cn1ncnn1)[OH:20]
Reactions: R5.1_0 nH-SN alkylation ||| Participating synthons: c1nn[nH:20]n1.NC(=O)OC(c1ccccc1Cl)[CH3:10]   

__________________________________________________________ 
                  Synthon Information       
__________________________________________________________ 

Synthon: c1nn[nH:20]n1
Synthon was not found in provided library of building blocks. 4 analog(s) has/have been found
BB analogues:  
c1nnc[nH:20]1 EN300-20608 
c1c[nH:20]nn1 EN300-27201  
Cc1nnn[nH:20]1 EN300-104530    
Nc1nnn[nH:20]1 EN300-33999 
Parent synthons: Clc1ccccc1C(Cn1ncnn1)[OH:20]
Children synthons: -     

__________________________________________________________      
                  Synthon Information   
__________________________________________________________         

Synthon: Clc1ccccc1C([CH3:10])[OH:20] 
Available. Corresponding BBs: EN300-43119  
Parent synthons: Clc1ccccc1C(Cn1ncnn1)[OH:20] 
Children synthons: -        

__________________________________________________________      
                  Synthon Information   
__________________________________________________________    
     
Synthon: N[CH:10]=O
Available. Corresponding BBs: EN300-50197 
Parent synthons: -   
Children synthons: -      
```

### Selecting customized list of reactions for fragmentation

[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

**Mode "include_only"**

Only the reactions, selected by user will be used for fragmentation. The list of RiDs should be specified using argument `reactionsToWorkWith`. The list should be provided inside " "; intervals separated via "-" and "," can be used (e.g. "R1-R10,R11.1-R11.4,R12.1"). Specification "R1" implicitly includes all (R1.1, R1.2, R1.3 and R1.4) subreactions in the group R1.
  
```python
>>> SynthIfragmentorIncludeOnlyR1_R9 = fragmentation(mode="include_only", reactionsToWorkWith = "R1-R9", 
...                                                   maxNumderOfReactionCentersPerFragment=3, MaxNumberOfStages = 5, 
...                                                   SynthLibrary=SynthLibrary, FragmentsToIgnore=FragmentsToIgnore, 
...                                                   FindAnaloguesOfMissingSynthons=True) 

>>> allSyntheticPathways, allSynthons = fragmentMolecule(smi, SynthIfragmentorIncludeOnlyR1_R9) 

>>> for key in allSyntheticPathways:  
...     allSyntheticPathways[key].printShortReagentSetInfo() 
```
```text
R2.2_0 N[CH:10]=O.Clc1ccccc1C(Cn1ncnn1)[OH:20] Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17
R5.1_0 c1nn[nH:20]n1.NC(=O)OC(c1ccccc1Cl)[CH3:10] Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.0 
R5.2_0 c1nn[nH:20]n1.NC(=O)OC(c1ccccc1Cl)[CH3:21] Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.0
R2.2_0|R5.1_0 c1nn[nH:20]n1.Clc1ccccc1C([CH3:10])[OH:20].N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.72 
R2.2_0|R5.2_0 c1nn[nH:20]n1.Clc1ccccc1C([OH:20])[CH3:21].N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17 
```
**Mode "exclude_some"**

The list of RiDs of reactions that need to be excluded should be specified using argument `reactionsToWorkWith`.

In the example below all reactions except R5.1 (nH-SN alkylation of NH-heterocycles) will be used for fragmentaion. 

```python
>>> SynthIfragmentorExcludeSomeR5_1 = fragmentation(mode="exclude_some", reactionsToWorkWith = "R5.1",
...                                                 maxNumderOfReactionCentersPerFragment=3, MaxNumberOfStages = 5, 
...                                                 SynthLibrary=SynthLibrary, FragmentsToIgnore=FragmentsToIgnore, 
...                                                 FindAnaloguesOfMissingSynthons=True) 

>>> allSyntheticPathways, allSynthons = fragmentMolecule(smi, SynthIfragmentorExcludeSomeR5_1) 

>>> for key in allSyntheticPathways:  
...     allSyntheticPathways[key].printShortReagentSetInfo() 
```
```text
R2.2_0 N[CH:10]=O.Clc1ccccc1C(Cn1ncnn1)[OH:20] Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17 
R5.2_0 c1nn[nH:20]n1.NC(=O)OC(c1ccccc1Cl)[CH3:21] Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.0 
R2.2_0|R5.2_0 c1nn[nH:20]n1.Clc1ccccc1C([OH:20])[CH3:21].N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17
R10.1_0|R2.2_0 c1nnn(C[CH2:10][OH:20])n1.Clc1cccc[cH:20]1.N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17
R10.1_1|R2.2_0 Clc1ccccc1[CH2:10][OH:20].c1nnn([CH3:20])n1.N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17 
```
**Mode "one_by_one"**

In this mode user-provided reactions are applied in the specified order; each bond can be cut using onlly one reaction rule. If `reactionsToWorkWith` is not specified, than all reactions in the default order will be applied.

```python
>>> SynthIfragmentorOneByOneR2_R10_R5 = fragmentation(mode="one_by_one", reactionsToWorkWith = "R2,R10,R5",
...                                                   maxNumderOfReactionCentersPerFragment=3, MaxNumberOfStages = 5, 
...                                                   SynthLibrary=SynthLibrary, FragmentsToIgnore=FragmentsToIgnore, 
...                                                   FindAnaloguesOfMissingSynthons=True) 

>>> allSyntheticPathways, allSynthons = fragmentMolecule(smi, SynthIfragmentorOneByOneR2_R10_R5 ) 

>>> for key in allSyntheticPathways:  
...     allSyntheticPathways[key].printShortReagentSetInfo() 
```
```text
R2.2_0 N[CH:10]=O.Clc1ccccc1C(Cn1ncnn1)[OH:20] Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17 
R10.1_0|R2.2_0 c1nnn(C[CH2:10][OH:20])n1.Clc1cccc[cH:20]1.N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17
R10.1_1|R2.2_0 Clc1ccccc1[CH2:10][OH:20].c1nnn([CH3:20])n1.N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17
```
If there are several ways to cut particular chemical bond, it will be cut only according to the rule that comes first in the customized ordered list of reactions to use. 
```python
>>> SynthIfragmentorOneByOneR2_R5_R10  = fragmentation(mode="one_by_one", reactionsToWorkWith = "R2,R5,R10",
...                                                    maxNumderOfReactionCentersPerFragment=3, MaxNumberOfStages = 5, 
...                                                    SynthLibrary=SynthLibrary, FragmentsToIgnore=FragmentsToIgnore,
...                                                    FindAnaloguesOfMissingSynthons=True) 

>>> allSyntheticPathways, allSynthons = fragmentMolecule(smi, SynthIfragmentorOneByOneR2_R5_R10  ) 

>>> for key in allSyntheticPathways:  
...     allSyntheticPathways[key].printShortReagentSetInfo() 
```
```text
R2.2_0 N[CH:10]=O.Clc1ccccc1C(Cn1ncnn1)[OH:20] Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.17
R2.2_0|R5.1_0 c1nn[nH:20]n1.Clc1ccccc1C([CH3:10])[OH:20].N[CH:10]=O Availability rate (% of atoms of fragmented molecule coming from available synthons): 0.72 
```

### Bulk compounds fragmentation

[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

Library of compounds can be fragmented using *SynthI_BulkFragmentationEnumerationAndAnaloguesDesign.py* script.

```shell script
$ python3 SynthI/SynthI_BulkFragmentationEnumerationAndAnaloguesDesign.py -h

usage: SynthI_BulkFragmentationEnumerationAndAnaloguesDesign [-h] [-i INPUT] [-oD OUTDIR] [--SynthLibrary SYNTHLIBRARY] 
                                                     [--nCores NCORES]  [--analoguesLibGen] [--strictAvailabilityMode]
                                                     [--simBBselection] [--Ro2Filtration] [--mode MODE] [--simTh SIMTH]
                                                     [--reactionsToWorkWith REACTIONSTOWORKWITH] [--MaxNumberOfStages MAXNUMBEROFSTAGES]
                                                     [--maxNumberOfReactionCentersPerFragment MAXNUMBEROFREACTIONCENTERSPERFRAGMENT]

Compound fragmentaitiona and analogues generation. 

optional arguments:
  -h, --help            show this help message and exit       
  -i INPUT, --input INPUT        
                        input file        
  -oD OUTDIR, --outDir OUTDIR   
                        Output directory to write analogues. 
  --SynthLibrary SYNTHLIBRARY                        
                        Library of available synthons. Generated from avaialable BBs using SynthI_BBsBulkClassificationAndSynthonization.py  
  --nCores NCORES       Number of CPUs available for parallelization.          
  --simTh SIMTH         Similarity threshold for BB analogues search. If not specified, only positional variational approach will be used for BBs search             
  --analoguesLibGen     Generate library of analogues from input mol                                                       
  --strictAvailabilityMode         
                        Only fully synthesizable analogues are generated. Alternatively, unavailable synthons resulted from compound fragmentation will still be used for its analogues generation.
  --simBBselection       Used always with analoguesLibGen.For library generation will be used not only synthons from input molecules but also their closest analogues                                              
  --Ro2Filtration       Filter input synthons library by Ro2 (MW <= 200, logP <= 2, H-bond donors count <= 2 and H-bond acceptors count <= 4) 
  --mode MODE           Mode of fragmentation (defines how the reaction list is specified)          
                        Possible options: use_all, include_only, exclude_some, one_by_one      
                        (default: use_all)        
  --reactionsToWorkWith REACTIONSTOWORKWITH    
                        List of RiDs to be used.        
                        (default: R1-R13 (all reactions)                
  --MaxNumberOfStages MAXNUMBEROFSTAGES           
                        Maximal number of stages during fragmentation.        
                        (default: 5)           
  --maxNumberOfReactionCentersPerFragment MAXNUMBEROFREACTIONCENTERSPERFRAGMENT  
                        Maximal number of reaction centers per fragment.    
                        (default: 3)          

_________________________________________________________________________________________________________________________ 

Code implementation:                Yuliana Zabolotna, Alexandre Varnek             
                                    Laboratoire de Chémoinformatique, Université de Strasbourg.        
Knowledge base (SMARTS library):    Dmitriy M.Volochnyuk, Sergey V.Ryabukhin, Kostiantyn Gavrylenko, Olexandre Oksiuta   
                                    Institute of Organic Chemistry, National Academy of Sciences of Ukraine    
                                    Kyiv National Taras Shevchenko University      
2021 Strasbourg, Kiev   

```
Example of launch:
```shell script
python3 SynthI/SynthI_FragmentationEnumerationUsedInArticle.py -i FDA_small_drugs_SMILES_ForFragmentation.smiles  
-oD testDir --SynthLibrary outENSynthmode.smi --MaxNumberOfStages 5 --maxNumberOfReactionCentersPerFragment 3 
--mode use_all --reactionsToWorkWith R1-R13 --nCores 15

```
It produces 2 files: 
* *InputName_out* - contains only one synthetic pathway per compound (selected by the availability rate)

    ```text
    1-InitialCompound 2-AllSynthons 3-ReactionsUsedInFragmentation 4-NumberOfSynthons 5-AvailabilityRate 6-AvailableSynthons 7-NotAvailableSynthons
    
    C1=CC=C(C(=C1)C(CN2N=CN=N2)OC(=O)N)Cl c1nn[nH:20]n1.Clc1ccccc1C([CH3:10])[OH:20].N[CH:10]=O R2.2_0|R5.1_0 3 0.72 AvailableSynthons:Clc1ccccc1C([CH3:10])[OH:20]->EN300-43119|N[CH:10]=O->EN300-50197 NotAvailableSynthons:c1nn[nH:20]n1
    ``` 
 
 * *allSythons_InputName_out* - contains all synthons from the "leaf" synthetic pathways  (largest number of stages).
    ```test
   1-InitialCompound 2-AvailableSynthons 3-NotAvailableSynthons
   
   C1=CC=C(C(=C1)C(CN2N=CN=N2)OC(=O)N)Cl AvailableSynthons:Clc1ccccc1C([CH3:10])[OH:20].N[CH:10]=O notAvaialableSynthons:Clc1ccccc1C([OH:20])[CH3:21].Clc1cccc[cH:20]1.c1nnn(C[CH2:10][OH:20])n1.c1nn[nH:20]n1.Clc1ccccc1[CH2:10][OH:20].c1nnn([CH3:20])n1
    ```
 

## SynthI-Enumeration

[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

### Generate analogues of a compound

[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

### Enumerate library of all possible compounds using given set of synthons

[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

### Detailed classes description

[Back to Table of Contents](https://github.com/Laboratoire-de-Chemoinformatique/SynthI#table-of-contents)

_**CLASS `synthon (smiles, cutLevel=1, directParent=None, directChildren=None, syntheticPathway=None, BBlibProvided=False)`**_

The class to store Synthons. The instances of this class are created during compound fragmentaion.

**Available attributes**
* ``synthonInstance.smiles`` - SMILES of the synthon
* ``synthonInstance.functionalityCount`` - number of reactive centers in the synthon
* ``synthonInstance.marks`` - labels defining reactivity of the reactive centers of the synthon
* ``synthonInstance.directParents`` - parent synthons (their fragmentation lead to the current synthon)
* ``synthonInstance.directChildren`` - synthons obtained via fragmentation of the current synthon
* ``synthonInstance.syntheticPathway`` - list containig all synthetic pathways, that include current synthon 
* ``synthonInstance.correspondingBB`` - IDs of the BBs that produce current synthon
* ``synthonInstance.bbAnalogues`` - dictionary containing BBs that produce analogues synthons

**Methods**
* ``synthonInstance.printSynthonInfo()`` - print information about synthon
* ``synthonInstance.searchForSynthonAnalogues(synthLib: dict, simTh=-1)`` - search of the analogues of the current synthon in the provided library of avaialble syntons. 

Dictionary *synthLib* can be obtained using `UsefulFunctions.readSyntonLib(synthLibFile, Ro2Filtration=False, FindAnaloguesOfMissingBBs=False)` or retrieved as an attribute of the instant of the class fragmentation if the library was provided during class initiation (`fragmentation.SynthLib`). 

_**CLASS `syntheticPathway (name , synthPathwayReactions,  reagentsNumber, cutLevel, directParentsSynthPathways=None, SynthLibProvided=False)`**_

The class to store information about possible syntheticPathway for compound synthesis - a set of synthons obtatined via fragmentation and list of reaction rules used for this fragmentation. The instances of this class are created during compound fragmentaion.

**Available attributes**
* ``syntheticPathwayInstance.name`` - the name of syntheticPathway contains RiDs of the reaction rules used for compound fragmentation in this particular pathway
* ``syntheticPathwayInstance.participatingSynthon`` - list of synthons obtained via fragmentation according to this pathway. This list contain instances of *class synthon*
* ``syntheticPathwayInstance.directParentsSynthPathways`` - syntheticPathways having one stage less than current pathway and thus placed upper in the disconnection hierarchy 
* ``syntheticPathwayInstance.directChildrenSynthPathways`` - syntheticPathways having one stage more than current pathway and thus placed lower in the disconnection hierarchy
* ``syntheticPathwayInstance.synthPathwayReactions`` - list of reations of the current pathway
* ``syntheticPathwayInstance.reagentsNumber`` - number of synthons resulted from molecule fragmentation according to the current synthetic pathway 
* ``syntheticPathwayInstance.availabilityRate`` - Availability rate (% of atoms of fragmented molecule coming from available synthons) for the current synthetic pathway


**Methods**

* ``syntheticPathwayInstance.printShortReagentSetInfo()`` - print short information about syntheticPathway (RiDs, participating synthons and avaialability rate)
* ``syntheticPathwayInstance.printDetailedReagentsSetInfo()`` - print detailed information about syntheticPathway
* ``syntheticPathwayInstance.checkAvailability(self, SynthLib: dict, simTh=-1, FindAnaloguesOfMissingSynthons=True)`` - if  SynthLib containing available synthons is provided, participating synthons will be looked up there and avaialability rate for the pathway will be calculated. Changes will be made directly in the *synthon* and *syntheticPathway* objects

Dictionary *synthLib* can be obtained using `UsefulFunctions.readSyntonLib(synthLibFile, Ro2Filtration=False, FindAnaloguesOfMissingBBs=False)` or retrieved as an attribute of the instant of the class fragmentation if the library was provided during class initiation (`fragmentation.SynthLib`).

_**CLASS `fragmentation (ode="use_all", reactionsToWorkWith = "R1-R13", maxNumderOfReactionCentersPerFragment = 3, MaxNumberOfStages = 5, FragmentsToIgnore = None,
                 FindAnaloguesOfMissingSynthons = False, parsedSynthLib = False, SynthLibrary=None, Ro2SynthonsFiltration = False)`**_

The class to store setup for the fragmentation (including parsed synthons library, if providded).

**Available attributes**
* ``fragmentationInstance.SynthLib`` - dictionary containig parsed synthons library (if it was provided during class initialization)

**Methods**

* ``fragmentationInstance.cutWithHierarchyStorred(mol)`` - fragment provided molecule (should be RdKit molecule object and not smiles). Return two dictionary - *allSyntheticPathways* and *allSynthons* obtained during fragmentation. 
* ``fragmentationInstance.getReactionForReconstruction()`` - get reconstruction reaction rules based on what was used for fragmentation. 

_**CLASS `enumeration (outDir, Synthons=None, reactionSMARTS=None, maxNumberOfReactedSnthons=6, MWupperTh=None, MWlowerTh=None,
                  minNumberOfNewMols = 1000, nCores=1, analoguesEnumeration=False)`**_
    
The class to store setup for the compounds enumeration.

*reactionSMARTS* should be obtained from fragmentationInstance (fragmentationInstance.getReactionForReconstruction()).

*MWupperTh* and *MWlowerTh* if present allow to filter enumerated compounds in order to store only molecules of specified size.

*nCores* specify number of CPUs that can be used for parallelized compound enumeration.

Generation of analogues of the specified molecule and unbiased library enumeration based on the provided list of synthons differ slightly in the setup. Therefore, `analoguesEnumeration=True` or `analoguesEnumeration=False` respectively should be specified.  

**Methods**

* ``enumerationInstance.getReconstructedMols(mol)`` - fragment provided molecule (should be RdKit molecule object and not smiles). Return two dictionary - *allSyntheticPathways* and *allSynthons* obtained during fragmentation.

