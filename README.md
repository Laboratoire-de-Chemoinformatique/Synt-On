# SynthI
Open-source tool for synthons-based library design.


All BBs structures need to be sanitized and standartized idependently by user prior to SynthI usage. Solvents and contriones should be deleted. There is no need to generate major tautomer form as soon as SynthI will do it for each generated synthons separately.    

# Table of Contents
* [Building blocks Classification]()
* [SynthI-BBs]()
    * [Scaffold analysis]()
    * [BBs synthonization]()
    * [Bulk synthons generation for the large BBs library]()


## Building blocks Classification
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

### Bulk BBs classification 

Also *SynthI_Classifier.py* can be launched as a comand line tool for the BBs library separation into several sublibraries according to the BBs classes:

```shell script
$ python3 SynthI/src/SynthI_Classifier.py -h

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
SynthI-BBs allows to perform scaffold analysis of BBs and generate exhaustively all possible synthons from a given BB. The position of the functional groups as well as type of the resulting intermediate product (cation, anion, radical etc.) is encoded in synthon’s SMILES by introducing special system of labels. 

|Synthon lable | Synthon example | Nature of the reaction center | Example of corresponding reagent classes |
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

SynthI-BBs allow to generate meaningful scaffolds from BBs by removing any ring-containing moieties that will not be kept in the reaction product and thus are irrelevant in BB analysis (e.g. protective (Bnz, Cbz, Fmoc) and leaving groups (boronics, oxiranes, etc.))

```python
>>> from SynthI.src.SynthI_BBs import generateScaffoldForBB
>>> generateScaffoldForBB("OC(=O)C=1C=CC=C(NC(=O)OCC2C=3C=CC=CC3C=4C=CC=CC24)C1")
```
```text
# here some RdKit messages can appear, they can be ignored

'c1ccccc1'

```


SynthI-BBs allow to generate meaningful scaffolds from BBs by removing any ring-containing moieties that are not parts that will not be kept in the reaction product and thus are irrelevant in BB analysis (e.g. protective (Bnz, Cbz, Fmoc) and leaving groups (boronics, oxiranes, etc.)):

### BBs synthonization

### Bulk synthons generation for the large BBs library
In case of large BB library, classification and synthonization of BBs can be performed using command line tool *SynthI/SynthI_BBsBulkClassificationAndSynthonization.py*:

## SynthI-Fragmentation

### Selecting customized list of reactions for fragmentation

**Mode "include only"**

The list of RiDs of selected reactions should be specified using argument `reactionsToWorkWith`. The list should be provided inside " "; intervals separated via "-" and "," can be used (e.g. "R1-R10,R11.1-R11.4,R12.1"). Specification "R1" implicitly includes all (R1.1, R1.2, R1.3 and R1.4) subreactions in the group R1.
  
**Mode "exclude some"**

### Detailed classes description
