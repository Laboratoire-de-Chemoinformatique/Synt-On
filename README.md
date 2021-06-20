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

### Scaffold analysis

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
