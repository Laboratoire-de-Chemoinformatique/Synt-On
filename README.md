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

from SynthI.src.SynthI_Classifier import BBClassifier
BBClassifier(molSmiles="CCOC(=O)C1=C(N)SC=C1C2CC2")
```
```
['Bifunctional_Amine_Ester', 'PrimaryAmines_PriAmines_Het-Anilines']
```

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
