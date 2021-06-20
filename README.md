# SynthI
Open-source tool for synthons-based library design.

```python
>>> from SynthI.src.SynthI import * 

>>> SynthLibrary = "/data/outENSynthmode.smi"
>>> FragmentsToIgnore = ["*C(C)C", "*C(=O)C", "*C=O", "*[V]C=O", "*[V]C(C)C", "*[V]C(=O)C"]

>>> SynthIfragmentor = fragmentation(mode="use all", maxNumderOfReactionCentersPerFragment=3, MaxNumberOfStages = 5, 
...                                     SynthLibrary=SynthLibrary, FragmentsToIgnore=FragmentsToIgnore, 
...                                     FindAnaloguesOfMissingSynthons=True)
```
```text
Processing BB library. It may take a few minutes, depending on the library size

Lib BB reading time:
0:01:33.865876
```
