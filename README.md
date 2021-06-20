# SynthI
Open-source tool for synthons-based library design.


```python
>>> synthonsDictionary = mainSynthonsGenerator("CCOC(=O)C1=C(N)SC=C1C2CC2", returnDict=True)

# here some RdKit messages can appear, they can be ignored

>>> synthonsDictionary

{'O=C(O)c1c(C2CC2)csc1[NH2:20]': {'Bifunctional_Amine_Ester'}, 
'O=C(O)c1c(C2CC2)csc1[NH2:40]': {'Bifunctional_Amine_Ester'},
'O=[CH:10]c1c(C2CC2)csc1[NH2:20]': {'Bifunctional_Amine_Ester'},
'O=[CH:10]c1c(C2CC2)csc1[NH2:40]': {'Bifunctional_Amine_Ester'}}
```

```python
>>> mainSynthonsGenerator("OC(=O)C(F)(F)F.O=C1CCCC2(CCCC2)C1")
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
