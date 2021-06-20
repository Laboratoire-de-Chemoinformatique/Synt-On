# SynthI
Open-source tool for synthons-based library design.

This synthetic pathways are organized in a disconnection hierarchically, that can be navigated with the help of functions: 
* *getLongestSyntheticPathways()* - creates a list of the synthetic pathways containing the largest number of stages (leafs of the hierarchy)
* ```python
    >>> LongestSyntheticPathways = getLongestSyntheticPathways(allSyntheticPathways)
    >>> for ind,reagentSet in enumerate(LongestSyntheticPathways):
    ...     print("reagentSet " + str(ind) + " :")
    ...     reagentSet.printShortReagentSetInfo()
    ```  
    ```text
    reagentSet 0 :
    R2.2_0|R5.1_0 c1nn[nH:20]n1.Clc1ccccc1C([CH3:10])[OH:20].N[CH:10]=O Availability rate (% of mol. atoms coming from available synthons): 0.72
    reagentSet 1 :
    R2.2_0|R5.2_0 c1nn[nH:20]n1.Clc1ccccc1C([OH:20])[CH3:21].N[CH:10]=O Availability rate (% of mol. atoms coming from available synthons): 0.17
    reagentSet 2 :
    R10.1_0|R2.2_0 c1nnn(C[CH2:10][OH:20])n1.Clc1cccc[cH:20]1.N[CH:10]=O Availability rate (% of mol. atoms coming from available synthons): 0.17 
    reagentSet 3 : 
    R10.1_1|R2.2_0 Clc1ccccc1[CH2:10][OH:20].c1nnn([CH3:20])n1.N[CH:10]=O Availability rate (% of mol. atoms coming from available synthons): 0.17
    ```
