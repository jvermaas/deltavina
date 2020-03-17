DeltaVina
======

A scoring function for rescoring protein-ligand binding affinity.

Project Organization
------------
    ├── LICENSE
    ├── README.md          <- README file.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── bin                <- Master script to run deltavina
    │   └── dvrf20.py
    |   |__ dvrf20-prep.py
    │
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   ├── pharma.py  
    │   │   ├── featureSASA.py 
    │   │   ├── featureVina.py 
    │   │   └── tests
    │   │
    │   └── models         <- Trained model rffit.rda and scripts to run model
    │       ├── rffit.rda
    │       ├── modelDV.py
    │       └── tests
    │
    └── examples           <- Examples


How to use
-------
Please see the DeltaVinaTutorial.pdf or visit http://www.nyu.edu/projects/yzhang/DeltaVina

How to really use
-------
If you are given a list of proteins and their ligands, the featurization of those protein/ligand pairs can all happen in parallel (`dvrf20-prep.py`), and the subsequent step of making a random forest model can just read in the precomputed features (`dvrf20.py`). Thus, through bash or some other means, you can run all of the featurization steps independently, and when these are all finished, then you'd run the second step. For a single desktop in bash, something like this would be useful:

```bash
echo "" > groupfile
for pdb in $pdblist
do
    for lig in $liglist
    do
        #"&" symbol means run the job in the background and proceed.
        dvrf20-prep.py -r $pdb -l $lig &
        echo "$pdb $lig" >> groupfile
    done
done
#The wait command will stall until all the backgrounded jobs finish.
wait
dvrf20.py -g groupfile
```

Reference
---------
Cheng Wang, and Yingkai Zhang, J. Comput. Chem. 2017, 38, 169-177.

[Improving Scoring-Docking-Screening Powers of Protein–Ligand Scoring Functions using Random Forest](http://onlinelibrary.wiley.com/doi/10.1002/jcc.24667/abstract)
    
