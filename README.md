PIDGIN : Prediction IncluDinG INactivity
===========

Author : Lewis Mervin, lhm30@cam.ac.uk

Supervisor : Dr. A. Bender

Protein Target Prediction Tool trained on SARs from PubChem (Mined 08/04/14) and ChEMBL18

Molecular Descriptors : 2048bit Morgan Binary Fingerprints (Rdkit) - ECFP4

Dependencies : rdkit, sklearn, numpy


All rights reserved 2014
==========================================================================================

Dependencies: 

Requires Python 2.7, Scikit-learn [1], Numpy [2] and Rdkit [3] to be installed on system.

Follow these steps on Linux:
 
1. ```sudo apt-get install python-numpy```
2. ```sudo apt-get install python-sklearn```
3. ```sudo apt-get install python-rdkit librdkit1 rdkit-data```

==========================================================================================

Instructions:

IMPORTANT:
*	The program currently currently recognises line separated SMILES in .csv format
*	Molecules Should be standardized before running models
*	ChemAxon is recommended and free for academic use at (http://www.chemaxon.com)
*	Protocol used to standardise these molecules is provided: StandMoleProt.xml
*	Do not modify the 'models' name or directory 
*	Cytotoxicity_library.csv is included for use as example dataset for test data

Two different python scripts can be used when performing target prediction.
Both utilise the Naive Bayes models created using Scikit-learn [1]. 


1. ```predict.py filename.csv```
    This script outputs only the raw probabilities for the compounds in a matrix. 
    
    Example of how to run the code:

    ```
    python predict.py input.csv
    ```

2. ```predict_binary.py theshold filename.csv```
    This script generates binary predictions for the models after application of Class-Specific Activity Cut-off Thresholds.
    
    This script requires an argument for the choice of Class-Specific Cut-off Thresholds.
    
    The options available are Precision (p), Recall (r), F1-score (f), Accuracy (a) and the Default 0.5 global threshold (0.5)
    
    Example of how to run the code:
    ```
    python predict_binary.py r input.csv
    ```
    
    where r would apply Thresholds calculated using the Recall metric
 
 [1] http://scikit-learn.org/stable/
 [2] http://www.numpy.org
 [3] http://www.rdkit.org
