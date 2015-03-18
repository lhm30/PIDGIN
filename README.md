PIDGIN : Prediction IncluDinG INactivity
===========
[![DOI](https://zenodo.org/badge/10824/lhm30/PIDGIN.svg)](http://dx.doi.org/10.5281/zenodo.15984)


Author : Lewis Mervin, lhm30@cam.ac.uk

Supervisor : Dr. A. Bender

Protein Target Prediction Tool trained on SARs from PubChem (Mined 08/04/14) and ChEMBL18
![](https://pubchem.ncbi.nlm.nih.gov/images/pubchemlogob.gif) ![](http://upload.wikimedia.org/wikipedia/commons/a/a1/Chembl_logo.png)

Molecular Descriptors : 2048bit Morgan Binary Fingerprints (Rdkit) - ECFP4

* Targets: 1,080 Human Targets (activities over 10)
* Bioactivities: 194,563,469
* Distinct compounds: 826,354 (Based on binary fingerprints)
* Distinct Actives:	295,973
* Distinct Inactives:	648,930
* Number of classes requiring sphere exclusion:	480 (Ratio of less than 100:1 Inactive:Active)
* Number of classes requiring under-sampling:	600 (Ratio of less than 100:1 Inactive:Active)


Dependencies : rdkit, sklearn, numpy

![](http://www.rdkit.org/Images/logo.png) ![](http://scikit-learn.org/stable/_static/scikit-learn-logo-small.png) ![](http://upload.wikimedia.org/wikipedia/ru/c/cc/Numpylogo.png)

ChemAxon Standardizer was used for structure canonicalization and transformation, JChem 6.0.2.

![](http://www.chemaxon.com/images/powered_100px.gif)  http://www.chemaxon.com

![](http://www.ebi.ac.uk/sites/ebi.ac.uk/files/field/image/logo_5.gif)


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
*	The program currently recognises line separated SMILES in .csv format
*	Molecules Should be standardized before running models
*	ChemAxon Standardizer should be used for structure canonicalization and is free for academic use at (http://www.chemaxon.com)
*	Protocol used to standardise these molecules is provided: StandMoleProt.xml
*	Do not modify the 'models' name or directory 
*	Cytotoxicity_library.csv is included for use as example dataset for test data

Two different python scripts can be used when performing target prediction.
Both utilise the Naive Bayes models created using Scikit-learn [1]. 


1. ```predict_raw.py filename.csv```
    This script outputs the raw Naive Bayes score probabilities for the compounds in a matrix. 
    
    Example of how to run the code:

    ```
    python predict_raw.py input.csv
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

3. ```predict_binary_heat.py threshold filename.csv```
    This script generates binary predictions for the models using the thresholds, but does not output predictions for targets with no hits (this reduces the number of classes that one attempts to plot on a heat map)
    
    Example of how to run the code:

    ```
    python predict_binary_heat.py r input.csv
    ```

4. ```predict_single.py filename.csv```
    This script generates outputs for the raw Naive Bayes score probabilities and the binary output for all the thresholds for a single compound (used to give more in depth analysis for single compounds)
    
    Example of how to run the code:

    ```
    python predict_single.py input.csv
    ```
    

5. ```predict_ranked.py filename.csv```
    This script generates ranked targets (based on the raw Naive Bayes score probabilities) for compounds. The final column includes the average ranking position for target classes given the input compounds.
    
    Example of how to run the code:

    ```
    python predict_ranked.py input.csv
    ```
    
6. ```predict_enriched_targets.py threshold filename.csv```
    This script enriched targets for a library of compounds, when compared to a background sample from PubChem. The script will ask for login details for MySQL on Calculon, and the number of background samples to compare the library.
    
    Example of how to run the code:

    ```
    python predict_enriched_targets.py threshold input.csv
    ```
    
7. ```predict_enriched_two_libraries.py threshold input_active_library.csv input_inactive_library.csv```
    This script enriched targets for a library of phenotypically active compounds compared to phenotypically inactive compounds.
    
    Example of how to run the code:

    ```
    python predict_enriched_two_libraries.py threshold filename_1.csv filename_2.csv
    ```
    
    
===========
PIDGIN now offers predictions excluding inactivity.
===========

* The predictions produced for this activity-only model reflect the probability that a compound is active for a given target class, when considering the probability of activity for the other activity classes. 

* The model and scripts to use this single model are found in the 'singlemodel' folder.

5. ```predict_singlemodel.py filename.csv```
    This script generates raw Naive Bayes score probabilities for compounds (this only considers activity information).
    
    Example of how to run the code:

    ```
    python predict_singlemodel.py input.csv
    ```

6. ```predict_singlemodel_ranked_number.py filename.csv```
    This script outputs the ranking for targets from compound libraries based on the raw Naive Bayes score probabilities. The final column includes the average ranking position for target classes given the input compounds (this only considers activity information).
    
    Example of how to run the code:

    ```
    python predict_singlemodel_ranked_number.py input.csv
    ```

==========================================================================================

 [1] http://scikit-learn.org/stable/
 [2] http://www.numpy.org
 [3] http://www.rdkit.org
