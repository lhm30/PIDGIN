PIDGIN : Prediction IncluDinG INactivity
===========
[![DOI](https://zenodo.org/badge/10824/lhm30/PIDGIN.svg)](http://dx.doi.org/10.5281/zenodo.15984)


This is a legacy version of PIDGIN, version 2 found here: https://github.com/lhm30/PIDGINv2
==========================================================================================


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

Pathway information from NCBI BioSystems

![](http://www.ncbi.nlm.nih.gov/Structure/IMG/banner_graphics/biosystems_entrez3.png) ![](http://www.genome.jp/Fig/kegg128.gif) ![](http://biocyc.org/BioCyc.gif) ![](http://blog.openhelix.eu/wp-content/uploads/2011/01/Reactome_logo.jpg) ![](http://i.picresize.com/images/2015/04/29/oAE7h.png) ![](https://s-media-cache-ak0.pinimg.com/216x146/e3/71/2d/e3712dd81b80c17e24d4fb529f6bafab.jpg) ![](http://www.wikipathways.org/skins/common/images/earth-or-pathway_text3_beta.png)

Dependencies : rdkit, sklearn, numpy

![](http://www.rdkit.org/Images/logo.png) ![](http://scikit-learn.org/stable/_static/scikit-learn-logo-small.png) ![](http://upload.wikimedia.org/wikipedia/ru/c/cc/Numpylogo.png)

ChemAxon Standardizer was used for structure canonicalization and transformation, JChem 6.0.2.

![](http://www.chemaxon.com/images/powered_100px.gif)  http://www.chemaxon.com

![](https://dnasu.org/DNASU/image/Uniprot300.jpg)


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
    
6. ```predict_enriched.py threshold filename.csv```
    This script enriched targets and pathways for a library of compounds, when compared to a background sample from PubChem. The script will ask for login details for MySQL on Calculon, and the number of background samples to compare the library. The script will live access the NCBI BioSystems repository for up-to-date pathways information.
    
    Example of how to run the code:

    ```
    python predict_enriched.py threshold input.csv
    ```
    
7. ```predict_enriched_two_libraries.py threshold input_active_library.csv input_inactive_library.csv```
    This script calculates enriched targets and pathways for a library of phenotypically active compounds compared to phenotypically inactive compounds.
    
    Example of how to run the code:

    ```
    python predict_enriched_two_libraries.py threshold filename_1.csv filename_2.csv
    ```
    
8. ```predict_fingerprints.py threshold filename.csv```
    This script calculates target and pathway hits and represents them as binary fingerprints in a matrix.
    
    Example of how to run the code:

    ```
    python predict_fingerprints.py threshold input.csv
    ```
    
9. ```ad_analysis_all.py number_of_neighbours input.csv```
    This script takes an input of SMILES and calculates the average nearest-neighbour Tanimoto similarity of the n nearest compounds in the active training set, for each target. From here you can average the distance to a model for your input compounds, or perhaps all the models in PIDGIN etc.
    
    Example of how to run the code:

    ```
    python ad_analysis_all.py 1 input.csv
    ```
    
10. ```ad_analysis.py number_of_neighbours input.csv```
    Similar to ad_analysis_all.py, but can be used for more in-depth analysis for near-neighbours. This script takes the list of smiles and retrieves the top n nearest active compound SMILES and Tanimoto's, giving a more detailed view of the distribution of compounds and which compounds are active.
    
    Example of how to run the code:

    ```
    python ad_analysis.py 1 input.csv
    ```
    
    
===========
PIDGIN also offers predictions considering only activity.
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
