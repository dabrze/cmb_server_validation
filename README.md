# CheckMyBlob server validation

This repository hosts the validation results of the CheckMyBlob web server, as described in "Recognizing and validating ligands with CheckMyBlob" by Brzezinski et al. 

The code for the CheckMyBlob algorithm is hosted at https://github.com/dabrze/CheckMyBlob, and that code was used to train a classification model. The model was created (and cross-validated) on a dataset of 696,887 ligands and then tested on a holdout dataset of 17,150 ligands. These datasets are available at Zenodo ([10.5281/zenodo.4554473](https://doi.org/10.5281/zenodo.4554473)). The CheckMyBlob predictions for these to datasets are in the `classification_results` folder. The `validation_results.ipynb` notebook contains the code that was used to calculate evaluation metrics based on these predictions.

The `ligand_clustering` folder contains the code that was used to cluster ligands according to their chemical properties (as described in the paper). The results of this clustering that where used to group ligands into classes for the CheckMyBlob server are in the `ligand_clustering/clustering/results/` folder.
