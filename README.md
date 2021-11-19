
![alt text](./logo_finsurf.png?raw=true "FINSURF")
## Introduction

FINSURF (Functional Identification of Non-coding Sequences Using Random Forests) is a tool designed to analyse lists of sequences variants in the human genome. 

It assigns a score to each variant, reflecting its functional importance and therefore its likelihood to disrupt the physiology of its carrier. FINSURF scores Single Nucleotide Variants (SNV), insertions and deletions. Among SNVs, transitions and transversions are treated separately. 
Insertions are characterised by a score given to each base flanking the insertion point. Deletions are characterised by a score at every deleted base. FINSURF can (optionally) use a list of known or suspected disease genes, in order to restrict results to variants overlapping cis-regulatory elements linked to these genes. 

For a variant of interest, users can generate a graphical representation of "feature contributions », showing the relative contributions of genomic, functional or evolutionary information to its score.



FINSURF is implemented as python3 scripts.

## License

This code may be freely distributed and modified under the terms of the GNU General Public License version 3 (GPL v3)
and the CeCILL licence version 2 of the CNRS. These licences are contained in the files:

1. [LICENSE-GPL.txt](LICENSE-GPL.txt) (or on [www.gnu.org](https://www.gnu.org/licenses/gpl-3.0-standalone.html))
2. [LICENCE-CeCILL.txt](LICENCE-CeCILL.txt) (or on [www.cecill.info](https://cecill.info/licences/Licence_CeCILL_V2-en.html))

Copyright for this code is held by the [Dyogen](http://www.ibens.ens.fr/?rubrique43) (DYnamic and Organisation of GENomes) team
of the Institut de Biologie de l'Ecole Normale Supérieure ([IBENS](http://www.ibens.ens.fr)) 46 rue d'Ulm Paris and the individual authors.

- Copyright © 2020 IBENS/Dyogen : **Lambert MOYON**, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Camille Berthelot and Hugues ROEST CROLLIUS

## Contact

Email finsurf {at} bio {dot} ens {dot} psl {dot} eu

*If you use FINSURF, please cite:*

Classification of non-coding variants with high pathogenic impact
Lambert Moyon, Camille Berthelot, Alexandra Louis, Nga Thi Thuy Nguyen, Hugues Roest Crollius
bioRxiv 2021.05.03.442347; doi: https://doi.org/10.1101/2021.05.03.442347


# Content

This repository contains all the code and notebooks related to the training, evaluation, and application
of the FINSURF-adjusted prediction model.

In `code_train/FINSURF_train/` you will find Python scripts and libraries that were developped in the
context of this project.
For illustration, configuration files required by some of these scripts are present in the
`configs_train` directory.

In `notebooks_train` you will find a set of Jupyter notebooks which show some of the experiments and
associated visualization.
Actually most of them contain figures that are presented in the papers :

- `FINSURF_model-creation_and_cross-performance.ipynb` contains the figures related to the training of
  FINSURF models following different sampling-schemes for negative controls, and how this sampling
  affect cross-performance of the models (supplementary figure 1 and main figure 2a,b,c)

- `FINSURF-adjusted_comparison-methods_kfold_and_genomiser-vs-clinvar.ipynb` contains the figures
  related to the evaluation of FINSURF-adjusted model on independent datasets of variants, as
  presented in main figure 2d and supplementary figure 3

- `models_analysis_HGMD_feature-contribs-profiles-training-set.ipynb` contains the exploration of the
  positive controls from the training-set of FINSURF-adjusted using the feature-contribution method,
  that allowed to identified different functional profiles of regulatory variants, as presented in
  main figure 03 and in supplementary figure 5

- `non-HGMD-Genomizer_ranks_diseases_visualization.ipynb` shows the results of the application of
  FINSURF-adjusted in simulated "real-case scenarios", as presented in main figure 4 and supplementary
  figure 6


In addition to these notebooks, a notebook named `FINSURF-adjusted_DEMO-COMPLETE_CV_TRAINING_and_functional-profiles-examples` can be found.

This notebook is a fully executable notebook that enables one to explore the FINSURF-adjusted trained model, its cross-validation results, as well as the method to extract feature contributions from the trained model.

Besides, the model contains the code and plots associated to the main figure 3c and d where we present two examples of non-coding variants, scored by FINSURF, and for which the functional profiles (combination of feature contributions and scaled features) were extracted.

To allow the execution of this notebook, we added a `data/` directory, which contain the following folders:

- `example_functional_profiles/` : all data related to the two example variants.
- `regulatory_regions/` : the table of FINSURF Regulatory Regions.
- `trained_model/` : all data related to the training and cross-validation of the model.
- `variants/` : metadata related to variant annotation.

Please note that the training data for the FINSURF models is under licensed access.
We cannot publish the genomic coordinates or identifiers of variants in the database, as this violates the terms of the licensing agreement.
Hence the data table available only contains an internal ID `row_id` and numeric annotations used by the model.
