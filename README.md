# InvivoBrainHSI-Benchmark

## Overview
This GitHub repository contains the MATLAB code of the research paper titled "Hyperspectral Imaging Benchmark based on Machine Learning for Intraoperative Brain Tumour Detection" published in npj Precision Oncology in 2023. The paper presents a robust k-fold cross-validation approach, that hyperspectral imaging combined with the proposed processing framework is a promising intraoperative tool for in-vivo identification and delineation of brain tumours, including both primary (high-grade and low-grade) and secondary tumours. Analysis of the in-vivo brain database, consisting of 61 HS images from 34 different patients, achieve a highest median macro F1-Score result of 70.2±7.9% on the test set using both spectral and spatial information. Here, we provide a benchmark based on machine learning for further developments in the field of in-vivo brain tumour detection and delineation using hyperspectral imaging to be used as a real-time decision support tool during neurosurgical workflows.

The code provided here is intended to replicate the experiments and results presented in the paper. This README file will guide you through the repository's structure, how to use the code, and how to reproduce the results. Please refer to the original paper for a detailed explanation of the research.

## Installation
```git clone https://github.com/HIRIS-Lab/InvivoBrainHSI-Benchmark.git```

MATLAB script requires:
   - Statistics and Machine Learning Toolbox

## Usage

Before using the code, make sure you have set up the project as described in the Installation section. After that, you can run the code with the following command:

``` matlab [classification, pca_r, knn, clusterMap, mv, tmdMap ] = spectralSpatialMV(model, hsCube)```

This command will execute the main code, which replicates the experiments presented in the paper. You may need to configure certain parameters or settings within the code to match the specifics of your experiments or dataset.

## Dataset
The dataset used in the paper can be found at [HSI Brain Database](https://hsibraindatabase.iuma.ulpgc.es). Follow the instructions provided there to download and prepare the dataset for use with this code.

## Code Structure
The project structure is organized as follows:

* data/: This directory contains input data.
* models/: Contains an example of pre-trained models. User-trained models must be inserted here. 
* utils/: Utility functions and helper scripts.
* spectralSpatialMV.m: Code to run the spectral-spatial framework.
* main.m: The main script to run the experiments and plot the results.

## Citation
If you use this code or the results obtained from it in your research, please cite the original paper:

Leon, R., Fabelo, H., Ortega, S. et al. Hyperspectral imaging benchmark based on machine learning for intraoperative brain tumour detection. npj Precis. Onc. 7, 119 (2023). https://doi.org/10.1038/s41698-023-00475-9

## License

Copyright 2023 Himar Fabelo, Raquel Leon and Samuel Ortega

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
