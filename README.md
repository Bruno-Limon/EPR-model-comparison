<strong>Exploration and Preferential Return (EPR)</strong>, a generative human mobility model and comparison of some of its variants.

Conducted and developed by:
- Bruno Limon

Under the guidance of professors: 
- Luca Pappalardo 
- Mirco Nanni
@ Università di Pisa, Italy
Geospatial Analytics Course

<!-- --------------------------------------------------------------------------------------- -->
<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#Introduction">Introduction</a></li>
    <li><a href="#Method">Method</a></li>
    <li><a href="#Results">Results</a></li>
    <li><a href="#Contact">Contact</a></li>
    <li><a href="#Acknowledgements">Acknowledgements</a></li>
  </ol>
</details>

<!-- --------------------------------------------------------------------------------------- -->
## Introduction

This project aims to present a systematic comparison of EPR models to produce synthetic trajectories of human mobility, reproducing the movements of an individual in a realistic manner. 3 variants of this model are presented:

- Spatial EPR: Variant where individuals are constrained to move in a defined geographical space, given by a certain tessellation.
- Density EPR: Improving on the basic EPR model by using collective information and the gravity model to influence the decisions of an individual. The model thus exploits information about the relevance of locations within the given space, individuals are inclined to explore new places with a force depending on the relevance of such places at a collective level.
- Ditras: EPR model that acts in 2 steps, first, it produces a mobility diary based on a Markov Chains Model, this captures the temporal aspects of an individual's trajectories, and then using the explorations mechanisms of the density EPR to produce a trajectory that is now able to capture both the temporal as well as the spatial aspects of human mobility.

The goal is to observe and compare the statistical distributions of measures of realness for the synthetic trajectories obtained with these models. These measures are:

M1. Travel Distance
M2. Radius of Gyration
M3. Uncorrelated Entropy
M4. Distinct Visited Locations
M5. Visits per Location
M6. Location Frequency

<!-- --------------------------------------------------------------------------------------- -->
## Method

To test the capabilities of the previously presented models, a systematic comparison is made by using them to build synthetic trajectories in 4 different geographical areas, each with 4 different spatial tessellations. The areas to analyze are:

1. New York State, United States
2. San Francisco, California, United States
3. Houston, Texas, United States
4. Mexico City, Mexico

While the tessellations are as follows:

1. Squared
2. Hexagonal
3. Official, representing administrative divisions such as neighborhoods or counties
4. Voronoi, build from points of interest within the geo area, such as health care facilities

Furthermore, a set of real trajectory datasets has been obtained for each area, with social media check-ins for area 1 and 3 and taxi cab gps trajectories for areas 2 and 4.

Then, 3 models, namely S-EPR, D-EPR and ditras are built and used to generate a set of trajectories for each pair of areas and tessellations, which will then be compared both quantitatively by the distributions of their measures and their root mean squared error as well as qualitatively through the use of plotting. 
Throughout this project, the notation 'a1[]_t[]' is used to identify the area and its tessellation, for example 'a1_t1' would correspond to New York State with squared tessellation.

<!-- --------------------------------------------------------------------------------------- -->
## Results

Overall, all three models perform relative well, with measures such as M1 and M6 with good results all-around, regardless of area, tessellation or model, while for some other measures, such as M4 and M5, the synthetic distributions greatly differ from those observed on the real trajectories. However, Ditras and D-EPR seem to perform better, as their distributions at least somewhat resemble the power-law distribution expected of measures such as M3 for example, whereas S-EPR produces a completely different type of distribution.
All in all, as expected, the S-EPR model and its trajectories seem to be the worst performing due to their lack of spatio-temporal awareness, which D-EPR and Ditras seem to solve.

<!-- --------------------------------------------------------------------------------------- -->
## Contact

Bruno Limon 
- bruno.limon.avila@gmail.com

<!-- --------------------------------------------------------------------------------------- -->
## Acknowledgements

- <a href="https://scikit-mobility.github.io/scikit-mobility/"> Scikit-mobility</a>

<img alt="License" src="https://img.shields.io/badge/license-MIT-brightgreen"> <img alt="Python" src="https://img.shields.io/badge/python->=3.8-blue">

