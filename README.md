# <Your-Project-Title>

## Description

Provide a short description explaining the what, why, and how of your project. Use the following questions as a guide:

- What was your motivation?
- Why did you build this project? (Note: the answer is not "Because it was a homework assignment.")
- What problem does it solve?
- What did you learn?

## Table of Contents (Optional)

If your README is long, add a table of contents to make it easy for users to find what they need.

- [Installation](#installation)
- [Usage](#usage)
- [Credits](#credits)
- [License](#license)

## Installation

What are the steps required to install your project? Provide a step-by-step description of how to get the development environment running.

## Usage

Provide instructions and examples for use. Include screenshots as needed.

To add a screenshot, create an `assets/images` folder in your repository and upload your screenshot to it. Then, using the relative filepath, add it to your README using the following syntax:

    ```md
    ![alt text](assets/images/screenshot.png)
    ```

## Credits

List your collaborators, if any, with links to their GitHub profiles.

If you used any third-party assets that require attribution, list the creators with links to their primary web presence in this section.

If you followed tutorials, include links to those here as well.

## License

The last section of a high-quality README file is the license. This lets other developers know what they can and cannot do with your project. If you need help choosing a license, refer to [https://choosealicense.com/](https://choosealicense.com/).

---

🏆 The previous sections are the bare minimum, and your project will ultimately determine the content of this document. You might also want to consider adding the following sections.

## Badges

![badmath](https://img.shields.io/github/languages/top/lernantino/badmath)

Badges aren't necessary, per se, but they demonstrate street cred. Badges let other developers know that you know what you're doing. Check out the badges hosted by [shields.io](https://shields.io/). You may not understand what they all represent now, but you will in time.

## Features

If your project has a lot of features, list them here.

## How to Contribute

If you created an application or package and would like other developers to contribute it, you can include guidelines for how to do so. The [Contributor Covenant](https://www.contributor-covenant.org/) is an industry standard, but you can always write your own if you'd prefer.

## Tests

Go the extra mile and write tests for your application. The


================================================================================================================
Exploration and Preferential Return (EPR), a generative human mobility model and comparison of some of its variants.

<!-- LATESTCOMMIT:START -->

<!-- LATESTCOMMIT:END -->

Conducted and developed by:
- Bruno Limón

Under the guidance of professors: 
- Luca Pappalardo 
- Mirco Nanni
@ Università di Pisa, Italy
Geospatial Analytics Course

================================================================================================================
INTRODUCTION
================================================================================================================

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

================================================================================================================
METHOD
================================================================================================================

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

================================================================================================================
RESULTS
================================================================================================================
Overall, all three models perform relative well, with measures such as M1 and M6 with good results all-around, regardless of area, tessellation or model, while for some other measures, such as M4 and M5, the synthetic distributions greatly differ from those observed on the real trajectories. However, Ditras and D-EPR seem to perform better, as their distributions at least somewhat resemble the power-law distribution expected of measures such as M3 for example, whereas S-EPR produces a completely different type of distribution.
All in all, as expected, the S-EPR model and its trajectories seem to be the worst performing due to their lack of spatio-temporal awareness, which D-EPR and Ditras seem to solve.
