# Python Implementation of Binary Tree Partition Method for Joint Probabilistic Data Association (JPDA)

## Introduction

Rezatofighi et al. introduced a method to generate the *M* best solutions to the binary linear program used for associating detections to tracks. Their work is proposed in *Joint Probabilistic Data Association Revisited* (https://ieeexplore.ieee.org/document/7410706). 

## MATLAB Implementation

The authors provide a MATLAB implementation of the algorithm at http://www.milanton.de/#publications under *2015*. 

## Python Version 

This project aims to replicate the code hosted at the site above in Python. 

### Hyperparameters

The linear objective features costs for linking established tracks to new detections as well as a cost for linking existing tracks to no detection. The original authors use a Normal probabilistic model to determine these costs. The proposed implementation uses a hyperparameter *rad* to specify this cost, although in practice in be variable for each track. 

### Functions 

Added functions *build_initial_constr* and *build_C* construct the inequality and equality constraints as well as the linear objective, respectively. Then *guobi_ilp*, *CalcNextBestSolution* and *MBest* are verbatim Python equivalents of the MATLAB functions. 

### Data Input

The motivating goal is to track nuclei in microscopy image volumes. The input are two NumPy arrays. The first contains the ground truth for the intial frame containing the established tracks. The second NumPy array contains detections from the image volumes. The BTP algorithm is used to associate tracks to detections with the goal of accurately updating tracks and establishing which nuclei are undetected in each volume. Each run of the main function *MBest* requires two arrays to determine the linear program objective and constraints as well as *M*, the number of best solutions to return. 



