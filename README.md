# Kinetic Monte Carlo Coarse Graining Library 

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/142f5448ab0243acabe198e6632b8e84)](https://app.codacy.com/app/JoshuaSBrown/CoarseGrainSites?utm_source=github.com&utm_medium=referral&utm_content=JoshuaSBrown/CoarseGrainSites&utm_campaign=Badge_Grade_Dashboard)
[![Build Status](https://travis-ci.com/JoshuaSBrown/CoarseGrainSites.svg?branch=master)](https://travis-ci.com/JoshuaSBrown/CoarseGrainSites)

Library is designed to coarse grain kinetic Monte Carlo Simulations. Primary motivation for this library is to reduce the number of unnecessary compute cycles.

This library is still under development. 

## The Problem

Though kinetic Monte Carlo simulations are a very powerful tool they can at times suffer from significant performance hits. This can occur in particular when a relationship between two events dominates the dynamics. As an illustration consider an electron moving between molecules. If the rate of transfer from (molecule A -> molecule B) and (molecule B -> molecule A) is much greater than the rates to any other molecule C, D, ...etc a signifcant number of compute cycles will be wasted simply moving the charge back and forth between these two molecules. This can become such a resource intensive situation that 99% of the compute cycles are spent moving the charge between molecule A and B. Thus making simulations of relevant time scales infeasible. This library attempts to mitigate this problem by coarse graining out these rapid oscillations during the simulation.  

## Example Applications

The library is designed to work with rates and could for instance be used with:
 * Monte Carlo Charge Transport Simulations
 * Tracking the spread of disease
 * Modeling the circadian cycle
 
## Dependencies

The library makes use of c++14 features so requires gcc 6 or a compiler with similar support. 

## Download
    
    git clone --recursive https://github.com/JoshuaSBrown/CoarseGrainSites

## Installation 

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ../
    make 
    make install

## Tutorials and Further Documentation

[Wiki](https://github.com/JoshuaSBrown/CoarseGrainSites/wiki)

## Developers

[Developers Guide](CoarseGrainSites/doc/DEVELOPERS_GUIDE.md)

## Author List

* Joshua Brown
* Riley Hadjis
