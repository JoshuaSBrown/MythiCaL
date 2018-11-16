# Kinetic Monte Carlo Course Graining Library 
Library is designed to course grain kinetic Monte Carlo Simulations. Primary motivation for this library is to reduce the number of unnecessary compute cycles. 

## The Problem

Though kinetic Monte Carlo simulations are a very powerful tool they can at times suffer from significant performance hits. This can occur in particular when a relationship between two events dominates the dynamics. As an illustration consider an electron moving between molecules. If the rate of transfer from (molecule A -> molecule B) and (molecule B -> molecule A) is much greater than the rates to any other molecule C, D, ...etc a signifcant number of compute cycles will be used simply moving the charge back and forth between these two molecules. This can become such a resource intensive situation that 99% of the compute cycles are spent moving the charge between molecule A and B. Thus making simulations of relevant time scales infeasible. This library attempts to mitigate this problem by course graining out these rapid oscillations during the course of the simulation.  

## Example Applications

The library is designed to work with rates and could for instance be used with:
 * Monte Carlo Charge Transport Simulations
 * Tracking the spread of disease
 
## Dependencies

The library makes use of c++14 features so requires gcc 6 or a compiler with similar support. 

## Download
    
    git clone --recursive https://github.com/JoshuaSBrown/CourseGrainSites

## Installation 

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ../
    make 
    make install

## Developers

[Developers Guide](doc/DEVELOPERS_GUIDE.md)

## Author List

* Joshua Brown
* Riley Hadjis
