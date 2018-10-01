
# KMC Coursegraining Library Developers Guide

Welcome! Here is some basic documentation to help guide potential developers on the standards used to guide the development of this library. 

## Header files

 * Header files that are placed in `include/kmccoursegrain/` if the contents are meant to be publicly accessible. A guiding principle would be to put as little content as possible in these files. As the simpler the public interface is the easier it will be for someone to take advantage of
 * All other header files are placed in the `src/libkmccoursegrain` folder and should therefore not be available to the public, this allows a little more flexibility durin the development, as external users should not be reliant on the interface. 
 * Every header file should have a guard of the following form, where NAMOVEOFFILE is replaced with the appropriate name of the file:
```
#ifndef KMCCOURSEGRAIN_NAMEOFTHEFILE_HPP
#define KMCCOURSEGRAIN_NAMEOFTHEFILE_HPP
:
code
:
#endif // KMCCOURSEGRAIN_NAMEOFTHEFILE_HPP

```

## Functions

 * Function names should start with small letters and follow Camelcase

## Classes

 * Classes should begin with a capital letter
 * Private variables and functions should end with an underscore
