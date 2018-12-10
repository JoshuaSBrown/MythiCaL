
# KMC Coursegraining Library Developers Guide

Welcome! Here is some basic documentation to help guide potential developers on the standards used to guide the development of this library. 

1. [Compiling](#compiling)
2. [Header Files](#header-files)
3. [Functions](#functions)
4. [Classes](#classes)

## Compiling

To compile with unit testing enable you will need to use the `-DENABLE_TESTING=ON` flag when calling cmake. Once compiled you can run the tests with `make test`. If you want to continue and install the library if it is working you can call `suod make install`. The build directory should be made in the source directory which should have the name `CourseGrainSites` if cloned from github. 

```
cd CourseGrainSites
mkdir build
cd build
cmake -DENABLE_TESTING=ON ../
make
make test
sudo make install
```

If one of the tests fails you can run the unit test by itself to see more detailed output. For instance lets assume that `unit_test_site` has failed. 

```
cd build/src/tests
./unit_test_site
```

## Header Files

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

## Variables
 
 * Unlike function names, the variable names should all be lower case and if 
   they are composed of more than one word they should be separated by an
   underscore.
