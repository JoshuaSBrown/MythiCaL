# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/joshbro42867/Documents/OPV_CourseGrainSites

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/joshbro42867/Documents/OPV_CourseGrainSites/build

# Include any dependencies generated for this target.
include CMakeFiles/kmc_course_grain.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/kmc_course_grain.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/kmc_course_grain.dir/flags.make

# Object files for target kmc_course_grain
kmc_course_grain_OBJECTS =

# External object files for target kmc_course_grain
kmc_course_grain_EXTERNAL_OBJECTS =

libkmc_course_grain.so: CMakeFiles/kmc_course_grain.dir/build.make
libkmc_course_grain.so: CMakeFiles/kmc_course_grain.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/joshbro42867/Documents/OPV_CourseGrainSites/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX shared library libkmc_course_grain.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/kmc_course_grain.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/kmc_course_grain.dir/build: libkmc_course_grain.so

.PHONY : CMakeFiles/kmc_course_grain.dir/build

CMakeFiles/kmc_course_grain.dir/requires:

.PHONY : CMakeFiles/kmc_course_grain.dir/requires

CMakeFiles/kmc_course_grain.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/kmc_course_grain.dir/cmake_clean.cmake
.PHONY : CMakeFiles/kmc_course_grain.dir/clean

CMakeFiles/kmc_course_grain.dir/depend:
	cd /home/joshbro42867/Documents/OPV_CourseGrainSites/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/joshbro42867/Documents/OPV_CourseGrainSites /home/joshbro42867/Documents/OPV_CourseGrainSites /home/joshbro42867/Documents/OPV_CourseGrainSites/build /home/joshbro42867/Documents/OPV_CourseGrainSites/build /home/joshbro42867/Documents/OPV_CourseGrainSites/build/CMakeFiles/kmc_course_grain.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/kmc_course_grain.dir/depend

