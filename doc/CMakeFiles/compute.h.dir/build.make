# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.18.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.18.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ngoccuongnguyen/lammps/cmake

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ngoccuongnguyen/lammps/doc

# Utility rule file for compute.h.

# Include the progress variables for this target.
include CMakeFiles/compute.h.dir/progress.make

CMakeFiles/compute.h: includes/lammps/compute.h


includes/lammps/compute.h: /Users/ngoccuongnguyen/lammps/src/compute.h
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/ngoccuongnguyen/lammps/doc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating includes/lammps/compute.h"
	/usr/local/Cellar/cmake/3.18.2/bin/cmake -E copy_if_different /Users/ngoccuongnguyen/lammps/src/compute.h /Users/ngoccuongnguyen/lammps/doc/includes/lammps/compute.h

compute.h: CMakeFiles/compute.h
compute.h: includes/lammps/compute.h
compute.h: CMakeFiles/compute.h.dir/build.make

.PHONY : compute.h

# Rule to build all files generated by this target.
CMakeFiles/compute.h.dir/build: compute.h

.PHONY : CMakeFiles/compute.h.dir/build

CMakeFiles/compute.h.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/compute.h.dir/cmake_clean.cmake
.PHONY : CMakeFiles/compute.h.dir/clean

CMakeFiles/compute.h.dir/depend:
	cd /Users/ngoccuongnguyen/lammps/doc && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ngoccuongnguyen/lammps/cmake /Users/ngoccuongnguyen/lammps/cmake /Users/ngoccuongnguyen/lammps/doc /Users/ngoccuongnguyen/lammps/doc /Users/ngoccuongnguyen/lammps/doc/CMakeFiles/compute.h.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/compute.h.dir/depend

