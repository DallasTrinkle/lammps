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

# Utility rule file for fix-homepage.

# Include the progress variables for this target.
include CMakeFiles/fix-homepage.dir/progress.make

CMakeFiles/fix-homepage:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/ngoccuongnguyen/lammps/doc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Fix homepage URL errors"
	cd /Users/ngoccuongnguyen/lammps && /usr/local/Frameworks/Python.framework/Versions/3.9/bin/python3.9 /Users/ngoccuongnguyen/lammps/tools/coding_standard/homepage.py -f .

fix-homepage: CMakeFiles/fix-homepage
fix-homepage: CMakeFiles/fix-homepage.dir/build.make

.PHONY : fix-homepage

# Rule to build all files generated by this target.
CMakeFiles/fix-homepage.dir/build: fix-homepage

.PHONY : CMakeFiles/fix-homepage.dir/build

CMakeFiles/fix-homepage.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fix-homepage.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fix-homepage.dir/clean

CMakeFiles/fix-homepage.dir/depend:
	cd /Users/ngoccuongnguyen/lammps/doc && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ngoccuongnguyen/lammps/cmake /Users/ngoccuongnguyen/lammps/cmake /Users/ngoccuongnguyen/lammps/doc /Users/ngoccuongnguyen/lammps/doc /Users/ngoccuongnguyen/lammps/doc/CMakeFiles/fix-homepage.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fix-homepage.dir/depend

