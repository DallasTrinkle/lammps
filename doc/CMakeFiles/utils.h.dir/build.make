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

# Utility rule file for utils.h.

# Include the progress variables for this target.
include CMakeFiles/utils.h.dir/progress.make

CMakeFiles/utils.h: includes/lammps/utils.h


includes/lammps/utils.h: /Users/ngoccuongnguyen/lammps/src/utils.h
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/ngoccuongnguyen/lammps/doc/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating includes/lammps/utils.h"
	/usr/local/Cellar/cmake/3.18.2/bin/cmake -E copy_if_different /Users/ngoccuongnguyen/lammps/src/utils.h /Users/ngoccuongnguyen/lammps/doc/includes/lammps/utils.h

utils.h: CMakeFiles/utils.h
utils.h: includes/lammps/utils.h
utils.h: CMakeFiles/utils.h.dir/build.make

.PHONY : utils.h

# Rule to build all files generated by this target.
CMakeFiles/utils.h.dir/build: utils.h

.PHONY : CMakeFiles/utils.h.dir/build

CMakeFiles/utils.h.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/utils.h.dir/cmake_clean.cmake
.PHONY : CMakeFiles/utils.h.dir/clean

CMakeFiles/utils.h.dir/depend:
	cd /Users/ngoccuongnguyen/lammps/doc && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ngoccuongnguyen/lammps/cmake /Users/ngoccuongnguyen/lammps/cmake /Users/ngoccuongnguyen/lammps/doc /Users/ngoccuongnguyen/lammps/doc /Users/ngoccuongnguyen/lammps/doc/CMakeFiles/utils.h.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/utils.h.dir/depend

