# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/meduzen/Documents/vtm codes/cpp/6sem/svt_f"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/meduzen/Documents/vtm codes/cpp/6sem/svt_f/build"

# Include any dependencies generated for this target.
include CMakeFiles/sparskit.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/sparskit.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/sparskit.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/sparskit.dir/flags.make

CMakeFiles/sparskit.dir/main.f90.o: CMakeFiles/sparskit.dir/flags.make
CMakeFiles/sparskit.dir/main.f90.o: ../main.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/meduzen/Documents/vtm codes/cpp/6sem/svt_f/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/sparskit.dir/main.f90.o"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c "/home/meduzen/Documents/vtm codes/cpp/6sem/svt_f/main.f90" -o CMakeFiles/sparskit.dir/main.f90.o

CMakeFiles/sparskit.dir/main.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/sparskit.dir/main.f90.i"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E "/home/meduzen/Documents/vtm codes/cpp/6sem/svt_f/main.f90" > CMakeFiles/sparskit.dir/main.f90.i

CMakeFiles/sparskit.dir/main.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/sparskit.dir/main.f90.s"
	/usr/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S "/home/meduzen/Documents/vtm codes/cpp/6sem/svt_f/main.f90" -o CMakeFiles/sparskit.dir/main.f90.s

# Object files for target sparskit
sparskit_OBJECTS = \
"CMakeFiles/sparskit.dir/main.f90.o"

# External object files for target sparskit
sparskit_EXTERNAL_OBJECTS =

sparskit: CMakeFiles/sparskit.dir/main.f90.o
sparskit: CMakeFiles/sparskit.dir/build.make
sparskit: ../libskit.a
sparskit: CMakeFiles/sparskit.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/meduzen/Documents/vtm codes/cpp/6sem/svt_f/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable sparskit"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sparskit.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/sparskit.dir/build: sparskit
.PHONY : CMakeFiles/sparskit.dir/build

CMakeFiles/sparskit.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/sparskit.dir/cmake_clean.cmake
.PHONY : CMakeFiles/sparskit.dir/clean

CMakeFiles/sparskit.dir/depend:
	cd "/home/meduzen/Documents/vtm codes/cpp/6sem/svt_f/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/meduzen/Documents/vtm codes/cpp/6sem/svt_f" "/home/meduzen/Documents/vtm codes/cpp/6sem/svt_f" "/home/meduzen/Documents/vtm codes/cpp/6sem/svt_f/build" "/home/meduzen/Documents/vtm codes/cpp/6sem/svt_f/build" "/home/meduzen/Documents/vtm codes/cpp/6sem/svt_f/build/CMakeFiles/sparskit.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/sparskit.dir/depend

