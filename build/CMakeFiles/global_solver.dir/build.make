# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

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
CMAKE_COMMAND = /u/sw/toolchains/gcc-glibc/11.2.0/base/bin/cmake

# The command to remove a file.
RM = /u/sw/toolchains/gcc-glibc/11.2.0/base/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/nico/21-eikonal-21-eikonal

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/nico/21-eikonal-21-eikonal/build

# Include any dependencies generated for this target.
include CMakeFiles/global_solver.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/global_solver.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/global_solver.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/global_solver.dir/flags.make

CMakeFiles/global_solver.dir/src/main_global_solver.cpp.o: CMakeFiles/global_solver.dir/flags.make
CMakeFiles/global_solver.dir/src/main_global_solver.cpp.o: /home/nico/21-eikonal-21-eikonal/src/main_global_solver.cpp
CMakeFiles/global_solver.dir/src/main_global_solver.cpp.o: CMakeFiles/global_solver.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/nico/21-eikonal-21-eikonal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/global_solver.dir/src/main_global_solver.cpp.o"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/global_solver.dir/src/main_global_solver.cpp.o -MF CMakeFiles/global_solver.dir/src/main_global_solver.cpp.o.d -o CMakeFiles/global_solver.dir/src/main_global_solver.cpp.o -c /home/nico/21-eikonal-21-eikonal/src/main_global_solver.cpp

CMakeFiles/global_solver.dir/src/main_global_solver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/global_solver.dir/src/main_global_solver.cpp.i"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nico/21-eikonal-21-eikonal/src/main_global_solver.cpp > CMakeFiles/global_solver.dir/src/main_global_solver.cpp.i

CMakeFiles/global_solver.dir/src/main_global_solver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/global_solver.dir/src/main_global_solver.cpp.s"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nico/21-eikonal-21-eikonal/src/main_global_solver.cpp -o CMakeFiles/global_solver.dir/src/main_global_solver.cpp.s

CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.o: CMakeFiles/global_solver.dir/flags.make
CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.o: /home/nico/21-eikonal-21-eikonal/LocalProblem/DescentDirectionFactory.cpp
CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.o: CMakeFiles/global_solver.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/nico/21-eikonal-21-eikonal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.o"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.o -MF CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.o.d -o CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.o -c /home/nico/21-eikonal-21-eikonal/LocalProblem/DescentDirectionFactory.cpp

CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.i"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nico/21-eikonal-21-eikonal/LocalProblem/DescentDirectionFactory.cpp > CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.i

CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.s"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nico/21-eikonal-21-eikonal/LocalProblem/DescentDirectionFactory.cpp -o CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.s

CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.o: CMakeFiles/global_solver.dir/flags.make
CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.o: /home/nico/21-eikonal-21-eikonal/LocalProblem/DescentDirections.cpp
CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.o: CMakeFiles/global_solver.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/nico/21-eikonal-21-eikonal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.o"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.o -MF CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.o.d -o CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.o -c /home/nico/21-eikonal-21-eikonal/LocalProblem/DescentDirections.cpp

CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.i"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nico/21-eikonal-21-eikonal/LocalProblem/DescentDirections.cpp > CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.i

CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.s"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nico/21-eikonal-21-eikonal/LocalProblem/DescentDirections.cpp -o CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.s

CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.o: CMakeFiles/global_solver.dir/flags.make
CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.o: /home/nico/21-eikonal-21-eikonal/LocalProblem/LineSearchSolver.cpp
CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.o: CMakeFiles/global_solver.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/nico/21-eikonal-21-eikonal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.o"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.o -MF CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.o.d -o CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.o -c /home/nico/21-eikonal-21-eikonal/LocalProblem/LineSearchSolver.cpp

CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.i"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nico/21-eikonal-21-eikonal/LocalProblem/LineSearchSolver.cpp > CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.i

CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.s"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nico/21-eikonal-21-eikonal/LocalProblem/LineSearchSolver.cpp -o CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.s

CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.o: CMakeFiles/global_solver.dir/flags.make
CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.o: /home/nico/21-eikonal-21-eikonal/LocalProblem/solveEikonalLocalProblem.cpp
CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.o: CMakeFiles/global_solver.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/nico/21-eikonal-21-eikonal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.o"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.o -MF CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.o.d -o CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.o -c /home/nico/21-eikonal-21-eikonal/LocalProblem/solveEikonalLocalProblem.cpp

CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.i"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nico/21-eikonal-21-eikonal/LocalProblem/solveEikonalLocalProblem.cpp > CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.i

CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.s"
	/u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nico/21-eikonal-21-eikonal/LocalProblem/solveEikonalLocalProblem.cpp -o CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.s

# Object files for target global_solver
global_solver_OBJECTS = \
"CMakeFiles/global_solver.dir/src/main_global_solver.cpp.o" \
"CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.o" \
"CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.o" \
"CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.o" \
"CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.o"

# External object files for target global_solver
global_solver_EXTERNAL_OBJECTS =

global_solver: CMakeFiles/global_solver.dir/src/main_global_solver.cpp.o
global_solver: CMakeFiles/global_solver.dir/LocalProblem/DescentDirectionFactory.cpp.o
global_solver: CMakeFiles/global_solver.dir/LocalProblem/DescentDirections.cpp.o
global_solver: CMakeFiles/global_solver.dir/LocalProblem/LineSearchSolver.cpp.o
global_solver: CMakeFiles/global_solver.dir/LocalProblem/solveEikonalLocalProblem.cpp.o
global_solver: CMakeFiles/global_solver.dir/build.make
global_solver: CMakeFiles/global_solver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/nico/21-eikonal-21-eikonal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable global_solver"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/global_solver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/global_solver.dir/build: global_solver
.PHONY : CMakeFiles/global_solver.dir/build

CMakeFiles/global_solver.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/global_solver.dir/cmake_clean.cmake
.PHONY : CMakeFiles/global_solver.dir/clean

CMakeFiles/global_solver.dir/depend:
	cd /home/nico/21-eikonal-21-eikonal/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/nico/21-eikonal-21-eikonal /home/nico/21-eikonal-21-eikonal /home/nico/21-eikonal-21-eikonal/build /home/nico/21-eikonal-21-eikonal/build /home/nico/21-eikonal-21-eikonal/build/CMakeFiles/global_solver.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/global_solver.dir/depend
