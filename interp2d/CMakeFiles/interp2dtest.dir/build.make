# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_SOURCE_DIR = /home/unkokusei/prog/nfit12.16/interp2d

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/unkokusei/prog/nfit12.16/interp2d

# Include any dependencies generated for this target.
include CMakeFiles/interp2dtest.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/interp2dtest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/interp2dtest.dir/flags.make

CMakeFiles/interp2dtest.dir/test.c.o: CMakeFiles/interp2dtest.dir/flags.make
CMakeFiles/interp2dtest.dir/test.c.o: test.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/unkokusei/prog/nfit12.16/interp2d/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/interp2dtest.dir/test.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/interp2dtest.dir/test.c.o   -c /home/unkokusei/prog/nfit12.16/interp2d/test.c

CMakeFiles/interp2dtest.dir/test.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/interp2dtest.dir/test.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/unkokusei/prog/nfit12.16/interp2d/test.c > CMakeFiles/interp2dtest.dir/test.c.i

CMakeFiles/interp2dtest.dir/test.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/interp2dtest.dir/test.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/unkokusei/prog/nfit12.16/interp2d/test.c -o CMakeFiles/interp2dtest.dir/test.c.s

CMakeFiles/interp2dtest.dir/test.c.o.requires:
.PHONY : CMakeFiles/interp2dtest.dir/test.c.o.requires

CMakeFiles/interp2dtest.dir/test.c.o.provides: CMakeFiles/interp2dtest.dir/test.c.o.requires
	$(MAKE) -f CMakeFiles/interp2dtest.dir/build.make CMakeFiles/interp2dtest.dir/test.c.o.provides.build
.PHONY : CMakeFiles/interp2dtest.dir/test.c.o.provides

CMakeFiles/interp2dtest.dir/test.c.o.provides.build: CMakeFiles/interp2dtest.dir/test.c.o

# Object files for target interp2dtest
interp2dtest_OBJECTS = \
"CMakeFiles/interp2dtest.dir/test.c.o"

# External object files for target interp2dtest
interp2dtest_EXTERNAL_OBJECTS =

interp2dtest: CMakeFiles/interp2dtest.dir/test.c.o
interp2dtest: CMakeFiles/interp2dtest.dir/build.make
interp2dtest: libinterp2d.a
interp2dtest: CMakeFiles/interp2dtest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable interp2dtest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/interp2dtest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/interp2dtest.dir/build: interp2dtest
.PHONY : CMakeFiles/interp2dtest.dir/build

CMakeFiles/interp2dtest.dir/requires: CMakeFiles/interp2dtest.dir/test.c.o.requires
.PHONY : CMakeFiles/interp2dtest.dir/requires

CMakeFiles/interp2dtest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/interp2dtest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/interp2dtest.dir/clean

CMakeFiles/interp2dtest.dir/depend:
	cd /home/unkokusei/prog/nfit12.16/interp2d && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/unkokusei/prog/nfit12.16/interp2d /home/unkokusei/prog/nfit12.16/interp2d /home/unkokusei/prog/nfit12.16/interp2d /home/unkokusei/prog/nfit12.16/interp2d /home/unkokusei/prog/nfit12.16/interp2d/CMakeFiles/interp2dtest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/interp2dtest.dir/depend

