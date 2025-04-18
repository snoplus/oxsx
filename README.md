<h1> OXSX </h1>
Signal Extraction framework for the SNO+ experiment

<h2> Dependencies </h2>
1. GCC compiler capable of compiling C++ code to the C++17 standard

2. [Armadillo](http://arma.sourceforge.net/) is a linear algebra package used for quick matrix multiplication

3. [GSL](https://gcc.gnu.org/libstdc++/) - likely you already have this installed, especially if you are running RAT

4. [SCons](http://www.scons.org/) Is used for the build, also a dependency for RAT

5. [HDF5](https://www.hdfgroup.org/HDF5/release/obtain5.html) Should be configured to install with c++ support `./configure --enable-cxx && make && make install`

6. [ROOT](https://root.cern.ch/downloading-root) Should be installed with Minuit2 enabled `./configure --enable-minuit2`

7. [Catch2](https://github.com/catchorg/Catch2) Version 3.0.0+ is needed. Install with:
   ```
   git clone git@github.com:catchorg/Catch2.git
   cd Catch2/
   git checkout v3.7.0 // Or whichever tag you want (must be >= 3.0.0)
   cmake -B build -S . -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=/path/to/Catch2/
   cmake --build build/ --target install
   ```
   (Of course, if you have permissions to write to default locations you don't necessarily need to invoke the `DCMAKE_INSTALL_PREFIX` option)

<h2>Installation Instructions </h2>
Follow the installation instructions for each of the above using either the default install location or a different directory if you would prefer. Be careful to start the install with a clean environment.

1. Clone this repository with `git clone https://github.com/snoplus/oxsx.git`

2. If your dependencies are somewhere the compiler can't find them, copy `config/userconfig.ini.template` to `config/userconfig.ini` and add the relevant paths. Missing entries are assumed to be in standard locations. e.g.

   ```
   [root]
   header_path : <path/to/headers>
   lib_path    : <path/to/libraries>
   ```

3. Run `scons && scons units`: this will compile the OXSX library and subsequently the unit tests.

4. Test the build was sucessful with `./test/RunUnits`

<h2> Alternative Installation Instructions Using Singularity/Apptainer & CMake </h2>
OXSX now comes with the ability to compile the repository via the build system CMake, 
and the definition file needed to create a container which will contain all of the necessary 
external repositories.

1. Clone this repository with `git clone https://github.com/snoplus/oxsx.git`

2. Have either [Docker](https://www.docker.com/), [Singularity](https://sylabs.io/singularity/), or [Apptainer](https://apptainer.org/) installed on your system. These are programs which allow you to build Containers on your system.

3. Navigate into the `oxsx` repository, and create an OXSX container `oxsx_container.sif` with the following command (this is for Apptainer; very similar commands for what follows are used for Docker/Singularity):
```
apptainer build oxsx_container.sif oxsx_container.def
```

4. Open the container:
```
apptainer shell oxsx_container.sif
```

5. Build the repository, using CMake:
```
cmake -S . -B cmake-build-debug
cmake --build cmake-build-debug
```
This will create a new build directory, `cmake-build-debug`, as part of the build process. Feel free to use a different name, such as `build` - though that may clash with any existing `build` directory if you've also compiled OXSX with Sconscript. This build procedure generates the `oxsx` library, compiles all of the code in `example/` and builds all of the unit tests within `test/`.

6. Test the build was successful with `./cmake-build-debug/test/unit/RunUnits`

<h2> Compiling Your Own Scripts</h2>

scons auto-generates a script that compiles and links your c++ against the source code and dependencies just run `. <oxsx root>/bin/compile.sh <your cpp file>` to produce an executible of the same name

Look in `<oxsx root>/examples` for help getting started

<h2> Creating ntuples </h2>
One way to read in data is using a ROOT ntuple. If you are looking at SNO+ data, you probably have a flat root tree that can be easily pruned down into this format with only the branches you are interested in.

To create ntuples for oxsx run `./util/prune_flat_tree <path_to_file_to_convert> -treename <name_of_the_tree> <branch1> <branch2> <branch3> -newtreename <name_of_your_tree> -outfile <name_of_your_file> -nevents <number_of_events>`

- The name of the tree in an input file is optional, as a default it is "output"
- The name of the output file is optional, as a default is is <the_name_of_input_file>+"\_oxsx.root"
- The name of the tree in an output file is optional, as a default it is "oxsx"
- The number of events of an output file is optional
