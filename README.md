<h1> OXSX </h1>

Signal Extraction framework for the SNO+ experiment

<h2> Dependencies </h2>

1. GCC compiler capable of compiling C++ code to the C++17 standard

2. [Armadillo](http://arma.sourceforge.net/) is a linear algebra package used for quick matrix multiplication

3. [GSL](https://gcc.gnu.org/libstdc++/) - likely you already have this installed, especially if you are running RAT

4. [SCons](http://www.scons.org/) Is used for the build, also a dependency for RAT

5. [HDF5](https://www.hdfgroup.org/HDF5/release/obtain5.html) Should be configured to install with C++ support `./configure --enable-cxx && make && make install`

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

Currently there are two build systems supported: CMake and SCons. It is recommended you use CMake as we will likely stop supporting SCons in the future.

<h3>CMake</h3>

1. Clone this repository with `git clone https://github.com/snoplus/oxsx.git`

2. If any dependencies are installed in non-standard locations, add their install prefixes to ```CMAKE_PREFIX_PATH``` so that CMake can locate their ```<Package>Config.cmake``` files. For example:

   ```
   export CMAKE_PREFIX_PATH=/path/to/deps/install:$CMAKE_PREFIX_PATH
   ```

3. Configure and build:
   ```
   cmake -S . -B cmake-build
   cmake --build cmake-build
   ```   

4. You can check the build was successful by running the unit tests with the command:

   ```
   ./cmake-build/test/unit/RunUnits
   ```

<h3>SCons</h3>

1. Clone this repository with ```git clone https://github.com/snoplus/oxsx.git```

2. If your dependencies are somewhere the compiler can't find them, copy `config/userconfig.ini.template` to `config/userconfig.ini` and add the relevant paths. Missing entries are assumed to be in standard locations. e.g.

   ```
   [root]
   header_path : <path/to/headers>
   lib_path    : <path/to/libraries>
   ```

3. Run `scons && scons units`: this will compile the OXSX library and subsequently the unit tests.

4. Test the build was sucessful with `./test/RunUnits`

<h2> Alternative Installation Instructions using the Container </h2>

There is a container that includes OXSX and all the required dependencies. It is stored on the GitHub Container Registry (GHCR) attached to the [GitHub OXSX repository](https://github.com/snoplus/oxsx).

<h3>Using the Prebuilt Image</h3>
The following uses Apptainer, though you can use similar commands for using Docker, which is sometimes more convenient when running locally.
Firstly, setup a Github Personal Access Token (PAT) to be able to read packages. You can allow access to the registry with:

```
apptainer registry login -u <GitHub username> docker://ghcr.io
```

and enter your PAT as the password when prompted. Alternatively, your username and password (PAT) can be set as environment variables:

```
export APPTAINER_DOCKER_USERNAME=<GitHub username>
export APPTAINER_DOCKER_PASSWORD=<GitHub PAT>
```

You can then pull the container with:

```
apptainer pull oras://ghcr.io/snoplus/oxsx:<tag>
```

Replace <tag> with the desired OXSX tag, eg. 1.4.0 . Not all tags are available as containers; see [here](https://github.com/snoplus/oxsx/pkgs/container/oxsx) for the available container tags. If you do not have permission to view that page, you will need to be granted access on the GitHub repo by one of the admins. Note: the container is stored as an OCI artifact and must be pulled using the ```oras://``` prefix.

You can open the container with a command like:

```
apptainer shell oxsx_container.sif
```

The default build of OXO inside the container is fine if you plan on simply running code that uses OXO. However, if you would like to make active changes to OXO, you should have a clone of the repo outside the container, and when opening the container bind your oxsx directory to `/oxsx/` in the container. You can compile your own modified version within the container by using CMake (similarly to building outside the container):

```
cmake -S . -B cmake-build
cmake --build cmake-build
```

You can run the unit tests (which get compiled as part of the CMake build process) with the command:

```
./cmake-build/test/unit/RunUnits
```

Note: the container filesystem is read-only. To write outputs (including running unit tests), you must bind a writable directory when launching the container.

<h3>Building the Image</h3>

OXSX also comes with a definition file needed to create the container which will contain OXO and all of the necessary external repositories. Currently, you need a system you have sudo rights to for building the container yourself; you can always build the container locally and then copy the SIF file to the remote machine you are likely working on (presuming you don't have sudo rights there).

1. Clone this repository with `git clone https://github.com/snoplus/oxsx.git`

2. Have either [Docker](https://www.docker.com/), [Singularity](https://sylabs.io/singularity/), or [Apptainer](https://apptainer.org/) installed on your system. These are programs which allow you to build Containers on your system.

3. Navigate into the `oxsx` repository, and create an OXSX container `oxsx_container.sif` with the following command (this is for Apptainer; very similar commands for what follows are used for Docker/Singularity):

```
sudo apptainer build oxsx_container.sif oxsx_container.def
```

This will build OXO through CMake, and also build all of the required external repositories; this repo will be located inside the container at `/oxsx/`. A build sub-directory, `cmake-build`, is also made as part of the build process. This build procedure generates the `oxsx` library, compiles all of the code in `example/` and builds all of the unit tests within `test/`. The unit tests are run at the end of the container build process, to demonstrate the container's validity.

You can then use the newly built container in the same way as for the pre-built version.

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
