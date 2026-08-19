============
Installation
============

------------------
Obtaining the code
------------------

The latest development version of MRCPP can be found on the ``master``
branch on GitHub::

    $ git clone https://github.com/MRChemSoft/mrcpp.git

The released versions can be found from Git tags ``vX.Y.Z`` under the
``release/X.Y`` branches in the same repository, or a zip file can be
downloaded from `Zenodo <https://doi.org/10.5281/zenodo.3606670>`_.

By default, all dependencies will be **fetched** at configure time if they are
not already available.


-------------------
Build prerequisites
-------------------

- CMake-3.22 (or later)
- GNU-11.2, Clang-14.0 or IntelLLVM-2022.1 (or later) compilers (C++17 standard)

.. hint::
    We have collected the recommended modules for the different Norwegian HPC
    systems under ``tools/<machine>.env``. These files can be sourced in order
    to get a working environment on the respective machines, and may also serve
    as a guide for other HPC systems.
  

C++ dependencies
----------------

- Linear algebra: `Eigen-3.4  <https://gitlab.com/libeigen/eigen>`_

Eigen will be downloaded automatically at configure time by CMake,
but can also be linked manually by setting the variable::

    EIGEN3_DIR=<path_to_eigen3>/share/eigen3/cmake


-----------------
Building the code
-----------------

Configure
---------

The ``setup`` script will create a directory called ``<build-dir>`` and run
CMake. There are several options available for the setup, the most
important being:

``--cxx=<CXX>``
  C++ compiler [default: g++]
``--omp``
  Enable OpenMP parallelization [default: False]
``--mpi``
  Enable MPI parallelization [default: False]
``--enable-tests``
  Enable tests [default: True]
``--enable-examples``
  Enable tests [default: False]
``--type=<TYPE>``
  Set the CMake build type (debug, release, relwithdebinfo, minsizerel) [default: release]
``--prefix=<PATH>``
  Set the install path for make install
``-h --help``
  List all options

The code can be built with four levels of parallelization:

 - no parallelization
 - only shared memory (OpenMP)
 - only distributed memory (MPI)
 - hybrid OpenMP + MPI

.. note::
    In practice we recommend the **shared memory version** for running on your
    personal laptop/workstation, and the **hybrid version** for running on a
    HPC cluster. The serial and pure MPI versions are only useful for debugging.

The default build is *without* parallelization and using GNU compilers::

    $ ./setup --prefix=<install-dir> <build-dir>

To use clang compilers you need to specify the ``--cxx`` option::

    $ ./setup --prefix=<install-dir> --cxx=clang++ <build-dir>

To build the code with shared memory (OpenMP) parallelization,
add the ``--omp`` option::

    $ ./setup --prefix=<install-dir> --omp <build-dir>

To build the code with distributed memory (MPI) parallelization, add the
``--mpi`` option *and* change to the respective MPI compilers (``--cxx=mpicxx``
for GNU)::

    $ ./setup --prefix=<install-dir> --omp --mpi --cxx=mpicxx <build-dir>

.. note::
    If you compile the MRCPP library manually as a separate project, the level
    of parallelization **must be the same** for MRCPP and MRChem. Similar
    options apply for the MRCPP setup, see
    `mrcpp.readthedocs.io <https://mrcpp.readthedocs.io/en/latest/>`_.


Build
-----

If the CMake configuration is successful, the code is compiled with::

    $ cd <build-dir>
    $ make


Test
----

A test suite is provided to make sure that everything compiled properly.
To run a collection of small tests::

    $ cd <build-dir>
    $ ctest


Install
-------

After the build has been verified with the test suite, it can be installed with
the following command::

    $ cd <build-dir>
    $ make install

Now libraries, headers and CMake configuration files can be found under the
given prefix::

    mrcpp/
    ├── include/
    │   └── MRCPP/
    ├── lib64/
    │   ├── libmrcpp.a
    │   ├── libmrcpp.so -> libmrcpp.so.1*
    │   └── libmrcpp.so.1*
    └── share/
        └── cmake/

Please refer to the :ref:`User's Manual` for instructions for how to run the program.

.. hint::
    We have collected scripts for configure and build of the hybrid OpenMP + MPI
    version on the different Norwegian HPC systems under ``tools/<machine>.sh``.
    These scripts will build the current version under ``build-${version}``,
    run the unit tests and install under ``install-${version}``, e.g. to build
    version v1.5.0 on Olivia::

        $ cd mrcpp
        $ git checkout v1.5.0
        $ tools/olivia.sh

    The configure step requires internet access, so the scripts must be run on
    the login nodes, and it will run on a single core, so it might take some
    minutes to complete. 


----------------
Running examples
----------------

In addition to the test suite, the code comes with a number of small code
snippets that demonstrate the features and the API of the library. These are
located in the ``examples`` directory. To compile the example codes, add the
``enable-examples`` option to setup, and the example executables can be found
under ``<build-dir>/bin/``. E.g. to compile and run the MW projection example::

    $ ./setup --enable-examples build-serial
    $ cd build-serial
    $ make
    $ bin/projection

The shared memory parallelization (OpenMP) is controlled by the environment
variable ``OMP_NUM_THREADS`` (make sure you have compiled with the ``--omp``
option to setup). E.g. to compile and run the Poisson solver example using 10
CPU cores::

    $ ./setup --enable-examples --omp build-omp
    $ cd build-omp
    $ make
    $ OMP_NUM_THREADS=10 bin/poisson

To run in MPI parallel, use the ``mpirun`` (or equivalent) command (make sure
you have compiled with the ``--mpi`` option to setup, and used MPI compatible
compilers, e.g. ``--cxx=mpicxx``). Only examples with an `mpi` prefix will be
affected by running in MPI::

    $ ./setup --cxx=mpicxx --enable-examples --mpi build-mpi
    $ cd build-mpi
    $ make
    $ mpirun -np 4 bin/mpi_send_tree

To run in hybrid OpenMP/MPI parallel, simply combine the two above::

    $ ./setup --cxx=mpicxx --enable-examples --omp --mpi build-hybrid
    $ cd build-hybrid
    $ make
    $ export OMP_NUM_THREADS=5
    $ mpirun -np 4 bin/mpi_send_tree

Note that the core of MRCPP is *only* OpenMP parallelized. All MPI data or work
distribution must be done manually in the application program, using the tools
provided by MRCPP (see the Parallel section of the API).


----------
Pilot code
----------

Finally, MRCPP comes with a personal sandbox where you can experiment and test
new ideas, without messing around in the git repository. In the ``pilot/``
directory you will find a skeleton code called ``mrcpp.cpp.sample``. To trigger
a build, re-name (copy) this file to ``mrcpp.cpp``::

    $ cd pilot
    $ cp mrcpp.cpp.sample mrcpp.cpp

Now a corresponding executable will be build in ``<builddir>/bin/mrcpp-pilot/``.
Feel free to do whatever you like in your own pilot code, but please don't add
this file to git. Also, please don't commit any changes to the existing examples
(unless you know what you're doing).

As an example, the ``pilot`` sample can be built with the following ``CMakeLists.txt``:

.. literalinclude:: snippets/CMakeLists.txt

This will set up the include paths and library paths correctly.
During configuration you will have to specify *where* the CMake configuration
file for MRCPP is located::

   $ cmake -H. -Bbuild -DMRCPP_DIR=$HOME/Software/share/cmake/MRCPP
