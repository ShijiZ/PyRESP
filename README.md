# PyRESP
This is the repository for the PyRESP (Python Restrained Electrostatic Potential) program. This program fits the quantum mechanically calculated electrostatic potential at molecular surfaces using electrostatic models with atom-centered (1) permanent charges and (2) induced dipoles and (3) permanent dipoles. The molecular surfaces should be generated beyond Van der Waal surface in order to minimize other contributions such as exchange repulsion and charge transfer. 

## History Versions and Authors
    PyRESP version 1.0     Feb 2022 - Shiji Zhao

    RESP   version 2.4     Nov 2013 - q4md-forcefieldtools.org
    RESP   version 2.2     Jan 2011 - q4md-forcefieldtools.org
    RESP   version 2.1     Oct 1994 - Jim Caldwell
    RESP   version 2.0     Sep 1992 - Christopher Bayly

    ESPFIT version 1.0 (modified)   - Ian Gould
    ESPFIT version 1.0              - U.Chandra Singh and P.A.Kollman

## Library Dependencies
This program runs with Python 3. To run this program, make sure the following dependencies are installed:
- [numpy](https://numpy.org/): A Python library supporting for large, multi-dimensional arrays and matrices.
- [scipy](https://scipy.org/): A Python library used for scientific computing and technical computing.
- [f90nml](https://github.com/marshallward/f90nml): A Python module providing interface for reading, writing, and modifying Fortran namelists.

The following command should install all required libraries:

`$ pip install numpy scipy f90nml`

## Test Cases
Four test cases are provided which covers a range of posibilities. Run the following script in each subfolder under [test](https://github.com/ShijiZ/PyRESP/tree/master/test) for testing.

`$ ./py_resp.run` 

1. [water](https://github.com/ShijiZ/PyRESP/tree/master/test/water)
- **Models:** resp, resp-ind, resp-perm and resp-perm-v.
- **Test:** One-stage fitting on a single conformation.

2. [ethylene](https://github.com/ShijiZ/PyRESP/tree/master/test/ethylene)
- **Models:** resp, resp-ind, resp-perm and resp-perm-v.
- **Test:** One stage fitting on a single conformation. No total charge constraint applied.

3. [peptoid](https://github.com/ShijiZ/PyRESP/tree/master/test/peptoid)
- **Models:** resp, resp-ind and resp-perm.
- **Test:** Two stage fitting on a single conformation. Intra-molecular fractional charge constraint applied.
- **Note:** A peptide-type residue with consistent constraints with Cornell et al. '95 force fields.

4. [bis-naphthyl](https://github.com/ShijiZ/PyRESP/tree/master/test/bis-naphthyl)
- **Models:** resp, resp-ind and resp-perm.
- **Test:** Two stage fitting on a two conformations. Inter-molecular fractional charge constraint applied.
- **Note:** Two *2-Methyl-3-naphthylpropionic acid* molecules fitted together to obtain charges for the "super molecule" *bis-(naphthyl-1-methyl) acetic acid*.

## Citation
To cite PyRESP, see the following publication:

[Shiji Zhao, Haixin Wei, Piotr Cieplak, Yong Duan, and Ray Luo, "PyRESP: A Program for Electrostatic Parameterizations of Additive and Induced Dipole Polarizable Force Fields". J. Chem. Theory Comput. 2022, 18, 6, 3654-3670](https://pubs.acs.org/doi/full/10.1021/acs.jctc.2c00230).