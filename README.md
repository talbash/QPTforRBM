# QPTforRBM
Quantum-Inspired Tempering for Ground State Approximation using RBM's. Repository for code to reproduce the main results of arXiv:2210.11405 

## Compilation instructions
To compile, you will need an MPI enabled C++ compiler with c++14.  For example, gcc 11.3 and openmpi 4.1.4 should allow compilation.

There are two source codes available.  The first, XXXX, uses a wavefunction ansatze based on a standard RBM, as in Eqt.(1) of arXiv:2210.11405. The second, XXXX, uses a wavefunction ansatze based on a symmetric RBM, as in Eqt.(17) of arXiv:2210.11405.

To create an executable called QPTRBM, run from the command line:

`$ mpicxx -O3 -g -std=c++14 -o QPTRBM main.cpp  -I ./ -I ./tminres-0.2 -I ./other -I ./minresqlp `

where main.cpp should be replaced with the appropriate source file from above.

## Usage instructions for symmetric RBM

To execute the program:

`$ mpiexec -np $input1 QPTRBM $input2 $input3 $input4 $input5`

where there are a total of 5 input arguments:

- $input1: (integer) Number of MPI tasks.  This controls the number of parallel tempering replicas to be used. The distribution of temperatures of these replicas is a cubic function.
 
- $input2: (integer) Number of visible nodes of the RBM.  This should equal the number of qubits in the problem Hamiltonian.
 
- $input3: (integer) Number of hidden nodes of the RBM.

- $input4: (integer) Number of RBM parameter updates using Stochastic Reconfiguration to perform.

- $input5: (integer) Assign a value to the simulation.  The output files will be labeled with this value.


## Usage instructions for standard RBM

To execute the program:

`$ mpiexec -np $input1 QPTRBM $input2 $input3 $input4 $input5 $input6 $input7 $input8`

where there are a total of 8 input arguments:

- $input1: (integer) Number of MPI tasks.  This controls the number of parallel tempering replicas to be used. The distribution of temperatures of these replicas is a cubic function.
 
- $input2: (integer) Number of visible nodes of the RBM.  This should equal the number of qubits in the problem Hamiltonian.
 
- $input3: (integer) Number of hidden nodes of the RBM.
 
- $input4: (integer) Chooses whether to perform full sampling (0) or restrict the sampling to a fixed Hamming weight sector (1).
 
- $input5: (integer) If $input4 = 1, chooses the Hamming weight of the spin configurations to sample.

- $input6: (integer) Number of RBM parameter updates using Stochastic Reconfiguration to perform.

- $input7: (integer) Assign a value to the simulation.  The output files will be labeled with this value.

- $input8: (*char) The name of the input file without its extension; the extension is assumed to be .txt. The name of the file should not start with a number.


## Input file format for standard RBM
The input file must have a specific format.  Each line corresponds to a unique Pauli operator term that appears in the Hamiltonian.  The first argument in the line is the locality of the operator, which gives how many non-identity terms will appear in the operator.  After the locality, we have pairs of numbers: the first corresponding to the type of Pauli operator (1 = X, 2 = Y, 3 = Z) and the second corresponding to the qubit index that the Pauli operators acts.  The number of pairs should equal the locality. If the locality is 0, it is assumed that this is an identity term corresponding to an overall energy shift term.  The last two arguments in the line are the real part and imaginary part of the coefficient of the operator. For example, an input file of the form:

```
0   0.000
1   1   0   -1.00   0.00
2   3   0   3   1   -1.00   0.00
```

corresponds to the Hamiltonian H = - X_0 - Z_0 Z_1

The input files for the H4 rectangle used in arXiv:2210.11405 can be found in the InstanceFiles folder.  They are H4.80.txt, H4.85.txt, and H4.90.txt, corresponding to the three angles of 80, 85, and 90 degrees studied. 

## Output files
In an effort to be as complete as possible, the code outputs a lot of information that is organized in different files.  The filenames will have input16 as their prefix, in addition to some of the input parameters, such as the number of visible and hidden nodes, number of samples, number of sweeps between samples, and the number of ANN updates performed.  The following files will be output when the training algorithm is used ($input16 = 0):

- ($input17)_n=($input3)_m=($input4)_nS=($input5)_nSweepsPerSwap=($input6)_ANNParameterFile_r($input12).dat

    Each line gives the current state of the ANN for the zero temperature replica.  The first n data points is the spin configuration of the visible layer, followed by the 2 \times n elements of the ANN parameter a (remember that a is complex so a_1 has two parts), follewed by the 2 \times m elements of the ANN parameter b, followed by the 2 \times n \times m elements of W (indexed as W[i + j*m] with i = 0,...,m-1, j=0,...,n-1).

- ($input17)_n=($input3)_m=($input4)_nS=($input5)_nSweepsPerSwap=($input6)_DataFile_Replica($k)_r($input12).dat

    Identical to the previous, except it gives the ANN parametesr for the k-th temperature replica.  These files will only appear if $input1 is greater than 1.

- ($input17)_n=($input3)_m=($input4)_nS=($input5)_nSweepsPerSwap=($input6)_EnergyForceFile_r($input12).dat

    Five column file listing: the local energy, the wavefunction change, the expected energy change, the norm of the force vector, and the implemented learning rate.

- ($input17)_n=($input3)_m=($input4)_nS=($input5)_nSweepsPerSwap=($input6)_SimulationFile_r($input12).dat

    File that lists all the simulation parameters.  It also outputs the replica swap probabilities.

- ($input17)_n=($input3)_m=($input4)_nS=($input5)_nSweepsPerSwap=($input6)_TemperatureFile_r($input12).dat

    File that lists that trajectory of the temperature replicas.

If the sampling algorithm is used ($input16 = 1), then two files are output:

-  ($input17)_n=($input3)_m=($input4)_nS=($input5)_nSweepsPerSwap=($input6)_EnergySamples_r($input12).dat

    File that lists the energy samples for the averaged ANN parameters in the file ($input17)_n=($input3)_m=($input4)_nS=($input5)_nSweepsPerSwap=($input6)_ANNParameterFile_r($input12).dat.

-  ($input17)_n=($input3)_m=($input4)_nS=($input5)_nSweepsPerSwap=($input6)_SamplingANNParameters_r($input12).dat

    File that lists the averaged ANN parameters used for sampling
    
## MINRES
Implementation of a conjugate-gradient type method for solving sparse symmetric linear equations.  Returns a least-squares solution.

Cite:
-  C. C. Paige and M. A. Saunders (1975). Solution of sparse indefinite systems of linear equations, SIAM J. Numerical Analysis 12, 617-629.
- S.-C. Choi (2006). Iterative Methods for Singular Linear Equations and Least-Squares Problems, PhD thesis, Stanford University.

## MINRES-QLP
Implementation of a conjugate-gradient type method for solving sparse symmetric/Hermitian linear equations.  The method is based on Lanczos tridiagonalization. Returns minimum-length solution.

Cite:
- MINRES-QLP: A Krylov subspace method for indefinite or singular symmetric systems, SIAM J. Sci. Comput. 33:4, 1810-1836, published electronically Aug 4, 2011.  
- S.-C. T. Choi and M. A. Saunders. Algorithm 937: MINRES-QLP for symmetric and Hermitian linear equations and least-squares problems, ACM Trans. Math. Softw. 40:2, Article 16 (Feb 2014), 12 pp.
- S.-C. Choi (2006). Iterative Methods for Singular Linear Equations and Least-Squares Problems, PhD thesis, Stanford University.

MINRES-QLP algorithm C++ implementation is based on the MATLAB code available [here](https://web.stanford.edu/group/SOL/software/minresqlp/).
