# QPTforRBM
Quantum-Inspired Tempering for Ground State Approximation using RBM's. Repository for code to reproduce the main results of arXiv:2210.11405.  It implements quantum parallel tempering according to the training Hamiltonian in Eqts.(14-19) using the swap probability defined in Eqt.(11).

## Compilation instructions
To compile, you will need an MPI enabled C++ compiler with c++14.  For example, gcc 11.3 and openmpi 4.1.4 should allow compilation.

There are two source codes available.  The first, main_StandardRBM.cpp, uses a wavefunction ansatze based on a standard RBM, as in Eqt.(1) of arXiv:2210.11405. The second, main_SymmetricRBM.cpp, uses a wavefunction ansatze based on a symmetric RBM, as in Eqt.(17) of arXiv:2210.11405.

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
The filenames will have the problem name as their prefix.  For example for the symmetric RBM, this will be restricted to the Precipice problem.  For the standard RBM, it will inherit the name of the input file. The following files will be output when the training algorithm is used:

- ($ProblemName)_StandardRBM_n=($input2)_m=($input3)_nUpdates=($input6)_PTEnergyFile_r($input7).dat

This file will have (1 + $input1) columns.  The first column carries the update step, while the remaining columns give the expectation value of the target Hamiltonian for each replica.

- ($ProblemName)_StandardRBM_n=($input2)_m=($input3)_nUpdates=($input6)_PTReplicaFile_r($input7).dat

This file will have (1 + $input1) columns.  The first column carries the update step, while the remaining columns give the original replica index of the configuration currently at that replica location.

- ($ProblemName)_StandardRBM_n=($input2)_m=($input3)_nUpdates=($input6)_SimulationFile_r($input7).dat

This file includes a summary of the parameters used in the simulation.



## MINRES-QLP
The code uses MINRES-QLP, an implementation of a conjugate-gradient type method for solving sparse symmetric/Hermitian linear equations.  MINRES-QLP algorithm C++ implementation is based on the MATLAB code available [here](https://web.stanford.edu/group/SOL/software/minresqlp/).  To help our implementaiton, we make use of the SimpleVector implementation by Umberto Villa, Michael Saunders, Santiago Akle.
