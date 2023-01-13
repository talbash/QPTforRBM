# QPTforRBM
Quantum-Inspired Tempering for Ground State Approximation using RBM's. Repository for code to reproduce the main results of arXiv:2210.11405 

## Compilation instructions
To compile, you will need an MPI enabled C++ compiler with c++14. On CARC, you will need to load the appropriate module before compiling:

`$ module load openmpi-4.0.3-gcc-10.1.0-qm5tmqr`

To create an executable called NNet, run from the command line:

`$ mpicxx -O3 -g -std=c++14 -fopenmp -o NNet main.cpp  -I ./ -I ./tminres-0.2 -I ./other -I ./minresqlp `

or you can use the included Makefile:

`$ make mpi`

## Usage instructions

To exectute the program:

`$ mpiexec -np $input1 NNet $input2 $input3 $input4 $input5 $input6 $input7 $input8 $input9 $input10 $input11 $input12 $input13 $input14 $input15 $input16 $input17`

where there are a total of 17 input arguments:

- $input1: (integer) Number of MPI tasks.  This controls the number of parallel tempering replicas to be used.  At the moment, the distribution of temperatures of these replicas cannot be altered.  Something for a future update.
 
- $input2: (integer) Number of OMP threads.  This controls the number of parallel threads to be used for the MCMC sampling used to estimate the covariance matrix and force vector.  When running on CARC, make sure that the number of cores equals $input1 times $input2.
 
- $input3: (integer) Number of visible nodes of the RBM.  This should equal the number of qubits.
 
- $input4: (integer) Number of hidden nodes of the RBM.
 
- $input5: (integer) Number of MCMC samples to take in order to estimate the covariance matrix and force vector.
 
- $input6: (integer) Number of MCMC sweeps to perform before a sample measurement is taken.  A single sweep corresponds performing a single Metropolis update on all spins.

- $input7: (integer) Number of RBM parameter updates to perform.

- $input8: Two possible situations here.

    if $input16 = 0, then (double) Initial value for the learning rate \gamma. RBM parameters are updated according to -\gamma S^{-1} F.  A good choice seems to be 0.01.
    
    if $input16 = 1, then (integer), the line number in the ANNparameter file to start averaging (inclusive)

- $input9: Two possible situations here.

    if $input16 = 0, then (double) Final value for the learning rate \gamma. A good choice seems to be 0.001.
    
    if $input16 = 1,  then (integer), the line number in the ANNparameter file to end averaging (inclusive)

- $input10: Two possible situations here.

    if $input16 = 0, then (integer) log_10 of the maximum number of iterations used by the MINRES or MINRES-QLP algorithm. A good value is 3.
    
    if $input16 = 1, then (integer) log_10 of the number of energy samples to take.

- $input11: (double) -log_10 of the tolerance used by the MINRES or MINRES-QLP algorithm. A good value is 6.

- $input12: (integer) A label for the training run.  Give different values when wanting to run independent training runs.

- $input13: (integer) Chooses between which method to use to calculate the inverse of the covariance matrix S.  For MINRES, choose value 0; for MINRES-QLP choose value 1.

- $input14: (integer) Chooses what kind of MCMC simulation to use.  For single spin updates, choose 0.  For updates that preserve the Hamming weight of the spin configuration, choose 1.
 
- $input15: (integer) If $input14 = 1, chooses the Hamming weigt of the spin configuration.

- $input16: (integer) If $input16 = 0, it is the training algorithm.  If $input16 = 1, it is the sampling algorithm.

- $input17: (*char) The name of the input file without its extension; the extensin is assumed to be .txt. The name of the file should not start with a number.


## Input file format
The input file must have a specific format.  Each line corresponds to a unique Pauli operator term that appears in the Hamiltonian.  The first argument in the line is the locality of the operator, which gives how many non-identity terms will appear in the operator.  After the locality, we have pairs of numbers: the first corresponding to the type of Pauli operator (1 = X, 2 = Y, 3 = Z) and the second corresponding to the qubit index that the Pauli operators acts.  The number of pairs should equal the locality. If the locality is 0, it is assumed that this is an identity term corresponding to an overall energy shift term.  The last two arguments in the line are the real part and imaginary part of the coefficient of the operator. For example, an input file of the form:

```
0   0.000
1   1   0   -1.00   0.00
2   3   0   3   1   -1.00   0.00
```

corresponds to the Hamiltonian H = - X_0 - Z_0 Z_1

Some example input files can be found in the InstanceFiles folder.

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
