using namespace std;
#include <random>

#include <minresqlp.hpp>
#include "SimpleVector.hpp"
#include "SimpleVector.cpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <mpi.h>
#include <cstring>

#include "SR_StandardRBM.cpp"

static double epsM= 2.22*pow(10.0,-16);

void ReadSpinSystem(ifstream& Instance, spin &SpinSystem, SimulationParameters &param)
{
    unsigned int OpLength, OpType, OpIndex;
    double OpStrengthReal,OpStrengthImag,Eoffset,MaxStrength, MinStrength;
    std::complex <double> OpStrength;
    
    
    
    MaxStrength = 0;
    MinStrength = 100;
    Instance >> OpLength;
    while(!Instance.eof()){
        if(OpLength>0){
            SpinSystem.OpType.resize(SpinSystem.OpType.size()+1);
            SpinSystem.OpIndices.resize(SpinSystem.OpIndices.size()+1);
            SpinSystem.OpIx.resize(SpinSystem.OpIx.size()+1);
            SpinSystem.OpIy.resize(SpinSystem.OpIy.size()+1);
            SpinSystem.OpIz.resize(SpinSystem.OpIz.size()+1);
            for (unsigned int i = 0; i< OpLength; i++){
                Instance >> OpType;
                Instance >> OpIndex;
                
                SpinSystem.OpType[SpinSystem.OpType.size()-1].push_back(OpType);
                SpinSystem.OpIndices[SpinSystem.OpIndices.size()-1].push_back(OpIndex);
                
                if(OpType==1){
                    SpinSystem.OpIx[SpinSystem.OpIx.size()-1].push_back(OpIndex);
                }
                else if(OpType==2){
                    SpinSystem.OpIy[SpinSystem.OpIy.size()-1].push_back(OpIndex);
                }
                else if(OpType==3){
                    SpinSystem.OpIz[SpinSystem.OpIz.size()-1].push_back(OpIndex);
                }
            } //for
            Instance >> OpStrengthReal;
            Instance >> OpStrengthImag;
            SpinSystem.OpStrength.push_back(OpStrengthReal+1i*OpStrengthImag);
            if(abs(OpStrengthReal+1i*OpStrengthImag) > MaxStrength){
                MaxStrength = abs(OpStrengthReal+1i*OpStrengthImag);
            }
            
            if(abs(OpStrengthReal+1i*OpStrengthImag) < MinStrength){
                MinStrength = abs(OpStrengthReal+1i*OpStrengthImag);
            }
        }
        
        
        else{
            Instance >> OpStrengthReal;
            SpinSystem.Eoffset = OpStrengthReal;
        }
        
        Instance >> OpLength;
        
    } //while
    
    SpinSystem.MaxStrength = MaxStrength;
    SpinSystem.MinStrength = MinStrength;
    
    if(param.FixedHammingWeight){
        
        SpinSystem.DriverType.resize(3*param.n*(param.n-1)/2);
        SpinSystem.DriverIndices.resize(3*param.n*(param.n-1)/2);
        SpinSystem.DriverIx.resize(3*param.n*(param.n-1)/2);
        SpinSystem.DriverIy.resize(3*param.n*(param.n-1)/2);
        SpinSystem.DriverIz.resize(3*param.n*(param.n-1)/2);
        SpinSystem.DriverStrength.resize(3*param.n*(param.n-1)/2);
        
        OpIndex = 0;
        for (unsigned int i = 0; i< param.n; i++){
            for (unsigned int k=i+1; k<param.n; k++){
                
                for(unsigned int j = 0; j< 3; j++){
                    SpinSystem.DriverType[j + 3*OpIndex].push_back(j+1);
                    SpinSystem.DriverIndices[j + 3*OpIndex].push_back(i);
                    SpinSystem.DriverType[j + 3*OpIndex].push_back(j+1);
                    SpinSystem.DriverIndices[j + 3*OpIndex].push_back(k);
                    if(j==0){
                        SpinSystem.DriverIx[j + 3*OpIndex].push_back(i);
                        SpinSystem.DriverIx[j + 3*OpIndex].push_back(k);
                    }
                    else if(j==1){
                        SpinSystem.DriverIy[j + 3*OpIndex].push_back(i);
                        SpinSystem.DriverIy[j + 3*OpIndex].push_back(k);
                    }
                    else if(j==2){
                        SpinSystem.DriverIz[j + 3*OpIndex].push_back(i);
                        SpinSystem.DriverIz[j + 3*OpIndex].push_back(k);
                    }
                    SpinSystem.DriverStrength[j + 3*OpIndex]=(-1.0);
                }
                
                OpIndex++;
            } //for k
            
            
        }//for i
        
        
        
        
    }//if fixed hamming weight
    else{
        SpinSystem.DriverType.resize(param.n);
        SpinSystem.DriverIndices.resize(param.n);
        SpinSystem.DriverIx.resize(param.n);
        SpinSystem.DriverIy.resize(param.n);
        SpinSystem.DriverIz.resize(param.n);
        SpinSystem.DriverStrength.resize(param.n);
        
        for (unsigned int i = 0; i< param.n; i++){
            
            
            SpinSystem.DriverType[i].push_back(1);
            SpinSystem.DriverIndices[i].push_back(i);
            
            SpinSystem.DriverIx[i].push_back(i);
            
            SpinSystem.DriverStrength[i]=(-1.0);
            
        }
    }
    
    
} //ReadSpinSystem



int main(int argc, char **argv)
{
    
    /* Examples to run code
     mpiexec -n 4 NNet 8 48 0 4 20 0 H4.80
     */
    
    /* initialize MPI stuff */
    int rank,nNodes;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    
    SimulationParameters param;
    double rn,Tmax,ratio,pMet,T1,T2;
    int info, iter,count;
    unsigned long long int nBits;
    unsigned int SwapFlag,tempReplica,replica1,rank1,rank2;
    char buffer[100];
    ofstream OutputFile,OutputFile4,OutputFile5;
    ifstream InstanceFile;
    string Case,FileName;
    spin SpinSystem;
    
    param.nNodes = nNodes;
    
    /*! Read simulations parameters */
    param.n = atoi(argv[1]); //number of visible nodes; should be equal to the number of qubits
    param.m = atoi(argv[2]); //number of hidden nodes;
    param.FixedHammingWeight = atoi(argv[3]); //whether to sample only bit strings of a fixed hamming weight (1) or not (0)
    param.HammingWeight = atoi(argv[4]); //what hamming weight needs to be sampled if param.FixedHammingWeight=1
    param.nUpdates = atoi(argv[5]); //number of ANN updates to perform in total
    param.RepVersion = atoi(argv[6]); //Assigns a value to the training run
    strcpy(buffer, argv[7]); //Prefix for the output files
    
    
    //perform exact sampling
    if(param.FixedHammingWeight==0){
        //not fixed Hamming weight, so need to test a total of 2^n inputs
        param.nSperNode = pow(2,param.n);
        param.nS = param.nSperNode;
    }
    else{
        //Fixed hamming weight so only need n choose k inputs
        nBits = choose(param.n,param.HammingWeight);
        param.nSperNode = nBits;
        param.nS = nBits;
    }
    
    param.gammaInit = 0.01;  //Initial value of learning rate to use
    param.gammaFinal = 0.01; //Final value of learning rate to reach in training.  We won't need this since we will be using RK method
    param.gamma = param.gammaInit;
    param.max_iter = pow(10,3);//Maximum number of iterations used by MINRES or MINRES-QLP;
    param.tol = pow(10.0,-6.0); //Tolerance used by MINRES or MINRES-QLP;
    param.Delta = pow(10,-6); //controls how small the maximum wavefunction change can be
    param.updatesPerSwap = 20;//number of ANN updates to perform before a PT swap move is attempted.
    
    
    /*!Covariance regularization parameters */
    param.lambda0 = 100;
    param.b0 = 0.9;
    param.lambdamin = pow(10.0,-4);
    
    //Assign interpolation parameter according to a cubic function
    Tmax = 10; //sets a maximum temperature to use
    if(param.nNodes>1){
        ratio = 1.0/(param.nNodes-1.0);
        param.InterpolationParameter =  Tmax*pow(rank*ratio,3);
    }
    else{
        param.InterpolationParameter = 0.0;
    }
    
    /*! Read in Hamiltonian here */
    
    FileName.append(buffer);
    FileName.append(".txt");
    InstanceFile.open(FileName.c_str());
    if(!InstanceFile){
        cout << "Failed to open " << FileName.c_str() << endl;
        MPI_Finalize();
        exit(1);
    } //if
    else{
        ReadSpinSystem(InstanceFile, SpinSystem, param);
        InstanceFile.close();
    } //else
    
    
    if(param.InterpolationParameter>0){
        param.beta = 1/(SpinSystem.MaxStrength*param.InterpolationParameter);
    }
    else{
        param.beta = -1;
    }
    

    /*! initalize random number seeds*/
    std::random_device r;
    param.seed = rank+r(); //2574165069
    param.seedThreads = rank+r(); //1909953757
    
    
    /*!Initialize random number generator */
    std::mt19937_64 rng(param.seed);
    vector <std::mt19937_64> rngThreads(1);
    for(int i = 0; i< 1; i++){
            rngThreads[i].seed(param.seedThreads+i);
    }
    std::uniform_real_distribution<> rand01(0.0,1.0);
    
    /*! Initialize network parameters */
    param.nNetworkParam = param.n + param.m + param.n*param.m;
    vector <std::complex <double>> a(param.n);
    vector <std::complex <double>> b(param.m);
    vector <std::complex <double>> W(param.n*param.m);
    vector <std::complex <double>> a0(param.n);
    vector <std::complex <double>> b0(param.m);
    vector <std::complex <double>> W0(param.n*param.m);
    vector <int> Nodev(param.n);
    
    /*! Vectors for storing matrices and vectors associated with the local energy and covariance */
    vector <std::complex <double> > NodeAverage;
    std::complex <double> NodeAverageElocTarget;
    vector <std::complex <double> > NodeEloc;
    vector <std::complex <double> > NodeElocTarget;
    vector <std::complex <double> > NodeX;
    vector <std::complex <double> > GlobalAverage;
    std::complex <double> GlobalAverageElocTarget;
    vector <std::complex <double> > GlobalEloc;
    vector <std::complex <double> > GlobalX;
    
    vector <double> GatherTargetE(param.nNodes);
    vector <double> GatherExactTargetE(param.nNodes);
    vector <double> GatherBeta(param.nNodes);

    
    vector <double> ReplicaLearningRate(param.nNodes);
    vector <double> ReplicaLearningRateInit(param.nNodes);
    //vector <double> ReplicaLearningRateDecay(param.nNodes);
    vector <double> ReplicaInterpParameter(param.nNodes);
    
    //vector <double> pMetVector(param.nNodes);
    //vector <double> pSwapVector(param.nNodes);
    vector <unsigned int> Rank2Replica(param.nNodes);
    vector <unsigned int> Replica2Rank(param.nNodes);
    vector <unsigned int> SwapFlagVector(param.nNodes);


    param.replica = rank;
    
    for (unsigned int i =0; i<Replica2Rank.size(); i++){
        Replica2Rank[i] = i;
    }
    for (unsigned int i =0; i<Rank2Replica.size(); i++){
        Rank2Replica[i] = i;
    }
    
    
    
    MPI_Gather(&param.beta,1,MPI::DOUBLE,&GatherBeta.front(),1,MPI::DOUBLE,0,MPI_COMM_WORLD);
    MPI_Allgather(&param.InterpolationParameter,1,MPI::DOUBLE,&ReplicaInterpParameter.front(),1,MPI::DOUBLE,MPI_COMM_WORLD);
    MPI_Allgather(&param.gammaInit,1,MPI::DOUBLE,&ReplicaLearningRateInit.front(),1,MPI::DOUBLE,MPI_COMM_WORLD);
    MPI_Allgather(&param.gamma,1,MPI::DOUBLE,&ReplicaLearningRate.front(),1,MPI::DOUBLE,MPI_COMM_WORLD);
    

    /*! Initialize the network weights. */
    InitializeRBMParameters2(0,0.01,a,b,W,rng);
    
    /*! Initialize the output files */
    Case.append(buffer);
    InitializeFiles4(rank, Case, param, GatherBeta, ReplicaInterpParameter, ReplicaLearningRate, OutputFile, OutputFile4, OutputFile5);
    
    /*! Initialize the states*/
    InitializeState2(param,Nodev,rngThreads);
    

    NodeAverage.resize(param.nNetworkParam + 1);
    NodeX.resize(param.nSperNode*param.nNetworkParam);
    NodeEloc.resize(param.nSperNode);
    NodeElocTarget.resize(param.nSperNode);
    GlobalAverage.resize(param.nNetworkParam + 1);
    GlobalX.resize(param.nS*param.nNetworkParam);
    GlobalEloc.resize(param.nS);
    
    /*! Create objects needed for linear system */
    Preconditioner * prec = NULL;
    Operator X(param.m,param.n,param.nS);
    SimpleVector F(2*param.nNetworkParam);
    SimpleVector x(2*param.nNetworkParam);
    SimpleVector x1(2*param.nNetworkParam);
    SimpleVector x2(2*param.nNetworkParam);
    SimpleVector x3(2*param.nNetworkParam);
    SimpleVector x4(2*param.nNetworkParam);
    SimpleVector x5(2*param.nNetworkParam);
    
    
    /*! Parameters used for PT*/
    param.RunningTargetE = 0.0;
    param.RunningCount = 0;
    param.EnergyTarget = 0;
    
    
    /*! Perform variational updates */
    for (int update = 0; update < param.nUpdates; update++){
        
        count = 0;
        a0=a;
        b0=b;
        W0=W;
        /* Runge-Kutte (RK) loop */
        while (count<5){
            /*! Perform SR sampling to construct covariance matrix */
            SRexact(param, SpinSystem, a0, b0, W0, Nodev, NodeAverage, NodeAverageElocTarget, NodeX, NodeEloc, NodeElocTarget,GlobalAverage, GlobalAverageElocTarget, GlobalX, GlobalEloc);
            
            /*! Set values of covariance matrix */
            X.SetR(max(param.lambda0*pow(param.b0,update),param.lambdamin));
            X.SetVals(GlobalX);
            X.Apply1(GlobalEloc,F); //calculate force vector
            
            
            /*! Solve the linear problem F = S x with MINRES-QLP */
            x = 0;
            info = MINRESQLP(X, x, F, param.tol,param.max_iter, iter);
            
            /*!   If linear system is successfully solved, then update ANN parameters.  */
            if(info==1){
                if(count==0){
                    x1 = x;
                    UpdateNetworkRK(update, param.gamma, param, F,  x1, X, GlobalAverage, GlobalAverageElocTarget,a0,b0,W0, rng);
                }
                else if(count==1){
                    a0=a;
                    b0=b;
                    W0=W;
                    x2 = x;
                    UpdateNetworkRK(update, param.gamma/2.0, param, F,  x2, X, GlobalAverage, GlobalAverageElocTarget,a0,b0,W0, rng);
                }
                else if(count==2){
                    a0=a;
                    b0=b;
                    W0=W;
                    x3 = x;
                    UpdateNetworkRK(update, param.gamma/4.0, param, F,  x1, X, GlobalAverage, GlobalAverageElocTarget,a0,b0,W0, rng);
                    UpdateNetworkRK(update, param.gamma/4.0, param, F,  x3, X, GlobalAverage, GlobalAverageElocTarget,a0,b0,W0, rng);
                }
                else if(count==3){
                    a0=a;
                    b0=b;
                    W0=W;
                    x4 = x;
                    UpdateNetworkRK(update, param.gamma/4.0, param, F,  x1, X, GlobalAverage, GlobalAverageElocTarget,a0,b0,W0, rng);
                    UpdateNetworkRK(update, param.gamma/4.0, param, F,  x3, X, GlobalAverage, GlobalAverageElocTarget,a0,b0,W0, rng);
                    UpdateNetworkRK(update, param.gamma/2.0, param, F,  x4, X, GlobalAverage, GlobalAverageElocTarget,a0,b0,W0, rng);
                    
                }
                else if(count==4){
                    x5 = x;
                    
                }
                count++;
                if(count==5){
                    
                    UpdateNetworkRK2exact(update, param, F,  x1, x2, x3, x4, x5,X, GlobalAverage, GlobalAverageElocTarget,a,b,W, rng);
                    
                }
                
                
            }//if successfully solved problem
            else{
                //when this happens it typically is some blow up in the calculation of the local energy
                InitializeRBMParameters2(0,0.01,a,b,W,rng);
                param.RunningCount=1;
                param.RunningTargetE = 0;
                param.gamma = param.gammaInit;
                count = 6;
                
            }

        } //RK while loop
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* PT swap loop */
        {
            if(modulo(update,param.updatesPerSwap)==param.updatesPerSwap/2-1 || (modulo(update,param.updatesPerSwap)==param.updatesPerSwap-1)){
                
                
                SRexact(param, SpinSystem, a, b, W, Nodev, NodeAverage, NodeAverageElocTarget, NodeX, NodeEloc, NodeElocTarget,GlobalAverage, GlobalAverageElocTarget, GlobalX, GlobalEloc);
                
                
                param.RunningCount++;
                param.RunningTargetE += GlobalAverageElocTarget.real();
                param.RunningTargetE = param.RunningTargetE/(param.RunningCount);
                
                MPI_Allgather(&param.RunningTargetE,1,MPI::DOUBLE,&GatherTargetE.front(),1,MPI::DOUBLE,MPI_COMM_WORLD);
                MPI_Allgather(&param.EnergyTarget,1,MPI::DOUBLE,&GatherExactTargetE.front(),1,MPI::DOUBLE,MPI_COMM_WORLD);
                
                param.RunningCount = 0;
                param.RunningTargetE = 0.0;
                
                
                if(rank==0){
                    OutputFile4 << update+1 << "\t";
                    for(unsigned int i= 0; i<GatherExactTargetE.size(); i++){
                        OutputFile4 << std::fixed << std::setprecision(8) << GatherExactTargetE[Replica2Rank[i]] << "\t";
                    }
                    OutputFile4 << "\n";
                    OutputFile4.flush();
                    
                    OutputFile5 << update << "\t";
                    for(unsigned int i= 0; i<GatherExactTargetE.size(); i++){
                        OutputFile5  << Replica2Rank[i] << "\t";
                    }
                    OutputFile5 << "\n";
                    OutputFile5.flush();
                }
                
            }
            //Swap pattern A
            if(modulo(update,param.updatesPerSwap)==param.updatesPerSwap/2-1){
                
                replica1 = param.replica;
                
                if(modulo(replica1,2)==0 && replica1<param.nNodes-1){
                    
                    rank2 = Replica2Rank[replica1+1];
                    if(ReplicaInterpParameter[replica1]==0){
                        if(GatherTargetE[rank] > GatherTargetE[rank2]){
                            pMet = 1.0;
                        }
                        else{
                            pMet = 0.0;
                        }
                    }
                    else{
                        T1 = (SpinSystem.MaxStrength)*ReplicaInterpParameter[replica1];
                        T2 = (SpinSystem.MaxStrength)*ReplicaInterpParameter[replica1+1];
                        pMet = min(1.0,exp((GatherTargetE[rank] - GatherTargetE[rank2])*(1.0/T1-1.0/T2)));
                    }
                    
                    rn = rand01(rng);
                    
                    if(rn<pMet){
                        SwapFlag=1;
                    }
                    else{
                        SwapFlag=0;
                    }
                }
                else if(modulo(replica1,2)==1){
                    rank1 = Replica2Rank[replica1-1];
                    SwapFlag=0;
                    pMet=0;
                }
                
                MPI_Allgather(&SwapFlag,1,MPI::UNSIGNED,&SwapFlagVector.front(),1,MPI::UNSIGNED,MPI_COMM_WORLD);
                
                
                if(modulo(replica1,2)==0 && replica1<param.nNodes-1){
                    if(SwapFlagVector[rank]==1){
                        
                        MPI_Send(&param.replica,1, MPI::UNSIGNED, rank2,0,MPI_COMM_WORLD);
                        MPI_Recv(&param.replica,1, MPI::UNSIGNED, rank2,0,MPI_COMM_WORLD,&status);
                        
                        param.InterpolationParameter = ReplicaInterpParameter[param.replica];
                        param.gamma = ReplicaLearningRate[param.replica];
                        param.gammaInit = ReplicaLearningRateInit[param.replica];
                        
                    }
                }
                else if(modulo(replica1,2)==1){
                    if(SwapFlagVector[rank1]==1){
                        tempReplica = param.replica;
                        
                        MPI_Recv(&param.replica,1, MPI::UNSIGNED, rank1,0,MPI_COMM_WORLD,&status);
                        MPI_Send(&tempReplica,1, MPI::UNSIGNED, rank1,0,MPI_COMM_WORLD);
                        
                        param.InterpolationParameter = ReplicaInterpParameter[param.replica];
                        param.gamma = ReplicaLearningRate[param.replica];
                        param.gammaInit = ReplicaLearningRateInit[param.replica];
                    }
                }
                
            }// if even update
            //Swap pattern B
            if(modulo(update,param.updatesPerSwap)==param.updatesPerSwap-1){
                
                replica1 = param.replica;
                
                if(modulo(replica1,2)==1 && replica1<param.nNodes-1){
                    
                    rank2 = Replica2Rank[replica1+1];
                    
                    T1 = (SpinSystem.MaxStrength)*ReplicaInterpParameter[replica1];
                    T2 = (SpinSystem.MaxStrength)*ReplicaInterpParameter[replica1+1];
                    pMet = min(1.0,exp((GatherTargetE[rank] - GatherTargetE[rank2])*(1.0/T1-1.0/T2)));
                    
                    rn = rand01(rng);
                    
                    
                    
                    if(rn<pMet){
                        SwapFlag=1;
                    }
                    else{
                        SwapFlag=0;
                    }
                }
                else if(modulo(replica1,2)==0 && replica1>0){
                    rank1 = Replica2Rank[replica1-1];
                    SwapFlag=0;
                    pMet = 0;
                }
                
                MPI_Allgather(&SwapFlag,1,MPI::UNSIGNED,&SwapFlagVector.front(),1,MPI::UNSIGNED,MPI_COMM_WORLD);
                

                if(modulo(replica1,2)==1 && replica1<param.nNodes-1){
                    if(SwapFlagVector[rank]==1){
                        
                        MPI_Send(&param.replica,1, MPI::UNSIGNED, rank2,0,MPI_COMM_WORLD);
                        MPI_Recv(&param.replica,1, MPI::UNSIGNED, rank2,0,MPI_COMM_WORLD,&status);
                        
                        param.InterpolationParameter = ReplicaInterpParameter[param.replica];
                        param.gamma = ReplicaLearningRate[param.replica];
                        param.gammaInit = ReplicaLearningRateInit[param.replica];
                    }
                }
                else if(modulo(replica1,2)==0 && replica1>0){
                    if(SwapFlagVector[rank1]==1){
                        tempReplica = param.replica;
                        
                        MPI_Recv(&param.replica,1, MPI::UNSIGNED, rank1,0,MPI_COMM_WORLD,&status);
                        MPI_Send(&tempReplica,1, MPI::UNSIGNED, rank1,0,MPI_COMM_WORLD);
                        
                        param.InterpolationParameter = ReplicaInterpParameter[param.replica];
                        param.gamma = ReplicaLearningRate[param.replica];
                        param.gammaInit = ReplicaLearningRateInit[param.replica];
                    }
                }
                
            } //odd update
            
            if(modulo(update,param.updatesPerSwap)==param.updatesPerSwap/2-1 || (modulo(update,param.updatesPerSwap)==param.updatesPerSwap-1)){
                
                MPI_Gather(&param.replica,1,MPI::UNSIGNED,&Rank2Replica.front(),1,MPI::UNSIGNED,0,MPI_COMM_WORLD);
                
                if(rank==0){
                    for (unsigned int i = 0; i< Replica2Rank.size(); i++){
                        Replica2Rank[Rank2Replica[i]] = i;
                    }
                }
                
                MPI_Bcast(&Replica2Rank.front(),Replica2Rank.size(),MPI::UNSIGNED,0,MPI_COMM_WORLD);
            } //output data to file
            
            
        }
    }//for update
    
    
    if(rank==0){
        /*! Close files */
        OutputFile4.close();
        OutputFile5.close();
    }
    
    MPI_Finalize();
    
    return 0;
}
