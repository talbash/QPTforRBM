struct SimulationParameters {
    int updatesPerSwap,n,m,nNetworkParam, alpha,nS,nSperNode,nSweepsPerSwap,nUpdates,RepVersion,max_iter,nThreads,nNodes,RunningCount,ExactSampling,PrecipicePosition;
    double beta,gammaInit,gammaFinal,tol,Delta,lambda0,b0,lambdamin,Energy,EnergyTarget,Gamma,DeltaE,fnorm,gamma,gammaDecay,RunningE,RunningTargetE,InterpolationParameter, AnnealingParameter;
    unsigned long long int seed,seedThreads;
    unsigned int FixedHammingWeight, HammingWeight, replica;
    vector <double> BinomialCoefficient;
};


struct spin {
    unsigned int Labelphy;
    unsigned int Labellog;
    double Eoffset,MaxStrength, MinStrength;
    std::complex <double> E2offset;
    vector <vector <unsigned int>> OpType,Op2Type;
    vector <vector <unsigned int>> OpIndices,Op2Indices;
    vector <vector <unsigned int>> OpIx,OpIy,OpIz,Op2Ix,Op2Iy,Op2Iz,ObIx,ObIy,ObIz;
    vector <std::complex <double>> OpStrength,Op2Strength;
    
    vector <vector <unsigned int>> DriverType;
    vector <vector <unsigned int>> DriverIndices;
    vector <vector <unsigned int>> DriverIx,DriverIy,DriverIz;
    vector <std::complex <double>> DriverStrength;
    
};

struct qubit {
    int label;
    double hx, hy, hz;
    
    vector <std::tuple<int, double>> Jxx;
    vector <std::tuple<int, double>> Jyy;
    vector <std::tuple<int, double>> Jzz;
    vector <std::tuple<int, int, double>> Jzzz;
    
};





//!@class Preconditioner
/*
 * @brief An abstract preconditioner class.
 * It should be used to call minres without preconditioner
 */
class Preconditioner
{
public:
    //! Y = M\X
    virtual void Apply(const SimpleVector & X, SimpleVector & Y) const = 0;
};

//!@class Operator
/*
 * @brief The operator X that will be passwed to minres-QLP.
 */
class Operator
{
public:
    
    //!Constructor
    /*
     * @param m_ number of hidden nodes
     * @param n_ number of visible nodes
     * @param nS_ number of samples
     */
    Operator(int m_, int n_, int nS_):
    m(m_), n(n_), nS(nS_)
    
    {
        nRows = m + n + m*n;
        nColumns = nS;
        size = nRows*nS;
        //vals = new double[size];
        /*
         int i(0);
         for( ; i < size/2; ++i)
         vals[i] = i+1;
         
         for( ; i < size; ++i)
         vals[i] = -i;
         */
    }
    
    ~Operator()
    {
        //delete[] vals;
    }
    void SetVals(vector <std::complex <double> > &x)
    {
        //cout << "hello\n";
        vals = &x.front();
    }//set
    
    void Set(int i, std::complex <double> x)
    {
        vals[i] = x;
    }//set
    
    void SetR(double x)
    {
        lambda = x;
    }//set
    
    std::complex <double> Return(int i)
    {
        return vals[i];
    }//set
    
    //! a = X*b; used to calculate the force vector
    void Apply1(vector <std::complex <double> > &X, SimpleVector &Y) const
    {
        //conjugate(X)*vec(e)
        int index;
        
        for(int i = 0; i< nRows; i++){
            Y[i] = 0;
            Y[i+nRows] = 0;
            
            for(int s = 0; s < nColumns; s++){
                index = i + s*nRows;
                Y[i] += (vals[index].real() * X[s].real() + vals[index].imag() * X[s].imag());
                Y[i+nRows] += (vals[index].real() * X[s].imag() -vals[index].imag() * X[s].real());
            }
            
        }
    }
    
    void Apply2(vector <std::complex <double> > &V, vector <std::complex <double> > &X, SimpleVector &Y) const
    {
        std::complex <double>  Z;
        
        for(int i = 0; i< nRows; i++){
            Y[i] = 0;
            Y[i+nRows] = 0;
            
            for(int s = 0; s < nColumns; s++){
                Z = 0;
                for(int j = 0; j< nRows; j++){
                    Z += std::conj(V[j + i*nRows])*std::conj(vals[j + s*nRows]);
                }// for j
                Z = Z * X[s];
                Y[i] += Z.real()*Z.real() + Z.imag()*Z.imag(); //taking into account that X and val both have been normalized by 1/(nColumns-1)^0.5
            }
            Y[i] = Y[i]*(nColumns-1);
        }
    }
    
    void Apply(const SimpleVector & X, SimpleVector & Y) const
    {
        SimpleVector Zrr(nColumns);
        SimpleVector Zri(nColumns);
        SimpleVector Zir(nColumns);
        SimpleVector Zii(nColumns);
        double temp;
        
        {
            for(int s = 0; s < nColumns; ++s){
                Zrr[s] = 0;
                Zri[s] = 0;
                Zir[s] = 0;
                Zii[s] = 0;
                for(int i = 0; i< nRows; i++){
                    Zrr[s] += vals[i + s*nRows].real() * X[i];
                    Zri[s] += vals[i + s*nRows].imag() * X[i];
                    Zir[s] += vals[i + s*nRows].real() * X[i + nRows];
                    Zii[s] += vals[i + s*nRows].imag() * X[i + nRows];
                }// for i
            }// for s
            for(int i = 0; i< nRows; i++){
                Y[i] = 0;
                Y[i + nRows] = 0;
                temp = 0;
                for(int s = 0; s < nColumns; ++s){
                    temp += (vals[i + s*nRows].real()*vals[i + s*nRows].real() + vals[i + s*nRows].imag()*vals[i + s*nRows].imag());
                    Y[i] += vals[i + s*nRows].real() * Zrr[s] + vals[i + s*nRows].imag() * Zri[s] ;
                    Y[i] += -vals[i + s*nRows].real() * Zii[s] + vals[i + s*nRows].imag() * Zir[s] ;
                    Y[i + nRows] += vals[i + s*nRows].real() * Zri[s] - vals[i + s*nRows].imag() * Zrr[s] ;
                    Y[i + nRows] += vals[i + s*nRows].real() * Zir[s] + vals[i + s*nRows].imag() * Zii[s] ;
                } //for s
                Y[i] += lambda*temp*X[i];
                Y[i+nRows] += lambda*temp*X[i + nRows];
            } //for i
        }
        
        
    }//Apply
    
    
    
private:
    int size;
    int m, n, nS, nRows, nColumns;
    std::complex <double> *vals;
    //vector <double> Smatrix;
    double lambda;
};


const unsigned long long int MAX_VALUE = std::numeric_limits< unsigned long long >::max();
const int MAX_NUMBER_OF_BITS = sizeof(unsigned long long int)*8;


int modulo(int a, int q)
{
    if(q == 0)
        return q;
    
    int result = fmod( a,  q);
    return ((result >= 0 && q > 0) || (a <= 0 && q  < 0)) ? result : (result + q);
    
}

void InitializeRBMParameters2(double mean, double stdev, vector <std::complex <double>> &a,vector <std::complex <double>> &b,vector <std::complex <double>> &W,std::mt19937_64 &rng)
{
    //! Initalize RBM network weights with gaussian
    
    std::normal_distribution<> normal_dist(mean, stdev);
    
    
    for (int i = 0; i<a.size(); i++){
        a[i] = normal_dist(rng) + normal_dist(rng)*1i;
    }
    
    for (int i = 0; i<b.size(); i++){
        b[i] = normal_dist(rng) + normal_dist(rng)*1i;
    }
    
    for (int i = 0; i<W.size(); i++){
        W[i] = normal_dist(rng) + normal_dist(rng)*1i;
    }
    
    
    
} //InitializeRBMParameters2



void InitializeState2(SimulationParameters &param, vector <int> &v, vector <std::mt19937_64> &rng)
{
    int n = param.n;
    int nThreads = 1;
    int offsetv;
    std::vector<int>::iterator it1,it2;
    unsigned long long int mask,rn;
    std::uniform_int_distribution<uint64_t> randInt64(0, MAX_VALUE);
    
    
    for (int i = 0; i<nThreads; i++){
        offsetv = i*n;
        it1 = v.begin() + offsetv;
        it2 = v.begin() + offsetv + n;
        for (int k =0; k<n; k++){
            if(modulo(k,64)==0){
                rn = randInt64(rng[i]);
                mask = 1;
            }
            else{
                mask = mask << 1;
            }
            if( (mask & rn) !=0){
                v[k + offsetv] = -1;
            }
            else{
                v[k + offsetv] = 1;
            }
        }//for k
        
    }//for i
    
}//InitializeState2


void UpdateTheta2(SimulationParameters &param, vector <int> &v, vector <std::complex <double>> &b, vector <std::complex <double>> &W, vector <std::complex <double>> &Theta)
{
    //! Recalculate the value of \theta with the current network weights \theta_i = b_i + \sum_j W_{ij} v_j
    int m = param.m;
    int n = param.n;
    
    for (int i = 0; i < m; i++){
        Theta[i] = b[i];
        for (int j = 0; j< n; j++){
            Theta[i]  += ((double) v[j])*W[i + m*j];
        }//for j
    }//for i
} //UpdateTheta2

void UpdateTheta2HW(SimulationParameters &param,int &Magnetization, vector <std::complex <double>> &b, vector <std::complex <double>> &W, vector <std::complex <double>> &Theta)
{
    //! Recalculate the value of \theta with the current network weights \theta_i = b_i + \sum_j W_{ij} v_j
    int m = param.m;
    int n = param.n;
    for (int i = 0; i < m; i++){
        Theta[i] = b[i] + W[i]*((double)Magnetization);
        /*for (int j = 0; j< n; j++){
         Theta[i]  += ((double) v[j])*W[i + m*j];
         }//for j
         */
    }//for i
} //UpdateTheta2HW





void MeasureExact4HW(int count, int update, double sqrtProb, unsigned int s, SimulationParameters &param, vector <std::complex <double>> &Theta, vector <int> &v,   vector <std::complex <double>>  &a, vector <std::complex <double>>  &W, vector <std::complex <double>> &X, vector <std::complex <double>> &Eloc, vector <std::complex <double>> &ElocTarget, spin &SpinSystem)
{
    
    //identical to MeasureParallel2 except multiplies by the unnormalized probability
    
    int m = param.m;
    int n = param.n;
    std::complex <double>  piCoshTheta,temp,temp2,dtemp, atemp, Wtemp;
    vector <int> v2;
    int j1,j2,shift1,shift2,index, index3, index2,index1,offsetW,offsetW1,offsetW2,offseta,offsetb;
    double Ed,Edp, Edm, A,B,AnnealParameter,tempI, tempR,temp2I, temp2R, Jzzz,Jzz, Jxx, Jyy, hx, hy, hz,spin,spin1,spin2,spin3,phase1,phase2,eigT;
    
    //cout << "beta = " << Beta[BetaPosition[0]] << "\n";
    //offsetv = ThreadIndex*n;
    //offsetTheta = 0*ThreadIndex*m;
    
    int offsetParameters = (1+m + m);
    
    tempR = 0;
    tempI = 0;
    temp2R = 0;
    temp2I = 0;
    
    
    AnnealParameter = param.InterpolationParameter;
    A = AnnealParameter;
    
    B=1;
    
    
    //diagonal piece first; this describes the spike problem
    temp=0;
    
    //Precipice use an annealing parameter of 0.8
    if(s==param.PrecipicePosition){
        temp += -1.0;
        temp2 += 1.0;
        Ed = -1.0;
        Edp = (double) s+1.0;
        Edm = (double) s-1.0;
    }
    else{
        Ed =(double) s;
        temp += (double) s;
        temp2 += (double) s*s;
        if(s==param.PrecipicePosition-1){
            Edm=(double) s-1.0;
            Edp= -1.0;
        }
        else if (s==param.PrecipicePosition+1){
            Edm= -1.0;
            Edp=(double) s+1.0;
        }
        else{
            Edp = (double) s+1.0;
            Edm = (double) s-1.0;
            
        }
    }
    
    tempR += param.AnnealingParameter*temp.real() ;
    tempI += param.AnnealingParameter*temp.imag();
    
    
    //transverse field contribution
    temp=0;
    if(s==0){
        //only connects to HW+1 states
        atemp = exp(-2.0*a[0]);
        piCoshTheta=1;
        for(unsigned int j = 0; j < param.m; j++){
            piCoshTheta = piCoshTheta*(cosh(Theta[j] - 2.0*W[j]))/(cosh(Theta[j]));
        }
        temp += ((double)param.n)*atemp*piCoshTheta;
    }
    else if(s==param.n){
        //only connects to HW-1 states
        atemp = exp(2.0*a[0]);
        piCoshTheta=1;
        for(unsigned int j = 0; j < param.m; j++){
            piCoshTheta = piCoshTheta*(cosh(Theta[j] + 2.0*W[j]))/(cosh(Theta[j]));
        }
        temp += ((double)param.n)*atemp*piCoshTheta;
    }
    else{
        //connects to HW+1 states
        atemp = exp(-2.0*a[0]);
        piCoshTheta=1;
        for(unsigned int j = 0; j < param.m; j++){
            piCoshTheta = piCoshTheta*(cosh(Theta[j] - 2.0*W[j]))/(cosh(Theta[j]));
        }
        temp += ((double)(param.n-s))*atemp*piCoshTheta;
        
        //connects to HW-1 states
        atemp = exp(2.0*a[0]);
        piCoshTheta=1;
        for(unsigned int j = 0; j < param.m; j++){
            piCoshTheta = piCoshTheta*(cosh(Theta[j] + 2.0*W[j]))/(cosh(Theta[j]));
        }
        temp += ((double)s)*atemp*piCoshTheta;
    }//else
    
    tempR += (1.0-param.AnnealingParameter)*(param.n/2.0 + 0.5*(-temp.real())) ;
    tempI += (1.0-param.AnnealingParameter)*( 0.5*(-temp.imag()));
    
    ElocTarget[s] = SpinSystem.Eoffset + tempR + 1i*tempI;
    
    
    tempR = 0;
    tempI = 0;
    //Training Hamiltonian contribution; same as transverse field
    /*
     //transverse field contribution
     temp=0;
     if(s==0){
     //only connects to HW+1 states
     atemp = exp(-2.0*a[0]);
     piCoshTheta=1;
     for(unsigned int j = 0; j < param.m; j++){
     piCoshTheta = piCoshTheta*(cosh(Theta[j] - 2.0*W[j]))/(cosh(Theta[j]));
     }
     temp += ((double)param.n)*atemp*piCoshTheta;
     }
     else if(s==param.n){
     //only connects to HW-1 states
     atemp = exp(2.0*a[0]);
     piCoshTheta=1;
     for(unsigned int j = 0; j < param.m; j++){
     piCoshTheta = piCoshTheta*(cosh(Theta[j] + 2.0*W[j]))/(cosh(Theta[j]));
     }
     temp += ((double)param.n)*atemp*piCoshTheta;
     }
     else{
     //connects to HW+1 states
     atemp = exp(-2.0*a[0]);
     piCoshTheta=1;
     for(unsigned int j = 0; j < param.m; j++){
     piCoshTheta = piCoshTheta*(cosh(Theta[j] - 2.0*W[j]))/(cosh(Theta[j]));
     }
     temp += ((double)(param.n-s))*atemp*piCoshTheta;
     
     //connects to HW-1 states
     atemp = exp(2.0*a[0]);
     piCoshTheta=1;
     for(unsigned int j = 0; j < param.m; j++){
     piCoshTheta = piCoshTheta*(cosh(Theta[j] + 2.0*W[j]))/(cosh(Theta[j]));
     }
     temp += ((double)s)*atemp*piCoshTheta;
     }//else
     */
    tempR += -temp.real() ;
    tempI += -temp.imag();
    
    Eloc[s] = B*ElocTarget[s] + (SpinSystem.MaxStrength)*A*(0.5*param.n + 0.5*(tempR + 1i*tempI));
    
    
    double Magnetization = ((double) (param.n-2.0*s));
    offseta = 0 + s*offsetParameters;
    offsetb = 1 + s*offsetParameters;
    offsetW = 1+m + s*offsetParameters;
    
    for (int k=0; k<offsetParameters; k++){
        
        if(k<1){
            index = k;
            index1 = index;
            index += offseta;
            X[index] = Magnetization;
        }
        else if(k<1+m){
            
            index = k-1;
            index1 = index;
            index += offsetb;
            X[index] = tanh(Theta[index1]);
            
        }
        else{
            index = k-m-1;
            index2 = index/m; //index j
            index1 = index - index2*m; //index i
            index += offsetW;
            
            X[index] = Magnetization*tanh(Theta[index1]);
        }
        
        
    } //for k
    
} //MeasureExact4HW



void NormalizeMeasureExact(vector <double> UnnormalizedP, double Normalization, int s, SimulationParameters &param, vector <std::complex <double> > &Average, vector <std::complex <double> > &X,  vector <std::complex <double> > &Eloc)
{
    int nS = param.nS;
    int nNetworkParam = Average.size()-1;
    int offset = s*(nNetworkParam);
    for (int i=0; i<nNetworkParam; i++){
        X[i + offset] = ((double)pow(UnnormalizedP[s]/Normalization,0.5))*(X[i + offset] - Average[i]);
    }
    Eloc[s] = ((double)pow(UnnormalizedP[s]/Normalization,0.5))*(Eloc[s] - Average[nNetworkParam]);
    
}//NormalizeMeasureExact


void InitializeFiles4(int rank, string &Case2, SimulationParameters &param, vector <double> &Beta, vector <double> &IP, vector <double> &LearningRate, ofstream &OutputFile, ofstream &OutputFile4, ofstream &OutputFile5)
{
    string Case,FileName;
    char buffer[100];
    
    Case = Case2;
    
    Case.append("_SymmetricRBM_n=");
    sprintf(buffer,"%d",param.n);
    Case.append(buffer);
    Case.append("_m=");
    sprintf(buffer,"%d",param.m);
    Case.append(buffer);
    Case.append("_nUpdates=");
    sprintf(buffer,"%d",param.nUpdates);
    Case.append(buffer);
    
    if(rank==0){
        FileName=Case;
        FileName.append("_SimulationFile_r");
        sprintf(buffer,"%d",param.RepVersion);
        FileName.append(buffer);
        FileName.append(".dat");
        OutputFile.open(FileName.c_str());
        
        OutputFile << "Seed for main RNG:" << param.seed << "\n";
        OutputFile << "Seed for threads RNG:" << param.seedThreads << "\n";
        OutputFile << "Number of hidden nodes in network:\t" << param.m << "\n";
        OutputFile << "Total number of network updates performed:\t" << param.nUpdates << "\n";
        
        OutputFile << "Number of measurements per network update:\t" << param.nS << "\n";
        
        OutputFile << "(Note that the total number of single site MC updates per measurement is " << param.n*param.nSweepsPerSwap << ")\n";
        OutputFile << "Perform ANN Replica swap attempts at every " << param.updatesPerSwap/2 << "\n";
        OutputFile << "Covariance matrix regularization max(lambda_0 b_0^p,lambda_min)\n";
        OutputFile << "lambda_0:\t" << param.lambda0 << "\n";
        OutputFile << "b_0:\t" << param.b0 << "\n";
        OutputFile << "lambda_min:\t" << param.lambdamin << "\n";
        
        OutputFile << "Wavefunction change:\t" << param.Delta << "\n";
        OutputFile << "Parallel tempering inverse temperature:\t";
        for (int i = 0; i< Beta.size(); i++){
            OutputFile << Beta[i] << "\t";
        }
        OutputFile << "\n";
        OutputFile << "Interpolating Hamiltonian (1-s) H_TF + (s) H_Target: \t";
        for (int i = 0; i< IP.size(); i++){
            OutputFile << IP[i] << "\t";
        }
        OutputFile << "\n";
        OutputFile << "Learning rate at each replica: \t";
        for (int i = 0; i< IP.size(); i++){
            OutputFile << LearningRate[i] << "\t";
        }
        OutputFile << "\n";
        
        
        
        OutputFile << "Using MINRES-QLP.\n";
        
        OutputFile << "MinRes tolerance parameter:\t" << param.tol << "\n";
        OutputFile << "MinRes maximum number of iterations:\t" << param.max_iter << "\n";
        
        OutputFile << "Number of nodes:\t" << param.nNodes << "\n";
        
        OutputFile.flush();
        OutputFile.close();
    }
    
    
    if(rank==0){
        FileName=Case;
        FileName.append("_PTEnergyFile_r");
        sprintf(buffer,"%d",param.RepVersion);
        FileName.append(buffer);
        FileName.append(".dat");
        OutputFile4.open(FileName.c_str());
        
        FileName=Case;
        FileName.append("_PTReplicaFile_r");
        sprintf(buffer,"%d",param.RepVersion);
        FileName.append(buffer);
        FileName.append(".dat");
        OutputFile5.open(FileName.c_str());
        
        
    }
    
} //InitializeFiles4




double CalculateUnnormalizedProbabilityHW(SimulationParameters &param, unsigned int &HW, int &Magnetization, vector <std::complex <double>> a, vector <std::complex <double>> Theta)
{
    std::complex <double> temp,p;
    
    temp = ((double)Magnetization)*a[0];
    
    p = exp(temp);
    //cout << "p = " << p << "\n";
    for (unsigned int i = 0; i< Theta.size(); i++){
        p = p*cosh(Theta[i]);
    }
    //should include combintorial factor here
    return param.BinomialCoefficient[HW]*pow(abs(p),2.0);
    
} //CalculateUnnormalizedProbabilityHW





double SRexact5HW(int count, int update, unsigned int rank, SimulationParameters &param, spin &SpinSystem, vector <std::complex <double>> &a0, vector <std::complex <double>> &b0,vector <std::complex <double>> &W0, vector <int> &Nodev, vector <std::complex <double> > &NodeAverage,  std::complex <double> &NodeAverageElocTarget, vector <std::complex <double> > &NodeX, vector <std::complex <double> > &NodeEloc, vector <std::complex <double> > &NodeElocTarget, vector <std::complex <double> > &GlobalAverage, std::complex <double> &GlobalAverageElocTarget, vector <std::complex <double> > &GlobalX, vector <std::complex <double> > &GlobalEloc, vector <std::mt19937_64> &rngThreads)
{
    GlobalAverageElocTarget=0;
    NodeAverageElocTarget=0;
    double NodeNormalization = 0;
    double Normalization = 0;
    vector <double> UnnormalizedP(param.nSperNode);
    
    
    bool Proceed;
    double ThreadNormalization;
    unsigned int offset;
    int ThreadMagnetization;
    unsigned long long int Mask, spinConfig;
    std::complex <double> ThreadAverageTargetEnergy;
    vector <std::complex <double> > ThreadAverage(param.nNetworkParam + 1);
    vector <std::complex <double>> Theta(param.m);
    vector <int> v(param.n);
    vector <std::complex <double>> a(a0.size());
    vector <std::complex <double>> b(b0.size());
    vector <std::complex <double>> W(W0.size());
    
    a = a0;
    b = b0;
    W = W0;
    
    ThreadNormalization =0;
    
    for(unsigned int i = 0; i< NodeAverage.size(); i++){
        NodeAverage[i] = 0;
        GlobalAverage[i] = 0;
    }
    ThreadAverageTargetEnergy = 0;
    for (int i = 0; i< ThreadAverage.size(); i++){
        ThreadAverage[i] = 0;
    }
    
    
    for(unsigned int s=0; s<param.nSperNode; s++){ //s is the Hamming weight
        
        spinConfig = s;
        
        if(s < param.n+1){
            Proceed = true;
        }
        
        ThreadMagnetization = param.n - 2*s;
        
        if(Proceed){
            
            UpdateTheta2HW(param,ThreadMagnetization,b,W,Theta);
            UnnormalizedP[s] = CalculateUnnormalizedProbabilityHW(param,s,ThreadMagnetization,a,Theta);
            
            
        }
        else{
            UnnormalizedP[s] = 0;
        }
        
        ThreadNormalization += UnnormalizedP[s];
        
        MeasureExact4HW(count, update,UnnormalizedP[s],s, param, Theta, v,  a, W, NodeX, NodeEloc, NodeElocTarget, SpinSystem);
        
        offset = s*param.nNetworkParam;
        for (int k=0; k<param.nNetworkParam; k++){
            ThreadAverage[k] += UnnormalizedP[s]*NodeX[k + offset];
        }
        ThreadAverage[param.nNetworkParam] += UnnormalizedP[s]*NodeEloc[s];
        ThreadAverageTargetEnergy += UnnormalizedP[s]*NodeElocTarget[s]; //be careful need square here
    } //for samples
    
    
    
    {
        NodeNormalization += ThreadNormalization;
        for (int k=0; k<ThreadAverage.size(); k++){
            NodeAverage[k] += ThreadAverage[k];
        }
        NodeAverageElocTarget += ThreadAverageTargetEnergy;
    }
    
    
    
    Normalization = NodeNormalization;
    GlobalAverage = NodeAverage;
    GlobalAverageElocTarget = NodeAverageElocTarget;
    GlobalAverageElocTarget = GlobalAverageElocTarget/Normalization;
    
    
    for(int s=0; s<GlobalAverage.size(); s++){
        GlobalAverage[s] = GlobalAverage[s]/Normalization;
    }
    
    //Calculate X and Eloc minus mean
    for(int s=0; s<param.nS; s++){
        NormalizeMeasureExact(UnnormalizedP,Normalization, s, param, GlobalAverage, NodeX, NodeEloc);
    } // for
    
    
    
    GlobalX = NodeX;
    GlobalEloc = NodeEloc;
    param.EnergyTarget = GlobalAverageElocTarget.real();
    
    return GlobalAverageElocTarget.real();
    
    
    
} //SRexact5HW



void UpdateNetworkRKHW(int update, double gamma, SimulationParameters &param, SimpleVector &F, SimpleVector &x, Operator &X, vector <std::complex <double> > &GlobalAverage, std::complex <double> &GlobalAverageElocTarget, vector <std::complex <double>> &a,vector <std::complex <double>> &b,vector <std::complex <double>> &W, std::mt19937_64 &rng)
{
    int count;
    std::complex <double> temp;
    
    
    
    count = 0;
    for (int i = 0; i< 1; i++){
        temp =(x[i + count] + 1i*x[i+count + param.nNetworkParam]) ;
        a[i] += -gamma*temp ;
        
        
    }
    
    count = 1;
    for (int i = 0; i< param.m; i++){
        temp =(x[i + count] + 1i*x[i+count + param.nNetworkParam]) ;
        b[i] += -gamma*temp ;
        
        
    }
    
    count = param.m + 1;
    
    for (int i = 0; i< param.m; i++){
        temp =(x[i + count] + 1i*x[i + count + param.nNetworkParam]);
        
        W[i] += -gamma*temp;
        
    }
    
    
} //UpdateNetworkRKHW



void UpdateNetworkRK2exactHW(int rank, int update, SimulationParameters &param, SimpleVector &F, SimpleVector &x1, SimpleVector &x2, SimpleVector &x3, SimpleVector &x4, SimpleVector &x5, Operator &X, vector <std::complex <double> > &GlobalAverage, std::complex <double> &GlobalAverageElocTarget, vector <std::complex <double>> &a,vector <std::complex <double>> &b,vector <std::complex <double>> &W, std::mt19937_64 &rng)
{
    int count;
    double magnitude,gamma;
    std::complex <double> temp;
    
    SimpleVector y(2*param.nNetworkParam);
    SimpleVector y2(2*param.nNetworkParam);
    
    
    for(unsigned int i = 0; i< 2*param.nNetworkParam; i++){
        y[i] = param.gamma*(0.5*(x1[i]+x2[i])-0.25*(x1[i]+x3[i]+x4[i]+x5[i]));
    }
    
    //Calculate error delta
    y2=0;
    X.Apply(y,y2); //calculate y2 = S*y
    param.Gamma = 0;
    for (int i = 0; i< param.nNetworkParam; i++){
        param.Gamma += (y[i])*(y2[i]) + (y[i+param.nNetworkParam])*(y2[i+param.nNetworkParam]);
    }
    param.Gamma = pow(param.Gamma,0.5)/param.nNetworkParam;
    
    if(param.Gamma<pow(10.0,-12)){
        InitializeRBMParameters2(0,0.01,a,b,W,rng);
        param.RunningCount=1;
        param.RunningTargetE = 0;
        param.gamma = param.gammaInit;
        
    }
    else{
        
        param.RunningCount++;
        param.RunningTargetE += GlobalAverageElocTarget.real();
        
        gamma = param.gamma/2.0;
        count = 0;
        for (int i = 0; i< 1; i++){
            temp =(x1[i + count] + 1i*x1[i+count + param.nNetworkParam]) ;
            a[i] += -gamma*temp ;
            
            temp =(x2[i + count] + 1i*x2[i+count + param.nNetworkParam]) ;
            a[i] += -gamma*temp ;
            
        }
        
        count = 1;
        for (int i = 0; i< param.m; i++){
            temp =(x1[i + count] + 1i*x1[i+count + param.nNetworkParam]) ;
            b[i] += -gamma*temp ;
            
            temp =(x2[i + count] + 1i*x2[i+count + param.nNetworkParam]) ;
            b[i] += -gamma*temp ;
        }
        
        count = param.m + 1;
        for (int i = 0; i< param.m; i++){
            
            temp =(x1[i + count] + 1i*x1[i +  count + param.nNetworkParam]);
            W[i] += -gamma*temp;
            
            temp =(x2[i + count] + 1i*x2[i +  count + param.nNetworkParam]);
            W[i] += -gamma*temp;
            
        }
        
        
        param.gamma = min(param.gamma*pow(6.0*param.Delta/param.Gamma,1/3.0),10.0);
        
        
    }//else
    
} //UpdateNetworkRK2exactHW

