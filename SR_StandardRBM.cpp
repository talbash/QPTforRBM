struct SimulationParameters {
    int updatesPerSwap,n,m,nNetworkParam, alpha,nS,nSperNode,nSweepsPerSwap,nUpdates,RepVersion,max_iter,nThreads,nNodes,RunningCount,ExactSampling;
    double beta,gammaInit,gammaFinal,tol,Delta,lambda0,b0,lambdamin,minReplicaT,Energy,EnergyTarget,Gamma,DeltaE,fnorm,gamma,gammaDecay,RunningE,RunningTargetE, RunningExactTargetE, InterpolationParameter;
    unsigned long long int seed,seedThreads;
    unsigned int FixedHammingWeight, HammingWeight, replica;
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
    
    if(param.FixedHammingWeight){
        for (int i = 0; i<nThreads; i++){
            offsetv = i*n;
            it1 = v.begin() + offsetv;
            it2 = v.begin() + offsetv + n;
            for (int k =0; k<n; k++){

                if(k<param.HammingWeight){
                    v[k + offsetv] = -1;
                }
                else{
                    v[k + offsetv] = 1;
                }
            }//for k
            

            std::shuffle (it1, it2,rng[i]);
        }//for i
    }
    else{
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
    }
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


void MeasureExact(double sqrtProb, int s, SimulationParameters &param, vector <std::complex <double>> &Theta, vector <int> &v,   vector <std::complex <double>>  &a, vector <std::complex <double>>  &W, vector <std::complex <double>> &X, vector <std::complex <double>> &Eloc, vector <std::complex <double>> &ElocTarget, spin &SpinSystem)
{
    
    //identical to MeasureParallel2 except multiplies by the unnormalized probability
    
    int m = param.m;
    int n = param.n;
    std::complex <double>  piCoshTheta,temp,dtemp, atemp, Wtemp;
    vector <int> v2;
    int j1,j2,shift1,shift2,index, index3, index2,index1,offsetW,offsetW1,offsetW2,offseta,offsetb;
    double A,B,AnnealParameter,tempI, tempR, Jzzz,Jzz, Jxx, Jyy, hx, hy, hz,spin,spin1,spin2,spin3,phase1,phase2,eigT;
    
    
    int offsetParameters = (n+m + n*m);
    
    tempR = 0;
    tempI = 0;
    
    
    AnnealParameter = param.InterpolationParameter;
    A = AnnealParameter;
    B=1;
    
    /* Calculate expectation value of target Hamiltonian */
    for (int k=0; k<SpinSystem.OpStrength.size(); k++){
        
        temp = 1;
        
        for(unsigned int i = 0; i< SpinSystem.OpIz[k].size(); i++){
            index = SpinSystem.OpIz[k][i];
            spin = v[index];
            temp = temp*pow(-1.0,(1.0-spin)/2.0);
        }
        
        if(SpinSystem.OpIy[k].size() > 0){
            temp = temp*pow(1i,SpinSystem.OpIy[k].size());
            for(unsigned int i = 0; i< SpinSystem.OpIy[k].size(); i++){
                index = SpinSystem.OpIy[k][i];
                spin = v[index];
                temp = temp*pow(-1.0,(1.0+spin)/2.0);
            }
        }
        
        atemp = 0;
        
        for(unsigned int i = 0; i< SpinSystem.OpIx[k].size(); i++){
            index = SpinSystem.OpIx[k][i];
            spin = v[index];
            atemp += a[index]*spin;
            
        }
        for(unsigned int i = 0; i< SpinSystem.OpIy[k].size(); i++){
            index = SpinSystem.OpIy[k][i];
            spin = v[index];
            atemp += a[index]*spin;
        }
        temp = temp*exp(-2.0*atemp);
        
        piCoshTheta = 1;
        for(unsigned int j = 0; j< m; j++){
            Wtemp = 0;
            for(unsigned int i = 0; i< SpinSystem.OpIx[k].size(); i++){
                index = SpinSystem.OpIx[k][i];
                
                spin = v[index];
                offsetW = index*m;
                
                Wtemp += W[j + offsetW]*spin;
            } //for i
            for(unsigned int i = 0; i< SpinSystem.OpIy[k].size(); i++){
                index = SpinSystem.OpIy[k][i];
                //cout << index << "\n";
                spin = v[index];
                offsetW = index*m;
                Wtemp += W[j + offsetW]*spin;
            } //for i
            piCoshTheta = piCoshTheta*cosh(Theta[j] - 2.0*Wtemp)/cosh(Theta[j]);
        } //for j
        
        temp = temp*piCoshTheta;
        temp = SpinSystem.OpStrength[k]*temp;
        tempR += temp.real() ;
        tempI += temp.imag();
        
        
    } //for q
    
    ElocTarget[s] = SpinSystem.Eoffset + tempR + 1i*tempI;
    
    
    /* Calculate expectation value of mixer Hamiltonian */
    tempR = 0;
    tempI = 0;
    for (int k=0; k<SpinSystem.DriverStrength.size(); k++){
        
        temp = 1;
        
        for(unsigned int i = 0; i< SpinSystem.DriverIz[k].size(); i++){
            index = SpinSystem.DriverIz[k][i];
            spin = v[index];
            temp = temp*pow(-1.0,(1.0-spin)/2.0);
        }
        
        if(SpinSystem.DriverIy[k].size() > 0){
            temp = temp*pow(1i,SpinSystem.DriverIy[k].size());
            for(unsigned int i = 0; i< SpinSystem.DriverIy[k].size(); i++){
                index = SpinSystem.DriverIy[k][i];
                spin = v[index];
                temp = temp*pow(-1.0,(1.0+spin)/2.0);
            }
        }
        
        atemp = 0;
        
        for(unsigned int i = 0; i< SpinSystem.DriverIx[k].size(); i++){
            index = SpinSystem.DriverIx[k][i];
            spin = v[index];
            atemp += a[index]*spin;
            
        }
        for(unsigned int i = 0; i< SpinSystem.DriverIy[k].size(); i++){
            index = SpinSystem.DriverIy[k][i];
            spin = v[index];
            atemp += a[index]*spin;
        }
        temp = temp*exp(-2.0*atemp);
        
        piCoshTheta = 1;
        for(unsigned int j = 0; j< m; j++){
            Wtemp = 0;
            for(unsigned int i = 0; i< SpinSystem.DriverIx[k].size(); i++){
                index = SpinSystem.DriverIx[k][i];
                spin = v[index];
                offsetW = index*m;
                Wtemp += W[j + offsetW]*spin;
            } //for i
            for(unsigned int i = 0; i< SpinSystem.DriverIy[k].size(); i++){
                index = SpinSystem.DriverIy[k][i];
                spin = v[index];
                offsetW = index*m;
                Wtemp += W[j + offsetW]*spin;
            } //for i
            piCoshTheta = piCoshTheta*cosh(Theta[j] - 2.0*Wtemp)/cosh(Theta[j]);
        } //for j
        
        temp = temp*piCoshTheta;
        temp = SpinSystem.DriverStrength[k]*temp;
        tempR += temp.real() ;
        tempI += temp.imag();
        
        // cout << temp << "  ";
        
    } //for q
    
    
    //remember that minus sign is included in the definition of the mixer Hamiltonian
    //the appropriate term for the heisenberg model should be \sum_{i<j} (0.5*Id - 0.5*(XX + YY + ZZ))/(n).  This way the eigenvalues of each term are always 0 and 1.
    if(param.FixedHammingWeight){
        Eloc[s] = B*ElocTarget[s] + (SpinSystem.MaxStrength)*A*(0.25*(param.n-1) + (0.5/param.n)*(tempR + 1i*tempI));
    }
    else{
        Eloc[s] = B*ElocTarget[s] + (SpinSystem.MaxStrength)*A*(0.5*param.n + 0.5*(tempR + 1i*tempI));
    }
    
    /* Calculate gradients with respect to RBM parameters */
    offseta = 0 + s*offsetParameters;
    offsetb = n + s*offsetParameters;
    offsetW = n+m + s*offsetParameters;
    
    
    for (int k=0; k<offsetParameters; k++){
        
        if(k<n){
            index = k;
            index1 = index;
            
            index += offseta;
            
            X[index] = ((double) v[index1]);
        }
        else if(k<n+m){
            
            index = k-n;
            index1 = index;
            
            index += offsetb;
            
            X[index] = tanh(Theta[index1]);
            
        }
        else{
            index = k-m-n;
            index2 = index/m;
            index1 = index - index2*m;
            index += offsetW;
            
            X[index] = ((double) v[index2])*tanh(Theta[index1]);
            
        }
        
        
    } //for k
    
    
} //MeasureExact




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
    
    Case.append("_StandardRBM_n=");
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
        OutputFile << "Perform ANN Replica swap attempts at every " << param.updatesPerSwap/2 << "\n";
        OutputFile << "Covariance matrix regularization max(lambda_0 b_0^p,lambda_min)\n";
        OutputFile << "lambda_0:\t" << param.lambda0 << "\n";
        OutputFile << "b_0:\t" << param.b0 << "\n";
        OutputFile << "lambda_min:\t" << param.lambdamin << "\n";
        
        OutputFile << "Wavefunction change:\t" << param.Delta << "\n";
        OutputFile << "Parallel tempering inverse temperate:\t";
        for (int i = 0; i< Beta.size(); i++){
            OutputFile << Beta[i] << "\t";
        }
        OutputFile << "\n";
        OutputFile << "Interpolating Hamiltonian Gamma H_TF +  H_Target: \t";
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
        
        if (param.FixedHammingWeight){
            OutputFile << "Sampling restricted to fixed Hamming weight: " << param.HammingWeight << "\n";
        }
        
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



double CalculateUnnormalizedProbability(vector <int> v, vector <std::complex <double>> a, vector <std::complex <double>> Theta)
{
    std::complex <double> temp,p;
    
    temp = 0.0;
    for (unsigned int i = 0; i<v.size(); i++){
        temp += ((double)v[i])*a[i];
    }
    
    p = exp(temp);
    
    for (unsigned int i = 0; i< Theta.size(); i++){
        p = p*cosh(Theta[i]);
    }
    
    return pow(abs(p),2.0);
    
} //CalculateUnnormalizedProbability


double SRexact(SimulationParameters &param, spin &SpinSystem, vector <std::complex <double>> &a0, vector <std::complex <double>> &b0,vector <std::complex <double>> &W0, vector <int> &Nodev, vector <std::complex <double> > &NodeAverage,  std::complex <double> &NodeAverageElocTarget, vector <std::complex <double> > &NodeX, vector <std::complex <double> > &NodeEloc, vector <std::complex <double> > &NodeElocTarget, vector <std::complex <double> > &GlobalAverage, std::complex <double> &GlobalAverageElocTarget, vector <std::complex <double> > &GlobalX, vector <std::complex <double> > &GlobalEloc)
{
    GlobalAverageElocTarget=0;
    NodeAverageElocTarget=0;
    double NodeNormalization = 0;
    double Normalization = 0;
    vector <double> UnnormalizedP(pow(2,param.n));
    
    
    bool Proceed;
    double ThreadNormalization;
    unsigned int offset,sPrevious;
    unsigned long long int Mask, spinConfig, c,r;
    std::complex <double> ThreadAverageTargetEnergy;//,ThreadAverageTargetEnergy2;
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
    
    
    sPrevious = 0;
    spinConfig = 0;
    if(param.FixedHammingWeight){
        for (unsigned int k = 0; k<param.HammingWeight; k++){
            spinConfig = spinConfig << 1;
            spinConfig += 1;
            
        }
    }
    
    
    
    for(unsigned int s=0; s<param.nS; s++){
        
        
        if(param.FixedHammingWeight){
            for (unsigned int k = sPrevious; k<s; k++ ){
                c = spinConfig & (~spinConfig +1);
                r = spinConfig + c;
                spinConfig = (((r ^ spinConfig) >> 2)/c) | r;
            }
            sPrevious = s;
            
            if(__builtin_popcount(spinConfig) == param.HammingWeight && spinConfig<pow(2,param.n)){
                Proceed = true;
            }
            else{
                Proceed = false;
                cout << "This should never happen now\n";
            }
        }
        else{
            spinConfig = s;
            if(spinConfig<pow(2,param.n)){
                Proceed = true;
            }
            else{
                Proceed = false;
            }
        }
        
        
        
        
        for (int k= 0; k< param.n; k++){
            Mask = 1;
            Mask = Mask << (param.n-1-k);
            if ( (spinConfig & Mask) !=0){
                v[k]=-1;
            }
            else{
                v[k]=1;
            }
        }
        
        if(Proceed){
            
            UpdateTheta2(param,v,b,W,Theta);
            UnnormalizedP[s] = CalculateUnnormalizedProbability(v,a,Theta);
            
            
            
        }
        else{
            UnnormalizedP[s] = 0;
        }
        
        
        ThreadNormalization += UnnormalizedP[s];
        
        
        MeasureExact(UnnormalizedP[s],s, param, Theta, v,  a, W, NodeX, NodeEloc, NodeElocTarget,SpinSystem);
        
        offset = s*param.nNetworkParam;
        for (int k=0; k<param.nNetworkParam; k++){
            ThreadAverage[k] += UnnormalizedP[s]*NodeX[k + offset];
        }
        
        
        ThreadAverage[param.nNetworkParam] += UnnormalizedP[s]*NodeEloc[s];
        ThreadAverageTargetEnergy += UnnormalizedP[s]*NodeElocTarget[s]; //be careful need square here
    } //for samples
    
    
    NodeNormalization += ThreadNormalization;
    for (int k=0; k<ThreadAverage.size(); k++){
        NodeAverage[k] += ThreadAverage[k];
    }
    NodeAverageElocTarget += ThreadAverageTargetEnergy;
    
    
    
    
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
    
    
} //SRexact5


unsigned long long int choose(unsigned int n, unsigned int k)
{
    if (k == 0) {
        return 1;
    }
    else{
        return (n * choose(n - 1, k - 1)) / k;
    }
}






void UpdateNetworkRK(int update, double gamma, SimulationParameters &param, SimpleVector &F, SimpleVector &x, Operator &X, vector <std::complex <double> > &GlobalAverage, std::complex <double> &GlobalAverageElocTarget, vector <std::complex <double>> &a,vector <std::complex <double>> &b,vector <std::complex <double>> &W, std::mt19937_64 &rng)
{
    int count;
    std::complex <double> temp;
    
    
    
    count = 0;
    for (int i = 0; i< param.n; i++){
        temp =(x[i + count] + 1i*x[i+count + param.nNetworkParam]) ;
        a[i] += -gamma*temp ;
        
        
    }
    
    count = param.n;
    for (int i = 0; i< param.m; i++){
        temp =(x[i + count] + 1i*x[i+count + param.nNetworkParam]) ;
        b[i] += -gamma*temp ;
        
        
    }
    
    count = param.m + param.n;
    
    for (int i = 0; i< param.m; i++){
        for (int j = 0; j< param.n; j++){
            temp =(x[i + j*param.m + count] + 1i*x[i + j*param.m + count + param.nNetworkParam]);
            
            W[i + j*param.m] += -gamma*temp;
            
            
        }
    }
    
    
    
    
} //UpdateNetworkRK


void UpdateNetworkRK2exact(int update, SimulationParameters &param, SimpleVector &F, SimpleVector &x1, SimpleVector &x2, SimpleVector &x3, SimpleVector &x4, SimpleVector &x5, Operator &X, vector <std::complex <double> > &GlobalAverage, std::complex <double> &GlobalAverageElocTarget, vector <std::complex <double>> &a,vector <std::complex <double>> &b,vector <std::complex <double>> &W, std::mt19937_64 &rng)
{
    int count;
    double magnitude,gamma;
    std::complex <double> temp;
    
    SimpleVector y(2*param.nNetworkParam);
    SimpleVector y2(2*param.nNetworkParam);
    
    //Calculate norm of force vector.
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
        param.RunningE += GlobalAverage[param.nNetworkParam].real();
        param.RunningTargetE += GlobalAverageElocTarget.real();
        
        gamma = param.gamma/2.0;
        count = 0;
        for (int i = 0; i< param.n; i++){
            temp =(x1[i + count] + 1i*x1[i+count + param.nNetworkParam]) ;
            a[i] += -gamma*temp ;
            
            temp =(x2[i + count] + 1i*x2[i+count + param.nNetworkParam]) ;
            a[i] += -gamma*temp ;
            
        }
        
        count = param.n;
        for (int i = 0; i< param.m; i++){
            temp =(x1[i + count] + 1i*x1[i+count + param.nNetworkParam]) ;
            b[i] += -gamma*temp ;
            
            temp =(x2[i + count] + 1i*x2[i+count + param.nNetworkParam]) ;
            b[i] += -gamma*temp ;
        }
        
        count = param.m + param.n;
        for (int i = 0; i< param.m; i++){
            for (int j = 0; j< param.n; j++){
                temp =(x1[i + j*param.m + count] + 1i*x1[i + j*param.m + count + param.nNetworkParam]);
                W[i + j*param.m] += -gamma*temp;
                
                temp =(x2[i + j*param.m + count] + 1i*x2[i + j*param.m + count + param.nNetworkParam]);
                W[i + j*param.m] += -gamma*temp;
            }
        }
        
        
        param.gamma = min(param.gamma*pow(6.0*param.Delta/param.Gamma,1/3.0),10.0);
        
        
    }//else
    
} //UpdateNetworkRK2exact

