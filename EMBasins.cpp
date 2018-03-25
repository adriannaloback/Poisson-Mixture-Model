//_________________________________________________________________________
//  main.cpp
//  EM_PoissonMM
//  Description: Fits Poisson Mixture Model to input population
//               spike time data via batch EM algorithm.
//  Input Arguments:
//      (1) CV = 1 : Perform k-fold cross-validation & return CV-LL
//       or CV = 0 : Instead, train on all data & return fit params
//
//  Created by adrianna on 4/20/17.
//  Copyright © 2017 adrianna. All rights reserved.
//
//  04/27/2017 - Verified that the cross-validation is fully functional.
//_________________________________________________________________________

#include "EMBasins.h"
#include "BasinModel.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include <exception>

using namespace std;

// Selects which basin model to use
typedef PoissonBasin BasinType;

int main(int argc, char *argv[]) {
    
    // ** Command-Line Argument & Initializations: **
    bool CV        = argv[2];  //user-specified (in Xcode, use argv[1] for true, argv[2] for false)
    int N          = 152;      //total population size (# of neurons)
    double binsize = 200;      //time bin width (units: 10 kHz)
    int nbasins    = 2;        //total # of latent states in the model
    int niter      = 100;      //specified # of iterations of EM until termination
    
    string datafile = "ST_NRN_N152.txt"; //File w/ all population spike time data (input)
    
    // ** Poisson Mixture Model **
    cout << "Initializing EM object..." << endl;
    EMBasins<BasinType> basin_obj(datafile, binsize, nbasins, N, CV); //Instantiate
    
    if (CV==1) { //Perform cross-validation:
        cout << "Performing k-fold cross-validation..." << endl;
        int kfolds = 2;
        vector<double> cv_logli = basin_obj.crossval(niter, kfolds, CV);
    
        // ** Write CV-LL Results to Output File **
        string outfile = "CVLL_NRN_N152_nb5_2foldcv.txt";  //File will have the cv LL results (output)
        ofstream cv_output(outfile);
        for (int i=0; i<kfolds; i++) {
            for (int j=0; j<niter; j++) {
                cv_output << cv_logli[(i*niter) + j] << " ";
            }
            cv_output << endl;
        }
        cv_output.close(); //Note: Will just need to take the transpose of this matrix in Matlab
    }
    else { //Train on *all* data & return fit model params:
        cout << "Training mixture model on all data..." << endl;
        vector<double> train_logli  = basin_obj.train(niter, CV); //train_logli = log(L(D;\theta^*)) for all data D
        vector<paramsStruct> params = basin_obj.basin_params();   //get fit model params, \theta^*
        
        // ** Write Fit Params \Lambda^* to Output File **
    }

    return 0;
}

// ***** RNG Methods *****
//Constructor:
RNG::RNG() {
    // Initialize mersenne twister RNG
    rng_pr = gsl_rng_alloc(gsl_rng_mt19937);
}
//Destructor:
RNG::~RNG() {
    gsl_rng_free(rng_pr);
}

int RNG::poisson(double rate) { //NEEDS TO BE FINISHED!
    return rate;
}

vector<int> RNG::randperm(int nmax) {
//Description: This fn is analogous to Matlab's randperm(nmax)
    vector<int> nvals (nmax);
    for (int i=0; i<nmax; i++) {
        nvals[i] = i;
    }
    for (int i=0; i<nmax; i++) {
        // select random integer ix between i and nmax-1
        // swap i with ix
        unsigned long ix = i + gsl_rng_uniform_int(rng_pr, nmax-i);
        int tmp   = nvals[i];
        nvals[i]  = nvals[ix];
        nvals[ix] = tmp;
    }
    return nvals;
}

// ***** EMBasins Methods *****
//Constructor #1:
template <class BasinT>
EMBasins<BasinT>::EMBasins(int N, int nbasins) : N(N), nbasins(nbasins), w(nbasins) {
    rng = new RNG();
    srand(0);
}
//Constructor #2:
template <class BasinT>
EMBasins<BasinT>::EMBasins(const string& filename, double binsize, int nbasins, int N, bool CV) : nbasins(nbasins),nsamples(0),w(nbasins),N(N) {
    
    rng = new RNG();
    srand(0);
    
    //Load the population spike time data from the input filename:
    cout << "Loading in spike time data..." << endl;
    vector<Spike> all_spikes; //Init cache
    ifstream infile;
    infile.open(filename);
    double st;
    int nidx;
    
    while( infile >> st >> nidx ) {
        Spike s;
        s.bin = floor(st/binsize);
        s.neuron_ind = nidx;
        all_spikes.push_back(s);
    }
    infile.close();
    
    //Now sort spikes to be in chronological order:
    all_spikes = sort_spikes(all_spikes); //I verified on 4/26/17 that this works correctly
    cout << "Sorted spikes..." << endl;
    
    //Now add the all-silent population response with frequency of 0:
    string silent_str(N,'0');
    State this_state;
    string this_str = silent_str;
    this_state.freq = 0;
    this_state.P.assign(nbasins, 0);
    this_state.weight.assign(nbasins, 0);
    this_state.word.assign(N,0);
    
    all_states.insert(pair<string,State> (silent_str,this_state));  //protected data of EMBasins obj
    test_states.insert(pair<string,State> (silent_str,this_state)); //protected data of EMBasins obj
    if (CV==1) {
        train_states.insert(pair<string,State> (silent_str,this_state)); //protected data of EMBasins obj
    }
    
    int curr_bin = 0;
    for (vector<Spike>::iterator it=all_spikes.begin(); it!=all_spikes.end(); ++it) {
        int next_bin  = it-> bin;
        int next_cell = it->neuron_ind;
        
        //Go through this following portion when have completed adding all
        //active neurons constituting this response pattern:
        if (next_bin > curr_bin){
            //Add new response pattern; if it's already been discovered, increment its frequency:
            this_state.active_constraints = BasinT::get_active_constraints(this_state);
            
            raster.push_back(this_str); //raster is protected data of EMBasins (type vector<string>)
            
            pair<state_iter, bool> ins = all_states.insert(pair<string,State> (this_str,this_state));
            if (!ins.second) {
                (((ins.first)->second).freq)++;
            } //map object has *unique* keys, and map::insert() checks that the key of the pair to insert is unique
            
            //All responses between curr_bin and next_bin (exclusive) are silent;
            //update frequency of all-silent response accordingly:
            for (int i=0; i<(next_bin-curr_bin-1); i++) {
                raster.push_back(silent_str);
            }
            all_states[silent_str].freq += (next_bin - curr_bin - 1);
            
            // Reset response \vec{y}(t) and jump to next time bin:
            this_str = silent_str;
            this_state.freq = 1;
            this_state.on_neurons.clear();
            this_state.P.assign(nbasins,0);
            this_state.weight.assign(nbasins,0);
            this_state.word.assign(N,0);
            
            curr_bin = next_bin;
        }//end outer if
        
        //Add next_cell to this population response (enable Poisson var representation, i.e. \in \Z_\ge0):
        //(Still in same time bin, so update the response for this time bin accordingly)
        if (this_state.word[next_cell] == 0) {
            this_state.on_neurons.push_back(next_cell); //unique set of active neurons defining this pop. response
        }
        this_state.word[next_cell]++;
        this_str[next_cell] = this_state.word[next_cell];
    }//end for
    
    //Remove the all-silent response if it never occurred (VERY unlikely for our data, but just in case):
    if (all_states[silent_str].freq == 0) {
        all_states.erase(silent_str);
    }
    
    //Note: Now all_states contains all pop. responses found in the data together with their frequencies.
    //      And raster contains the full chronological data, D \in \Z^{N x nsamples}. (Verified correct 4/26/17).

    if (CV==0) {
        //Compute the total # of samples & set nsamples:
        for (state_iter it=all_states.begin(); it!=all_states.end(); ++it) {
            nsamples += (it->second).freq;
        }
        train_states = all_states; //protected data of EMBasins
    }
}//end 2nd EMBasins constructor

//Destructor:
template <class BasinT>
EMBasins<BasinT>::~EMBasins() {
    delete rng;
}

//EMBasins::sort_spikes
//Sorts input spike times in ascending chronological order (a_spikes is map object)
template <class BasinT>
vector<Spike> EMBasins<BasinT>::sort_spikes(vector<Spike>& a_spikes) { //use call-by-reference to modify input vector
    sort(a_spikes.begin(), a_spikes.end(), SpikeComparison());
    return a_spikes;
}

//EMBasins::crossval
//Note: vector<string> raster = D (all observed population responses in chronological order)
//Note: raster.size = total # of times bins (T)
template <class BasinT>
vector<double> EMBasins<BasinT>::crossval(int niter, int k, bool CV) {
    int blocksize = floor(raster.size() / k);
    vector<string> check = raster; //debug
    // -- Generate random permutation of time bins: --
    vector<int> tperm = rng->randperm((int)raster.size()); //tperm = random permutation vec of 1:raster.size()
    vector<double> all_logli (k*niter);
    
    for (int cv_i=0; cv_i<k; cv_i++) {
        // -- Populate train_states and test_states: --
        train_states.clear();
        test_states.clear();
        for (int t=0; t<tperm.size(); t++) {
            string this_str  = raster[tperm[t]];
            State this_state = all_states[this_str]; //returns corresponding State obj for key this_str
            this_state.freq  = 1;                    //will only be saved if this_state is not already in curr_map
            map<string, State>& curr_map = (t < cv_i*blocksize || t >= (cv_i+1)*blocksize)
            ? train_states : test_states;
            pair<state_iter, bool> ins = curr_map.insert(pair<string,State> (this_str,this_state));
            if (!ins.second) {
                (((ins.first)->second).freq)++;
            }
        }//inner for
        
        // -- Compute the size of the training data set: --
        nsamples = 0; //clear from constructor value
        for (state_iter it=train_states.begin(); it!=train_states.end(); ++it) {
            nsamples += (it->second).freq;
        }
        
        // -- Train (returns logli of *test* set): --
        vector<double> cv_logli = train(niter, CV);
        for (int j=0; j<niter; j++) {
            all_logli[cv_i*niter + j] = cv_logli[j];
        }
    }//outer for (over cv folds)
    return all_logli;
}

//EMBasins::train
//Note: This function is called by EMBasins::crossval
//Learns the model parameters (that minimize KL divergence) via batch EM algorithm
template <class BasinT>
vector<double> EMBasins<BasinT>::train(int niter, bool CV) {
    
    cout << "Initializing the mixing weight parameters..." << endl;
    w.assign(nbasins,1/(double)nbasins); //Initialize mixing weights to uniform
    
    basins.clear();
    
    // Initialize each basin model
    for (int i=0; i<nbasins; i++) {
        basins.push_back(BasinT(N,i,rng));
    }
    cout << "Initialized each basin model..." << endl;
    
    update_P(); //Note: update_P() updates the posterior this_state.P[a] = p(a|y(t);\theta^(p)) for each pop response y(t)
                //in the training data & each latent mode a, using the pth (i.e. initialized at first) param values.
                //It also updates this_state.weight[a] = this_state.freq * this_state.P[a].
                //It can also return the training LL if the output argument is requested.
    cout << "Updated the posterior for each training response & latent mode a..." << endl;
    
    test_logli.assign(niter,0); //will cache the CV-LL values, i.e. log(L(D_test;\theta_train)) for each iteration
    vector<double> logli(niter);
    for (int p=1; p<=niter; p++) {
        cout << "Iteration " << p << endl;
        
        // ** E-Step: **
        //---------------------------------------------------
        for (int a=0; a<nbasins; a++) {
            basins[a].reset_stats();
            for (state_iter it = train_states.begin(); it!=train_states.end(); ++it) {
                basins[a].increment_stats(it->second);
            }
            basins[a].normalize_stats();
        }
        
        // ** M-Step **
        //---------------------------------------------------
        double alpha = 0;
        //        if (p >= niter/2) {
        //            alpha = 0.002 + (1-0.002)*exp(-(double) (i-niter/2) * (10.0/(((double)(niter/2)-1))));
        //        }
        //        cout << alpha << endl;
        for (int a=0; a<nbasins; a++) {
            basins[a].doMLE(alpha); //updates the \lambda_{a,i} values; note: alpha not actually used
        }
        update_w();                 //verified correct on 4/26/2017
        
        logli[p-1]      = update_P();      //update the posterior using the new param values, & return training LL
        test_logli[p-1] = update_P_test(); //=log L(D_test; \theta^(p))
    }
    if (CV==1) { return test_logli; }
    else { return logli; }
}

//EMBasins:set_state_P
//Computes the posterior p(\alpha|y(t);\theta^(p)) & updates this_state.weight[\alpha]
template <class BasinT>
double EMBasins<BasinT>::set_state_P(State& this_state) {
    double Z = 0;
    for (int a=0; a<nbasins; a++) {
        this_state.P[a] = w[a] * basins[a].P_state(this_state); //P_state() computes the observation model value
        Z += this_state.P[a];
    }
    for (int a=0; a<nbasins; a++) {
        this_state.P[a] /= Z; //now = the value of the posterior, p(a|y(t);\theta^(p))
        this_state.weight[a] = this_state.freq * this_state.P[a]; //will use these values to compute the norm, S(a)
    }
    this_state.pred_prob = Z;
    return Z; //Z is the denominator of the posterior p(a|y(t);\theta^(p)), aka marginal p(y(t);\theta^(p))
}

//EMBasins::update_P
//Computes and returns the *training* log-likelihood (i.e. log L(D_train; \theta) )
template <class BasinT>
double EMBasins<BasinT>::update_P() {
    double logli = 0;
    for (state_iter it=train_states.begin(); it != train_states.end(); ++it) {
        State& this_state = it->second;     //acquire the State object from the map iterator
        double Z = set_state_P(this_state); //Z = marginal p(y(t);\theta^(p)) = \sum_{a=1}^nbasins w[a] * p(y(t)|a;theta^(p))
                                            //set_state_P() also updates this_state.weight[a]
        logli += log(Z);
    }
    return logli; //TRAIN LL = l(D_train;\theta^(p)) := log( L(D_train;\theta^(p)) )
}

//EMBasins::update_w
template <class BasinT>
void EMBasins<BasinT>::update_w() {
    //cout << "nsamples = " << nsamples << endl; //for debugging
    for (int a=0; a<nbasins; a++) {
        w[a] = ( basins[a].get_norm() / nsamples );
    }
    //_________________________________________________
    //Debugging: Check that mixing weights sum to 1:
    //double check = 0;
    //for (int a=0; a<nbasins; a++) {
    //    check += w[a];
    //}
    //cout << "sum_a w[a] = " << check << endl;
    //_________________________________________________
    return;
}

//EMBasins::update_P_test
//Computes and returns the *test* log-likelihood
template <class BasinT>
double EMBasins<BasinT>::update_P_test() {
    double cv_logli = 0;
    for (state_iter it=test_states.begin(); it != test_states.end(); ++it) {
        State& this_state = it->second;
        double Z = set_state_P(this_state); //Z = marginal p(y(t);\theta^(p)) =
                                            //\sum_{a=1}^nbasins w[a] * p(y(t)|a;theta^(p))
        cv_logli += log(Z);
    }
    return cv_logli; //TEST LL = l(D_test;\theta^(p)) := log( L(D_test;\theta^(p)) )
}

//EMBasins::basins_params
template <class BasinT>
vector<paramsStruct> EMBasins<BasinT>::basin_params() {
    //    int nparams = basins[0].nparams();
    vector<paramsStruct> params (nbasins);
    
    for (int i=0; i<nbasins; i++) {
        params[i] = basins[i].get_params();
    }
    
    return params;
}
