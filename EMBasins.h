//                                                            |
//  EMBasins.h                                                |
//                                                            |
//  Modified by Adrianna Loback on 4/20/2017 for Poisson MM.  |
//                                                            |
//____________________________________________________________|

#ifndef ____EMBasins__
#define ____EMBasins__

#include <gsl/gsl_rng.h>

#include <vector>
#include <string>
#include <map>

using namespace std;


// ************ RNG ***************
class RNG
{
public:
    RNG();
    ~RNG();
    int discrete(const vector<double>&);
    //bool bernoulli(double);
    int poisson(double);
    vector<int> randperm(int);
private:
    gsl_rng* rng_pr;
};

struct State
{
    vector<int> active_constraints;
    vector<int> on_neurons;
    vector<double> P;               //vector of posterior probabilities
    vector<double> weight;
    double freq;
    double pred_prob;

    vector<char> word;
};
// *********************************
typedef map<string,State>::iterator state_iter;
typedef map<string,State>::const_iterator const_state_iter;

// ************ Spike ***************
struct Spike
{
    int bin;
    int neuron_ind;
};
// *********************************

// ************ SpikeComparison ***************
class SpikeComparison
{
public:
    bool operator() (const Spike& lhs, const Spike& rhs) const
    {
        return (lhs.bin < rhs.bin);
    }
};
// *********************************

//typedef priority_queue<Spike, vector<Spike>, SpikeComparison> SpikeHeap;

// ************ EMBasins ***************
class paramsStruct;

template <class BasinT>
class EMBasins
{
public:
    EMBasins(int N, int nbasins);
    EMBasins(const string& filename, double binsize, int nbasins, int N, bool CV);
    ~EMBasins();

    vector<double> crossval(int niter, int k, bool cv); //implements k-fold cross-validation
    vector<double> train(int niter, bool cv);           //will be called by crossval or stand-alone by main()
    
    int nstates() const {return all_states.size();};
    vector<unsigned long> state_hist() const;
    vector<unsigned long> test_hist() const;
    vector<double> all_prob() const;
    vector<double> test_prob() const;
    vector<double> P() const;
    vector<paramsStruct> basin_params();
    vector<char> sample(int);
    vector<char> word_list();
    
    vector<double> w;       //1 x nbasins
    vector<double> Lambda;  //N x nbasins
    
    vector<double> test_logli;
    
protected:
    int nbasins;
    int N;
    double nsamples;
    
    map<string, State> all_states;
    map<string, State> train_states;
    map<string, State> test_states;
    
    vector<BasinT> basins;
    
    RNG* rng;
    
    vector<string> raster;
    
    void update_w();
    double update_P();
    double update_P_test();
    
    double set_state_P(State&); //computes posterior p(a|y(t); \theta^(p)) for input y(t)
    vector<Spike> sort_spikes(vector<Spike>&);
};

// *********************************

#endif /* defined(____EMBasins__) */
