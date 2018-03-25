//                                                            |
//  BasinModel.cpp                                            |
//                                                            |
//  Modified by Adrianna Loback on 4/21/2017 for Poisson MM.  |
//                                                            |
//____________________________________________________________|


#include "BasinModel.h"
#include "EMBasins.h"

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <math.h>


// ***** paramsStruct *****
paramsStruct::paramsStruct() : nfields(0) {}

int paramsStruct::get_nfields() {return nfields;}

void paramsStruct::addField(string name, myMatrix<double>& value) {
    nfields++;
    fields[name] = &value;
    
    fieldNames.clear();
    for (map<string, myMatrix<double>* >::iterator it = fields.begin(); it!=fields.end(); ++it) {
        fieldNames.push_back((it->first).data());
    }
    
    return;
}

const char** paramsStruct::fieldNamesArray() {
    return fieldNames.data();
}

vector<double>* paramsStruct::getFieldData(int index) {
    return fields[fieldNames[index]]->data();
}

int paramsStruct::getFieldN(int index) {
    if (index >= nfields || index < 0) {
        cerr << "Index out of range." << endl;
        return 0;
    } else {
        return fields[fieldNames[index]]->get_N();
    }
}

int paramsStruct::getFieldM(int index) {
    if (index >= nfields || index < 0) {
        cerr << "Index out of range." << endl;
        return 0;
    } else {
        return fields[fieldNames[index]]->get_M();
    }
}

const char* paramsStruct::getFieldName(int index) {
    if (index >= nfields || index < 0) {
        cerr << "Index out of range." << endl;
        return 0;
    } else {
        return fieldNames[index];
    }
}


// ***** BasinModel *****
void BasinModel::reset_stats() {
    norm = 0;
    for (vector<double>::iterator it=stats.begin(); it!=stats.end(); ++it) {
        *it = 0;
    }
    return;
}

void BasinModel::increment_stats(const State& this_state) {
    double wt = this_state.weight[basin_num]; //Acquires this_state.freq * p(\alpha|y(t);\theta(p)), which
                                              //was computed by EMBasins::set_state_P.
    norm += wt;                               //This is equivalent to summing over all time bins in main() fn.

    //Only need to iterate though the non-silent neurons for this given state, i.e. i s.t. y_i(t) > 0:
    for (vector<int>::const_iterator it=this_state.active_constraints.begin(); it!=this_state.active_constraints.end(); ++it)
    {
        int y_i     = this_state.word[*it]; //=y_i(t) > 0
        if (y_i<0) {cerr << "y_i < 0" << endl; ;}
        stats[*it] += (y_i * wt);
    }
    return;
}

void BasinModel::normalize_stats() {
    for (vector<double>::iterator it=stats.begin(); it!=stats.end(); ++it) {
        *it /= norm;
    }
    return;
}

int BasinModel::factorial(int n) const {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


// ***** PoissonBasin *****
//Constructor:
PoissonBasin::PoissonBasin(int N, int basin_num, RNG* rng) : BasinModel(N,basin_num,rng), prefactor(1)
{
    stats.assign(N, 0);
    for (vector<double>::iterator it = stats.begin(); it != stats.end(); ++it) {
        double u = 0.1*((double) rand() / (double) RAND_MAX) + 0.45;
        (*it) = u;
    }
    lambda.assign(stats,N,1); //initializes the lambda_{a,i} values 
}

vector<int> PoissonBasin::get_active_constraints(const State& this_state) {
    return this_state.on_neurons; //returns the *non-silent* neurons for this population response
}

void PoissonBasin::doMLE(double alpha) {
    lambda.assign(stats, N, 1);
    return;
}

double PoissonBasin::P_state(const State& this_state) const {
    //Compute the value of the observation model, p(this_state | \alpha; \theta^(p)):
    double P = prefactor; //initialized in constructor to = 1
    double logP = log(P);

    for (int i=0; i<N; ++i) {
        int y_i = this_state.word[i]; //=y_i(t)
        logP += log( ( exp(-1 * lambda.at(i)) * pow(lambda.at(i),y_i) ) / factorial(y_i) ); //to prevent underflow
        //P *= ( exp(-lambda.at(i)) * pow(lambda.at(i),y_i) ) / factorial(y_i);
    }
    
    P = exp(logP);  //P = p(\vec{y}(t) | \alpha; \theta^{(p)}), i.e. value of observation model
    return P;       //for the given response \vec{y}(t) for this latent state
}

vector<char> PoissonBasin::sample() { //NEED TO FINISH THIS!
    vector<char> this_sample (N);
    //for (int i=0; i<N; i++) {s
    //    this_sample[i] = (rng->bernoulli(m.at(i))) ? 1 : 0;
    //}
    return this_sample;
}

paramsStruct PoissonBasin::get_params() {
    paramsStruct params;
    params.addField("lambda", lambda);
    return params;
}
