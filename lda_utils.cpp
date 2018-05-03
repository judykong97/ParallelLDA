#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <ctime>
#include <string.h>


using namespace std;

/* Reference: taken from stackoverflow: 
 * https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
 */
template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
    return idx;
}

inline double logDirichlet_const(double alpha, int k) {
    return k * lgamma(alpha) - lgamma(k * alpha);
}

double logDirichlet_vector(vector<double> alpha) {
    double sumLogGamma = 0.0;
    double logSumGamma = 0.0;
    for(int i = 0; i < alpha.size(); i++) {
        sumLogGamma += lgamma(alpha[i]);
        logSumGamma += alpha[i];
    }
    return sumLogGamma - lgamma(logSumGamma);
}

double getLogLikelihood(int* wordTopicTable, int* docTopicTable, double alpha, double beta, int numWords, int numDocs, int numTopics, int process_id, int process_count) {
    double lik = 0.0;
    // int numWords = wordTopicTable.size();
    // int numDocs = docTopicTable.size();
    // int numTopics = docTopicTable[0].size();
    vector<double> temp(numWords, 0.0);
    for (int k = 0; k < numTopics; k++) {
        for (int w = 0; w < numWords; w++) {
            temp[w] = beta + wordTopicTable[w * numTopics + k];
        }
        lik += logDirichlet_vector(temp);
        lik -= logDirichlet_const(beta, numWords);
    }
    vector<double> temp2(numTopics, 0.0);
    for (int d = process_id; d < numDocs; d+=process_count) {
        int offset = d * numTopics;
        for (int k = 0; k < numTopics; k++) {
            temp2[k] = alpha + docTopicTable[offset + k];
        }
        lik += logDirichlet_vector(temp2);
        lik -= logDirichlet_const(alpha, numTopics);
    }
    return lik;
}

