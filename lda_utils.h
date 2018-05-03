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

inline double logDirichlet_const(double alpha, int k);

double logDirichlet_vector(vector<double> alpha);

double getLogLikelihood(int* wordTopicTable, int* docTopicTable, double alpha, double beta, int numWords, int numDocs, int numTopics, int process_id, int process_count); 
