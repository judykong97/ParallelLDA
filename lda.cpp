#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <ctime>


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

double getLogLikelihood(vector<vector<int>> wordTopicTable, vector<vector<int>> docTopicTable, double alpha, double beta) {
    double lik = 0.0;
    int numWords = wordTopicTable.size();
    int numDocs = docTopicTable.size();
    int numTopics = docTopicTable[0].size();
    vector<double> temp(numWords, 0.0);
    for (int k = 0; k < numTopics; k++) {
        for (int w = 0; w < numWords; w++) {
            temp[w] = beta + wordTopicTable[w][k];
        }
        lik += logDirichlet_vector(temp);
        lik -= logDirichlet_const(beta, numWords);
    }
    vector<double> temp2(numTopics, 0.0);
    for (int d = 0; d < numDocs; d++) {
        for (int k = 0; k < numTopics; k++) {
            temp2[k] = alpha + docTopicTable[d][k];
        }
        lik += logDirichlet_vector(temp2);
        lik -= logDirichlet_const(alpha, numTopics);
    }
    return lik;
}

void runLDA(vector<vector<int>> w, vector<vector<int>> &z, 
        int numDocs, int numWords, int numTopics, double alpha, double beta, int numIterations) {

    clock_t start;
    double duration;

    vector<vector<int>> docTopicTable(numDocs, vector<int>(numTopics, 0));
    vector<vector<int>> wordTopicTable(numWords, vector<int>(numTopics, 1));
    vector<int> topicTable(numTopics, 0);
    vector<double> p(numTopics, 0.0);

    for (int i = 0; i < numDocs; i++) { // numDocs = w.size();
        vector<int> z_line(w[i].size(), -1);
        z.push_back(z_line);
    }

    for (int d = 0; d < numDocs; d++) {
        for (int j = 0; j < w[d].size(); j++) {
            int word = w[d][j];
            int topic = rand() % numTopics;
            z[d][j] = topic;
            docTopicTable[d][topic]++;
            wordTopicTable[word][topic]++;
            topicTable[topic]++;
        }
    }

    for (int i = 0; i < numIterations; i++) {

        start = clock();

        for (int d = 0; d < numDocs; d++) {
            for (int j = 0; j < w[d].size(); j++) {
                int word = w[d][j];
                int topic = z[d][j];
                docTopicTable[d][topic]--;
                wordTopicTable[word][topic]--;
                topicTable[topic]--;

                double norm = 0.0;
                int newk = topic;
                for (int k = 0; k < numTopics; k++) {
                    int z_dj_equals_k = (k == topic);
                    double ak = docTopicTable[d][k] - z_dj_equals_k + alpha;
                    double bk = (wordTopicTable[word][k] - z_dj_equals_k + beta) / (topicTable[k] - z_dj_equals_k + numWords + beta);
                    norm += ak * bk;
                    p[k] = norm;
                }

                double r = ((double) rand()) / RAND_MAX * norm;
                int lo = 0;
                int hi = numTopics - 1;
                int mid;
                while (lo < hi) {
                    if (hi - lo < 10) {
                        int i;
                        for (i = lo; i <= hi; i++) {
                            if (p[i] > r) {
                                lo = i;
                                hi = lo;
                                break;
                            }
                        }
                        break;
                    } else {
                        mid = lo + (hi - lo)/2;
                        if (r <= p[mid]) {
                            hi = mid;
                        } else {
                            lo = mid + 1;
                        }
                    }
                }
                newk = lo;
                // double sum_p_up_to_k = 0.0;
                // double r = ((double) rand()) / RAND_MAX;
                // for(int k = 0; k < numTopics; k++) {
                //     sum_p_up_to_k += p[k] / norm;
                //     if(r < sum_p_up_to_k) {
                //         newk = k;
                //         break;
                //     }
                // }
                z[d][j] = newk;
                docTopicTable[d][newk]++;
                wordTopicTable[word][newk]++;
                topicTable[newk]++;
            }
        }

        duration += (clock() - start) / (double)CLOCKS_PER_SEC;

        // Output log likelihood at each iteration
        double lik = getLogLikelihood(wordTopicTable, docTopicTable, alpha, beta);
        cout << lik << endl;
    }

    // Output top words for each topic
    vector<string> vocab;
    ifstream infile("./data/vocab.csv");
    string line;
    while (getline(infile, line))
    {
        vocab.push_back(line);
    }
    vector<vector<int>> outputTable(wordTopicTable[0].size(), vector<int>(numWords));
    for (int i = 0; i < numWords; i++) {
        for (int j = 0; j < wordTopicTable[0].size(); j++) {
            outputTable[j][i] = wordTopicTable[i][j];
        }
    }
    for (int t = 0; t < numTopics; t++) {
        vector<size_t> v = sort_indexes(outputTable[t]);
        for (int w = 0; w < 5; w++) {
            cout << vocab[v[w]] << " ";
        }
        cout << endl;
    }

    // Output duration of training
    cout << "Duration of " << numIterations << " Iterations: " << duration << endl;
}
