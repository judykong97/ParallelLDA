#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>

using namespace std;

inline double logDirichlet_const(double alpha, int k) {
    return k * lgamma(alpha) - lgamma(k * alpha);
}

double logDirichlet_vector(vector<double> alpha) {
    double sumLogGamma = 0.0;
    double logSumGamma = 0.0;
    for(int i=0; i<alpha.size(); i++) {
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
    for(int k=0; k<numTopics; k++) {
        for(int w=0; w<numWords; w++) {
            temp[w] = beta + wordTopicTable[w][k];
        }
        lik += logDirichlet_vector(temp);
        lik -= logDirichlet_const(beta, numWords);
    }
    vector<double> temp2(numTopics, 0.0);
    for(int d=0; d<numDocs; d++) {
        for(int k=0; k<numTopics; k++) {
            temp2[k] = alpha + docTopicTable[d][k];
        }
        lik += logDirichlet_vector(temp2);
        lik -= logDirichlet_const(alpha, numTopics);
    }
    return lik;
}

void runLDA(vector<vector<int>> w, vector<vector<int>> &z, 
	int numDocs, int numWords, int numTopics, double alpha, double beta, int numIterations) {

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

  	// for (int d = 0; d < numDocs; d++) {
  	// 	for (int t = 0; t < numTopics; t++) {
  	// 		cout << docTopicTable[d][t] << " ";
  	// 	}
  	// 	cout << endl;
  	// }

  	// for (int d = 0; d < numWords; d++) {
  	// 	for (int t = 0; t < numTopics; t++) {
  	// 		cout << wordTopicTable[d][t] << " ";
  	// 	}
  	// 	cout << endl;
  	// }



  	// for (int t = 0; t < numTopics; t++) {
  	// 	cout << topicTable[t] << " ";
  	// }
  	// cout << endl;

  	// cout << docTopicTable.size() << endl;
  	// cout << wordTopicTable.size() << endl;
  	// cout << topicTable.size() << endl;
	// cout << "LDA" << endl;

    for(int i=0; i<numIterations; i++) {
        for(int d=0; d<numDocs; d++) {
            for(int j=0; j<w[d].size(); j++) {
                int word = w[d][j];
                int topic = z[d][j];
                docTopicTable[d][topic]--;
                wordTopicTable[word][topic]--;
                topicTable[topic]--;

                double norm = 0.0;
                int newk = topic;
                for(int k = 0; k<numTopics; k++) {
                    int z_dj_equals_k = (k == topic);
                    double ak = docTopicTable[d][k] - z_dj_equals_k + alpha;
                    double bk = (wordTopicTable[word][k] - z_dj_equals_k + beta) / (topicTable[k] - z_dj_equals_k + numWords + beta);
                    p[k] = ak * bk;
                    norm += p[k];
                }

                double sum_p_up_to_k = 0.0;
                double r = ((double) rand()) / RAND_MAX;
                for(int k = 0; k<numTopics; k++) {
                    sum_p_up_to_k += p[k] / norm;
                    if(r < sum_p_up_to_k) {
                        newk = k;
                        break;
                    }
                }
                z[d][j] = newk;
                docTopicTable[d][newk]++;
                wordTopicTable[word][newk]++;
                topicTable[newk]++;
            }
        }
        double lik = getLogLikelihood(wordTopicTable, docTopicTable, alpha, beta);
        cout << lik << endl;
    }

    // for (int d = 0; d < numDocs; d++) {
    //  for (int j = 0; j < w[d].size(); j++) {
    //      cout << z[d][j] << " ";
    //  }
    //  cout << endl;
    // }
}