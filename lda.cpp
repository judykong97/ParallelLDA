#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

void runLDA(vector<vector<int>> w, vector<vector<int>> &z, 
	int numDocs, int numWords, int numTopics, double alpha, double beta, int numIterations) {

	vector<vector<int>> docTopicTable(numDocs, vector<int>(numTopics, 0));
	vector<vector<int>> wordTopicTable(numWords, vector<int>(numTopics, 1));
  	vector<int> topicTable(numTopics, 0);

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

  	for (int d = 0; d < numDocs; d++) {
  		for (int j = 0; j < w[d].size(); j++) {
  			cout << z[d][j] << " ";
  		}
  		cout << endl;
  	}

  	for (int t = 0; t < numTopics; t++) {
  		cout << topicTable[t] << " ";
  	}
  	cout << endl;

  	// cout << docTopicTable.size() << endl;
  	// cout << wordTopicTable.size() << endl;
  	// cout << topicTable.size() << endl;
	cout << "LDA" << endl;
}