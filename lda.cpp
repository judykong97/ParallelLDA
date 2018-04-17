#include <iostream>
#include <string>
#include <vector>

using namespace std;

void runLDA(vector<vector<int>> w, vector<vector<int>> &z, 
	int numDocs, int numWords, int numTopics, double alpha, double beta, int numIterations) {

	vector<vector<int>> docTopicTable(numDocs, vector<int>(numTopics, 0));
	vector<vector<int>> wordTopicTable(numWords, vector<int>(numTopics, 1));
  	vector<int> topicTable(numTopics, 0);

  	

  	cout << docTopicTable.size() << endl;
  	cout << wordTopicTable.size() << endl;
  	cout << topicTable.size() << endl;
	cout << "LDA" << endl;
}