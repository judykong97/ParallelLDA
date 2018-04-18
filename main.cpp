#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "lda.cpp"

using namespace std;

void readCorpus(vector<vector<int>> &w, string filename) {
	ifstream infile(filename);
	string line;
	while (getline(infile, line))
	{
    	stringstream ss(line);
		vector<int> result;
		int i;
		while (ss >> i)
    	{
        	result.push_back(i);
        	if (ss.peek() == ',')
            	ss.ignore();
    	}
    	w.push_back(result);
	}
}

int main() {

	vector<vector<int>> w, z;

	string filename = "20news.csv";
	readCorpus(w, "./data/" + filename);

	int numDocs = w.size();
	int numWords = 60057; // numWords for 20news.csv
	int numTopics = 20; // numTopics
	double alpha = 0.1; // alpha
	double beta = 0.1; // beta
	int numIterations = 10; // numIterations
	// int numClocksPerIteration = 25; // numClocksPerIteration
	// int staleness = 0; // staleness

	runLDA(w, z, numDocs, numWords, numTopics, alpha, beta, numIterations);

	return 0;
}