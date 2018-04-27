#ifndef MPI
#define MPI 0
#endif

#if MPI
#include <mpi.h>
#endif

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "lda.cpp"

using namespace std;

void readCorpus(vector<vector<int>> &w, string filename, int part, int numParts) {
	ifstream infile(filename);
	string line;
    int count = 0;
	while (getline(infile, line))
	{
        if (count == part) {
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
        count++;
        if (count == numParts) count = 0;
	}
}

int main(int argc, char *argv[]) {

    int process_count = 1;
    int process_id = 0;

#if MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
#endif

	vector<vector<int>> w, z;

	string filename = "20news.csv";
	readCorpus(w, "./data/" + filename, process_id, process_count);

	int numDocs = w.size();
	int numWords = 60057; // numWords for 20news.csv
	int numTopics = 20; // numTopics
	double alpha = 0.1; // alpha
	double beta = 0.1; // beta
	int numIterations = 10; // numIterations
	// int numClocksPerIteration = 25; // numClocksPerIteration
	// int staleness = 0; // staleness
    
    cout << "numDocs = " << numDocs << endl;

    bool mpi_master = process_id == 0;
    if (mpi_master) {
        cout << "process count: " << process_count << endl;
    }

	runLDA(w, z, numDocs, numWords, numTopics, alpha, beta, numIterations);

#if MPI
    MPI_Finalize();
#endif

	return 0;
}
