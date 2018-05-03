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

#define TOTAL_WORDS 1414350
#define	NUM_DOCS 18774
#define NUM_WORDS 60057
#define NUM_TOPICS 20

#include "lda.h"

using namespace std;

void readCorpus(int *w, int *w_start, string filename, int part, int numParts) {

	ifstream infile(filename);
	string line;
    // int count = 0;
    int num = 0;
    int doc = 0;
	while (getline(infile, line))
	{
        // if (count == part) {
            w_start[doc] = num;
    	    stringstream ss(line);
		    // vector<int> result;
		    int i;
		    while (ss >> i)
    	    {
        	    // result.push_back(i);
                w[num] = i;
                num++;
        	    if (ss.peek() == ',')
            	    ss.ignore();
    	    }
    	    // w.push_back(result);
            doc++;
        // }
        // count++;
        // if (count == numParts) count = 0;
	}
    w_start[doc] = num;
    cout << "Total number of words: " << num << endl;
}

int main(int argc, char *argv[]) {

    int process_count = 1;
    int process_id = 0;

#if MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
#endif
    
    bool mpi_master = process_id == 0;
    if (mpi_master) {
        cout << "Number of Processes: " << process_count << endl;
    }

	// vector<vector<int>> w, z;
    int* w = (int*)calloc(TOTAL_WORDS, sizeof(int));
    // int z[TOTAL_WORDS];
    int* w_start = (int*)calloc(NUM_DOCS + 1, sizeof(int));
    int totalWords = TOTAL_WORDS;
	int numDocs = NUM_DOCS;
	int numWords = NUM_WORDS;
	int numTopics = NUM_TOPICS;
	double alpha = 0.1; // alpha
	double beta = 0.1; // beta
	int numIterations = atoi(argv[1]); // = 1000; // numIterations
	// int numClocksPerIteration = 25; // numClocksPerIteration
	int staleness = atoi(argv[2]); // = 2; // staleness

    if (mpi_master) {
        cout << "Number of Iterations: " << numIterations << endl;
        cout << "Staleness: " << staleness << endl;
	    string filename = "20news.csv";
	    readCorpus(w, w_start, "./data/" + filename, process_id, process_count);
    }
    
#if MPI
    MPI_Bcast(w, TOTAL_WORDS, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(w_start, NUM_DOCS + 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

	runLDA(w, w_start, totalWords, numDocs, numWords, numTopics, alpha, beta, numIterations, staleness, process_id, process_count);

#if MPI
    MPI_Finalize();
#endif

    free(w);
    free(w_start);

	return 0;
}
