#ifndef MPI
#define MPI 0
#endif

#if MPI
#include <mpi.h>
#endif

void runLDA(int *w, int *w_start, 
        int totalWords, int numDocs, int numWords, int numTopics, double alpha, double beta, int numIterations, int staleness, int process_id, int process_count);

void runLDASync(int *w, int *w_start, 
        int totalWords, int numDocs, int numWords, int numTopics, double alpha, double beta, int numIterations, int staleness, int process_id, int process_count);

void runLDAAsync(int *w, int *w_start, 
        int totalWords, int numDocs, int numWords, int numTopics, double alpha, double beta, int numIterations, int staleness, int process_id, int process_count);

void runLDAAsyncBlocking(int *w, int *w_start, 
        int totalWords, int numDocs, int numWords, int numTopics, double alpha, double beta, int numIterations, int staleness, int process_id, int process_count);
