#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <ctime>
#include <string.h>
#include "lda_utils.h"
#include "lda.h"

using namespace std;

void runLDAAsyncBlocking(int *w, int *w_start, 
        int totalWords, int numDocs, int numWords, int numTopics, double alpha, double beta, int numIterations, int staleness, int process_id, int process_count) {

    clock_t start;
    double duration;

    int* z = (int*) calloc(totalWords, sizeof(int));
    int* docTopicTable = (int*) calloc(numDocs * numTopics, sizeof(int));
    int* wordTopicTable = (int*) calloc(numWords * numTopics, sizeof(int));    
    int* topicTable = (int*) calloc(numTopics, sizeof(int));
    double* p = (double*) calloc(numTopics, sizeof(double));

    int* updateW = (int*) calloc(numWords * numTopics, sizeof(int));
    int* updateT = (int*) calloc(numTopics, sizeof(int));
    int* globalW;
    int* globalT;

    int iter = 0;

#if MPI
    globalW = (int*) calloc(numWords * numTopics, sizeof(int));
    globalT = (int*) calloc(numTopics, sizeof(int));

    int* num_rec = (int*) calloc(process_count, sizeof(int));
#endif

    int mpi_master = process_id == 0;

    for (int d = process_id; d < numDocs; d += process_count) {
        int doffset = d * numTopics;
        for (int j = w_start[d]; j < w_start[d + 1]; j++) {
            int word = w[j];
            int topic = rand() % numTopics;
            z[j] = topic;
            docTopicTable[doffset + topic]++;
            wordTopicTable[word * numTopics + topic]++;
            topicTable[topic]++;
        }
    }

    int curr_updated = 0;

    for (int i = 0; i < numIterations; i++) {

        start = clock();

        if (staleness == 1 || (staleness > 1 && iter % staleness == 0)) {
            memset(updateW, 0, sizeof(int) * numWords * numTopics);
            memset(updateT, 0, sizeof(int) * numTopics);
        }

        for (int d = process_id; d < numDocs; d += process_count) {
            int doffset = d * numTopics;
            for (int j = w_start[d]; j < w_start[d + 1]; j++) {
                int word = w[j];
                int topic = z[j];
                int woffset = word * numTopics;
                docTopicTable[doffset + topic]--;
                updateW[woffset + topic]--; 
                updateT[topic]--; 

                double norm = 0.0;
                int newk = topic;
                for (int k = 0; k < numTopics; k++) {
                    int z_dj_equals_k = (k == topic);
                    double ak = docTopicTable[doffset + k] - z_dj_equals_k + alpha;
                    double bk = (wordTopicTable[woffset + k] + updateW[woffset + k] - z_dj_equals_k + beta) / (topicTable[k] + updateT[k] - z_dj_equals_k + numWords + beta);
                    norm += ak * bk;
                    p[k] = norm;
                }
                
                double sum_p_up_to_k = 0.0;
                double r = ((double) rand()) / RAND_MAX;
                for(int k = 0; k < numTopics; k++) {
                    sum_p_up_to_k += p[k] / norm;
                    if(r < sum_p_up_to_k) {
                        newk = k;
                        break;
                    }
                }
                z[j] = newk;
                docTopicTable[doffset + newk]++;
                updateW[woffset + newk]++; 
                updateT[newk]++;
            }
        }

#if MPI
        iter++;
        if (iter % staleness > 0) {
            duration += (clock() - start) / (double)CLOCKS_PER_SEC;
            continue;
        }

        MPI_Request wreqs[process_count - 1];
        MPI_Request treqs[process_count - 1];
        MPI_Request reqs[2 * process_count - 2];
        MPI_Status stats[2 * process_count - 2];
        if (!mpi_master) {
            MPI_Isend(updateW, numWords * numTopics, MPI_INT, 0, 1, MPI_COMM_WORLD, &wreqs[process_id - 1]);
            MPI_Isend(updateT, numTopics, MPI_INT, 0, 2, MPI_COMM_WORLD, &treqs[process_id - 1]);
            MPI_Irecv(wordTopicTable, numWords * numTopics, MPI_INT, 0, 1, MPI_COMM_WORLD, &reqs[2 * process_id - 2]);
            MPI_Irecv(topicTable, numTopics, MPI_INT, 0, 2, MPI_COMM_WORLD, &reqs[2 * process_id - 1]);
            MPI_Waitall(2, reqs + 2 * process_id - 2, stats + 2 * process_id - 2);
        } else {

            for (int j = 0; j < numWords * numTopics; j++) {
                wordTopicTable[i] += updateW[j];
            }
            for (int j = 0; j < numTopics; j++) {
                topicTable[j] += updateT[j];
            }

            memset(num_rec, 0, process_count * sizeof(int));
            MPI_Request update_reqs[2 * process_count - 2];
            MPI_Status stat;
            int source_process, message_length;

            for (int i = 0; i < 2 * process_count - 2; i++) {
                MPI_Status status;
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                source_process = status.MPI_SOURCE;
                MPI_Get_count(&status, MPI_INT, &message_length);
                if (message_length > numTopics) { // Received word topic table
                    MPI_Recv(updateW, numWords * numTopics, MPI_INT, source_process, 1, MPI_COMM_WORLD, &status);
                    for (int j = 0; j < numWords * numTopics; j++) {
                        wordTopicTable[j] += updateW[j];
                    }
                    MPI_Isend(wordTopicTable, numWords * numTopics, MPI_INT, source_process, 1, MPI_COMM_WORLD, &update_reqs[source_process * 2 - 2]);
                } else { // Received topic table
                    MPI_Recv(updateT, numTopics, MPI_INT, source_process, 2, MPI_COMM_WORLD, &status);
                    for (int j = 0; j < numTopics; j++) {
                        topicTable[j] += updateT[j];
                    }
                    MPI_Isend(topicTable, numTopics, MPI_INT, source_process, 2, MPI_COMM_WORLD, &update_reqs[source_process * 2 - 1]);
                }
            }
            
            /*
            int received_count = 0;
            while (received_count < 2 * process_count - 2) {
                // cout << received_count << endl;
                for (int i = 1; i < process_count; i++) {
                    if (num_rec[i] == 2) continue;
                    // if (num_rec[i] < 2) {
                        int flag = 0;
                        MPI_Status status;
                        MPI_Request request;
                        MPI_Iprobe(i, 1, MPI_COMM_WORLD, &flag, &status);
                        if (flag) {
                            received_count++;
                            num_rec[i]++;
                            MPI_Recv(updateW, numWords * numTopics, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                            for (int j = 0; j < numWords * numTopics; j++) {
                                wordTopicTable[j] += updateW[j];
                            }
                        }
                    // }
                    if (num_rec[i] < 2) {
                        int flag = 0;
                        MPI_Status status;
                        MPI_Request request;
                        MPI_Iprobe(i, 2, MPI_COMM_WORLD, &flag, &status);
                        if (flag) {
                            received_count++;
                            num_rec[i]++;
                            MPI_Recv(updateT, numTopics, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
                            for (int j = 0; j < numTopics; j++) {
                                topicTable[j] += updateT[j];
                            }
                        }
                    }
                    if (num_rec[i] == 2) {
                        MPI_Isend(wordTopicTable, numWords * numTopics, MPI_INT, i, 1, MPI_COMM_WORLD, &update_reqs[i * 2 - 2]);
                        MPI_Isend(topicTable, numTopics, MPI_INT, i, 2, MPI_COMM_WORLD, &update_reqs[i * 2 - 1]);
                    }
                    if (received_count == 2 * process_count - 2) break; 
                }
            }
            // bar_count = process_count - 1;
            */
        }
        
#else 
        globalW = updateW;
        globalT = updateT;
        for (int i = 0; i < numWords * numTopics; i++) {
            wordTopicTable[i] += globalW[i];
        }
        for (int i = 0; i < numTopics; i++) {
            topicTable[i] += globalT[i];
        }
#endif

        duration += (clock() - start) / (double)CLOCKS_PER_SEC;

        double lik = getLogLikelihood(wordTopicTable, docTopicTable, alpha, beta, numWords, numDocs, numTopics, process_id, process_count);
        double global_lik = lik;

#if MPI
        global_lik = 0;
        MPI_Reduce(&lik, &global_lik, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

        if (mpi_master) {
            cout << global_lik << endl;
        }

    }

    if (mpi_master) {
        // Output top words for each topic
        vector<string> vocab;
        ifstream infile("./data/vocab.csv");
        string line;
        while (getline(infile, line))
        {
            vocab.push_back(line);
        }
        vector<vector<int>> outputTable(numTopics, vector<int>(numWords));
        for (int i = 0; i < numWords; i++) {
            for (int j = 0; j < numTopics; j++) {
                outputTable[j][i] = wordTopicTable[i * numTopics + j];
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

    free(z);
    free(p);
    free(docTopicTable);
    free(wordTopicTable);
    free(topicTable);
    free(updateW);
    free(updateT);
#if MPI
    free(globalW);
    free(globalT);
    free(num_rec);
#endif

}
