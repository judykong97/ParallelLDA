#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

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
            cout << i << " ";
    	}
    	w.push_back(result);
    	cout << endl;
	}
}

int main() {

	vector<vector<int>> w;

	string filename = "20news.csv";
	readCorpus(w, "./data/" + filename);

	cout << endl << endl << w.size() << endl;
	return 0;
}