// DataMiningProject2.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include<iostream>
#include<string>
#include<vector>
#include <fstream>
#include <sstream>
#include <random>
#include<set>
#include<numeric>
using namespace std;

vector<vector<double>> getData(string fileName, vector<string>& geneName);
vector<vector<double>> KMeansClustering(vector<vector<double>>& data, int numOfClusters, int maxNumOfIterations);
double squaredDistance(vector<double>& x, vector<double>& y);
vector<vector<double>> FuzzyLogic(vector<vector<double>>& data, int numOfClusters, int maxIterations);
void analyzeGenes(string geneFile, vector<vector<double>>& data, vector<vector<double>>& clusters, vector<string> geneNames);

int main(int argc, char* argv[]) {
	if (argc != 3) {
		cout << "INVALID ENTRY: PASS FILE AS ARGUMENT" << endl;
        cout << "./DataMiningProject2.exe dataFile importantGenes" << endl;
		return 1;
	}
	string fileName = argv[1];
    string geneFile = argv[2];
    vector<string> geneNames;
	vector<vector<double>> data = getData(fileName, geneNames);
    //Work on the clustering methods
    //1. K Means Clustering
    cout << endl << "K Means Clustering" << endl;
    vector<vector<double>> cluster = KMeansClustering(data,10,100);
    for (int i = 0; i < cluster.size(); i++) {
        cout << "Cluster " << i << ": ";
        for (int j = 0; j < cluster[0].size(); j++) {
            cout << cluster[i][j] << " ";
        }
        cout << endl;
    }
    //analyze cluster information
    analyzeGenes(geneFile, data, cluster, geneNames);

    //2. Fuzzy Logic
    cout << endl << endl << "Fuzzy Logic" << endl;
    cluster = FuzzyLogic(data, 10, 100);
    for (int i = 0; i < cluster.size(); i++) {
        cout << "Cluster " << i << ": ";
        for (int j = 0; j < cluster[0].size(); j++) {
            cout << cluster[i][j] << " ";
        }
        cout << endl;
    }
    //analyze cluster information
    analyzeGenes(geneFile, data, cluster, geneNames);
}

vector<vector<double>> getData(string fileName, vector<string>& geneName) {
    ifstream file(fileName);
    string line;

    // Skip header
    getline(file, line);

    vector<vector<double>> data;
    int dataIndex = 0;

    while (getline(file, line)) {
        stringstream point(line);
        string newDouble;
        data.emplace_back();  // Add new row
        bool gotName = false;
        while (getline(point, newDouble, ',')) {
            if (!gotName) {
                geneName.push_back(newDouble);
                gotName = true;
            }
            try {
                double newValue = stod(newDouble);
                if (isnan(newValue))
                    throw runtime_error("NaN Value Detected");
                data[dataIndex].push_back(newValue);
            }
            catch (const exception& e) {
                // skip bad values
            }
        }
        dataIndex++;
    }

    if (data.empty()) return data; // avoid divide by zero

    size_t cols = data[0].size();
    vector<double> Min(cols, INFINITY);
    vector<double> Max(cols, -INFINITY);

    for (const auto& row : data) {
        for (size_t j = 0; j < cols; ++j) {
            Min[j] = min(Min[j], row[j]);
            Max[j] = max(Max[j], row[j]);
        }
    }

    for (auto& row : data) {
        for (size_t j = 0; j < cols; ++j) {
            if (Max[j] == Min[j]) {
                row[j] = 0;
            }
            else {
                row[j] = (row[j] - Min[j]) / (Max[j] - Min[j]);
            }
        }
    }

    file.close();
    return data;
}

vector<vector<double>> KMeansClustering(vector<vector<double>>& data, int numOfClusters, int maxNumOfIterations) {
    //start by selecting new clusters
    int numOfInstances = data.size();
    int sizeOfInstance = data[0].size();
    vector<vector<double>> clusters(numOfClusters, vector<double>(sizeOfInstance, 0.0));
    
    //Select the Initial Clusters Randomly
    // Use random_device for a better random seed
    random_device rd;
    mt19937 gen(rd());  // Initialize random number generator with seed
    uniform_int_distribution<int> dist(0, numOfInstances - 1);

    // Make a set for the data positions
    set<vector<double>> clusterData;

    for (int i = 0; i < numOfClusters; i++) {
        int randomNum;
        do {
            randomNum = dist(gen);
        } while (clusterData.find(data[randomNum]) != clusterData.end());  // Ensure no duplicates
        clusterData.insert(data[randomNum]);

        // Create a new cluster and copy the data from the chosen instance
        clusters[i] = data[randomNum];
    }
    
    //Now run it for a certain number of iterations or it stops
    for (int i = 0; i < maxNumOfIterations; i++) {
        vector<vector<double>> newClusters(numOfClusters, vector<double>(sizeOfInstance, 0.0));
        vector<int> count(numOfClusters, 0);
        //find each data points closest cluster and add it to the new clusters to find the mean of the cluster points
        for (vector<double> point : data) {
            int closestCluster = 0;
            double closestDist = squaredDistance(point, clusters[closestCluster]);
            for (int h = 1; h < numOfClusters; h++) {
                double dist = squaredDistance(point, clusters[h]);
                if (dist < closestDist) {
                    closestCluster = h;
                    closestDist = dist;
                }
            }

            // Add the data point to the closest cluster, and check for NaN
            for (int h = 0; h < sizeOfInstance; h++) {
                newClusters[closestCluster][h] += point[h];
                if (isnan(newClusters[closestCluster][h])) {
                    cout << "NaN detected in newClusters after addition at cluster " << closestCluster << " position [" << h << "]" << endl;
                }
            }

            count[closestCluster]++;
        }
        //now take the mean of each cluster point
        for (int h = 0; h < numOfClusters; h++) {
            
            for (int j = 0; j < sizeOfInstance; j++) {
                newClusters[h][j] /= count[h];
            }
        }
        
        //replace the old clusters
        /*for (int i = 0; i < clusters.size(); i++) {
            cout << "Cluster " << i << ": ";
            for (int j = 0; j < clusters[i].size(); j++) {
                cout << clusters[i][j] << " ";
            }
            cout << endl;
            //print the counts
            cout << "Count: " << count[i] << endl;
        }*/

        
        clusters = newClusters;
    }

    return clusters;
}
double squaredDistance(vector<double>& x, vector<double>& y) {
    double dist = 0;
    int sizeOfInstance = min(x.size(), y.size());
    for (int i = 0; i < sizeOfInstance; i++) {
        double diff = x[i] - y[i];
        dist += diff * diff;
    }
    if (isnan(dist)) {
        cout << "NaN detected in squaredDistance " << endl;
        for (int i = 0; i < sizeOfInstance; i++) {
            cout << "X: " << x[i] << endl;
            cout << "Y: " << y[i] << endl;
        }
        return NULL;
    }
    return dist;
}
vector<vector<double>> FuzzyLogic(vector<vector<double>>& data, int numOfClusters, int maxIterations) {
    int numOfInstances = data.size();
    int sizeOfInstance = data[0].size();
    double m = numOfClusters;
    // Initialize fuzzy membership matrix
    vector<vector<double>> membership(numOfInstances, vector<double>(numOfClusters));
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < numOfInstances; i++) {
        double sum = 0.0;
        for (int j = 0; j < numOfClusters; j++) {
            membership[i][j] = dis(gen);
            sum += membership[i][j];
        }
        for (int j = 0; j < numOfClusters; j++) {
            membership[i][j] /= sum;
        }
    }

    vector<vector<double>> centers(numOfClusters, vector<double>(sizeOfInstance, 0.0));

    for (int iter = 0; iter < maxIterations; iter++) {
        // Update cluster centers
        for (int j = 0; j < numOfClusters; j++) {
            vector<double> numerator(sizeOfInstance, 0.0);
            double denominator = 0.0;

            for (int i = 0; i < numOfInstances; i++) {
                double u_ij_m = pow(membership[i][j], m);
                for (int k = 0; k < sizeOfInstance; k++) {
                    numerator[k] += u_ij_m * data[i][k];
                }
                denominator += u_ij_m;
            }

            for (int k = 0; k < sizeOfInstance; k++) {
                centers[j][k] = numerator[k] / denominator;
            }
        }

        // Update membership matrix
        for (int i = 0; i < numOfInstances; i++) {
            for (int j = 0; j < numOfClusters; j++) {
                double sum = 0.0;
                double dist_ij = 0.0;
                for (int k = 0; k < sizeOfInstance; k++) {
                    dist_ij += pow(data[i][k] - centers[j][k], 2);
                }
                dist_ij = sqrt(dist_ij) + 1e-6; // Avoid division by zero

                for (int c = 0; c < numOfClusters; c++) {
                    double dist_ic = 0.0;
                    for (int k = 0; k < sizeOfInstance; k++) {
                        dist_ic += pow(data[i][k] - centers[c][k], 2);
                    }
                    dist_ic = sqrt(dist_ic) + 1e-6;

                    sum += pow(dist_ij / dist_ic, 2.0 / (m - 1));
                }

                membership[i][j] = 1.0 / sum;
            }
        }
    }

    return centers;
}

void analyzeGenes(string geneFile, vector<vector<double>>& data, vector<vector<double>>& clusters, vector<string> geneNames) {
    vector<vector<string>> geneClusters(clusters.size(), vector<string>());

    //locate which geneNames go into what
    for (int i = 0; i < data.size(); i++) {
        double minDist = INFINITY;
        int clusterId = 0;
        for (int j = 0; j < clusters.size(); j++) {
            double dist = squaredDistance(data[i], clusters[j]);
            if (dist < minDist) {
                clusterId = j;
                minDist = dist;
            }
        }
        //save the file
        geneClusters[clusterId].push_back(geneNames[i]);
    }
    //make a set that contains the genes we are looking for
    set<string> importantGenes;
    fstream file(geneFile);
    string line;
    while (getline(file, line)) {
        importantGenes.insert(line);
    }

    //Display the names inside of each cluster
    for (int i = 0; i < geneClusters.size(); i++) {
        cout << "Cluster: " << i << endl;
        for (int j = 0; j < geneClusters[i].size(); j++) {
            if(importantGenes.find(geneClusters[i][j]) != importantGenes.end())
                cout << geneClusters[i][j] << " ";
        }
        cout << endl;
    }
}
