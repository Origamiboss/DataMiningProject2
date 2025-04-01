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

vector<vector<double>> getData(string fileName);
vector<vector<double>> KMeansClustering(vector<vector<double>>& data, int numOfClusters, int maxNumOfIterations);
double squaredDistance(vector<double>& x, vector<double>& y);

int main(int argc, char* argv[]) {
	if (argc != 2) {
		cout << "INVALID ENTRY: PASS FILE AS ARGUMENT" << endl;
		return 1;
	}
	string fileName = argv[1];
	vector<vector<double>> data = getData(fileName);
    //work on the clustering methods
    //1. K Means Clustering
    vector<vector<double>> cluster = KMeansClustering(data,2,100);
    //2. Fuzzy Logic
}

vector<vector<double>> getData(string fileName) {
    ifstream file(fileName);
    string line;

    // Read the first and throw away the header
    getline(file, line);

    // Create a 2D array with numOfInstances rows and sizeOfInstances columns
    vector<vector<double>> data;
    int dataIndex = 0;

    // Go through the rest of the file and store the data that are doubles, throw everything else out
    while (getline(file, line)) {
        stringstream point(line);
        string newDouble;
        int index = 0;
        data.push_back(vector<double>());
        while (getline(point, newDouble, ',')) {
            try {
                double newValue = stod(newDouble);
                if (isnan(newValue))
                    throw runtime_error("NaN Value Detected");
                data[dataIndex].push_back(newValue);
                index++;
            }
            catch (exception e) {
                //do nothing, throw out bad numbers
            }
            
        }
        dataIndex++;
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
        //check if they changed
        bool converged = true;
        for (int h = 0; h < numOfClusters; h++) {
            if (squaredDistance(newClusters[h], clusters[h]) > 1e-6) {
                converged = false;
                break;
            }
        }
        if (converged) {
            cout << "Clusters converged on iteration: " << i << endl;
            
            break;
        }
        //replace the old clusters
        for (int i = 0; i < clusters.size(); i++) {
            cout << "Cluster " << i << ": ";
            for (int j = 0; j < clusters[i].size(); j++) {
                cout << clusters[i][j] << " ";
            }
            cout << endl;
            //print the counts
            cout << "Count: " << count[i] << endl;
        }

        
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

