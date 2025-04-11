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
vector<vector<double>> FuzzyLogic(vector<vector<double>>& data, int numOfClusters, int maxIterations);

int main(int argc, char* argv[]) {
	if (argc != 2) {
		cout << "INVALID ENTRY: PASS FILE AS ARGUMENT" << endl;
		return 1;
	}
	string fileName = argv[1];
	vector<vector<double>> data = getData(fileName);
    //work on the clustering methods
    //1. K Means Clustering
    cout << endl << "K Means Clustering" << endl;
    vector<vector<double>> cluster = KMeansClustering(data,5,100);
    for (int i = 0; i < cluster.size(); i++) {
        cout << "Cluster " << i << ": ";
        for (int j = 0; j < cluster[0].size(); j++) {
            cout << cluster[i][j] << " ";
        }
        cout << endl;
    }
    //2. Fuzzy Logic
    cout << endl << endl << "Fuzzy Logic" << endl;
    cluster = FuzzyLogic(data, 5, 100);
    for (int i = 0; i < cluster.size(); i++) {
        cout << "Cluster " << i << ": ";
        for (int j = 0; j < cluster[0].size(); j++) {
            cout << cluster[i][j] << " ";
        }
        cout << endl;
    }
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
            if (squaredDistance(newClusters[h], clusters[h]) > .001) {
                converged = false;
                break;
            }
        }
        if (converged) {
            cout << "Clusters converged on iteration: " << i << endl;
            
            break;
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

