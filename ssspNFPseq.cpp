#include <iostream>
#include <vector>
#include <limits>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <string>
#include <chrono>
#include <queue>
#include <map>
#include <unordered_map>

struct Edge {
    int src;
    int dest;
    int weight;
};

class Graph {
public:
    // graph struct represented as a vector of vertices, and each vertex has a vector of associated edges
    int V;
    int E = 0;
    std::vector<std::vector<Edge>> adjList;

    Graph(int V) {
        this->V = V;
        adjList.resize(V);
    }

    void addEdge(int src, int dest, int weight) {
        Edge edge = {src, dest, weight};
        adjList[src].push_back(edge);
        E++;
    }

    // determines the "ideal" delta for the given graph, according to the paper
    int getDelta() {
        float avgWeight = 0;
        float avgDegree = 0;
        for (int i = 0; i < V; i++) {
            for (const Edge& edge : adjList[i]) {
                avgWeight += edge.weight;
                avgDegree += 1;
            }
        }
        avgWeight = avgWeight / E;
        avgDegree = avgDegree / V;
        std::cout << "Average Weight: " << avgWeight << ":\n";
        std::cout << "Average Degree: " << avgDegree << ":\n";
        return std::ceil(20 * (avgWeight / avgDegree));
    }
};

class Djikstra { 
public:
    static std::vector<int> dijkstra(Graph graph, int src) {
        std::vector<std::vector<Edge>> adjList = graph.adjList;

        // init distance vector to be all infinity, then set src distance to 0
        std::vector<int> dist(graph.V, std::numeric_limits<int>::max());
        dist[src] = 0;


        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;
        pq.push({0, src});

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();

            for (const Edge& edge : adjList[u]) {
                int v = edge.dest;
                int weight = edge.weight;

                if (dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    pq.push({dist[v], v});
                }
            }
        }
        // std::cout << "Shortest paths from vertex " << src << ":\n";
        // for (int i = 0; i < graph.V; ++i) {
        //     std::cout << "Vertex " << i << ": " << dist[i] << "\n";
        //}
        return dist;
    }
};

class NearFarPiling {
public:
    static std::vector<int> nearFarPiling(Graph graph, int src) {
        // delta determines cutoff for near and far piles
        int delta = graph.getDelta();
        std::cout << "Delta: " << delta << "\n";
        int iterationNumber = 1;
        std::vector<std::vector<Edge>> adjList = graph.adjList;

        // initialize dist vector s.t. src vertex dist is 0 and all other dist are max_int
        std::vector<int> dist(graph.V, std::numeric_limits<int>::max());
        dist[src] = 0;

        // begin with just src vertex in near pile
        std::vector<int> nearPile;
        nearPile.reserve(graph.V);
        std::unordered_map<int,bool> nearMap;

        std::vector<int> farPile;
        nearPile.reserve(graph.V);
        std::unordered_map<int,bool> farMap;

        nearPile.push_back(src);

        // Process vertices in the near and far piles
        std::vector<Edge> edgeList;
        while (!nearPile.empty() || !farPile.empty()) {
            while (!nearPile.empty()) {
                // Extract all edges from nearPile
                for (int x : nearPile) {
                    edgeList.insert(edgeList.end(), adjList[x].begin(), adjList[x].end());
                }
                nearPile.clear();
                nearPile.reserve(graph.V);
                farMap.clear();

                // Iterate through edges of the current vertex
                for (const Edge& edge : edgeList) {
                    int u = edge.src;
                    int v = edge.dest;
                    int weight = edge.weight;

                    // Relax the edge and add to either near or far pile
                    if (dist[u] != std::numeric_limits<int>::max() && dist[u] + weight < dist[v]) {
                        dist[v] = dist[u] + weight;
                        if (nearMap.count(v) == 0 && dist[v] <= delta * iterationNumber) {
                            nearPile.push_back(v);
                            nearMap.insert({v,true});
                        }
                        else if (farMap.count(v) == 0 && dist[v] > delta * iterationNumber) {
                            farPile.push_back(v);
                            farMap.insert({v,true});
                        }
                    }
                }
                edgeList.clear();
                edgeList.reserve(graph.V);
            }
            // Once nearPile has been emptied, update cutoff, and remove vertices in the far pile that have been processed already in the near pile
            iterationNumber++;
            if (!farPile.empty()) {
                // add to near pile if its distance is less than new cut off, can be parallelized
                for(int vertex : farPile) {
                    if (dist[vertex] <= delta * iterationNumber) {
                        nearPile.push_back(vertex);
                        nearMap.insert({vertex,true});
                    }
                }
                // remove items already processed in the near pile in previous round num and items just added to the near pile
                std::vector<int>::iterator new_end = std::remove_if(farPile.begin(), farPile.end(), 
                                                                [&](int vertex){ return dist[vertex] < delta * (iterationNumber - 1); });
                
                farPile.erase(new_end, farPile.end());
                for (std::vector<int>::iterator new_end1 = new_end; new_end1 < farPile.end(); new_end1++) {
                    farMap.erase(*new_end1);
                }
            }
            // std::cout << "Iteration Number: " << iterationNumber << "\n";
        }

        // Print the shortest paths
        // std::cout << "Shortest paths from vertex " << src << ":\n";
        // for (int i = 0; i < graph.V; ++i) {
        //     std::cout << "Vertex " << i << ": " << dist[i] << "\n";
        // }
        return dist;
    }
};

using namespace std;

int main() {

    int V;
    int fileNumber;
    std::string fileName = "10kGraph.csv";

    std::map<std::string, int> fileMap;
    fileMap["1mGraph.csv"] = 1;
    fileMap["10kGraph.csv"] = 2;
    fileMap["1kGraph.csv"] = 3;
    fileMap["250Graph.csv"] = 4;
    fileMap["8Graph.csv"] = 5;

    fileNumber = fileMap[fileName];
    switch(fileNumber) {
        case 1:
            V = 1000000;
            break;
        case 2:
            V = 10000;
            break;
        case 3:
            V = 1000;
            break;  
        case 4:
            V = 250;
            break;
        case 5:
            V = 8;
            break;
        default:
            V = 8;
            fileName = "8Graph.csv";
        }

    auto startGraph = chrono::high_resolution_clock::now();
    Graph graph(V);
    std::string line;
    std::ifstream myfile("data/" + fileName);

    while (std::getline(myfile, line)) {
        std::stringstream s(line);
        std::string word; 

        int pos = 0;
        int args[3];
        while (std::getline(s, word, ',')) {
            std::stringstream(word) >> args[pos++];
        }
        graph.addEdge(args[0], args[1], args[2]);
    }
    
    // Make the graph
    auto endGraph = chrono::high_resolution_clock::now();
	double elapsedGraph = chrono::duration_cast<chrono::nanoseconds>(endGraph - startGraph).count();
    std::cout << "Time to build graph: " << elapsedGraph << "\n";

    // Djikstra's
    auto startDjikstra = chrono::high_resolution_clock::now();
    std::vector<int> distancesDjikstra = Djikstra::dijkstra(graph, 0);
    auto endDjikstra = chrono::high_resolution_clock::now();
	double elapsedDjikstra = chrono::duration_cast<chrono::nanoseconds>(endDjikstra- startDjikstra).count();
    std::cout << "Time to run Djikstra: " << elapsedDjikstra << "\n";

    // Near-Far Piling
    auto startNFP = chrono::high_resolution_clock::now();
    std::vector<int> distancesNFP = NearFarPiling::nearFarPiling(graph, 0);
    auto endNFP = chrono::high_resolution_clock::now();
	double elapsedNFP = chrono::duration_cast<chrono::nanoseconds>(endNFP - startNFP).count();
    std::cout << "Time to run NFP: " << elapsedNFP << "\n";

    // Error checking NFP against Djikstra's
    int count = 0;
    for (int l = 0; l < distancesDjikstra.size(); l++) {
        if (distancesDjikstra[l] != distancesNFP[l]) {
            // std::cout << l << "\n";
            count++;
        }
    }
    std::cout << "Number of mismatches: " << count << "\n";
}
