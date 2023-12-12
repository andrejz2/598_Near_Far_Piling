#include <iostream>
#include <vector>
#include <limits>
#include <list>
#include <algorithm>

struct Edge {
    int dest;
    int weight;
};

class Graph {
public:
    int V;
    std::vector<std::vector<Edge>> adjList;

    int numProc = 1;

    Graph(int V) {
        this->V = V;
        adjList.resize(V);
    }

    void addEdge(int src, int dest, int weight) {
        Edge edge = {dest, weight};
        adjList[src].push_back(edge);
    }

    int getDelta() {
        float avgWeight, avgDegree = 0;
        for (int i = 0; i < V; i++) {
            for (const Edge& edge : adjList[i]) {
                avgWeight += edge.weight;
                avgDegree += 1;
            }
        }
        std::cout << "Average Weight: " << avgWeight / V << ":\n";
        std::cout << "Average Degree: " << avgDegree / V << ":\n";
        return numProc * (avgWeight / avgDegree);
    }

    int delta;
    int roundNum;

    bool IsLessThanCutOff(int i) { 
        return (i < (delta * roundNum)); 
    }

    void nearFarPiling(int src, int numProc) {
        delta = getDelta();
        roundNum = 1;

        std::vector<int> dist(V, std::numeric_limits<int>::max());
        dist[src] = 0;

        std::vector<int> nearPile;
        std::vector<int> farPile;
        nearPile.push_back(src);

        // Process vertices in the queue
        while (!nearPile.empty() || !farPile.empty()) {
            #pragma omp parallel
            {
                while (!nearPile.empty()) {
                    int u;
                    #pragma omp critical
                    {
                        u = nearPile.back();
                        nearPile.pop_back();
                    }

                    // Iterate through edges of the current vertex
                    #pragma omp for
                    for (int i = 0; i < adjList[u].size(); ++i) {
                        const Edge& edge = adjList[u][i];
                        int v = edge.dest;
                        int weight = edge.weight;

                        // Relax the edge and add to either near or far pile
                        #pragma omp critical
                        {
                            if (dist[u] != std::numeric_limits<int>::max() && dist[u] + weight < dist[v]) {
                                dist[v] = dist[u] + weight;
                                // delta*roundNum is the criteria for adding to near or far pile
                                if (dist[v] <= delta*roundNum) {
                                    nearPile.push_back(v);
                                }
                                else {
                                    farPile.push_back(v);
                                }
                            }
                        }
                    }
                // remove duplicates
                std::sort(nearPile.begin(), nearPile.end());
                nearPile.erase(std::unique(nearPile.begin(), nearPile.end()), nearPile.end());
                }
            }
            // Once the nearQueue has been emptied, update cutoff, and remove vertices in the far pile that have been processed already in the near pile
            roundNum++;
            if (!farPile.empty()) {
                std::cout << "if far pile not empty " << (delta*roundNum) << ":\n";
                std::sort(farPile.begin(), farPile.end());
                farPile.erase(std::unique(farPile.begin(), farPile.end()), farPile.end());
                // add to near pile if its distance is less than new cut off
                for(int vertex : farPile) {
                    if (dist[vertex] <= delta*roundNum) {
                        nearPile.push_back(vertex);
                    }
                // remove items already processed in the near pile in previous round num and items just added to the near pile
                farPile.erase(std::remove_if(farPile.begin(), farPile.end(), 
                    [&](int vertex){ return dist[vertex] < delta * roundNum; }), farPile.end());
                }
            }
        }

        // Print the shortest paths
        std::cout << "Shortest paths from vertex " << src << ":\n";
        for (int i = 0; i < V; ++i) {
            std::cout << "Vertex " << i << ": " << dist[i] << "\n";
        }
    }
};

int main() {
    // Create a graph with 5 vertices
    std::cout << "gabba gool" << ":\n";
    Graph graph(5);
    std::cout << "gabba gool 2" << ":\n";

    // Add edges to the graph
    graph.addEdge(0, 1, 2);
    graph.addEdge(0, 3, 1);
    graph.addEdge(1, 2, 3);
    graph.addEdge(1, 4, 7);
    graph.addEdge(2, 3, 4);
    graph.addEdge(3, 4, 5);

    // Find and print the single-source shortest paths from node 0 using Bellman-Ford with a queue
    graph.nearFarPiling(0, 1);

    return 0;
}
