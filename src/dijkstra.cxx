// Copyright (C) 2020 Joseph R. Quinn
// SPDX-License-Identifier: MIT

#include <limits>
#include <iostream>
#include <vector>

using namespace std;

template<typename T>
using Matrix = vector<vector<T>>;

template<typename T>
class Dijkstra {
public:
    Dijkstra(const Matrix<T> g) {
        graph = g;
        size = g.size();
    }
    vector<T> from_vertex(const int source);
private:
    int size;
    Matrix<T> graph;
    int min_distance(const vector<T> dist, const vector<bool> is_processed);
};

template<typename T>
int Dijkstra<T>::min_distance(const vector<T> dist,
    const vector<bool> is_processed)
{
    float min = std::numeric_limits<float>::infinity();
    int min_index;

    for (int v = 0; v < size; v++) {
        if (is_processed[v] == false && dist[v] <= min) {
            min = dist[v];
            min_index = v;
        }
    }

    return min_index;
}

template<typename T>
vector<T> Dijkstra<T>::from_vertex(const int source)
{
    vector<T> dist;
    vector<bool> is_processed;

    for (int i = 0; i < size; i++) {
        dist.push_back(numeric_limits<float>::infinity());
        is_processed.push_back(false);
    }

    dist[source] = 0;

    for (int count = 0; count < size - 1; count++) {
        int u = min_distance(dist, is_processed);

        is_processed[u] = true;

        for (int v = 0; v < size; v++) {
            if (!is_processed[v] && graph[u][v]
                && dist[u] != numeric_limits<float>::infinity()
                && dist[u] + graph[u][v] < dist[v]) {
                dist[v] = dist[u] + graph[u][v];
            }
        }
    }

    return dist;
}

int main(int _argc, char* _argv[])
{
    Matrix<float> output;

    const Matrix<float> graph = {
        {0.0, numeric_limits<float>::infinity(), -2.0, numeric_limits<float>::infinity()},
        {4.0, 0.0, 3.0, numeric_limits<float>::infinity()},
        {numeric_limits<float>::infinity(), numeric_limits<float>::infinity(), 0.0, 2.0},
        {numeric_limits<float>::infinity(), -1.0, numeric_limits<float>::infinity(), 0.0}
    };

    Dijkstra<float> dijkstra = Dijkstra<float>(graph);

    for (int i = 0; i < graph.size(); i++) {
        vector<float> result = dijkstra.from_vertex(i);
        output.push_back(result);
    }

    cout << "Vertex\t\tDistance" << endl;
    for (int i = 0; i < graph.size(); i++) {
        for (int j = 0; j < graph.size(); j++) {
            cout << j << "\t\t" << output[i][j] << endl;
        }
    }

    return 0;
}
