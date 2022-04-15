/**
 * @file pro.h
 * @author Natália Holková (xholko02@stud.fit.vutbr.cz)
 * @brief PRL - project 2 - Priradenie poradia preorder vrcholom
 * @date 2022
 */

#ifndef PRO_H
#define PRO_H

#include <iostream>
#include <string>
#include <vector>
#include <mpi.h>

#define ROOT 0

using namespace std;

class Edge {
public:
    int id;
    char from;
    char to;
    static int idCounter;

    Edge();
    Edge(char from, char to);
    Edge(int id, char from, char to);
    void printEdge();
};

// Strom T = (V,E), V - vrcholy, E - hrany
class Tree {
public:
    string nodes;
    vector<Edge> edges;
};

class Neighbour {
public:
    Edge edgeToNode;
    Edge edgeFromNode;
    Neighbour();
    Neighbour(int edgeIds[2], string edgesDirs);
    string convertEdgeNodesToString();
    void convertEdgeIdToArr(int ids[]);
    void printNeighbour();
};

Tree createTree(string nodes);
void printTree(Tree tree); // pomocná funkcia
void printEdges(vector<Edge> edges); // pomocná funkcia
void convertTreeToOrientedGraph(Tree* tree);
Edge createReverseEdge(Edge edge);
vector<vector<Neighbour>> createNeighboursVector(Tree og); // og - orientovaný graf
vector<Edge> getEdgesWithFrom(vector<Edge> edges, char from);
Edge findReverseEdge(vector<Edge> edges, Edge edge);
int sendNeighboursToProcessor(vector<Neighbour> neighbours, int receiverRank);
int receiveNeighbours(vector<Neighbour> *neighbours);
Neighbour findNeighbourWithFromId(vector<vector<Neighbour>> neighboursVec, int fromId);
void getReverseEdgeIndexes(vector<vector<Neighbour>> neighboursVec, int revEdgeId, int *vecIdx, int *insideIdx);
int getEtourElement(int edgeIdToProcess, vector<vector<Neighbour>> neighboursVec);

#endif