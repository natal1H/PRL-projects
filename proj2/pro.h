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
#include <map>
#include <mpi.h>

using namespace std;

class Edge {
public:
    int id;
    char from;
    char to;

    static int idCounter;
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
};

Tree createTree(string nodes);
void printTree(Tree tree); // pomocná funkcia
void printEdges(vector<Edge> edges); // pomocná funkcia
void convertTreeToOrientedGraph(Tree* tree);
Edge createReverseEdge(Edge edge);
map<char, vector<Neighbour>> createNeighboursMap(Tree og); // og - orientovaný graf
vector<Edge> getEdgesWithFrom(vector<Edge> edges, char from);
Edge findReverseEdge(vector<Edge> edges, Edge edge);

#endif