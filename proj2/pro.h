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
#include <math.h>
#include <mpi.h>

#define ROOT 0

using namespace std;

class Edge {
public:
    int id; // ID hrany
    char from; // Hodnota uzla odkiaľ ide hrana
    char to; // Hodnota uzla kam smeruje hrana
    static int idCounter; // Počítadlo ako identifikovať nové hrany, začína od 0 pre jednoduchosť

    Edge();
    Edge(char from, char to);
    Edge(int id, char from, char to);
    void printEdge();
};

// Strom T = (V,E), V - vrcholy, E - hrany
class Tree {
public:
    string nodes; // Reťazec uzlov
    vector<Edge> edges; // vektor hrán
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

/**
 * Vytvorí binárny strom na základe reťazca hodnot uzlov
 * @param nodes string Reťazec s hodnotami uzlov
 * @return Tree Binárny strom
 */
Tree createTree(string nodes);
void printTree(Tree tree); // pomocná funkcia
void printEdges(vector<Edge> edges); // pomocná funkcia

/**
 * Konvertuje binárny strom na orientovaný graf pridaním reverzných hrán
 * @param tree Ukazovateľ na binárny strom, ktorý má byť menený
 */
void convertTreeToOrientedGraph(Tree* tree);

/**
 * Vytvorí novú hranu, ktorý bude reverzná k poskytnutej hrany. Automaticky zvyšuje counter na id nového uzlu.
 * @param edge Hrana ku ktorej sa má vytvoriť opačná hrana
 * @return edge Reverzná hrana
 */
Edge createReverseEdge(Edge edge);

/**
 * Vytvorí zoznam susedov pre všetky uzly. Reprezentované ako vektor vektorov susedov pre ľahšie narábanie.
 * Prvý index určuje o ktorý uzol sa zaujímame, druhý index už indexuje jeho susedov.
 * @param og Orientovaný graf
 * @return Vektor vektorov typu Neighbour - reprezentuje zozname susedov
 */
vector<vector<Neighbour>> createNeighboursVector(Tree og); // og - orientovaný graf

/**
 * Vytvorí vektor hrán, ktoré všetky vychádzajú z uzla so zvolenou hodnotou
 * @param edges Vektor všetkých hrán
 * @param from Hodnota uzla, z ktorého má hrana vychádzať
 * @return Vektor hrán splňujúcich podmienku
 */
vector<Edge> getEdgesWithFrom(vector<Edge> edges, char from);

/**
 * Nájde opačnú hranu k zvolenej hrane. Predpokladá sa, že určite existuje.
 * @param edges Vektor všetkých hrán
 * @param edge Vzorová hrana ku ktorej sa hľadá opačná hrana
 * @return Objekt opačnej hrany
 */
Edge findReverseEdge(vector<Edge> edges, Edge edge);
int sendNeighboursToProcessor(vector<Neighbour> neighbours, int receiverRank);
int receiveNeighbours(vector<Neighbour> *neighbours);
Neighbour findNeighbourWithFromId(vector<vector<Neighbour>> neighboursVec, int fromId);
void getReverseEdgeIndexes(vector<vector<Neighbour>> neighboursVec, int revEdgeId, int *vecIdx, int *insideIdx);
int getEtourElement(int edgeIdToProcess, vector<vector<Neighbour>> neighboursVec);
int determineEdgeWeight(int edgeId, vector<vector<Neighbour>> neighboursVec);
void sendEtourElement(int eTour, int receiverRank);
void receiveEtourElement(int senderRank, int *eTourVal);
void correctEtour(int eTourLen, int eTours[], Tree tree);
int suffixSum(int numProcessors, int edgeId, int eTours[], int weights[]);
Edge findEdgeWithId(vector<Edge> edges, int id);

#endif
