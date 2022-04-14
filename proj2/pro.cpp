/**
 * @file pro.cpp
 * @author Natália Holková (xholko02@stud.fit.vutbr.cz)
 * @brief PRL - project 2 - Priradenie poradia preorder vrcholom
 * @date 2022
 */

#include "pro.h"

using namespace std;

Tree createTree(string nodes) {
    Tree t;
    t.nodes = nodes;

    vector<Edge> edges;

    int nodesNum = nodes.length();
    for (int i = 0; i < nodesNum; i++) {
        char parentChar = nodes.at(i);

        // Vypočítaj indexy ľavého a pravého potomka v strome
        int leftChildIdx = (i * 2) + 1;
        int rightChildIdx = (i * 2) + 2; // because indexing starts from 0

        if (leftChildIdx < nodesNum) {
            // Má ľavého potomka - vytvor hranu z rodiča do ľavého potomka
            char leftChildChar = nodes.at(leftChildIdx);
            Edge parentToLeftChild;
            parentToLeftChild.id = Edge::idCounter++; // prirať hrane ID a zároveň zvýš counter
            parentToLeftChild.from = parentChar;
            parentToLeftChild.to = leftChildChar;

            // Pridaj hrany do zoznamu (vektoru) hrán
            edges.push_back(parentToLeftChild);
        }
        if (rightChildIdx < nodesNum) {
            char rightChildChar = nodes.at(rightChildIdx);
            Edge parentToRightChild;
            parentToRightChild.id = Edge::idCounter++; // prirať hrane ID a zároveň zvýš counter
            parentToRightChild.from = parentChar;
            parentToRightChild.to = rightChildChar;

            // Pridaj hrany do zoznamu (vektoru) hrán
            edges.push_back(parentToRightChild);
        }
    }

    t.edges = edges;
    return t;
}

// Pomocná funkcia na výpis hrán
void printEdges(vector<Edge> edges) {
    cout << "Edges: " << endl;
    for(const auto& edge: edges) {
        cout << "ID " << edge.id << " from " << edge.from << " to " << edge.to << "\n";
    }
}

// Pomocná funkcia na výpis stromu
void printTree(Tree tree) {
    cout << "Nodes: " << tree.nodes << endl;
    printEdges(tree.edges);
}

Edge createReverseEdge(Edge edge) {
    Edge reverseEdge;
    reverseEdge.id = Edge::idCounter++;
    reverseEdge.from = edge.to;
    reverseEdge.to = edge.from;
    return reverseEdge;
}

void convertTreeToOrientedGraph(Tree* tree) {
    vector<Edge> newEdges;

    // Prechádzame všetkými hranami stromu
    // Každú hranu (u,v) nahradíme dvoma hranami (u,v) a (v,u)
    for(const auto& oldEdge: tree->edges) {
        newEdges.push_back(oldEdge); // Vlož pôvodnú hranu

        // Vytvor opačnú hranu
        Edge reverseEdge = createReverseEdge(oldEdge);
        newEdges.push_back(reverseEdge);
    }
    tree->edges = newEdges;
}

vector<Edge> getEdgesWithFrom(vector<Edge> edges, char from) {
    vector<Edge> foundEdges;
    for(const auto& edge: edges) {
        if (edge.from == from) {
            // Vrchol FROM sa zhoduje s hľadaným
            foundEdges.push_back(edge);
        }
    }

    return foundEdges;
}

Edge findReverseEdge(vector<Edge> edges, Edge edge) {
    // Pozor: predpokladá sa, že táto hrana určite existuje (poznáme hranu a chceme nájsť jej reverse)
    Edge reverseEdge;
    for(const auto& e: edges) {
        if (e.from == edge.to && e.to == edge.from) {
            reverseEdge = e;
            break;
        }
    }
    return reverseEdge;
}

map<char, vector<Neighbour>> createNeighboursMap(Tree og) {
    map<char, vector<Neighbour>> neighboursMap;
    // Prechádzaj všetkými vrcholmi orientovaného grafu
    for (int i = 0; i < og.nodes.length(); i++) {
        char node = og.nodes.at(i);
        // Nájdi všetky hrany vychádzajúce z tohto uzlu
        vector<Edge> edgesFromThisNode = getEdgesWithFrom(og.edges, node);
        vector<Neighbour> neighbours;

        for(const auto& edgeFrom: edgesFromThisNode) {
            // Nájdi reverse hranu k hrane
            Edge edgeTo = findReverseEdge(og.edges, edgeFrom);
            // Vytvor suseda
            Neighbour neighbour;
            neighbour.edgeFromNode = edgeFrom;
            neighbour.edgeToNode = edgeTo;
            neighbours.push_back(neighbour);
        }

        neighboursMap[node] = neighbours;
    }

    return neighboursMap;
}

int Edge::idCounter = 0; // TODO: od akého čísla číslovať hrany?

int main(int argc, char *argv[]) {
    int rank, num_proc; 

    string nodesString = argv[1];

    MPI_Init(&argc, &argv); // initialize the MPI environment
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc); // get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get the rank of the process

    cout << "Hello world from rank " << rank << endl;
    
    if (rank == 0) {
        // Vytvoríme strom zo zadaných vrcholov - treba vytvoriť hrany
        Tree tree = createTree(nodesString);
        // Konvertujeme strom na orientovaný graf
        convertTreeToOrientedGraph(&tree);
        // Vytvor zoznam susedov
        map<char, vector<Neighbour>> neighboursMap = createNeighboursMap(tree);
        cout << "Neighbours determined." << endl;
    }

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}