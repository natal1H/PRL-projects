/**
 * @file pro.cpp
 * @author Natália Holková (xholko02@stud.fit.vutbr.cz)
 * @brief PRL - project 2 - Priradenie poradia preorder vrcholom
 * @date 2022
 */

#include "pro.h"

using namespace std;

Edge::Edge() {}

Edge::Edge(char from, char to) {
    this->from = from;
    this->to = to;
    this->id = Edge::idCounter++;
}

Edge::Edge(int id, char from, char to) {
    this->from = from;
    this->to = to;
    this->id = id;
}

void Edge::printEdge() {
    cout << "ID: " << this->id << "; from: " << this->from << "; to: " << this->to << endl;
}

Neighbour::Neighbour() {}

Neighbour::Neighbour(int edgeIds[2], string edgesDirs) {
    Edge edgeFrom(edgeIds[0], edgesDirs.at(0), edgesDirs.at(1));
    Edge edgeTo(edgeIds[1], edgesDirs.at(2), edgesDirs.at(3));
    this->edgeFromNode = edgeFrom;
    this->edgeToNode = edgeTo;
}

string Neighbour::convertEdgeNodesToString() {
    char arr[] = {this->edgeFromNode.from, this->edgeFromNode.to,
                 this->edgeToNode.from, this->edgeToNode.to, '\0'};
    string res(arr);
    return res;
}

void Neighbour::convertEdgeIdToArr(int ids[]) {
    ids[0] = this->edgeFromNode.id;
    ids[1] = this->edgeToNode.id;
}

void Neighbour::printNeighbour() {
    cout << "Neighbour\nFROM: ";
    this->edgeFromNode.printEdge();
    cout << "TO: ";
    this->edgeToNode.printEdge();
}

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
            Edge parentToLeftChild(parentChar, leftChildChar);
            edges.push_back(parentToLeftChild); // Pridaj hrany do zoznamu (vektoru) hrán
        }
        if (rightChildIdx < nodesNum) {
            char rightChildChar = nodes.at(rightChildIdx);
            Edge parentToRightChild(parentChar, rightChildChar);
            edges.push_back(parentToRightChild); // Pridaj hrany do zoznamu (vektoru) hrán
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
    Edge reverseEdge(edge.to, edge.from);
    return reverseEdge;
}

void convertTreeToOrientedGraph(Tree* tree) {
    vector<Edge> newEdges;
    newEdges = tree->edges;

    // Prechádzame všetkými hranami stromu
    // Každú hranu (u,v) nahradíme dvoma hranami (u,v) a (v,u)
    for(const auto& oldEdge: tree->edges) {
        Edge reverseEdge = createReverseEdge(oldEdge); // Vytvor opačnú hranu
        newEdges.push_back(reverseEdge);
    }
    tree->edges = newEdges;
}

vector<Edge> getEdgesWithFrom(vector<Edge> edges, char from) {
    vector<Edge> foundEdges;
    for(const auto& edge: edges)
        if (edge.from == from) // Vrchol FROM sa zhoduje s hľadaným
            foundEdges.push_back(edge);
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

vector<vector<Neighbour>> createNeighboursVector(Tree og) {
    vector<vector<Neighbour>> neighboursVec;
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

        neighboursVec.push_back(neighbours);
    }

    return neighboursVec;
}

int sendNeighboursToProcessor(vector<Neighbour> neighbours, int receiverRank) {
    // Konvertovať vektor na pole
    int numNeighbours = neighbours.size();
    int *edgeIds = new int[numNeighbours * 2];
    string edgesStr = "";

    // V cykle každého suseda konvertovať na int[] ID hrán a string odkiaľ-kam idú hrany
    for (int i = 0; i < numNeighbours; i++) {
        int tmpIds[2] = {};
        neighbours[i].convertEdgeIdToArr(tmpIds);
        edgeIds[2*i] = tmpIds[0]; edgeIds[2*i+1] = tmpIds[1];
        edgesStr.append(neighbours[i].convertEdgeNodesToString());
    }

    // Poslať 3 správy: najprv počet susedov, potom int array s edge IDs, potom char[] s odkiaľ kam idú hrany
    MPI_Send(&numNeighbours, 1, MPI_INT, receiverRank, 0, MPI_COMM_WORLD);
    MPI_Send(edgeIds, numNeighbours * 2 ,MPI_INT, receiverRank, 0, MPI_COMM_WORLD);
    MPI_Send(edgesStr.c_str(), numNeighbours * 4 ,MPI_CHAR, receiverRank, 0, MPI_COMM_WORLD);

    return 0;
}

int receiveNeighbours(vector<Neighbour> *neighbours) {
    // Prijmi počet susedov
    int numNeighbours;
    MPI_Recv(&numNeighbours, 1, MPI_INT, ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Prijmi pole s ID hrán susedov
    int *edgeIds = new int[numNeighbours * 2];
    MPI_Recv(edgeIds, numNeighbours * 2, MPI_INT, ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Prijmi reťazec so smerovaniami hrán
    char *edgesCharArr = new char[numNeighbours * 4];
    MPI_Recv(edgesCharArr, numNeighbours * 4, MPI_CHAR, ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    string edgesStr(edgesCharArr);

    // Prekonvertovať do Neighbour objektov
    for (int i = 0; i < numNeighbours; i++) {
        int tmpIds[] = {edgeIds[i*2], edgeIds[i*2+1]};
        Neighbour tmpNeighbour(tmpIds, edgesStr.substr(i*4,4));
        neighbours->push_back(tmpNeighbour);
    }

    return 0;
}

Neighbour findNeighbourWithFromId(vector<vector<Neighbour>> neighboursVec, int fromId) {
    // Pozor: predpokladá že určite taký sused existuje

    Neighbour foundNeighbour;
    bool found = false;

    // Iterovať nad zoznamami susedov všetkých uzlov
    for (int vecIdx = 0; vecIdx < neighboursVec.size(); vecIdx++) {
        // Iterovať už v zozname susedov konkrétneho uzlu
        for (int insideIdx = 0; insideIdx < neighboursVec[vecIdx].size(); insideIdx++) {
            Neighbour tmpNeighbour = neighboursVec[vecIdx][insideIdx];
            if (tmpNeighbour.edgeFromNode.id == fromId) {
                foundNeighbour = tmpNeighbour;
                found = true;
                break;
            }
        }
        if (found) break;
    }

    return foundNeighbour;
}

void getReverseEdgeIndexes(vector<vector<Neighbour>> neighboursVec, int revEdgeId, int *vecIdx, int *insideIdx) {
    // Iterovať nad zoznamami susedov všetkých uzlov
    for (int i = 0; i < neighboursVec.size(); i++) {
        // Iterovať už v zozname susedov konkrétneho uzlu
        for (int j = 0; j < neighboursVec[i].size(); j++) {
            Neighbour tmpNeighbour = neighboursVec[i][j];
            if (tmpNeighbour.edgeFromNode.id == revEdgeId) {
                *vecIdx = i;
                *insideIdx = j;
                return;
            }
        }
    }
}

int getEtourElement(int edgeIdToProcess, vector<vector<Neighbour>> neighboursVec) {
    // Získať suseda, ktorý obsahuje hranu e
    Neighbour neighbourWithE = findNeighbourWithFromId(neighboursVec, edgeIdToProcess);

    // Vytiahnuť zo suseda ID hrany e_r (reverse hrana)
    int edgeReveseId = neighbourWithE.edgeToNode.id;

    // Vytiahnuť z neighboursVec indexy, kde sa nachádza e_r
    int vecIdx, insideIdx;
    getReverseEdgeIndexes(neighboursVec, edgeReveseId, &vecIdx, &insideIdx);

    // Pozri, či položka s reverseEdge má ešte nasledovníka (cez indexy vektorov)
    if (insideIdx != neighboursVec[vecIdx].size()-1) {
        // Nie je poslednou položkou vektora - vráť ID nasledovnej položky
        return neighboursVec[vecIdx][insideIdx+1].edgeFromNode.id;
    }
    else {
        // reverse hrana je poslednom položkou zoznamu
        return edgeReveseId;
    }
}

int determineEdgeWeight(int edgeId, vector<vector<Neighbour>> neighboursVec) {
    for (int i = 0; i < neighboursVec.size(); i++) {
        // Iterovať už v zozname susedov konkrétneho uzlu
        for (int j = 0; j < neighboursVec[i].size(); j++) {
            Neighbour tmpNeighbour = neighboursVec[i][j];
            if (tmpNeighbour.edgeFromNode.id == edgeId) {
                return 1; // Dopredná hrana
            }
            else if (tmpNeighbour.edgeToNode.id == edgeId) {
                return 0; // Spätná hrana
            }
        }
    }
    return -1; // Chyba, ale nemala by nastať
}

void sendEtourElement(int eTour) {
    // Pošli Etour element ROOT procesoru
    MPI_Send(&eTour, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD);
}

int receiveEtourElement(int senderRank) {
    // Prijmi eTour element a ulož na správne miesto v poli
    int eTourVal;
    MPI_Recv(&eTourVal, 1, MPI_INT, senderRank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return eTourVal;
}

void correctEtour(int eTourLen, int eTours[], Tree tree) {
    char rootNode = tree.nodes.at(0);
    int lastEdgeToRoot;

    for (int i = eTourLen-1; i >= 0; i--) {
        if (tree.edges[i].to == rootNode) {
            lastEdgeToRoot =  tree.edges[i].id;
            eTours[eTourLen-1] = lastEdgeToRoot;
            return;
        }
    }
}

int Edge::idCounter = 0;

int main(int argc, char *argv[]) {
    int rank, num_proc;
    string nodesString = argv[1];
    Tree tree;
    vector<vector<Neighbour>> neighboursVec;

    MPI_Init(&argc, &argv); // initialize the MPI environment
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc); // get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get the rank of the process

    if (rank == ROOT) {
        // Vytvoríme strom zo zadaných vrcholov - treba vytvoriť hrany
        tree = createTree(nodesString);
        // Konvertujeme strom na orientovaný graf
        convertTreeToOrientedGraph(&tree);
        //printTree(tree);

        // Vytvor zoznam susedov
        neighboursVec = createNeighboursVector(tree);

        // ROOT procesor rozošle príslušný zoznam (vektor) susedov každému procesoru
        for (int receiver = 0; receiver < num_proc; receiver++) {
            // Poslať ID hrany o ktorú sa má starať
            MPI_Send(&tree.edges[receiver].id, 1, MPI_INT, receiver, 0, MPI_COMM_WORLD);
            for (int vecIdx = 0; vecIdx < neighboursVec.size(); vecIdx++)
                sendNeighboursToProcessor(neighboursVec[vecIdx], receiver);
        }
    }

    // Všetky procesory príjmu info o susedoch svojho uzla
    int edgeIdToProcess;
    MPI_Recv(&edgeIdToProcess, 1, MPI_INT, ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int vecIdx = 0; vecIdx < nodesString.length(); vecIdx++) {
        int numNeighbours;
        vector<Neighbour> thisNodeNeighbours; // Vektor susedov pre uzol, ktorý spracováva tento rank
        receiveNeighbours(&thisNodeNeighbours);
        neighboursVec.push_back(thisNodeNeighbours);
    }

    // Paralelný výpočet poľa Etour, každý procesor jeden prvok
    int eTourElem = getEtourElement(edgeIdToProcess, neighboursVec);

    // Paralelný výpočet váhy hrany - 1 ak dopredná, 0 ak spätná
    int edgeWeight = determineEdgeWeight(edgeIdToProcess, neighboursVec);

    // Korekcia Eulerovej cesty procesorom ROOT, ROOT prijme
    if (rank == ROOT) {
        int eTours[num_proc];
        eTours[0] = eTourElem; // ROOT rovno uloží svoj vypočítaný element
        for (int proc = 1; proc < num_proc; proc++)
            eTours[proc] = receiveEtourElement(proc);

        cout << "ETOUR: ";
        for (int proc = 0; proc < num_proc; proc++)
            cout << eTours[proc] << " ";
        cout << endl;

        // Uprav Eulerovu cestu - hrana, ktorá posledná vedie do koreňového uzla sa dá na koniec
        correctEtour(num_proc, eTours, tree);

        cout << "ETOUR after correction: ";
        for (int proc = 0; proc < num_proc; proc++)
            cout << eTours[proc] << " ";
        cout << endl;
    }
    else {
        // non-ROOT procesor zasiela svoj eTourElem ROOT procesoru
        sendEtourElement(eTourElem);
    }

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}