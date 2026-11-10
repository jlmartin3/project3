
#include "utils.h"
#include <iostream>
#include <string>
#include <sstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

/* 
	James Martin
	UFID: 6121-1142
	11/2/2026
	COP3530
 
    Project 2: Simplified PageRank
    
    Description: A basic implimentation of Google's original PageRank Algorithim 
    using an adjacency list
*/
 
int main() {

    std::string smiles;
    int no_of_lines;

    std::cin >> no_of_lines; 

    for(int i = 0;i < no_of_lines; i++){

        


        std::cin >> smiles ;
        std::cout << "---------------------------------------------------------------------" << std::endl;

        std::cout << "Input: " << i << std::endl ;
        std::cout << "SMILES string: " << smiles << std::endl;
    
        //std::string smiles = "O=C(O)[C@H]1CCCNC1"; // Benzene ring example
        std::map<int, std::string> NodesMap;
        std::vector<std::pair<int, int>> Edges;
        std::map<int, int> ringOpenMap;

        int startIndex = 0;
        ExtractEdges(smiles, NodesMap, Edges, startIndex, ringOpenMap);

        std::cout << "Nodes:\n";
        for (size_t i = 1; i <= NodesMap.size(); i++) {
            std::cout << i << ": " << NodesMap[i] << "\n";
        }

        std::cout << "\nEdges:\n";
        for (const auto& edge : Edges) {
            std::cout << edge.first << " -- " << edge.second << "\n";
        }

        std::vector<std::vector<int>> TestMatrix = ConvertToAdjacencyMatrix(NodesMap, Edges);

        std::cout << "\nAdjacency Matrix:\n";
        PrintAdjacencyMatrix(NodesMap, TestMatrix);



    }

    


    return 0;

    
}

