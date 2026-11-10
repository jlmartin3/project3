#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <regex>
#include <map>
#include <iomanip>



// 1. Finds the matching closing parenthesis using depth tracking
size_t findMatchingParenthesis(const std::string& smiles, size_t openPos) {
    size_t j = openPos + 1;
    int depth = 1;
    while (j < smiles.size() && depth > 0) {
        if (smiles[j] == '(') depth++;
        else if (smiles[j] == ')') depth--;
        if (depth > 0) j++;
    }
    return j;
}

// 2. Extracts the clean element symbol from inside brackets (e.g., [C@H] -> C, [Na+] -> Na)
std::string extractElementFromBracket(const std::string& content) {
    size_t idx = 0;
    // Skip numbers/symbols at the start (like isotopes e.g., 13C)
    while (idx < content.size() && !std::isupper(content[idx])) {
        idx++;
    }
    if (idx >= content.size()) return "C"; // Fallback

    std::string elem = "";
    elem += content[idx];
    
    // Check if the element has a second lowercase letter (e.g., Na, Cl, Br)
    if (idx + 1 < content.size() && std::islower(content[idx + 1])) {
        elem += content[idx + 1];
    }
    return elem;
}

// 3. Recognizes standard atoms, including 2-character elements like Cl and Br
std::string parseNextAtom(const std::string& smiles, size_t& i) {
    char ch = smiles[i];

    // Lookahead for 2-character elements
    if (ch == 'C' && i + 1 < smiles.size() && smiles[i + 1] == 'l') {
        i++; 
        return "Cl";
    }
    if (ch == 'B' && i + 1 < smiles.size() && smiles[i + 1] == 'r') {
        i++; 
        return "Br";
    }

    // Single-character elements (aliphatic & aromatic)
    if (ch == 'C' || ch == 'c' || ch == 'O' || ch == 'o' || ch == 'N' || ch == 'n' || ch == 'F' || ch == 'I' || ch == 'P' || ch == 'S' || ch == 's') {
        return std::string(1, ch);
    }

    return ""; // Not an atom symbol
}

// 4. Handles ring closure digits following an atom or bracket
void processRingClosures(const std::string& smiles, size_t& i, int newNode, 
                         std::map<int, int>& ringOpenMap, 
                         std::vector<std::pair<int, int>>& Edges) {
    while (i + 1 < smiles.size() && std::isdigit(smiles[i + 1])) {
        i++;
        int ringNum = smiles[i] - '0';

        if (ringOpenMap.find(ringNum) == ringOpenMap.end()) {
            ringOpenMap[ringNum] = newNode; // Open ring
        } else {
            int startNode = ringOpenMap[ringNum];
            Edges.push_back({startNode, newNode}); // Close ring
            ringOpenMap.erase(ringNum);
        }
    }
}

// 5. The main coordinator function
void ExtractEdges(const std::string& smiles, 
                  std::map<int, std::string>& NodesMap, 
                  std::vector<std::pair<int, int>>& Edges,
                  int& currentIndex,
                  std::map<int, int>& ringOpenMap) {

    int prevNode = currentIndex;
    bool pendingDoubleBond = false;
    bool pendingTripleBond = false;

    for (size_t i = 0; i < smiles.size(); i++) {
        char ch = smiles[i];

        // Check for bond markers and component separators
        if (ch == '=') {
            pendingDoubleBond = true;
            continue;
        }
        if (ch == '#') {
            pendingTripleBond = true;
            continue;
        }
        if (ch == '.') {
            // Reset prevNode so the next part of the SMILES string 
            // doesn't form a covalent bond with the previous component
            prevNode = 0;
            continue;
        }

        if (ch == '(') {
            size_t j = findMatchingParenthesis(smiles, i);
            std::string sub = smiles.substr(i + 1, j - (i + 1));
            
            int branchParent = prevNode;
            ExtractEdges(sub, NodesMap, Edges, branchParent, ringOpenMap);
            
            i = j; // Skip past branch
        } 
        else if (ch == '[') {
            size_t j = i + 1;
            while (j < smiles.size() && smiles[j] != ']') {
                j++;
            }

            if (j < smiles.size()) {
                std::string bracketContent = smiles.substr(i + 1, j - (i + 1));
                std::string atomStr = extractElementFromBracket(bracketContent);

                if (!atomStr.empty()) {
                    int newNode = NodesMap.size() + 1;
                    NodesMap[newNode] = atomStr;

                    if (prevNode != 0) {
                        Edges.push_back({prevNode, newNode});
                        if (pendingDoubleBond) {
                            Edges.push_back({prevNode, newNode});
                            pendingDoubleBond = false;
                        }
                        if (pendingTripleBond) {
                            Edges.push_back({prevNode, newNode});
                            Edges.push_back({prevNode, newNode});
                            pendingTripleBond = false;
                        }
                    }
                    prevNode = newNode;
                    currentIndex = newNode;

                    processRingClosures(smiles, j, newNode, ringOpenMap, Edges);
                }
                i = j; // Skip past the entire bracket block
            }
        }
        else {
            std::string atomStr = parseNextAtom(smiles, i);

            if (!atomStr.empty()) {
                if (islower(atomStr[0])) {
                    atomStr[0] = std::toupper(atomStr[0]);
                }

                int newNode = NodesMap.size() + 1;
                NodesMap[newNode] = atomStr;

                if (prevNode != 0) {
                    Edges.push_back({prevNode, newNode});

                    if (pendingDoubleBond) {
                        Edges.push_back({prevNode, newNode});
                        pendingDoubleBond = false;
                    }
                    if (pendingTripleBond) {
                        Edges.push_back({prevNode, newNode});
                        Edges.push_back({prevNode, newNode});
                        pendingTripleBond = false;
                    }
                }
                prevNode = newNode;
                currentIndex = newNode;

                processRingClosures(smiles, i, newNode, ringOpenMap, Edges);
            }
        }
    }
}


// Converts NodesMap and Edges into a 2D adjacency matrix (1-indexed to match your nodes)
std::vector<std::vector<int>> ConvertToAdjacencyMatrix(
    const std::map<int, std::string>& NodesMap, 
    const std::vector<std::pair<int, int>>& Edges) {
    
    int n = NodesMap.size();
    
    // Initialize an (n + 1) x (n + 1) matrix with zeros
    // (Size n + 1 allows 1-based indexing matching your node IDs 1 to n)
    std::vector<std::vector<int>> adjMatrix(n + 1, std::vector<int>(n + 1, 0));

    for (const auto& edge : Edges) {
        int u = edge.first;
        int v = edge.second;

        // Ensure indices are within bounds
        if (u >= 1 && u <= n && v >= 1 && v <= n) {
            // Increment for undirected graph (and handles bond orders via edge duplication)
            adjMatrix[u][v]++;
            adjMatrix[v][u]++;
        }
    }

    return adjMatrix;
}


void PrintAdjacencyMatrix(const std::map<int, std::string>& NodesMap, 
                          const std::vector<std::vector<int>>& matrix) {
    int n = NodesMap.size();

    // Print header row (Node IDs & Elements)
    std::cout << "     ";
    for (int i = 1; i <= n; i++) {
        std::cout << i << "(" << NodesMap.at(i) << ") ";
    }
    std::cout << "\n";

    // Print rows
    for (int i = 1; i <= n; i++) {
        std::cout << i << "(" << NodesMap.at(i) << ") ";
        for (int j = 1; j <= n; j++) {
            std::cout << "  " << matrix[i][j] << "   ";
        }
        std::cout << "\n";
    }
}





#endif