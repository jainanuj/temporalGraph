//
//  main.cpp
//  XuantemporalGraph
//
//  Created by Anuj Jain on 11/20/20.
//

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include "graph.h"
#include "regHeap.hpp"
using namespace std;
#define VERTS_TEST 10

void testRegHeap() {
    regHeap regMinHeap(VERTS_TEST);
    for (int i = 0; i < VERTS_TEST-1; i++)
    {
        int r = rand() % 100000;
        regMinHeap.insert(i, r);
        cout << "(" << i << "," << r << ")   ";
    }
    
//    pair<int, int> minNK = fibMinHeap.removeMin();
//    cout << "Min node: " << get<0>(minNK) << " Val: " << get<1>(minNK) << "\n";
    regMinHeap.decrease_key(5, 2);
    pair<int, int> minNK = regMinHeap.getMin();
    cout << "Min node: " << get<0>(minNK) << " Val: " << get<1>(minNK) << "\n";
    minNK = regMinHeap.removeMin();
    cout << "Min node: " << get<0>(minNK) << " Val: " << get<1>(minNK) << "\n";
//    fibMinHeap.printHeap();
    int i = 0;
    while (!regMinHeap.empty())
    {
        pair<int, int> minNK = regMinHeap.removeMin();
        cout << i++ << ": Min node: " << get<0>(minNK) << " Val: " << get<1>(minNK) << "\n";
    }
    std::cout << "Hello, World!\n";
}

int main(int argc, const char * argv[]) {
    
    
//    testRegHeap();
//    return 0;
    
    const char* option = argv[1];
    
    if(!strcmp(option,"wu"))
    {
        Graph::wuGraph(argv[2], 2);     //Drop the element after u, v.
        return 0;
    }

    Timer t;
    t.start();

    Graph g(argv[2]);
    g.initial_query();

    t.stop();
    cout << "Reading time: " << t.GetRuntime() << "\n";


    if(!strcmp(option,"earliest"))
    {
        g.run_earliest_arrival();
    }
    if(!strcmp(option,"shortest"))
    {
        g.run_shortest();
    }
    std::cout << "Hello, World!\n";
    return 0;
}
