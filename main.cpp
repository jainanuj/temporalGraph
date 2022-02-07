//
//  main.cpp
//  XuantemporalGraph
//
//  Created by Anuj Jain on 11/20/20.
//
// ./XuantemporalGraph earliest|shortest|mwf <filePath>

// ./XuantemporalGraph wu <fileName> <1/2>(drop num Lines to drop)  <0/1>(normalize or not) //

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
//#include "graph.h"
#include "graphDualCriteria.hpp"
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
    int contactSeq = 0;
    
    if(!strcmp(option,"wu"))
    {
        int numDrop = 1, normalize=0;
        if (argc > 3)
            sscanf(argv[3],"%d",&numDrop);
        if (argc > 4)
            sscanf(argv[4],"%d",&normalize);

        Graph::wuGraph(argv[2], 2, numDrop, normalize);     //2=Drop the element after u, v. numDrop= num of lines to drop. normalize=normalize the timestamps.
        return 0;
    }
    
    if(!strcmp(option,"xuan"))
    {
        Graph::readWuFile(argv[2]);
        return 0;
    }

    Timer t;
    t.start();

    if (argc > 4)
        contactSeq = 1;
    Graph *g;
    if(!strcmp(option,"mwf") || (!strcmp(option,"mhf")) )
        g = new GraphDualCriteria(argv[2], contactSeq, option);
    else
        g = new Graph(argv[2], contactSeq, option);
    if (argc == 4)
        g->initial_query(argv[3]);
    else
        g->initial_query();

    t.stop();
    cout << "Reading time: " << t.GetRuntime() << "\n";


    if(!strcmp(option,"earliest"))
    {
        g->run_earliest_arrival();
    }
    else if(!strcmp(option,"shortest"))
    {
        g->run_shortest();
    }
    else if(!strcmp(option,"mwf"))
    {
        g->run_mwf();
    }
    else if(!strcmp(option,"mhf"))
    {
        g->run_mhf();
    }
    std::cout << "Hello, World!\n";
    return 0;
}
