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
using namespace std;

int main(int argc, const char * argv[]) {
    const char* option = argv[1];
    
    
    Graph g(argv[2]);
//    g.wuGraph("/Users/anujjain/research/temporalGraph/WuTemporalGraph/tempath/sx-mathoverflow-c2a.txt", 1);
    g.initial_query();


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
