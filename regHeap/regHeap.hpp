//
//  regHeap.hpp
//  XuantemporalGraph
//
//  Created by Anuj Jain on 3/8/21.
//

#ifndef regHeap_hpp
#define regHeap_hpp

#include <stdio.h>
#include <vector>
#include <tuple>
#include <cstring>

using namespace std;
#define VECFIB_NULL -1

class regHeap
{
private:
    vector<int> nodeLocInHeap;
    vector<pair<int, int>> regHeapQ;      //nodeId, key
    void bubbleUp(int itemLoc);
    void bubbleDown(int itemLoc);
    int getParent(int itemLoc);
    int getLChild(int itemLoc);
    int getRChild(int itemLoc);
    void swapItems(int item1, int item2);

public:
    int numNodes;

    void insert(int nodeID, int key);
    void decrease_key(int itemLoc, int newValue);
    pair<int, int> removeMin();
    pair<int, int> getMin();
    bool partOfHeap(int node);
    bool empty();
    regHeap(int V=0);    
};

#endif /* regHeap_hpp */
