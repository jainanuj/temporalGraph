//
//  regHeap.cpp
//  XuantemporalGraph
//
//  Created by Anuj Jain on 3/8/21.
//

#include "regHeap.hpp"
#define HEAP_ROOT 0

regHeap::regHeap(int V)
{
    numNodes = 0;
    nodeLocInHeap.resize(V);
    nodeLocInHeap.assign(V, -1);
    
    regHeapQ.resize(V);     //HArr
}

void regHeap::insert(int nodeID, int key)
{
    numNodes++;
//    regHeapQ.push_back(make_pair(nodeID, key));   //HArr
    regHeapQ[numNodes-1] = make_pair(nodeID, key);  //HArr
    nodeLocInHeap[nodeID] = numNodes-1;
    bubbleUp(numNodes-1);
}

pair<int, int> regHeap::removeMin()
{
    int nodeId;
    pair<int, int> retVal = regHeapQ[HEAP_ROOT];
    nodeLocInHeap[get<0>(retVal)] = -1;     //Older Heap top is not in heap anymore
    regHeapQ[HEAP_ROOT] = regHeapQ[numNodes-1];
    numNodes--;
//    regHeapQ.pop_back();      //HArr
    nodeId = get<0>(regHeapQ[HEAP_ROOT]);
    nodeLocInHeap[nodeId] = HEAP_ROOT;
    bubbleDown(HEAP_ROOT);
    return retVal;
}

bool regHeap::empty()
{
    if (numNodes <= 0)
        return true;
    else
        return false;
}

bool regHeap::partOfHeap(int node)
{
    if ( (node > nodeLocInHeap.size()) || (node < 0) )
        return false;
    else
        return (nodeLocInHeap[node] == -1) ? false : true;
}



pair<int, int> regHeap::getMin()
{
    return regHeapQ[HEAP_ROOT];
}

void regHeap::decrease_key(int nodeId, int newValue)
{
    int itemLoc = nodeLocInHeap[nodeId];
    if (itemLoc == -1)
        return;
    else if (get<1>(regHeapQ[itemLoc]) <= newValue)
        return;
    else
    {
        get<1>(regHeapQ[itemLoc]) = newValue;
        bubbleUp(itemLoc);
    }
}

void regHeap::bubbleUp(int itemLoc)
{
    if (itemLoc == HEAP_ROOT)       //at root.
        return;
    int p = getParent(itemLoc);
    while (get<1>(regHeapQ[itemLoc]) < get<1>(regHeapQ[p]))
    {
        swapItems(itemLoc, p);
        itemLoc = p;
        if (itemLoc == HEAP_ROOT)       //reached root.
            return;
        p = getParent(itemLoc);
    }
    return;
}

void regHeap::bubbleDown(int itemLoc)
{
    int child1 = getLChild(itemLoc);
    int child2 = getRChild(itemLoc);
    int smallChild;
    if ((child1 == -1) && (child2 == -1))       //at leaf.
        return;
    else if(child2 == -1)
        smallChild = child1;
    else
        smallChild = (get<1>(regHeapQ[child1]) < get<1>(regHeapQ[child2])) ? child1 : child2;
    while (get<1>(regHeapQ[itemLoc]) > get<1>(regHeapQ[smallChild]))
    {
        swapItems(itemLoc, smallChild);
        itemLoc = smallChild;
        child1 = getLChild(itemLoc);
        child2 = getRChild(itemLoc);
        if ((child1 == -1) && (child2 == -1))       //Reached a leaf.
            return;
        else if(child2 == -1)
            smallChild = child1;
        else
            smallChild = (get<1>(regHeapQ[child1]) < get<1>(regHeapQ[child2])) ? child1 : child2;
    }
    return;
    
}

int regHeap::getParent(int itemLoc)
{
    return (int)((itemLoc+1)/2) -1;
}

void regHeap::swapItems(int pos1, int pos2)
{
    pair<int, int> nodeKey1 = regHeapQ[pos1];
    regHeapQ[pos1] = regHeapQ[pos2];
    regHeapQ[pos2] = nodeKey1;
    nodeLocInHeap[get<0>(regHeapQ[pos1])] = pos1;       //Make sure nodeLoc is pointing to correct position in heap.
    nodeLocInHeap[get<0>(regHeapQ[pos2])] = pos2;
}

int regHeap::getLChild(int itemLoc)
{
    int childPos = 2*(itemLoc+1);
    return (childPos <= numNodes)? (childPos -1) : -1;
}

int regHeap::getRChild(int itemLoc)
{
    int childPos = 2*(itemLoc+1)+1;
    return (childPos <= numNodes)? (childPos -1) : -1;
}
