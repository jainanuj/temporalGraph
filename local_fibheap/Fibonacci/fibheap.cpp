//
//  fibheap.cpp
//  Fibonacci
//
//  Created by Anuj Jain on 2/27/21.
//

#include <stdio.h>
#include <iostream>
#include "fibheap.h"
#define phi (1 + sqrt(5))/2
fibHeap::fibHeap(int V)
{
    int maxDegree = (int)(log(V)/log(phi)) + 1;
    fibHeapQ.resize(V);
    treeTable.resize(maxDegree);
    fibNodeElement initialNode(VECFIB_NULL, false);
    fibHeapQ.assign(V, initialNode);
    treeTable.assign(maxDegree, make_tuple(VECFIB_NULL, VECFIB_NULL, VECFIB_NULL));
    startTreeList = -1;
    minKeyPtr = VECFIB_NULL;     //Equivalent to NULL
    numMinTrees = 0; numNodes = 0;
}

void fibHeap::initialize(int nodeID)
{
    fibHeapQ[nodeID].leftSibling = nodeID; fibHeapQ[nodeID].rtSibling = nodeID;
    fibHeapQ[nodeID].parent= VECFIB_NULL; fibHeapQ[nodeID].child = VECFIB_NULL;
    fibHeapQ[nodeID].degree = 0; fibHeapQ[nodeID].childCut = false;
}

void fibHeap::insert(int nodeID, int key)
{
    if (key < 0)
    {
        cout << "In fib Heap insert, asked to insert -ve key" << key << " for node:" << nodeID << ". Exiting \n";
        exit(1);
    }
    if (fibHeapQ[nodeID].key == VECFIB_NULL)
    {
        initialize(nodeID);
        fibHeapQ[nodeID].key = key;
        if (minKeyPtr == VECFIB_NULL)
            minKeyPtr = nodeID;
        else
        {
            joinDoubly(minKeyPtr, nodeID);
            updateMinKey(nodeID);
        }
        fibHeapQ[nodeID].inHeap = true;
        numMinTrees++;
        numNodes++;
    }
}

void fibHeap::meld(int fibTree1, int fibTree2)
{
    int temp = fibHeapQ[fibTree1].rtSibling;
    fibHeapQ[fibTree1].rtSibling =  fibHeapQ[fibTree2].rtSibling; fibHeapQ[fibHeapQ[fibTree1].rtSibling].leftSibling = fibTree1;
    fibHeapQ[fibTree2].rtSibling = temp; fibHeapQ[fibHeapQ[fibTree2].rtSibling].leftSibling = fibTree2;
}

void fibHeap::decrease_key(int fibNode, int newValue)
{
    if (newValue < 0)
    {
        cout << "In fib Heap dec key, asked to dec -ve key" << newValue << " for node:" << fibNode << ". Exiting \n";
        exit(1);
    }
    if (fibHeapQ[fibNode].key <= newValue)      //nothing to decrease as new value is higher.
        return;
    fibHeapQ[fibNode].key = newValue;
    if (fibHeapQ[fibNode].parent == VECFIB_NULL)        //Already at root level.
    {
        updateMinKey(fibNode);
        return;
    }
    else
    {
        int p = fibHeapQ[fibNode].parent;
        if (fibHeapQ[p].key <= fibHeapQ[fibNode].key)       //No need to move any nodes.
            return;
        cut(fibNode);       //removes from it's doubly and places in the root level doubly.
        if (fibHeapQ[p].childCut) //If Parent not at root and has childCut set, cascade cut.
            cascadeCut(p);
        else
            fibHeapQ[p].childCut = true;
        updateMinKey(fibNode);
    }
}

pair<int, int> fibHeap::getMin()
{
    return make_pair(minKeyPtr, fibHeapQ[minKeyPtr].key);
}

pair<int, int> fibHeap::removeMin()
{
    int list1 = -1, list2 = -1, numItemsList1 = 0, numItemsList2 = 0;
    if (numMinTrees == 0)
        return make_pair(-1, -1);
    pair<int, int> retVal = make_pair(minKeyPtr, fibHeapQ[minKeyPtr].key);
    fibHeapQ[minKeyPtr].inHeap = false;     //This is not in the heap anymore.

    numItemsList1 = numMinTrees - 1; numItemsList2 = fibHeapQ[minKeyPtr].degree;
    if (numItemsList1 > 0)
    {
        list1 = fibHeapQ[minKeyPtr].leftSibling;
        leaveDoubly(list1, minKeyPtr, numMinTrees);
    }
    else
        list1 = VECFIB_NULL;
    list2 = fibHeapQ[minKeyPtr].child;
    int list2Node = list2;
    for (int i = 0; i < numItemsList2; i++)
    {
        fibHeapQ[list2Node].parent = VECFIB_NULL;
        list2Node = fibHeapQ[list2Node].rtSibling;
    }
    numNodes--;
    fibHeapQ[minKeyPtr].key = VECFIB_NULL;
//    pairwiseCombine(list1, list2, numItemsList1, numItemsList2);
    pairwiseCombine_s(list1, list2, numItemsList1, numItemsList2);
    return retVal;
}

bool fibHeap::empty()
{
    if (numNodes)
        return false;
    else
        return true;
}

bool fibHeap::partOfHeap(int node)
{
    return fibHeapQ[node].inHeap;
}


void fibHeap::joinDoubly(int doubly, int nodeID)
{
    fibHeapQ[nodeID].rtSibling = fibHeapQ[doubly].rtSibling; fibHeapQ[fibHeapQ[nodeID].rtSibling].leftSibling = nodeID;
    fibHeapQ[doubly].rtSibling = nodeID;                     fibHeapQ[nodeID].leftSibling = doubly;
}

void fibHeap::leaveDoubly(int leftOfNode, int nodeID, int count)
{
    if (count == 1)
        return;
    fibHeapQ[leftOfNode].rtSibling = fibHeapQ[nodeID].rtSibling; fibHeapQ[fibHeapQ[leftOfNode].rtSibling].leftSibling = leftOfNode;
    fibHeapQ[nodeID].rtSibling = nodeID; fibHeapQ[nodeID].leftSibling = nodeID;
}

void fibHeap::updateMinKey(int fibNode)
{
    if (fibHeapQ[fibNode].key < fibHeapQ[minKeyPtr].key)
        minKeyPtr = fibNode;
}

void fibHeap::cut(int cutFibNode)
{
    int p = fibHeapQ[cutFibNode].parent;
    if (fibHeapQ[p].degree == 1)
    {
        fibHeapQ[p].child = VECFIB_NULL;
    }
    else
    {
        if (fibHeapQ[p].child == cutFibNode)
            fibHeapQ[p].child = fibHeapQ[cutFibNode].leftSibling;
        leaveDoubly(fibHeapQ[cutFibNode].leftSibling, cutFibNode, fibHeapQ[p].degree);
    }
    fibHeapQ[p].degree--;
    joinDoubly(minKeyPtr, cutFibNode);
    fibHeapQ[cutFibNode].parent = VECFIB_NULL;  //Cut node is at root level now.
    numMinTrees++;
    return;
}

void fibHeap::cascadeCut(int fibNode)
{
    while ( (fibHeapQ[fibNode].parent != VECFIB_NULL) && (fibHeapQ[fibNode].childCut) )     //fibNode is not at root and it's childcut is true.
    {
        int p = fibHeapQ[fibNode].parent;
        cut(fibNode);
        fibNode = p;
    }
    fibHeapQ[fibNode].childCut = true;
}

int fibHeap::removeFromTreeTableList(int d)
{
    if ( (get<1>(treeTable[d]) == VECFIB_NULL) && (get<2>(treeTable[d]) == VECFIB_NULL) )       //this was the only node.
        startTreeList = -1;
    else if (get<1>(treeTable[d]) == -1)        //start needs to be removed.
    {
        startTreeList = get<2>(treeTable[d]);               //start now points to right of start.
        get<1>(treeTable[startTreeList]) = VECFIB_NULL;     //lPtr of the start is null.
    }
    else
    {
        int lPtr = get<1>(treeTable[d]); int rPtr = get<2>(treeTable[d]);
        get<2>(treeTable[lPtr]) = get<2>(treeTable[d]);     //right of previous is right of d.
        get<1>(treeTable[rPtr]) = get<1>(treeTable[d]);     //left of next is left of d.
    }
    treeTable[d] = make_tuple(VECFIB_NULL, VECFIB_NULL, VECFIB_NULL);         //node at d got combined. So, clear that position in treetable.
    return 0;
}

void fibHeap::addToTreeList(int d)
{
    if (startTreeList != -1)
    {
        get<1>(treeTable[startTreeList]) = d;       //left of start is new.
        get<2>(treeTable[d]) = startTreeList;       //right of new is start.
        get<1>(treeTable[d]) = -1;                  //left of new is null.
    }
    startTreeList = d;                              //start is now new.
}


//Go through both lists.
//make use of treeTable to add elements.
//Keep track of the minKeyPtr
//Keep track of numMinTrees;
void fibHeap::pairwiseCombine(int list1, int list2, int numItemsList1, int numItemsList2)
{
    int currentMaxDegree = 0;
    int currentList1Node = list1, currentList2Node = list2;
    minKeyPtr = list1;
    for (int i=0; i < numItemsList1; i++)
    {
        int x = currentList1Node;
        int nextList1Node = fibHeapQ[currentList1Node].rtSibling;
        int d = fibHeapQ[currentList1Node].degree;
        while (get<0>(treeTable[d]) != VECFIB_NULL)
        {
            x = combine(get<0>(treeTable[d]), x);
            //remove treeTable[d] from list.
            removeFromTreeTableList(d);
            treeTable[d] = make_tuple(VECFIB_NULL, VECFIB_NULL, VECFIB_NULL);         //node at d got combined. So, clear that position in treetable.
            d = fibHeapQ[x].degree;
        }
        treeTable[d] = make_tuple(x, VECFIB_NULL, VECFIB_NULL);
        addToTreeList(d);
//        updateMinKey(x);
        currentList1Node = nextList1Node;
    }
    for (int i=0; i < numItemsList2; i++)
    {
        int x = currentList2Node;
        int nextList2Node = fibHeapQ[currentList2Node].rtSibling;
        int d = fibHeapQ[currentList2Node].degree;
        while (get<0>(treeTable[d]) != VECFIB_NULL)
        {
            x = combine(get<0>(treeTable[d]), x);
            removeFromTreeTableList(d);
            d = fibHeapQ[x].degree;
        }
        treeTable[d] = make_tuple(x, VECFIB_NULL, VECFIB_NULL);
        addToTreeList(d);
//        updateMinKey(x);
        currentList2Node = nextList2Node;
    }
    collectItemsInTreetable(currentMaxDegree);      //This will update numMinTrees and minKeyPtr.
}

void fibHeap::pairwiseCombine_s(int list1, int list2, int numItemsList1, int numItemsList2)
{
    int currentMaxDegree = 0;
    int currentList1Node = list1, currentList2Node = list2;
    minKeyPtr = list1;
    if (numNodes > 0)
        currentMaxDegree = (int)(log(numNodes)/log(phi)) + 1;
    else
    {
        numMinTrees = 0; minKeyPtr = -1;
        return;
    }
    for (int i=0; i < numItemsList1; i++)
    {
        int x = currentList1Node;
        int nextList1Node = fibHeapQ[currentList1Node].rtSibling;
        int d = fibHeapQ[currentList1Node].degree;
        while (get<0>(treeTable[d]) != VECFIB_NULL)
        {
            x = combine(get<0>(treeTable[d]), x);
            treeTable[d] = make_tuple(VECFIB_NULL, VECFIB_NULL, VECFIB_NULL);         //node at d got combined. So, clear that position in treetable.
            d = fibHeapQ[x].degree;
        }
        treeTable[d] = make_tuple(x, VECFIB_NULL, VECFIB_NULL);
        currentList1Node = nextList1Node;
    }
    for (int i=0; i < numItemsList2; i++)
    {
        int x = currentList2Node;
        int nextList2Node = fibHeapQ[currentList2Node].rtSibling;
        int d = fibHeapQ[currentList2Node].degree;
        while (get<0>(treeTable[d]) != VECFIB_NULL)
        {
            x = combine(get<0>(treeTable[d]), x);
            treeTable[d] = make_tuple(VECFIB_NULL, VECFIB_NULL, VECFIB_NULL);         //node at d got combined. So, clear that position in treetable.
            d = fibHeapQ[x].degree;
        }
        treeTable[d] = make_tuple(x, VECFIB_NULL, VECFIB_NULL);
        currentList2Node = nextList2Node;
    }
    collectItemsInTreetable_s(currentMaxDegree);      //This will update numMinTrees and minKeyPtr.
}


int fibHeap::combine(int node1, int node2)
{
    int lNode, gtNode;
    if (fibHeapQ[node1].key < fibHeapQ[node2].key)
    {
        lNode = node1; gtNode = node2;
    }
    else
    {
        lNode = node2; gtNode = node1;
    }
    if (fibHeapQ[lNode].degree == 0)
    {
        fibHeapQ[lNode].child = gtNode;
        fibHeapQ[gtNode].leftSibling = gtNode; fibHeapQ[gtNode].rtSibling = gtNode;
    }
    else
    {
        joinDoubly(fibHeapQ[lNode].child, gtNode);
    }
    fibHeapQ[gtNode].parent = lNode;
    fibHeapQ[lNode].degree += 1;
    fibHeapQ[gtNode].childCut = false;
    return lNode;
}

void fibHeap::collectItemsInTreetable_s(int maxDegree)
{
    numMinTrees = 0; minKeyPtr = -1;
    for (int i=0; i < maxDegree; i++)
    {
        if (get<0>(treeTable[i]) != VECFIB_NULL)
        {
            numMinTrees++;
            int itemToAdd = get<0>(treeTable[i]);
            if (minKeyPtr == -1)
            {
                fibHeapQ[itemToAdd].leftSibling = itemToAdd; fibHeapQ[itemToAdd].rtSibling = itemToAdd;
                minKeyPtr = itemToAdd;
            }
            else
            {
                joinDoubly(minKeyPtr, itemToAdd);
                updateMinKey(itemToAdd);
            }
            treeTable[i] = make_tuple(VECFIB_NULL, VECFIB_NULL, VECFIB_NULL);     //Reset the entry in treeTable after collecting it.
        }
    }
    return;
}

void fibHeap::collectItemsInTreetable(int maxDegree)
{
    int topList = -1; numMinTrees = 0; minKeyPtr = -1;
    int treeList = startTreeList, nextTreeListItem, itemToAdd;
    if (startTreeList == -1)
        return;
    while (get<2>(treeTable[treeList]) != -1)
    {
        numMinTrees++;
        itemToAdd = get<0>(treeTable[treeList]);
        nextTreeListItem = get<2>(treeTable[treeList]);
        if (topList == -1)
        {
            topList = itemToAdd;
            fibHeapQ[topList].leftSibling = topList; fibHeapQ[topList].rtSibling = topList;
            minKeyPtr = topList;
        }
        else
        {
            joinDoubly(topList, itemToAdd);
            updateMinKey(itemToAdd);
        }
        treeTable[treeList] = make_tuple(VECFIB_NULL, VECFIB_NULL, VECFIB_NULL);     //Reset the entry in treeTable after collecting it.
        treeList = nextTreeListItem;
    }
    //Add the last one.
    numMinTrees++;
    itemToAdd = get<0>(treeTable[treeList]);
    if (topList == -1)
    {
        topList = itemToAdd;
        fibHeapQ[topList].leftSibling = topList; fibHeapQ[topList].rtSibling = topList;
        minKeyPtr = topList;
    }
    else
    {
        joinDoubly(topList, itemToAdd);
        updateMinKey(itemToAdd);
    }
    treeTable[treeList] = make_tuple(VECFIB_NULL, VECFIB_NULL, VECFIB_NULL);     //Reset the entry in treeTable after collecting it.
    startTreeList = -1;
    return;
}


void fibHeap::printHeap()
{
    int node = minKeyPtr;
    for (int i = 0; i < numMinTrees; i++)
    {
        printTree(node);
        node = fibHeapQ[node].rtSibling;
    }
}
void fibHeap::printTree(int node)
{
    queue<int> q;
    q.push(node);
    cout << "\nMinTree with root: " << node << "\n";
    int p = node;
    int c = fibHeapQ[p].child;
    while (!q.empty())
    {
        p = q.front(); q.pop();
        cout << p << " child of: " << fibHeapQ[p].parent << "; ";
        c = fibHeapQ[p].child;
        for (int i = 0; i < fibHeapQ[p].degree; i++)
        {
            q.push(c); c = fibHeapQ[c].rtSibling;
        }
    }
}
