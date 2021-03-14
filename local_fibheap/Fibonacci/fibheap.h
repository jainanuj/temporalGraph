//
//  fibheap.h
//  Fibonacci
//
//  Created by Anuj Jain on 2/26/21.
//
#ifndef fibheap_h
#define fibheap_h

#include <vector>
#include <tuple>
#include <cstring>
#include <algorithm>
//#include <set>
#include <queue>

using namespace std;
#define VECFIB_NULL -1


class fibNodeElement {
public:
    int leftSibling;
    int rtSibling;
    int parent;
    int child;
    int degree;
    bool childCut;
    bool inHeap;
    int key;
    fibNodeElement(int x=VECFIB_NULL, bool cc=false)
    {
        key = x;
        leftSibling = x; rtSibling = x;
        parent= x; child = x;
        degree = 0; childCut = cc; inHeap = false;
    }
};

class fibHeap {

private:
    vector<fibNodeElement> fibHeapQ;
    vector<tuple<int, int, int>> treeTable;
    int minKeyPtr;
    int numMinTrees;
    int startTreeList;
    
    void cut(int cutFibNode);
    void cascadeCut(int fibNode);
    void joinDoubly(int doubly, int nodeID);
    void leaveDoubly(int leftOfNode, int nodeID, int count);
    void updateMinKey(int fibNode);
    void pairwiseCombine(int list1, int list2, int numItemsList1, int numItemsList2);
    void pairwiseCombine_s(int list1, int list2, int numItemsList1, int numItemsList2);
    int combine(int node1, int node2);
    void collectItemsInTreetable(int maxDegree);
    void collectItemsInTreetable_s(int maxDegree);
    void addToTreeList(int d);
    int removeFromTreeTableList(int d);
    
    
public:
    int numNodes;
    fibHeap(int V=0);
    void insert(int nodeID, int key);
    void decrease_key(int fibNode, int newValue);
    void initialize(int nodeID);
    void meld(int fibTree1, int fibTree2);
    pair<int, int> removeMin();
    pair<int, int> getMin();
    bool partOfHeap(int node);
    void printHeap();
    void printTree(int node);
    bool empty();
};

#endif /* fibheap_h */
