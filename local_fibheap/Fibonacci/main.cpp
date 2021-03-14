//
//  main.cpp
//  Fibonacci
//
//  Created by Anuj Jain on 2/26/21.
//

#include <iostream>
#include "fibheap.h"
#define VERTS 20

int main(int argc, const char * argv[]) {
    fibHeap fibMinHeap(VERTS);
    for (int i = 0; i < VERTS-1; i++)
    {
        int r = rand() % 100000;
        fibMinHeap.insert(i, r);
        cout << "(" << i << "," << r << ")   ";
    }
    
//    pair<int, int> minNK = fibMinHeap.removeMin();
//    cout << "Min node: " << get<0>(minNK) << " Val: " << get<1>(minNK) << "\n";
    fibMinHeap.decrease_key(5, 2);
    pair<int, int> minNK = fibMinHeap.getMin();
    cout << "Min node: " << get<0>(minNK) << " Val: " << get<1>(minNK) << "\n";
    minNK = fibMinHeap.removeMin();
    cout << "Min node: " << get<0>(minNK) << " Val: " << get<1>(minNK) << "\n";
    fibMinHeap.printHeap();
    int i = 0;
    while (!fibMinHeap.empty())
    {
        pair<int, int> minNK = fibMinHeap.removeMin();
        cout << i++ << ": Min node: " << get<0>(minNK) << " Val: " << get<1>(minNK) << "\n";
    }
    std::cout << "Hello, World!\n";
    return 0;
}
