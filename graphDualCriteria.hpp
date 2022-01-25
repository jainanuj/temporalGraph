//
//  graphDualCriteria.hpp
//  XuantemporalGraph
//
//  Created by Anuj Jain on 12/31/21.
//

#ifndef graphDualCriteria_hpp
#define graphDualCriteria_hpp

#include <stdio.h>
#include "graph.h"

class GraphDualCriteria  : public Graph
{
public:
    GraphDualCriteria(const char* filePath, int contactSeq) :  Graph(filePath, contactSeq) {} // input file
    bool lessCompArrHop(std::pair<int, int> newArrHop, std::pair<int,int> oldArrHop);
    void initial_ds_eha();
    int earliest_arrival_minHop_pair(int source);
    void run_mhf();
    void printmhfResultsTest(int source);
    void mhfHopByHop(int source);
    void build_mhf_Journeys(int source, vector<std::tuple<int, int, int>>& mhfJourneyPointer, vector<vector<incrementalJourney>>& allHopJourneys);
    void printmhfResultsTest2(int source, vector<std::tuple<int, int, int>>& mhfJourney);
public:
    vector <pair<int,int>> arr_hop_time, f_time;
};

class nodeComparisonArrHops {
public:
    bool operator () (std::tuple<int, int, int> node1, std::tuple<int, int, int> node2)     //(arr,hops,nodeId)
    {
        if (get<0>(node1) != get<0>(node2))
            return (get<0>(node1) > get<0>(node2));
        else
            return (get<1>(node1) > get<1>(node2));
    }
};


#endif /* graphDualCriteria_hpp */
