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

struct intervalInfo {
    int u;
    int v;
    int nbrIndexFor_v;
    int intvlStart;
    int intvlEnd;
    int lambda;
    int intvlId;            //Index of this intvl in original graph defn. for PreKnownIntvls
    int prevJourneyIndex;   //Index of last journey before intvl in the list of journeys arriving at u
};

struct mwfJourney {
    int arrivalTime;
    int wtTime;
    int prevNode;
    int prevJourneyIndex;
    int prevDepTime;
    vector<tuple<int, int, int>> expandedAt;     //(intvl.startTime,lambda) of intvl in which this journey was expanded on each nbr. Last element is flag whether it was expanded or not.
};

class GraphDualCriteria  : public Graph
{
public:
    GraphDualCriteria(const char* filePath, int contactSeq, const char* option);
    bool lessCompArrHop(std::pair<int, int> newArrHop, std::pair<int,int> oldArrHop);
    void initial_ds_eha();
    int earliest_arrival_minHop_pair(int source);
    void run_mhf();
    void run_mwf();
    void printmhfResultsTest(int source);
    void mhfHopByHop(int source);
    void mwfStreamingIntvls(int source);
    void build_mhf_Journeys(int source, vector<std::tuple<int, int, int>>& mhfJourneyPointer, vector<vector<incrementalJourney>>& allHopJourneys);
    void printmhfResultsTest2(int source, vector<std::tuple<int, int, int>>& mhfJourney);
    void printmwfResultsTest2(int source);
    int getMinIntvl(intervalInfo &newIntvl, int indexPreKnownIntvls);
    int searchPrevJourney(int node, int beforeTime);
    void setupNewJourney(int v, int arrivalTime);
public:
    vector <pair<int,int>> arr_hop_time, f_time;

    vector<vector<mwfJourney>> listJourneys;    //List of journeys at each node.
    vector<intervalInfo> listOfPreKnownIntvls;
    vector<intervalInfo> listOfAdHocIntvls;
    vector<mwfJourney> finalMWFJourneys;
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

class compareIntvlsMWF {      //sort by t.
public:
    bool operator () (intervalInfo intvl1, intervalInfo intvl2)
    {
        int arrIntvl1=intvl1.intvlStart+intvl1.lambda;
        int arrIntvl2=intvl2.intvlStart+intvl2.lambda;
        if ( arrIntvl1 != arrIntvl2)
            return (arrIntvl1 < arrIntvl2);
        else
            return (intvl1.intvlStart < intvl2.intvlStart);
    }
};

class compareIntvlsMWFMinHeap {      //sort by t.
public:
    bool operator () (intervalInfo intvl1, intervalInfo intvl2)
    {
        int arrIntvl1=intvl1.intvlStart+intvl1.lambda;
        int arrIntvl2=intvl2.intvlStart+intvl2.lambda;
        if ( arrIntvl1 != arrIntvl2)
            return (arrIntvl1 > arrIntvl2);
        else
            return (intvl1.intvlStart > intvl2.intvlStart);
    }
};



#endif /* graphDualCriteria_hpp */
