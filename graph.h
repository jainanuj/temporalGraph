#ifndef GRAPH_H_HHWU
#define GRAPH_H_HHWU

#include <cstdio>
#include <vector>
#include <list>
#include <tuple>
#include <cstring>
#include <algorithm>
#include <set>
#include <queue>
#include <iostream>
#include "Timer.h"

#include "bitfield.hpp"

using namespace std;

const int infinity = 2e9;

struct Interval
{
    int intvlStart;
    int intvlEnd;
    int adjustedEnd;
    int traveTime;
    int prevJourneyIndex;
    int divTime;
public:
    int isValidFor(int t)
    {
        if (adjustedEnd < intvlStart)       //Make sure interval itself is not invalid.
            return 0;
        return (t > adjustedEnd) ? 0 : 1;
    }
};

struct Nbrs
{
public:
    int nbrId;
    int nbr_lambda;
    int numIntvls;
    int adjustedNumIntvls;
    vector<Interval> edgeSchedules;
};

struct Journey
{
public:
    vector<tuple<int, int, int>> rPath;      //Pair consists of node Id & the index of edge of that node used in the journey from s to a node.
    vector<tuple<int, int, int>> sigmaSchedule;   //pair consists of the schedule/time of travel and the interval ID used.
};

struct incrementalJourney
{
    int indexPrevIncrement;
    int prevNodeID;
    int prevEdgeID;
    int prevIntvlID;
    int currentNodeID;
    int prevDepartureTime;
    int currentArrivalTime;
public:
    incrementalJourney(int cNId, int cArrT, int pNId, int pEdgeId, int pIntvlId, int pDep, int pIndex)
    {
        currentNodeID = cNId; currentArrivalTime = cArrT; prevNodeID = pNId; prevEdgeID = pEdgeId; prevIntvlID = pIntvlId; prevDepartureTime = pDep; indexPrevIncrement = pIndex;
    }
    incrementalJourney() {}
    void initializeJourneyIncrement(int cNId, int cArrT, int pNId, int pEdgeId, int pIntvlId, int pDep, int pIndex)
       {
           currentNodeID = cNId; currentArrivalTime = cArrT; prevNodeID = pNId; prevEdgeID = pEdgeId; prevIntvlID = pIntvlId; prevDepartureTime = pDep; indexPrevIncrement = pIndex;
       }
};

struct Node
{
    int nodeId;
    int numNbrs;
    vector<Nbrs> neighbors;
    int inDegree;
};

class compareMFSetElements {  //These elements are tuples (arrival time, wait time, prev node). Last int is not part of comparison.
public:
    bool operator () (const tuple<int, int, int, int, int>& tup1, const tuple<int, int, int, int, int>& tup2) const
    {
        if (get<0>(tup1) != get<0>(tup2))
            return (get<0>(tup1) < get<0>(tup2));     //Arrival time
        else if (get<1>(tup1) != get<1>(tup2))
            return (get<1>(tup1) < get<1>(tup2));     //Arrival time     //Wait time.
        else
            return false;                                       //prev node is used just for tracing path back.
    }
};

class compareMFSetElementsPrioritized {  //These elements are tuples (arrival time, wait time, prev node). Last int is not part of comparison.
public:
    bool operator () (const tuple<int, int, int, int, int, bool>& tup1, const tuple<int, int, int, int, int, bool>& tup2) const
    {
        if (get<0>(tup1) != get<0>(tup2))
            return (get<0>(tup1) < get<0>(tup2));     //Arrival time
        else if (get<1>(tup1) != get<1>(tup2))
            return (get<1>(tup1) < get<1>(tup2));     //Wait time.
        else
            return false;                                       //prev node is used just for tracing path back.
    }
};

class Graph
{
public:
    Graph() {}
    Graph(const char* filePath, int contactSeq, const char * option); // input file
    void initial_query(const char* filePath); // query file
    void initial_query();
    void initial_ds_ea();
    void initial_ds_ld();
    void initial_ds_f();
    void initial_ds_s();
    void run_earliest_arrival();
    void run_shortest();
    virtual void run_mwf();
    void earliest_arrival(int source);
    tuple<int,int> earliest_arrival_pair(int source,int retRchd=0);
    void earliest_arrival_fibo(int source);
    void earliest_arrival_fibo_external(int source);
    void shortest_path(int source);
    void shortest_path_xuan(int source);
    int edgeAndScheduleSel(vector<std::tuple<int, int, int>>& e_min, vector<int>& t_min, vector<int>& t_LBD);
    void minWaitForemost(int source);
    void extendMWFJourneys(int exploreNode);
    void extendJourney(tuple<int,int, int> journey, int node);
    bool checkDominatedAndInsert(int destNode,tuple<int, int, int, int, int> newJourney);
    void insertInJourneySets(int destNode, tuple<int, int, int, int, int> newJourney, set <tuple<int, int, int, int, int>, compareMFSetElements >::iterator itInsertPos);
    bool checkDominance(tuple<int, int> journey1,  tuple<int, int> journey2);

    
    void minWaitForemostPrioritized(int source);
    void extendPrioritizedJourney(tuple<int,int, int,int> inJourney, int node, int source);
    bool checkDominatedAndInsertPrioritized(int destNode,tuple<int, int, int, int, int, bool> newJourney);
    void insertPrioritizedInJourneySets(int destNode, tuple<int,int,int,int, int,bool> newJourney, set <tuple<int, int, int,int, int, bool>, compareMFSetElementsPrioritized >::iterator itInsertPos);
    
    void minWaitForemostPrioritizedNoSet(int source);
    void printMWFWalksPrioritizedNoSet(int source);
    
    void build_shortestJourneys(int source, vector<std::pair<int, int>>& shortestJourneyPointer, vector<vector<incrementalJourney>>& allHopJourneys);
    void print_shortest_paths(int source);
    void print_shortest_results_test(int source);

    
    int earliestUseEdgeAfterT(int u, Nbrs& v, int t, int &intvlID);
    
	void print_avg_time();
    
    static void wuGraph(const char* filePath,int noL = 0, int numDrop=1, int normalizeWrite=0);     //input file in wu Format.
    static void collapseIntervalsWriteOuput(const char* filePath, const char* opFile);
    static vector<tuple<int, int, int>> adjustSlowIntervals(std::vector<tuple<int, int, int>>& intervalVector);
    static void buildXuanGraph(const char* filePath, vector<tuple<int, int, int, int>>& inputRows, int vertices, int wuEdges);
    static void readWuFile(const char* filePath);
    
    // for testing the correctness
    void run_earliest_arrival(const char* filePath); // output the result
    void earliest_arrival(int source, FILE * file);  
    void print_result(const int source, const vector<int>& t_time, FILE * file);
    void print_result_ld(const int source, const vector<int>& t_time, FILE * file); 
    void print_avg_time(const char* filePath1, const char* filePath2);
    
    void printResults(int source);
    void printEarliestResultsTest(int source);
    void printMWFWalks(int source);
    void printMWFWalksPrioritized(int source);
    
    virtual void run_mhf() {}


public:
    string earliestResults;
    string shortestResults;
    string minHeapMonitor;
    string mhfResults;
    int mNumJourExtInst;
    vector<Node> vertices;
    
    
    vector< int > sources;
    int V, static_E, dynamic_E;
    int t_start, t_end;
    double time_sum;
    vector <int> arr_time, f_time;
    vector <tuple<int, int, int>> father;
    vector<Journey> shortestJourneys;
    vector<tuple<
                set<tuple<int, int, int>>,    //set indicating new journeys.
                set <tuple<int, int, int, int, int>, compareMFSetElements >  //set of (arr_tm,wt_time). 3rd is prev node, dep from prev node.,num hops in this journey reqd to trce back
                >
          > vecFullList;
    list<int> toExpandList;     //List of all nodes that received new incoming journeys.
    vector<
            set <tuple<int, int, int, int, int, bool>, compareMFSetElementsPrioritized >  //set of (arr_tm,wt_time). 3rd is prev node, dep from prev node.,num hops in this journey, is journey dominated.
          > vecFullListPriority;
    vector<tuple<int, int, int, int>> minHeap;    //arrTime,wtTime,numHops,nodeId.
    vector<tuple<int, int, int, int>> vecTuple;    //arrTime,wtTime,numHops,nodeId.
    vector<tuple<int, int, int, int>> mwfJourneys;    //arrTime,wtTime,numHops,nodeId.
};

class nodeComparison {
private:
    Graph *graphToCompare;
public:
    nodeComparison(Graph *nodesGraph)
    {
        graphToCompare = nodesGraph;
    }
    bool operator () (int node1, int node2)
    {
        return (graphToCompare->arr_time[node1] > graphToCompare->arr_time[node2]);
    }
};

class nodeComparison2 {

public:
    bool operator () (std::pair<int, int> node1, std::pair<int, int> node2)
    {
        return (node1.first > node2.first);
    }
};

class nodeCompareArrWt {
public:
    bool operator () (std::tuple<int, int, int,int> node1, std::tuple<int, int, int,int> node2)
    {
        if (std::get<0>(node1) != std::get<0>(node2))
            return (std::get<0>(node1) > std::get<0>(node2));       //Compare ArrTimes
        else if (std::get<1>(node1) != std::get<1>(node2))          //If ArrTimes are same, compare Wt times.
            return (std::get<1>(node1) > std::get<1>(node2));
        else
            return (std::get<3>(node1) > std::get<3>(node2));       //If arr & wt times are same, should be diff. node. If not its duplicate
    }
};

class compareTuple {
public:
    bool operator () (std::tuple<int, int, int, int> tup1, std::tuple<int, int, int, int> tup2)
    {
        if (std::get<0>(tup1) != std::get<0>(tup2))
            return (std::get<0>(tup1) < std::get<0>(tup2));
        else if (std::get<1>(tup1) != std::get<1>(tup2))
            return (std::get<1>(tup1) < std::get<1>(tup2));
        else if (std::get<2>(tup1) != std::get<2>(tup2))
            return (std::get<2>(tup1) < std::get<2>(tup2));
        else
            return (std::get<3>(tup1) < std::get<3>(tup2));
    }
};


class compareTupleWu {      //sort by t.
public:
    bool operator () (std::tuple<int, int, int, int> tup1, std::tuple<int, int, int, int> tup2)
    {
        if (get<2>(tup1) != get<2>(tup2))           //third element in the tuple is 't' (time instant)
            return (get<2>(tup1) < get<2>(tup2));
        else if (get<0>(tup1) != get<0>(tup2))
            return (get<0>(tup1) < get<0>(tup2));
        else if (get<1>(tup1) != get<1>(tup2))
            return (get<1>(tup1) < get<1>(tup2));
        else
            return (get<3>(tup1) < get<3>(tup2));
    }
};



#endif

