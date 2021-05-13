#ifndef GRAPH_H_HHWU
#define GRAPH_H_HHWU

#include <cstdio>
#include <vector>
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
};

class Graph
{
public:
    Graph() {}
    Graph(const char* filePath, int contactSeq); // input file
    void initial_query(const char* filePath); // query file
    void initial_query();
    void initial_ds_ea();
    void initial_ds_ld();
    void initial_ds_f();
    void initial_ds_s();
    void run_earliest_arrival();
    void run_shortest();
    void earliest_arrival(int source);
    void earliest_arrival_pair(int source);
    void earliest_arrival_fibo(int source);
    void earliest_arrival_fibo_external(int source);
    void shortest_path(int source);
    void shortest_path_xuan(int source);
    void edgeAndScheduleSel(vector<std::tuple<int, int, int>>& e_min, vector<int>& t_min, vector<int>& t_LBD);

    
    void build_shortestJourneys(int source, vector<std::pair<int, int>>& shortestJourneyPointer, vector<vector<incrementalJourney>>& allHopJourneys);
    void print_shortest_paths(int source);
    void print_shortest_results_test(int source);

    
    int earliestUseEdgeAfterT(int u, Nbrs& v, int t, int &intvlID);
    
	void print_avg_time();
    
    static void wuGraph(const char* filePath,int noL = 0, int numDrop=1);     //input file in wu Format.
    static void collapseIntervalsWriteOuput(const char* filePath, const char* opFile);
    static vector<tuple<int, int, int>> adjustSlowIntervals(std::vector<tuple<int, int, int>>& intervalVector);
    
    // for testing the correctness
    void run_earliest_arrival(const char* filePath); // output the result
    void earliest_arrival(int source, FILE * file);  
    void print_result(const int source, const vector<int>& t_time, FILE * file);
    void print_result_ld(const int source, const vector<int>& t_time, FILE * file); 
    void print_avg_time(const char* filePath1, const char* filePath2);
    
    void printResults(int source);
    void printEarliestResultsTest(int source);


public:
    string earliestResults;
    string shortestResults;
    string minHeapMonitor;
    vector<Node> vertices;
    
    vector< int > sources;
    int V, static_E, dynamic_E;
    int t_start, t_end;
    double time_sum;
    vector <int> arr_time, f_time;
    vector <tuple<int, int, int>> father;
    vector<Journey> shortestJourneys;
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

