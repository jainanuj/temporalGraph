//
//  graphDualCriteria.cpp
//  XuantemporalGraph
//
//  Created by Anuj Jain on 12/31/21.
//

#include "graphDualCriteria.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
//#include "cppfibonacci/fibonacci.hpp"
#include <fstream>
#include "fibheap.h"
#include "regHeap.hpp"
#ifndef __TEST__
#define __TEST__
#endif

int GraphDualCriteria::earliest_arrival_minHop_pair(int source)
{
    ofstream minHeapMonitorOutput(minHeapMonitor);
    int numNodesReached=0;
    Timer t;
#ifdef MEASUREHEAP_DET
    unsigned long maxHeapSize = 0, avgHeapSize = 0, count = 0, currentAvg = 0;
    clock_t ticks;
    clock_t insertTimer = 0; //clock() - clock();
    clock_t remMinTimer = 0; //clock() - clock();
#endif
    Node *u;
    int nodeID;
    int tDepart, intvlID = -1;
    vector<tuple<int, int, int>> minHeap;    //Fix the pq to make it a min heap instead of max heap.
    nodeComparisonArrHops heapCompFn;

    bit_queue closedNodes((int)vertices.size());
        
    arr_hop_time[source]=make_pair(t_start,0);
    father[source] = make_tuple(source,-1,0);
    minHeap.push_back(make_tuple(arr_hop_time[source].first, arr_hop_time[source].second, source));
    int i = 0, numRepeatedNodes = 0, numInserts = 0;
    
    t.start();

    while (!minHeap.empty())
    {
        i++;
        nodeID = get<2>(minHeap.front());
#ifdef MEASUREHEAP_DET
        ticks = clock();
#endif
        pop_heap<vector<tuple<int, int, int>>::iterator, nodeComparisonArrHops>(minHeap.begin(), minHeap.end(), heapCompFn);
        minHeap.pop_back();      //moves the min element to last and then removes it from heap.
#ifdef MEASUREHEAP_DET
        remMinTimer += (clock() - ticks);
#endif
        if (closedNodes.check_bit_obj_present(nodeID))      //This was already closed.
        {
            numRepeatedNodes++;
            continue;
        }
        closedNodes.queue_add_bit(nodeID);         //Add it to the closedNodes.
        numNodesReached++;
        u = &vertices[nodeID];
        for (int i = 0; i < u->numNbrs; i++)
        {
            if (! closedNodes.check_bit_obj_present( u->neighbors[i].nbrId))
            {
                tDepart = earliestUseEdgeAfterT(nodeID, u->neighbors[i], arr_hop_time[nodeID].first, intvlID );
                if (intvlID >= 0)
                {
                    int newArrTime=tDepart + u->neighbors[i].edgeSchedules[intvlID].traveTime;
                    int newHopCnt=arr_hop_time[nodeID].second+1;
                    if (lessCompArrHop(make_pair(newArrTime, newHopCnt), arr_hop_time[u->neighbors[i].nbrId]))
                    {
                        arr_hop_time[u->neighbors[i].nbrId].first = newArrTime;
                        arr_hop_time[u->neighbors[i].nbrId].second = newHopCnt;
                        father[u->neighbors[i].nbrId] = make_tuple(u->nodeId, intvlID, tDepart) ;
#ifdef MEASUREHEAP_DET
                        ticks = clock();
#endif
                        minHeap.push_back(make_tuple(newArrTime, newHopCnt, u->neighbors[i].nbrId));
                        push_heap(minHeap.begin(), minHeap.end(), heapCompFn);
#ifdef MEASUREHEAP_DET
                        insertTimer += (clock() - ticks);
#endif
                        numInserts++;
                    }
                }
            }
        }
/*        if (i % 50 == 0)
        {
            minHeapMonitorOutput << minHeap.size() << "; " ;
            if (minHeap.size() > maxHeapSize)
                maxHeapSize = (int)minHeap.size();
            avgHeapSize = (currentAvg*(count) + (unsigned long)minHeap.size() )/(count + 1);
            currentAvg = avgHeapSize; count++;
            minHeapMonitorOutput << avgHeapSize << "\n";
        }*/
    }

    t.stop();
/*    minHeapMonitorOutput << "Max Heap size: " << maxHeapSize << "\n";
    minHeapMonitorOutput << "Avg Heap size: " << avgHeapSize << "\n";
    minHeapMonitorOutput << "Num Repeated Nodes: " << numRepeatedNodes << "\n";
    minHeapMonitorOutput << "Num Inserts: " << numInserts << "\n";
#ifdef MEASUREHEAP_DET
    minHeapMonitorOutput << "Insert Time: " << insertTimer << "\n";
    minHeapMonitorOutput << "Rem Min Time: " << remMinTimer << "\n";
#endif*/
//    printResults(source);
    cout << numNodesReached << endl;
    time_sum += t.GetRuntime();
#ifdef __TEST__
    printmhfResultsTest(source);
#endif
    return numNodesReached;
}

bool GraphDualCriteria::lessCompArrHop(std::pair<int, int> newArrHop, std::pair<int,int> oldArrHop)
{
    if (get<0>(newArrHop) != get<0>(oldArrHop))
        return (get<0>(newArrHop) < get<0>(oldArrHop));
    else
        return (get<1>(newArrHop) < get<1>(oldArrHop));
}

void GraphDualCriteria::initial_ds_eha()
{
    arr_hop_time.resize(V);
    father.resize(V);
    for(int i=0; i<V; i++){
        arr_hop_time[i]= make_pair(infinity,infinity);
        father[i] = make_tuple(i,-1,0);
    }

}

void GraphDualCriteria::run_mhf()
{
    time_sum=0;
    
    for(int i = 0 ;i < sources.size() ;i ++)
    {
        initial_ds_eha();
        //earliest_arrival_minHop_pair(sources[i]);
        mhfHopByHop(sources[i]);
    }
    
    print_avg_time();
}

void GraphDualCriteria::printmhfResultsTest(int source)
{
    int rv = 0;
    ofstream earliestMHOut(mhfResults);
    earliestMHOut << V << "\n";
    for (int i = 0; i < arr_hop_time.size(); i++)
    {
//        if (arr_hop_time[i].first >= infinity)
//            continue;
        earliestMHOut << i << " " << arr_hop_time[i].first << "  "  <<  arr_hop_time[i].second <<"\n";
        if (arr_hop_time[i].first < infinity)
            rv++;
/*#ifdef __TEST__
        int x = 0;
        x = i;
        while (get<0>(father[x]) != x)
        {
            earliestOut << x <<"(" << get<1>(father[x]) << " " << get<2>(father[x]) << ") ";
            x = get<0>(father[x]);
        }
        earliestOut << x << "\n";
#endif*/
    }
    cout << "Source: " << source << "\n";
    cout << "Num reachable vertices: " << rv << "\n";
    cout << "Ratio of reachable vertices " << (float)rv/V << "\n";
}

void GraphDualCriteria::mhfHopByHop(int source)
{
    Timer t;

    vector<vector<incrementalJourney>> allHopJourneys; //At each hop there is a vector of foremost incremental journeys, discovered in that hop.
    vector<std::tuple<int, int, int>> earliestKnownTimeArrival;              // <earliestArrivalTime, hop in which this time achieved, index in allHopJourneys[thishop] vector.>

    allHopJourneys.resize(V);   //Max number of hops is the number of nodes.
    earliestKnownTimeArrival.resize(V);
    earliestKnownTimeArrival.assign(V, std::make_tuple(infinity, 0, -1));
    incrementalJourney journeyIncrement(source, 0, -1, -1, -1, -1, -1);   //current node, arrTime, prevNode, prevEdgeId, prevIntvlId, pDep, indexInPrevHopV

    t.start();
    int intvlID = -1, departTime, newArrivTime, hopCount = 0, nextNode;
    earliestKnownTimeArrival[source] = std::make_tuple(0, 0, 0);
    allHopJourneys[hopCount].push_back(journeyIncrement);

    int  numNodesSeen = 1, numNewNodesInCurrentHop = 1;     //number of new foremost nodes in current hop.
    while ( (hopCount < (V-1)) && (numNewNodesInCurrentHop > 0))// && (numNodesSeen < V))
    {
        hopCount++;
        numNewNodesInCurrentHop = 0;    //number of foremost nodes
        for (int i = 0; i < allHopJourneys[hopCount -1].size(); i++)        //extend all journeys from prev hop
        {
            incrementalJourney* incToExtend = &(allHopJourneys[hopCount-1][i]);
            int node = incToExtend->currentNodeID;
            int time = incToExtend->currentArrivalTime;
            
            for (int j = 0; j < vertices[node].numNbrs; j++)
            {
                departTime = earliestUseEdgeAfterT(node, vertices[node].neighbors[j], time, intvlID);
                if ( (departTime == -1) || (departTime >= infinity) || (intvlID == -1) )
                    continue;
                nextNode = vertices[node].neighbors[j].nbrId;
                newArrivTime = departTime + vertices[node].neighbors[j].edgeSchedules[intvlID].traveTime;
                if (newArrivTime < std::get<0>(earliestKnownTimeArrival[nextNode]))     //New foremost time to nextNode Discovered.
                {
                    if (std::get<1>(earliestKnownTimeArrival[nextNode]) < hopCount)     //This means the node is seen for first time in this hop with the lowest time.
                    {
                        incrementalJourney newIncrement(nextNode, newArrivTime, node, j, intvlID, departTime, i);
                        allHopJourneys[hopCount].push_back(newIncrement);
                        std::get<1>(earliestKnownTimeArrival[nextNode]) = hopCount;
                        std::get<2>(earliestKnownTimeArrival[nextNode]) = numNewNodesInCurrentHop++;
                    }
                    else        //It has new lowest time in this hop but it was already seen with a lowest time in this hop.
                    {
                        incrementalJourney *pNewIncrement = &(allHopJourneys[hopCount][std::get<2>(earliestKnownTimeArrival[nextNode])]);
                        pNewIncrement->currentArrivalTime = newArrivTime; pNewIncrement->prevNodeID = node; pNewIncrement->prevEdgeID = j;
                        pNewIncrement->prevIntvlID = intvlID; pNewIncrement->prevDepartureTime = departTime; pNewIncrement->indexPrevIncrement = i;
                    }
                    if (std::get<0>(earliestKnownTimeArrival[nextNode]) >= infinity)      //This is shortest journey to nextNode as no path found yet.
                        numNodesSeen++;
                    std::get<0>(earliestKnownTimeArrival[nextNode]) = newArrivTime;
                }       //found a better arrival time.
            }       //Looking at all the nbrs of a node to extend in current hop count.
        }       //for loop extending journey from hopCount-1.
    }       //While loop to cover all hop counts.

    t.stop();
    time_sum += t.GetRuntime();
    cout << "Num Nodes seen: " << numNodesSeen << endl;

#ifdef __TEST__
//    build_mhf_Journeys(source, earliestKnownTimeArrival, allHopJourneys);
    printmhfResultsTest2(source, earliestKnownTimeArrival);
#endif

//    print_shortest_paths(source);
}


void GraphDualCriteria::build_mhf_Journeys(int source, vector<std::tuple<int, int, int>>& mhfJourneyPointer, vector<vector<incrementalJourney>>& allHopJourneys)
{
    shortestJourneys.resize(V);                     //Shortest Journeys for all nodes needs to be found.
    shortestJourneys[source].rPath.push_back(std::make_tuple(source, -1, -1));       //start node, nbr ID.
    shortestJourneys[source].sigmaSchedule.push_back(std::make_tuple(0, -1, -1));        //depart time, IntvlId
    for (int node = 0; node < V; node++)
    {
        if (node == source)
            continue;
        int hopCount = std::get<1>(mhfJourneyPointer[node]);
        int incrIndex = std::get<2>(mhfJourneyPointer[node]);
        for (int mhfHopCount = hopCount; mhfHopCount > 0; mhfHopCount--)  //BuildShortestJourney for nextNode from Incremental Vectors.
        {
            incrementalJourney *pJourneyPath = &(allHopJourneys[mhfHopCount][incrIndex]);    //start node, it's nbr id.
            int currentNode = pJourneyPath->currentNodeID;int startNode = pJourneyPath->prevNodeID; int edgeId = pJourneyPath->prevEdgeID;
            int currentArrivalTime = pJourneyPath->currentArrivalTime;int depart = pJourneyPath->prevDepartureTime; int intvlUsed = pJourneyPath->prevIntvlID;
            shortestJourneys[node].rPath.push_back(std::make_tuple(currentNode, startNode, edgeId));
            shortestJourneys[node].sigmaSchedule.push_back(std::make_tuple(currentArrivalTime, depart, intvlUsed));    //dep. time on edge (node, nextnode), intvl used.
            incrIndex = pJourneyPath->indexPrevIncrement;
        }
    }
}

void GraphDualCriteria::printmhfResultsTest2(int source, vector<std::tuple<int, int, int>>& mhfJourney)
{
    int rv = 0; int sumHops = 0, maxHopCount=0; int avgHops;
    ofstream earliestMHFOut(mhfResults);
    earliestMHFOut << V << "\n";
    for (int i = 0; i < mhfJourney.size(); i++)
    {
//        if (get<0>(mhfJourney[i]) >= infinity)
//            continue;
        earliestMHFOut << i << " " << get<0>(mhfJourney[i]) << "  "  <<  get<1>(mhfJourney[i]) <<"\n";
        if (get<0>(mhfJourney[i]) < infinity)
        {
            rv++;
            sumHops += get<1>(mhfJourney[i]);
            if (get<1>(mhfJourney[i]) > maxHopCount)
                maxHopCount = get<1>(mhfJourney[i]);
        }
/*#ifdef __TEST__
        int x = 0;
        x = i;
        while (get<0>(father[x]) != x)
        {
            earliestOut << x <<"(" << get<1>(father[x]) << " " << get<2>(father[x]) << ") ";
            x = get<0>(father[x]);
        }
        earliestOut << x << "\n";
#endif*/
    }
    cout << "Source: " << source << "\n";
    cout << "Num reachable vertices: " << rv << "\n";
    cout << "Ratio of reachable vertices " << (float)rv/V << "\n";
    avgHops = sumHops/rv;
    cout << "Max Hops = " << maxHopCount << "\n";
    cout << "Avg # hops = " << avgHops << "\n";
}
