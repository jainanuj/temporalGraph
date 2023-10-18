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

#define __NO_CLASSES__

GraphDualCriteria::GraphDualCriteria(const char* filePath, int contactSeq, const char* option) :  Graph(filePath, contactSeq, option)
{
    struct intervalInfo intvl;
    compareIntvlsMWF intvlCompare;
    if (!strcmp(option, "mwf"))
    {
        //Traverse the whole graph, add intervals to the vector PreKnownIntvls
        for (int i=0; i < V; i++)
        {
            for (int j = 0; j < vertices[i].numNbrs; j++)
            {
                intvl.u=i;intvl.v=vertices[i].neighbors[j].nbrId;
                intvl.nbrIndexFor_v=j;
                for (int k=0; k< vertices[i].neighbors[j].numIntvls;k++)
                {
                    intvl.intvlId=k;
                    intvl.intvlStart=vertices[i].neighbors[j].edgeSchedules[k].intvlStart;
                    intvl.intvlEnd=vertices[i].neighbors[j].edgeSchedules[k].intvlEnd;
                    intvl.lambda=vertices[i].neighbors[j].edgeSchedules[k].traveTime;
                    intvl.prevJourneyIndex=-1;
                    vertices[i].neighbors[j].edgeSchedules[k].prevJourneyIndex=-1;      //Initialize prevJourneyIndex for all preKnown Intvls.
                    listOfPreKnownIntvls.push_back(intvl);
                }
            }
        }
        //Sort the vector w.r.t arrival time and start time.
        std::sort<vector<intervalInfo>::iterator, compareIntvlsMWF>(listOfPreKnownIntvls.begin(), listOfPreKnownIntvls.end(), intvlCompare);
#ifdef __NO_CLASSES__
        listJourneys.resize(V);
        finalMWFJourneys.resize(V);
#else
        listJourneyClasses.resize(V);
        finalMWFJourneyClass.resize(V);
#endif
    }
}

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
    arr_hop_time.clear();
    arr_hop_time.resize(V);
    father.resize(V);
    for(int i=0; i<V; i++){
        arr_hop_time[i]= make_pair(infinity,infinity);
        father[i] = make_tuple(i,-1,0);
    }

}

void GraphDualCriteria::initial_ds_hbhshrtst()
{
    incrementalShortestJourney shrtstJrny(infinity,infinity);     //shortest at source (arrTime, length).
    shrtstPathAllVertices.clear();
    shrtstPathAllVertices.resize(V);   //vector of incoming journeys at each vertex.
    shrtstPathAllVertices.assign(V, shrtstJrny);
}


void GraphDualCriteria::initial_ds_ewa()
{
    for (int i=0; i< V; i++)
    {
        finalMWFJourneys[i].arrivalTime = infinity;
        finalMWFJourneys[i].wtTime = infinity;
        finalMWFJourneys[i].prevNode=-1;
        listJourneys[i].clear();
    }
}

void GraphDualCriteria::initial_ds_ewa_classes()
{
    for (int i=0; i< V; i++)
    {
        finalMWFJourneyClass[i].arrivalTimeStart = infinity;
        finalMWFJourneyClass[i].arrivalTimeEnd = infinity;
        finalMWFJourneyClass[i].wtTime = infinity;
        finalMWFJourneyClass[i].prevNode=-1;
        listJourneyClasses[i].clear();
    }
}


void GraphDualCriteria::run_mhf()
{
    time_sum=0;
    max_runTime = 0;
    max_Src = -1;
    
    for(int i = 0 ;i < sources.size() ;i ++)
    {
        initial_ds_eha();
        //earliest_arrival_minHop_pair(sources[i]);
        mhfHopByHop(sources[i]);
    }
    
    print_avg_time();
    cout << "max Runtime: " << max_runTime << "Max Src: " << max_Src << endl;
}

void GraphDualCriteria::run_hbh_shortest()
{
    time_sum=0;
    max_runTime = 0;
    max_Src = -1;
    
    for(int i = 0 ;i < sources.size();i ++)
    {
        initial_ds_hbhshrtst();
        //earliest_arrival_minHop_pair(sources[i]);
        shortestHopByHop(sources[i]);
    }
    
    print_avg_time();
    cout << "max Runtime: " << max_runTime << "Max Src: " << max_Src << endl;
}


void GraphDualCriteria::run_mwf()
{
    time_sum=0;
    
    for(int i = 0 ;i < sources.size() ; i++)
    {
        initial_ds_ea();
#ifdef __NO_CLASSES__
        initial_ds_ewa();
#else
        initial_ds_ewa_classes();
#endif
        //earliest_arrival(sources[i]);
        //initial_ds_ea();earliest_arrival_pair(sources[i]);
        //initial_ds_ea();
//        shortest_path(sources[i]);
        //minWaitForemost(sources[i]);
 //       minWaitForemostPrioritized(sources[i]);
        //minWaitForemostPrioritizedNoSet(sources[i]);
#ifdef __NO_CLASSES__
        //mwfStreamingIntvls(sources[i]);
        mwfStreamingIntvlsV2(sources[i]);
#else
        mwfStreamingIntvlsWJourneyClasses(sources[i]);
#endif
    }
    print_avg_time();
}


void GraphDualCriteria::printmhfResultsTest(int source)
{
    int rv = 0;
    ofstream earliestMHOut(resultsFile);
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
    if (t.GetRuntime() > max_runTime)
    {
        max_runTime = t.GetRuntime();
        max_Src = source;
    }
    //cout << "Num Nodes seen: " << numNodesSeen << endl;

#ifdef __TEST__
//    build_mhf_Journeys(source, earliestKnownTimeArrival, allHopJourneys);
//    printmhfResultsTest2(source, earliestKnownTimeArrival);
#endif

    print_shortest_paths(source);
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
    ofstream earliestMHFOut(resultsFile);
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


void GraphDualCriteria::shortestHopByHop(int source)
{
    Timer t;

    vector<list<incrementalShortestJourney>> aggShrtstListAtAllVerts;
    vector<tuple<int, list<incrementalShortestJourney>>> *k_1arrAtAllVertices =
                                                        new vector<tuple<int, list<incrementalShortestJourney>>>;  //list of vectors of new paths seen at k-1
    vector<tuple<int, list<incrementalShortestJourney>>> *k_arrAtAllVertices =
                                                        new vector<tuple<int, list<incrementalShortestJourney>>>;//List of vectors of new paths seen at k.
    vector<tuple<int, list<incrementalShortestJourney>>> *temp_arrAtAllVertices = NULL;
    vector<tuple<int, int>> lastk_ShortestHop;  //for each vertex, (last hopcount vert seen, index of this vertex in k_arrAllVertices for that hopcount)

    lastk_ShortestHop.resize(V);
    lastk_ShortestHop.assign(V, make_tuple(-1,-1)); //last hop in which vert seen, index in k_arrAtAllVerts for that hop.
    aggShrtstListAtAllVerts.resize(V);
    incrementalShortestJourney shrtstJrny(infinity,infinity);     //shortest at source (arrTime, length).
    
    shrtstJrny.assign(0, 0);
    shrtstPathAllVertices[source] = shrtstJrny;
    aggShrtstListAtAllVerts[source].emplace_front(shrtstJrny);
    list<incrementalShortestJourney> tempListK; //list used for building shrtst journeys from a prev journey.
    tempListK.emplace_front(shrtstJrny);
    (*k_1arrAtAllVertices).push_back(make_tuple(source, tempListK));
    int hopCount = 0;

    t.start();
    int numNewNodesInCurrentHop = 1;     //number of new nodes discovered.
    while ( (hopCount < V-1) && (numNewNodesInCurrentHop > 0))// && (numNodesSeen < V))
    {
        hopCount++;
        numNewNodesInCurrentHop = 0;
        int intvlID = -1, nextNodeId;
        Nbrs nextNbr;
        for (int i = 0; i < (*k_1arrAtAllVertices).size(); i++)        //extend all journeys from prev hop
        {
            int node = get<0>((*k_1arrAtAllVertices)[i]);
            list<incrementalShortestJourney> *listIncToExtend = &get<1>((*k_1arrAtAllVertices)[i]);
            if (listIncToExtend->size() == 0)   //If paths in k hop were dominated by paths in k-1 hops.
                continue;
            for (int j = 0; j < vertices[node].numNbrs; j++)    //look at all outgoing nbrs of node
            {
                nextNbr = vertices[node].neighbors[j];
                nextNodeId = nextNbr.nbrId;
                if (nextNodeId == source)
                    continue;                   //no point going back to source.
                tempListK.clear();
                int prevEarliestIntvlDepTime = infinity;
                list<incrementalShortestJourney>::iterator extendListIt;
                for ( extendListIt = (*listIncToExtend).end(); extendListIt != (*listIncToExtend).begin(); extendListIt--)//For each nbr proc. in paths at node in k-1 hop
                {
                    int arrTime = prev(extendListIt)->currentArrivalTime;
                    int jrnyLngth = prev(extendListIt)->journeyLength;
                    int earliestIntvldepartTime = earliestUseEdgeAfterT(node, nextNbr, arrTime, intvlID);
                    if ( (earliestIntvldepartTime == -1) || (earliestIntvldepartTime >= infinity)
                        || (intvlID == -1) || (earliestIntvldepartTime >= prevEarliestIntvlDepTime))
                        continue;
                    int trvlTime=nextNbr.edgeSchedules[intvlID].traveTime;
                    int depTime = earliestIntvldepartTime;
                    incrementalShortestJourney tempShortestJourney(depTime+trvlTime ,jrnyLngth+trvlTime);
                    list<incrementalShortestJourney>::iterator posTempListK = tempListK.begin();
                    checkShrtstDominanceAndPush(tempListK, tempShortestJourney, posTempListK);
                    int prevTrvlTime = trvlTime;
                    for (int intvlIndex=intvlID+1; intvlIndex< nextNbr.numIntvls; intvlIndex++) //Coll. all ND paths in tempArr from node to next forAll edgeSchedules. This could be optimized using priority search tree instead of all intervals
                    {
                        trvlTime=nextNbr.edgeSchedules[intvlIndex].traveTime;
                        if (trvlTime >= prevTrvlTime)
                            continue;
                        prevTrvlTime=trvlTime;
                        depTime=nextNbr.edgeSchedules[intvlIndex].intvlStart;
                        if (depTime >= prevEarliestIntvlDepTime)
                            break;
                        tempShortestJourney.assign(depTime+trvlTime, jrnyLngth+trvlTime);
                        checkShrtstDominanceAndPush(tempListK, tempShortestJourney, posTempListK);
                    }
                    prevEarliestIntvlDepTime = earliestIntvldepartTime;
                } //All possible extensions of incoming journeys at this node to nextNode have been explored.
                int idxK_arr;
                if ( get<0>(lastk_ShortestHop[nextNodeId]) == hopCount) //check nextNodeId in k_ArrAtAllVerts if any.
                {
                    idxK_arr=get<1>(lastk_ShortestHop[nextNodeId]);
                    mergeJourneys(get<1>((*k_arrAtAllVertices)[idxK_arr]), tempListK); //merge with jrnys at nextNode in k_ArrAtAllVerts
                }
                else if (tempListK.size() > 0) //If none, tempListK is pushed at end and lastk vector updated.
                {
                    (*k_arrAtAllVertices).push_back(make_tuple(nextNodeId, tempListK));
                    get<0>(lastk_ShortestHop[nextNodeId]) = hopCount;
                    get<1>(lastk_ShortestHop[nextNodeId]) = (int)(*k_arrAtAllVertices).size()-1;
                    //numNewNodesInCurrentHop++;
                }
            }  //All the nbrs of a node to extend in current hop count have been explored.
        }  //All new journeys in k-1 hop to any new node have been extended to all their nbrs.
        for (int k_jrnysIndex = 0; k_jrnysIndex < (*k_arrAtAllVertices).size(); k_jrnysIndex++) //upd shrtst jrny
        {
            int vertId = get<0>((*k_arrAtAllVertices)[k_jrnysIndex]);
            mergeJourneys(aggShrtstListAtAllVerts[vertId], get<1>((*k_arrAtAllVertices)[k_jrnysIndex]));
            if (get<1>((*k_arrAtAllVertices)[k_jrnysIndex]).size() > 0)
            {
                numNewNodesInCurrentHop++;
                incrementalShortestJourney k_shrtstJrny = get<1>((*k_arrAtAllVertices)[k_jrnysIndex]).back();
                if (k_shrtstJrny.journeyLength < shrtstPathAllVertices[vertId].journeyLength)
                    shrtstPathAllVertices[vertId] = k_shrtstJrny;
            }
        }
        k_1arrAtAllVertices->clear();         //Delete k-1 as they have been all explored.
        temp_arrAtAllVertices = k_1arrAtAllVertices;
        k_1arrAtAllVertices = k_arrAtAllVertices; //Make k-1 point to k.
        k_arrAtAllVertices = temp_arrAtAllVertices;
    }       //While loop to cover all hop counts.
    t.stop();
    time_sum += t.GetRuntime();
    if (t.GetRuntime() > max_runTime)
    {
        max_runTime = t.GetRuntime();
        max_Src = source;
    }
    //cout << "Num Nodes seen: " << numNodesSeen << endl;

#ifdef __TEST__
//    build_mhf_Journeys(source, earliestKnownTimeArrival, allHopJourneys);
//    printmhfResultsTest2(source, earliestKnownTimeArrival);
#endif

    print_shortest_hbh_results(source);
}

void GraphDualCriteria::print_shortest_hbh_results(int source)
{
    string resultSource = resultsFile + "_" + to_string(source);
    ofstream hbh_out(resultSource);

    //cout << "Results as (t,l) at every vertex:" <<endl;
    hbh_out << V << endl;
    for (int i=0; i< shrtstPathAllVertices.size(); i++)
    {
        if (shrtstPathAllVertices[i].currentArrivalTime >= infinity)
            continue;
        hbh_out << i << " " << shrtstPathAllVertices[i].currentArrivalTime << " " << shrtstPathAllVertices[i].journeyLength << endl;
    }
}


bool GraphDualCriteria::mergeJourneys(list<incrementalShortestJourney>& toList, list<incrementalShortestJourney>& fromList)
{
    list<incrementalShortestJourney>::iterator fromIt = fromList.begin();
    list<incrementalShortestJourney>::iterator toIt = toList.begin();
    while ((fromIt != fromList.end()) && (toIt != toList.end()) ) {
        if ((*toIt).currentArrivalTime < (*fromIt).currentArrivalTime) {
            if ((*toIt).journeyLength <= (*fromIt).journeyLength) {
                fromIt = fromList.erase(fromIt); //journey in fromIt is dominated, so remove it.
            }
            else
                toIt++;
        }
        else if ((*fromIt).currentArrivalTime < (*toIt).currentArrivalTime) {
            toIt = toList.insert(toIt, (*fromIt));
            toIt++;
            while (((*fromIt).journeyLength <= (*toIt).journeyLength)
                   && (toIt != toList.end()))//eliminate dominated journeys.
                toIt = toList.erase(toIt);
            fromIt++;
        }
        else {      //arrival time is same for both from and to
            if ((*toIt).journeyLength <= (*fromIt).journeyLength) //journey in fromIt is dominated
                fromIt = fromList.erase(fromIt); //journey in fromIt is dominated, so remove it.
            else
            {
                toIt = toList.insert(toIt, (*fromIt));
                toIt++;
                while (((*fromIt).journeyLength <= (*toIt).journeyLength)
                       && (toIt != toList.end()))//eliminate dominated journeys.
                    toIt = toList.erase(toIt);
                fromIt++;
            }
        }
    }
    incrementalShortestJourney lastMergedJourney(infinity,infinity); // = *(prev(toIt));
    if ((toList.size() > 0) && (toIt == toList.end()))
        lastMergedJourney = *(prev(toIt));
    if ( (toIt == toList.end()) && (fromIt != fromList.end()))
    {
        while ( (lastMergedJourney.journeyLength <= (*fromIt).journeyLength)
               && (fromIt != fromList.end()) )
            fromIt = fromList.erase(fromIt);       //skip these as they are dominated.
        while (fromIt != fromList.end()) {
            toIt = toList.insert(toIt, *fromIt);
            toIt++; fromIt++;
        }
    }
    if (toList.size() > 0)
        return true;
    else
        return false;
}

bool GraphDualCriteria::checkShrtstDominanceAndPush(list<incrementalShortestJourney>& listShrtstJrnys, incrementalShortestJourney latestShrtstJrny,
                                                    list<incrementalShortestJourney>::iterator& posList)
{
    if (posList == listShrtstJrnys.begin())
    {
        posList= listShrtstJrnys.insert(posList, latestShrtstJrny);
    }
    else
    {
        incrementalShortestJourney prevListJrny = *(prev(posList));
        if (latestShrtstJrny.journeyLength < prevListJrny.journeyLength) //length is less of new journey.
        {
            posList= listShrtstJrnys.insert(posList, latestShrtstJrny);     //replace last journey in array with this one.
        }
        else //if (latestShrtstJrny.journeyLength >= prevListJrny.journeyLength).
            return false;   //lngth is same or bigger and arriv is bigger. So dominated, hence ignore.
    }
    
    if (posList==listShrtstJrnys.end())
        return false;
    //If journey was inserted, clean up any subsequent dominated journeys.
    bool cleaned = false;
    list<incrementalShortestJourney>::iterator cleanIt = next(posList);
    while((cleanIt != listShrtstJrnys.end()) && (cleaned==false))
    {
        if (posList->journeyLength <= cleanIt->journeyLength)
            cleanIt = listShrtstJrnys.erase(cleanIt);
        else
            cleaned=true;
    }
    posList++;
    return true;
}




/*************************
 Data structures used:
 1. Graph data structure for intvl temporal graphs (vector of nodes. Each Node has vector of nbrs. Each nbr has vector of intvls. ( Same as we used earlier.)
 2. Interval as {u,v,u_nbrIndx,(s,e),lambda, intvlId, prevJourneyIndex}. //u_nbrIndx represents index of v in nbrs list of u in graph of (1). IntvlId represents index of interval on the edge (u,v) as represented in (1) above.
 3. Sorted list listIntervals[] of all intervals present in the input graph. Sort key is (interval.s + lambda). For contact sequence graph, this is similar to sorted list of edges with sort key as arrival time of the edge.
 4. Journey is represented as (a,w,prevNode,lastExpAt)
 5. lastExpAt used in (4) is a vector with num elements equal to num neighbors of destination of the journey. For each neihbor, it stores info (startTime, lambda) about interval in which this journey was last expanded to this nbr.
 6. listJourneys is a vector of (list of journeys) for each node u. listJourneys[u] is a sorted list (by arrival time) of journeys arrived at node u
 7. PQ mergeIntvls = {}            //Priority Queue of new Intervals created during Algorithm with arrival time (s+lambda) as the key.

mwfWalks(source)
{
 for each u { listJ[u] = {} }
 newJourney = {t_start,0,-1,{(-1,-1)} }       //(arrival time, wait time, prevNode, lastExpandedAt is a vector with (-1,-1) for each nbr of source as this journey is not expanded yet to any nbr.
 listJ[source].push_back(newJourney)
 setupNewJourneyForNeighbors(source,newJourney.arrivalTime);
 numNodesRchd = 1;      //source has been rchd.
 newIntvl = min (top(Intvls, mergeIntvls));     //Comparison of arrTime= s+lambda and secondary start time
 while ( (listIntervals[] || mergeIntvls[]) && (numNodesRchd < reachable) )
 {
    currArrivalTime = newIntvl.arrTime
    numNewNodesRchd=0;
 
 //Keep processing intervals without checking numNodes reached until the intervals have same arrival time. This is because even though all nodes may have been reached, there may still be journey coming in with same arrival time but less wait time.
    while (newIntvl != NULL && currArrialTime == newIntvl.arrTime)
    {
        removeMin(top(listIntervals,mergeIntvls)
        u=newIntvl.u; v=newIntvl.v;
        prevJourneyIndex = getPreviousJourney(newIntvl)
        If (prevJourneyIndex==-1)
            newIntvl = min (top(Intvls, mergeIntvls));     //Comparison of arrTime= s+lambda and secondary start time
            continue;
        If (listJ[u][prevJourneyIndex].expandedIntvl[v].lambda <  newIntvl.lambda or listJ[u][prevJourneyIndex].expandedIntvl[v] == {-1,-1})  //Based on intvl dominance criteria. Need to expand only when(lambda > lastLambda )
        {
            newJourney.arrTime=newIntvl.s+newIntvl.lambda;            //This will be extension of the journey at prevJourneyIndex
            newJourney.wtTime=listJ[u][prevJourneyIndex].w+ newIntvl.s-listJ[u][prevJourneyIndex].arr;      //This will be extension of the journey at prevJourneyIndex
            listJ[u][prevJourneyIndex].expandedIntvl[v] = {newIntvl.s,newIntvl.lambda};
            If (newJourney not dominated by listJ[v][listJ[v].last()] || listJ[v] == null)     //Dominance based on mwf criteria.
            {
                If (mwf[v] == NULL)
                    mwf[v]=newJourney;
                    numNewNodesRchd++;
                else if ((newJourney.arrTime == mwf[v].arrTime) && (newJourney.wtTime < mwf[v].wtTime) )
                    mwf[v]=newJourney;
                
                if ((newJourney.arrTime == listJ[v][listJ[v].last].arrTime) && (newJourney.wtTime < listJ[v][listJ[v].last].wtTime) )
                    listJ[v][listJ[v].last()] = newJourney;         //Overwrite last journey in the list with newJourney as newJourney has same arr Time but less wait time.
                else
                    listJ[v].append(newJourney);
                setupNewJourneyForNeighbors(v, arrivalTime)
            }
        }       //else intvl is ignored for expansion.
        newIntvl = min (top(Intvls, mergeIntvls));     //fetch next intvl for consideration. (it is not removed yet).
    }
    numNodesRchd+= numNewNodesRchd;
 }
}

getPreviousJourney(interval newIntvl)
{
    If (newIntvl.Id != -1)      //this means this intvl came from static list of listIntervals.
        prevJourneyIndex=vert[u][v][newIntvl.Id].prevJourneyIndex;
    else        //This means this intvl came from mergeIntvls list and was created by breaking up an existing interval.
        prevJourneyIndex=newIntvl.prevJourneyIndex
    if (prevJuorneyIndex==-1)
    {
        If (listJ[u].last().arr <= newIntvl.strt)
            prevJourneyIndex=listJ[u].last()
        else
            search prevJourney in listJ[u] (search journey in sorted list listJ[u] with arrTm <= newIntvl.strtTime)     //This is binary search in list of sorted journeys by arrival time.
    }
 }

 setupNewJourneyForNeighbors(int v, int arrivalTime)
 {
 for each nbr w of v
     nextIntvl = nextFn(v,w,t_start)
     If (nextIntvl.start >= arrivalTime)
         vert[v][w][nextIntvl].prevJourneyIndex = listJ[v].last();     //Can be stored as Id in listJ[v]
     else If ( (nextIntvl.start < arrivalTime) and (arrivalTime < nextIntvl.e) )
         e = nextIntvl.e;
         newIntvlCreated=(v,w,arrivalTime, e, nextIntvl.lambda);
         newIntvlCreated.prevJourney=listJ[v].last();
         mergeIntvls.insert({newIntvlCreated})
 }

 ****************************/
/*
void GraphDualCriteria::mwfStreamingIntvls(int source)
{
    mwfJourney newJourney;
    intervalInfo newIntvl;
    int indexPreKnownIntvls = 0;
    int prevJourneyIndex = -1;
    Timer t;

    int newIntvlFrom = 0;
    int numNodesRchd = 0;
    tuple<int,int> nodesRchable_maxFmstTime = earliest_arrival_pair(source,1); //4092652
    int nodesReachable= get<0>(nodesRchable_maxFmstTime), maxFmstTime=get<1>(nodesRchable_maxFmstTime);
//    cout << "Num Nodes Reachable: " << nodesReachable << endl;
    
    newJourney.arrivalTime = t_start;
    newJourney.prevNode = -1; newJourney.wtTime=0;newJourney.prevJourneyIndex=-1;
    newJourney.lastExpandedAt.resize(vertices[source].numNbrs);
    listJourneys[source].push_back(newJourney);        //known journeys so far at source.
    finalMWFJourneys[source] = newJourney;
    
    numNodesRchd++;
    int totalStaticIntvls = (int)listOfPreKnownIntvls.size();
//    cout << "Total num Intvls: " << totalStaticIntvls << endl;
    
    t.start();
//    setupNewJourney(source, newJourney.arrivalTime);      //Commented as may not be reqd for csg graphs. TBD
    newIntvlFrom = removeMinIntvl(newIntvl, indexPreKnownIntvls);
    int u=0,v=0,nbrIndex=0,intvlId=0;
    while ( ((!listOfAdHocIntvls.empty()) || (indexPreKnownIntvls < totalStaticIntvls))
           && (numNodesRchd < nodesReachable) && (newIntvlFrom != -1))// && (newIntvl.intvlStart+newIntvl.lambda <= maxFmstTime))
    {
        int currArrivalTime = newIntvl.intvlStart+newIntvl.lambda;
        int numNewNodesRchd=0;
        while ((currArrivalTime == newIntvl.intvlStart+newIntvl.lambda) && (newIntvlFrom != -1))
        {
            u = newIntvl.u; v=newIntvl.v; nbrIndex=newIntvl.nbrIndexFor_v; intvlId=newIntvl.intvlId;
            if (v == source)        //no point going back to source. move on to next intvl.     //TBD
            {
                newIntvlFrom = removeMinIntvl(newIntvl, indexPreKnownIntvls);  //Move to the next interval.
                continue;
            }
            prevJourneyIndex = getPrevJourney(newIntvl);
            if (prevJourneyIndex == -1)
            {
                newIntvlFrom = removeMinIntvl(newIntvl, indexPreKnownIntvls);  //Move to the next interval.
                continue;
            }
            int prevJourneyLastExtLmbda = get<1>(listJourneys[u][prevJourneyIndex].lastExpandedAt[nbrIndex]);
            if ((prevJourneyLastExtLmbda < newIntvl.lambda) || (u == source))  //Always expand from source as there is no waiting at source. TBD
            {
                newJourney.arrivalTime = newIntvl.intvlStart+newIntvl.lambda;
                if (u== source) //Never any waiting at source.      //TBD
                    newJourney.wtTime = 0;
                else
                    newJourney.wtTime = listJourneys[u][prevJourneyIndex].wtTime + newIntvl.intvlStart-listJourneys[u][prevJourneyIndex].arrivalTime;
                
                newJourney.prevNode=u;newJourney.prevJourneyIndex=prevJourneyIndex;newJourney.prevDepTime=newIntvl.intvlStart;
                //newJourney.prevLambda=newIntvl.lambda;newJourney.prevNodeNbrIndex=newIntvl.nbrIndexFor_v;
                //if carryover, the intvl will be start=carryoverJ.st - lambda. end = carryoverjJ.e-lambda.
                listJourneys[u][prevJourneyIndex].lastExpandedAt[nbrIndex] = make_tuple(newIntvl.intvlStart,newIntvl.lambda,1);
                newJourney.lastExpandedAt.clear();
                newJourney.lastExpandedAt.resize(vertices[v].numNbrs);
                int inserted = checkNewJourneyAndInsert(newJourney, v);
                if (inserted == 0)      //This journey was dominated by previous journey at v.
                {
                    newIntvlFrom = removeMinIntvl(newIntvl, indexPreKnownIntvls);
                    continue;
                }
                if (inserted == 2)
                {
                    numNewNodesRchd++;
                }
//                setupNewJourney(v, newJourney.arrivalTime);
            }
            newIntvlFrom = removeMinIntvl(newIntvl, indexPreKnownIntvls);
        }
        numNodesRchd += numNewNodesRchd;
    }
    t.stop();
    time_sum += t.GetRuntime();
//    printmwfResultsTest2(source);
}
*/
void GraphDualCriteria::mwfStreamingIntvlsV2(int source)
{
    mwfJourney newJourney;
    intervalInfo newIntvl;
    int indexPreKnownIntvls = 0;
    int prevJourneyIndex = -1;
    Timer t;

    int newIntvlFrom = 0;
    int numNodesRchd = 0;
    tuple<int,int> nodesRchable_maxFmstTime = earliest_arrival_pair(source,1); //4092652
    int nodesReachable= get<0>(nodesRchable_maxFmstTime), maxFmstTime=get<1>(nodesRchable_maxFmstTime);
//    cout << "Num Nodes Reachable: " << nodesReachable << endl;
    
    newJourney.arrivalTime = t_start;
    newJourney.arrivalTimeEnd=infinity;
    newJourney.isClass=true;
    newJourney.prevNode = -1; newJourney.wtTime=0;newJourney.prevJourneyIndex=-1;newJourney.lastTravelTime=-1;
    newJourney.prevNodeNbrIndex=-1;
//    newJourney.lastExpandedAt.resize(vertices[source].numNbrs);
    listJourneys[source].push_back(newJourney);        //known journeys so far at source.
    finalMWFJourneys[source] = newJourney;
    
    numNodesRchd++;
    int totalStaticIntvls = (int)listOfPreKnownIntvls.size();
//    cout << "Total num Intvls: " << totalStaticIntvls << endl;
    
    t.start();
    setupNewJourney(source, newJourney.arrivalTime);      //Commented as may not be reqd for csg graphs. TBD
    newIntvlFrom = removeMinIntvl(newIntvl, indexPreKnownIntvls);
    int u=0,v=0,nbrIndex=0,intvlId=0;
    while ( ((!listOfAdHocIntvls.empty()) || (indexPreKnownIntvls < totalStaticIntvls))
           && (numNodesRchd < nodesReachable) && (newIntvlFrom != -1))// && (newIntvl.intvlStart+newIntvl.lambda <= maxFmstTime))
    {
        int currArrivalTime = newIntvl.intvlStart+newIntvl.lambda;
        int numNewNodesRchd=0;
        while ((currArrivalTime == newIntvl.intvlStart+newIntvl.lambda) && (newIntvlFrom != -1))
        {
            u = newIntvl.u; v=newIntvl.v; nbrIndex=newIntvl.nbrIndexFor_v; intvlId=newIntvl.intvlId;
            if (v == source)        //no point going back to source. move on to next intvl.     //TBD
            {
                newIntvlFrom = removeMinIntvl(newIntvl, indexPreKnownIntvls);  //Move to the next interval.
                continue;
            }
            prevJourneyIndex = getPrevJourney(newIntvl);
            if (prevJourneyIndex == -1)
            {
                newIntvlFrom = removeMinIntvl(newIntvl, indexPreKnownIntvls);  //Move to the next interval.
                continue;
            }
            newJourney.isClass=false;
            newJourney.arrivalTimeEnd=-1;
            if ((listJourneys[u][prevJourneyIndex].isClass) &&
                (newIntvl.intvlEnd > newIntvl.intvlStart))  //intvl is not instantaneous and prevJourney was a class
            {
                if (listJourneys[u][prevJourneyIndex].arrivalTimeEnd > newIntvl.intvlStart)
                {
                    //newJourney.arrivalTime=newIntvl.intvlStart+newIntvl.lambda;
                    if (listJourneys[u][prevJourneyIndex].arrivalTimeEnd > newIntvl.intvlEnd)
                        newJourney.arrivalTimeEnd=newIntvl.intvlEnd+newIntvl.lambda;
                    else
                        newJourney.arrivalTimeEnd=
                        listJourneys[u][prevJourneyIndex].arrivalTimeEnd+newIntvl.lambda;
//                        listJourneys[u][prevJourneyIndex].lastExpandedAt[nbrIndex] = make_tuple(newIntvl.intvlStart,newIntvl.lambda,1);
                    newJourney.isClass = true;
                }
            }
//            int prevJourneyLastExtLmbda = get<1>(listJourneys[u][prevJourneyIndex].lastExpandedAt[nbrIndex]);
//            if ((prevJourneyLastExtLmbda < newIntvl.lambda) ||
//                (u == source)
//                || (newJourney.isClass))  //Always expand from source as there is no waiting at source. TBD
//            {
                newJourney.arrivalTime = newIntvl.intvlStart+newIntvl.lambda;
                if ((u== source) || (newJourney.isClass)) //Never any waiting at source.      //TBD
                    newJourney.wtTime = 0;
                else
                    newJourney.wtTime = listJourneys[u][prevJourneyIndex].wtTime + newIntvl.intvlStart-listJourneys[u][prevJourneyIndex].arrivalTime;
                
                newJourney.prevNode=u;newJourney.prevJourneyIndex=prevJourneyIndex;newJourney.prevDepTime=newIntvl.intvlStart;
                newJourney.lastTravelTime=newIntvl.lambda;newJourney.prevNodeNbrIndex=newIntvl.nbrIndexFor_v;
//                if (!newJourney.isClass)
//                    listJourneys[u][prevJourneyIndex].lastExpandedAt[nbrIndex] = make_tuple(newIntvl.intvlStart,newIntvl.lambda,1);
//                newJourney.lastExpandedAt.clear();
//                newJourney.lastExpandedAt.resize(vertices[v].numNbrs);
                int inserted = checkNewJourneyAndInsertV2(newJourney, v); //include check for classes.
                if (inserted == 0)      //This journey was dominated by previous journey at v.
                {
                    newIntvlFrom = removeMinIntvl(newIntvl, indexPreKnownIntvls);
                    continue;
                }
                if (inserted == 2)
                {
                    numNewNodesRchd++;
                }
                setupNewJourney(v, newJourney.arrivalTime);
//            }
            newIntvlFrom = removeMinIntvl(newIntvl, indexPreKnownIntvls);
        }
        numNodesRchd += numNewNodesRchd;
    }
    t.stop();
    time_sum += t.GetRuntime();
//    printmwfResultsTest2(source);
}



int GraphDualCriteria::checkNewJourneyAndInsert(mwfJourney& newJourney, int v)
{
    int inserted = 0;
    if (!listJourneys[v].empty())
    {
        tuple<int,int> j1; tuple<int,int> j2;
        get<0>(j1) = listJourneys[v][listJourneys[v].size()-1].arrivalTime; get<1>(j1) = listJourneys[v][listJourneys[v].size()-1].wtTime;
        get<0>(j2) = newJourney.arrivalTime;get<1>(j2) = newJourney.wtTime;
        if (!checkDominance(j1, j2))
        {
            if (get<0>(j2) > get<0>(j1))        //If new journey has bigger arrival time, append it.
                listJourneys[v].push_back(newJourney);
            else                                 //If arrival time is same, just replace last journey.
                listJourneys[v][listJourneys[v].size()-1] = newJourney; //Replace last journey as arrTm is same and wt time more.
            inserted = 1;
        }
    }
    else
    {
        inserted = 2;
        listJourneys[v].push_back(newJourney);
        finalMWFJourneys[v] = newJourney;
    }
    if ((finalMWFJourneys[v].arrivalTime == newJourney.arrivalTime) &&
        (finalMWFJourneys[v].wtTime > newJourney.wtTime))
        finalMWFJourneys[v] = newJourney;
    return inserted;
}

int GraphDualCriteria::checkNewJourneyAndInsertV2(mwfJourney& newJourney, int v)
{
    int inserted = 0; int j1End=-1, j2End=-1;
    if (!listJourneys[v].empty())
    {
        tuple<int,int> j1; tuple<int,int> j2;
        get<0>(j1) = listJourneys[v][listJourneys[v].size()-1].arrivalTime; get<1>(j1) = listJourneys[v][listJourneys[v].size()-1].wtTime;
        get<0>(j2) = newJourney.arrivalTime;get<1>(j2) = newJourney.wtTime;

        bool j1Class=listJourneys[v][listJourneys[v].size()-1].isClass;
        bool j2Class=newJourney.isClass;
        j1End=listJourneys[v][listJourneys[v].size()-1].arrivalTimeEnd;
        j2End=newJourney.arrivalTimeEnd;
        if ((get<0>(j2) <= j1End) && (j1Class)) //j1 is class and there was overlap with j2.
        {
            if (j2Class)    //j2 is also a class
            {
                newJourney.arrivalTime=j1End+1;
                buildAndPushIntvlInHeapV2(newJourney, v);
            }
            return 0;   //there was overlap with an existing class so nothing inserted.
        }
        if (j1Class)        //There is no overlap and j1 was a class but j2 is not.
            get<0>(j1)=j1End;   //dominance will need to be tested with last instance in j1Class
        if (!checkDominance(j1, j2))    //both are regular walks.
        {
            if (get<0>(j2) > get<0>(j1))        //If new journey has bigger arrival time, append it.
                listJourneys[v].push_back(newJourney);
            else                                 //If arrival time is same, just replace last journey.
                listJourneys[v][listJourneys[v].size()-1] = newJourney; //Replace last journey as arrTm is same and wt time more.
            inserted = 1;
        }
    }
    else
    {
        inserted = 2;
        listJourneys[v].push_back(newJourney);
        finalMWFJourneys[v] = newJourney;
    }
    if ((finalMWFJourneys[v].arrivalTime == newJourney.arrivalTime) &&
        (finalMWFJourneys[v].wtTime > newJourney.wtTime))
        finalMWFJourneys[v] = newJourney;
    return inserted;
}


int GraphDualCriteria::getPrevJourney(intervalInfo& intvl)
{
    int prevJourneyIndex = -1;
    int u=intvl.u, nbrIndex=intvl.nbrIndexFor_v,intvlId=intvl.intvlId;
    
    if (intvlId != -1)
    {
        prevJourneyIndex = vertices[u].neighbors[nbrIndex].edgeSchedules[intvlId].prevJourneyIndex;
    }
    else
        prevJourneyIndex = intvl.prevJourneyIndex;
    
    if ((prevJourneyIndex == -1) && (!listJourneys[u].empty()) )
    {
        prevJourneyIndex = searchPrevJourney(u,intvl.intvlStart);
    }
    return prevJourneyIndex;
}

void GraphDualCriteria::setupNewJourney(int v, int arrivalTime)
{
    int nextTravelTime=0,nextIntvl=0;
    compareIntvlsMWFMinHeap intvlCompareObjForHeap;
    for (int i=0; i< vertices[v].numNbrs;i++)
    {
//        listJourneys[v][(int)listJourneys[v].size()-1].lastExpandedAt[i] = make_tuple(-1,-1,0); //starting journey isnt expanded yet
        nextTravelTime = earliestUseEdgeAfterT(v, vertices[v].neighbors[i], arrivalTime, nextIntvl);
        if ((nextIntvl == -1) || (nextTravelTime >= infinity))
            continue;
        if (arrivalTime <= vertices[v].neighbors[i].edgeSchedules[nextIntvl].intvlStart)
        {
            vertices[v].neighbors[i].edgeSchedules[nextIntvl].prevJourneyIndex=(int)listJourneys[v].size()-1;
        }
        else if ( (vertices[v].neighbors[i].edgeSchedules[nextIntvl].intvlStart < arrivalTime)
                 && (arrivalTime <= vertices[v].neighbors[i].edgeSchedules[nextIntvl].intvlEnd)
                 && arrivalTime > vertices[v].neighbors[i].edgeSchedules[nextIntvl].divTime) //If divTime is same,intvl has already be divided at this arrival time.
        {
            intervalInfo adHocIntvl;
            adHocIntvl.intvlEnd = vertices[v].neighbors[i].edgeSchedules[nextIntvl].intvlEnd;
            adHocIntvl.intvlStart = nextTravelTime; adHocIntvl.intvlId=-1;
            adHocIntvl.lambda= vertices[v].neighbors[i].edgeSchedules[nextIntvl].traveTime;
            adHocIntvl.u = v; adHocIntvl.nbrIndexFor_v=i; adHocIntvl.v = vertices[v].neighbors[i].nbrId;
            adHocIntvl.prevJourneyIndex = (int)listJourneys[v].size()-1;
            listOfAdHocIntvls.push_back(adHocIntvl); std::push_heap(listOfAdHocIntvls.begin(), listOfAdHocIntvls.end(), intvlCompareObjForHeap);
            vertices[v].neighbors[i].edgeSchedules[nextIntvl].divTime = arrivalTime;  //The intvl is subdivided at arrivalTime. It need not be divided again at same time. update divtime.
        }
    }
}

void GraphDualCriteria::setupNewJourneyClass(int v, int arrivalTime)
{
    int nextTravelTime=0,nextIntvl=0;
    compareIntvlsMWFMinHeap intvlCompareObjForHeap;
    for (int i=0; i< vertices[v].numNbrs;i++)       //This could be the problem. vector lastExpandedAt not been resized
    {
        listJourneyClasses[v][(int)listJourneyClasses[v].size()-1].lastExpandedAt[i] = make_tuple(-1,-1,0); //starting journey isnt expanded yet
        nextTravelTime = earliestUseEdgeAfterT(v, vertices[v].neighbors[i], arrivalTime, nextIntvl);
        if ((nextIntvl == -1) || (nextTravelTime >= infinity))
            continue;
        if (arrivalTime <= vertices[v].neighbors[i].edgeSchedules[nextIntvl].intvlStart)
        {
            vertices[v].neighbors[i].edgeSchedules[nextIntvl].prevJourneyIndex=(int)listJourneyClasses[v].size()-1;
        }
        else if ( (vertices[v].neighbors[i].edgeSchedules[nextIntvl].intvlStart < arrivalTime)
                 && (arrivalTime <= vertices[v].neighbors[i].edgeSchedules[nextIntvl].intvlEnd) )
        {
            intervalInfo adHocIntvl;
            adHocIntvl.intvlEnd = vertices[v].neighbors[i].edgeSchedules[nextIntvl].intvlEnd;
            adHocIntvl.intvlStart = nextTravelTime; adHocIntvl.intvlId=-1;
            adHocIntvl.lambda= vertices[v].neighbors[i].edgeSchedules[nextIntvl].traveTime;
            adHocIntvl.u = v; adHocIntvl.nbrIndexFor_v=i; adHocIntvl.v = vertices[v].neighbors[i].nbrId;
            adHocIntvl.prevJourneyIndex = (int)listJourneyClasses[v].size()-1;
            listOfAdHocIntvls.push_back(adHocIntvl);
            std::push_heap(listOfAdHocIntvls.begin(), listOfAdHocIntvls.end(), intvlCompareObjForHeap);
        }
    }
}

//CheckOverlap & dominance & carryover. Then insert in journeys list.
//put Carryover back in intvl heap.

//Datastructures:
//Graph
//intvlInfo
//intvlsList
//intvlHeap or AdHocIntvlsHeap
//mwfJourneyClass
//listClasses at each node
//finalMWFClass

void GraphDualCriteria::mwfStreamingIntvlsWJourneyClasses(int source)
{
    mwfJourneyClass newJourneyClass;
    intervalInfo newIntvl;
    int indexPreKnownIntvls = 0;
    int prevJourneyClassIndex = -1;
    Timer t;

    int newIntvlFrom = 0;
    int numNodesRchd = 0;
    tuple<int,int> nodesRchable_maxFmstTime = earliest_arrival_pair(source,1); //4092652
    int nodesReachable= get<0>(nodesRchable_maxFmstTime), maxFmstTime=get<1>(nodesRchable_maxFmstTime);
//    cout << "Num Nodes Reachable: " << nodesReachable << endl;
    
    newJourneyClass.arrivalTimeStart = t_start;
    newJourneyClass.arrivalTimeEnd = infinity;
    newJourneyClass.prevNode = -1; newJourneyClass.wtTime=0; newJourneyClass.prevJourneyIndex=-1;
    newJourneyClass.lastExpandedAt.resize(vertices[source].numNbrs);
    listJourneyClasses[source].push_back(newJourneyClass);        //known journeys so far at source.
    recordFinalJourney(source, newJourneyClass);
    //finalMWFJourneyClass[source] = newJourneyClass;
    
    numNodesRchd++;
    int totalStaticIntvls = (int)listOfPreKnownIntvls.size();
//    cout << "Total num Intvls: " << totalStaticIntvls << endl;
    
    t.start();
    setupNewJourneyClass(source, newJourneyClass.arrivalTimeStart);      //Commented as may not be reqd for csg graphs. TBD
    newIntvlFrom = removeMinIntvl(newIntvl, indexPreKnownIntvls);
    int u=0,v=0,nbrIndex=0,intvlId=0;
    while ( ((!listOfAdHocIntvls.empty()) || (indexPreKnownIntvls <= totalStaticIntvls))
           && (numNodesRchd < nodesReachable) && (newIntvlFrom != -1))// && (newIntvl.intvlStart+newIntvl.lambda <= maxFmstTime))
    {
        int currArrivalTime = newIntvl.intvlStart+newIntvl.lambda;
        int numNewNodesRchd=0;
        while ((currArrivalTime == newIntvl.intvlStart+newIntvl.lambda) && (newIntvlFrom != -1))
        {
            u = newIntvl.u; v=newIntvl.v; nbrIndex=newIntvl.nbrIndexFor_v; intvlId=newIntvl.intvlId;
            if (v != source)        //no point going back to source. move on to next intvl.
            {
                prevJourneyClassIndex = getPrevJourneyClass(newIntvl);
                if (prevJourneyClassIndex != -1)
                {
                    int newJourneyCreated = createNewJourneyClass(listJourneyClasses[u][prevJourneyClassIndex], newIntvl, newJourneyClass, prevJourneyClassIndex);
                    if (newJourneyCreated)
                    {
                        int inserted = checkNewJourneyClassAndInsert(newJourneyClass, v);
                        if (inserted != 0)      //This journey was dominated by previous journey at v.
                        {
                            if (inserted == 2)
                            {
                                numNewNodesRchd++;
                            }
                            setupNewJourneyClass(v, newJourneyClass.arrivalTimeStart);
                        }
                    }
                }
            }
            newIntvlFrom = removeMinIntvl(newIntvl, indexPreKnownIntvls);
        }
        numNodesRchd += numNewNodesRchd;
    }
    t.stop();
    time_sum += t.GetRuntime();
//    printmwfResultsTest2(source);
    printmwfClassesResultsTest2(source);
}

int GraphDualCriteria::getPrevJourneyClass(intervalInfo& intvl)
{
    int prevJourneyIndex = -1;
    int u=intvl.u, nbrIndex=intvl.nbrIndexFor_v,intvlId=intvl.intvlId;
    
    if (intvlId != -1)
    {
        prevJourneyIndex = vertices[u].neighbors[nbrIndex].edgeSchedules[intvlId].prevJourneyIndex;
    }
    else
        prevJourneyIndex = intvl.prevJourneyIndex;
    
    if ((prevJourneyIndex == -1) && (!listJourneyClasses[u].empty()) )
    {
        prevJourneyIndex = searchPrevJourneyClass(u,intvl.intvlStart);
    }
    return prevJourneyIndex;
}


int GraphDualCriteria::createNewJourneyClass(mwfJourneyClass& prevJourneyClass, intervalInfo& intervalToExpand, mwfJourneyClass& newJourneyClass, int prevJourenyClassIndex)
{
    int retVal = 0;
    int u = intervalToExpand.u, v_nbrIndex = intervalToExpand.nbrIndexFor_v, v = intervalToExpand.v;
    if (intervalToExpand.intvlStart < prevJourneyClass.arrivalTimeEnd)
    {
        retVal = 1;
        newJourneyClass.wtTime = prevJourneyClass.wtTime;
        if (intervalToExpand.intvlEnd < prevJourneyClass.arrivalTimeEnd)
            newJourneyClass.arrivalTimeEnd = intervalToExpand.intvlEnd + intervalToExpand.lambda;
        else
        {
            newJourneyClass.arrivalTimeEnd = prevJourneyClass.arrivalTimeEnd + intervalToExpand.lambda;
            prevJourneyClass.lastExpandedAt[v_nbrIndex] = make_tuple(intervalToExpand.intvlStart, intervalToExpand.lambda,1);       //TBD - check v_nbrIndex < numNbrs.
        }
    }
    else
    {
        int prevJourneyClassEndLastExtLmbda = get<1>(prevJourneyClass.lastExpandedAt[v_nbrIndex]);
        if (prevJourneyClassEndLastExtLmbda < intervalToExpand.lambda)// || (u == source))  //Always expand from source as there is no waiting at source. TBD
        {
            retVal = 1;
            newJourneyClass.arrivalTimeEnd = intervalToExpand.intvlStart+intervalToExpand.lambda;
            newJourneyClass.wtTime = prevJourneyClass.wtTime + intervalToExpand.intvlStart - prevJourneyClass.arrivalTimeEnd;
            prevJourneyClass.lastExpandedAt[v_nbrIndex] = make_tuple(intervalToExpand.intvlStart, intervalToExpand.lambda,1);       //TBD - check v_nbrIndex < numNbrs.
        }
    }
    
    if (retVal == 1)            //TBD lastExpAt has not been resized properly.
    {
        newJourneyClass.arrivalTimeStart = intervalToExpand.intvlStart+intervalToExpand.lambda;
        newJourneyClass.lastTravelTime = intervalToExpand.lambda;
        newJourneyClass.prevJourneyIndex = prevJourenyClassIndex;
        newJourneyClass.prevNode=u; newJourneyClass.prevDepTime=intervalToExpand.intvlStart;
        newJourneyClass.prevNodeNbrIndex= intervalToExpand.nbrIndexFor_v;
        newJourneyClass.lastExpandedAt.clear();
        newJourneyClass.lastExpandedAt.resize(vertices[v].numNbrs);
    }
    return retVal;
}

int GraphDualCriteria::checkNewJourneyClassAndInsert(mwfJourneyClass& newJourneyClass, int v)
{
    int inserted = 0;
    bool newJourneyDominated = false;
    mwfJourneyClass carryOverJClass;
    carryOverJClass.arrivalTimeEnd = -1;
    if (!listJourneyClasses[v].empty())
    {
        mwfJourneyClass lastJClass = listJourneyClasses[v][listJourneyClasses[v].size()-1];
        newJourneyDominated = resolveOverlap(lastJClass,newJourneyClass,carryOverJClass);     //This results in 1 or 2 non-overlapping journeys. Overlap portion taken care of.
        if (lastJClass.arrivalTimeStart <= lastJClass.arrivalTimeEnd)   //lastJClass survived.
        {
            if (!newJourneyDominated) //&& (listJourneyClasses[v][listJourneyClasses[v].size()-1].arrivalTimeEnd == lastJClass.arrivalTimeEnd))
                newJourneyDominated = checkJourneyClassDominance(lastJClass,newJourneyClass);
            if (!newJourneyDominated)     //New Class survived.
            {
//                listJourneyClasses[v][listJourneyClasses[v].size()-1].arrivalTimeEnd = lastJClass.arrivalTimeEnd;       //should comment out.
                listJourneyClasses[v].push_back(newJourneyClass);
                inserted = 1;
//                checkJourneyClassDominance(newJourneyClass,carryOverJClass);        //should comment out.
//                buildAndPushIntvlInHeap(carryOverJClass, v);                 //should comment out. //push whatever part carryOver survived back on intvl heap.
            }
            else        //new journey was dominated
            {
                buildAndPushIntvlInHeap(newJourneyClass, v);    //the new journey class was dominated, push back the portion of new journey class that survived.
                return 0;
            }
        }
        else
        {
            listJourneyClasses[v][listJourneyClasses[v].size()-1] = newJourneyClass;        //TBD record all values. Don't copy as lastExpAt may not be sized appropriately.
            inserted=1;
//            checkJourneyClassDominance(newJourneyClass,carryOverJClass);        //not needed comment out
//            buildAndPushIntvlInHeap(carryOverJClass, v);        //push back carryover that survived.        //not needed comment out.
        }
    }
    else
    {
        inserted = 2;
        listJourneyClasses[v].push_back(newJourneyClass);
//        finalMWFJourneyClass[v] = newJourneyClass;
        recordFinalJourney(v,newJourneyClass);
    }

    if (inserted == 1)
    {
        if ((finalMWFJourneyClass[v].arrivalTimeStart == newJourneyClass.arrivalTimeStart) &&
            (finalMWFJourneyClass[v].wtTime > newJourneyClass.wtTime))
            //finalMWFJourneyClass[v] = newJourneyClass;
            recordFinalJourney(v, newJourneyClass);
    }
    return inserted;
}

void GraphDualCriteria::recordFinalJourney(int v, mwfJourneyClass& newJourneyClass)
{
    finalMWFJourneyClass[v].arrivalTimeStart = newJourneyClass.arrivalTimeStart;
    finalMWFJourneyClass[v].arrivalTimeEnd = newJourneyClass.arrivalTimeEnd;
    finalMWFJourneyClass[v].wtTime = newJourneyClass.wtTime;
    finalMWFJourneyClass[v].lastTravelTime=newJourneyClass.lastTravelTime;
    finalMWFJourneyClass[v].prevNode=newJourneyClass.prevNode;
    finalMWFJourneyClass[v].prevDepTime=newJourneyClass.prevDepTime;
    finalMWFJourneyClass[v].prevJourneyIndex=newJourneyClass.prevJourneyIndex;
    finalMWFJourneyClass[v].prevNodeNbrIndex=newJourneyClass.prevNodeNbrIndex;
}
        
bool GraphDualCriteria::resolveOverlap(mwfJourneyClass &lastJClass, mwfJourneyClass &newJourneyClass, mwfJourneyClass &carryOverJClass)
{
    bool newJourneyDominated = false;
    if (newJourneyClass.arrivalTimeStart > lastJClass.arrivalTimeEnd)   //No overlap.
        return false;
    else if (newJourneyClass.wtTime >= lastJClass.wtTime)        //new is Dominated, so just trim it at the beginnig and return;
    {
        newJourneyClass.arrivalTimeStart = lastJClass.arrivalTimeEnd+1;
        newJourneyDominated = true;
        return newJourneyDominated;
    }
    else if (newJourneyClass.arrivalTimeStart == lastJClass.arrivalTimeStart) //new arriv at same time with less wt
    {
        //carryOverJClass.arrivalTimeStart = newJourneyClass.arrivalTimeEnd+1;
        lastJClass.arrivalTimeEnd=-1;
    }/*
    else //There is some overlap and newJourneyClass.wt < lastJClass.wt.
    {
        carryOverJClass = lastJClass;       //
        if (newJourneyClass.arrivalTimeStart == lastJClass.arrivalTimeStart)
        {
            carryOverJClass.arrivalTimeStart = newJourneyClass.arrivalTimeEnd+1;
            lastJClass.arrivalTimeEnd=-1;
        }
        else if (newJourneyClass.arrivalTimeStart > lastJClass.arrivalTimeStart)
        {
            lastJClass.arrivalTimeEnd = newJourneyClass.arrivalTimeStart-1;
            carryOverJClass.arrivalTimeStart=newJourneyClass.arrivalTimeEnd+1;
        }
    }*/
    return newJourneyDominated;
}

bool GraphDualCriteria::checkJourneyClassDominance(mwfJourneyClass& firstJClass, mwfJourneyClass& nextJClass) //the two journeys have to be non-overlapping
{
    int cutOverT;
    bool dominated = false;
    if (nextJClass.arrivalTimeStart > nextJClass.arrivalTimeEnd)    //nextJourney is not valid so it is dominated
        return true;
    if (firstJClass.arrivalTimeStart > firstJClass.arrivalTimeEnd)
        return  false;
    if ((nextJClass.wtTime <= firstJClass.wtTime) &&
        (nextJClass.arrivalTimeStart > firstJClass.arrivalTimeEnd))
        return false;
    
    cutOverT = nextJClass.wtTime - firstJClass.wtTime + firstJClass.arrivalTimeEnd;
    if (cutOverT < nextJClass.arrivalTimeStart)
        return false;
    else
    {
        //nextJClass.arrivalTimeStart=cutOverT+1;
        nextJClass.arrivalTimeEnd=-1;
        dominated = true;
    }
    return dominated;
}

void GraphDualCriteria::buildAndPushIntvlInHeap(mwfJourneyClass &carryOverJClass, int v)
{
    if (carryOverJClass.arrivalTimeStart > carryOverJClass.arrivalTimeEnd)      //not a valid journey remains to be pushed back.
        return;
    compareIntvlsMWFMinHeap intvlCompareObjForHeap;
    intervalInfo carryOverIntvl;
    carryOverIntvl.u = carryOverJClass.prevNode;
    carryOverIntvl.v = v;
    carryOverIntvl.intvlStart = carryOverJClass.arrivalTimeStart-carryOverJClass.lastTravelTime;
    carryOverIntvl.intvlEnd = carryOverJClass.arrivalTimeEnd - carryOverJClass.lastTravelTime;
    carryOverIntvl.lambda = carryOverJClass.lastTravelTime;
    carryOverIntvl.nbrIndexFor_v = carryOverJClass.prevNodeNbrIndex;
    carryOverIntvl.prevJourneyIndex = -1; carryOverIntvl.intvlId=-1;    //setting prev journeyIndex to -1 as there may have been some new prev journeys by now.
    listOfAdHocIntvls.push_back(carryOverIntvl);
    std::push_heap(listOfAdHocIntvls.begin(), listOfAdHocIntvls.end(), intvlCompareObjForHeap);
}

void GraphDualCriteria::buildAndPushIntvlInHeapV2(mwfJourney &carryOverJ, int v)
{
    if (carryOverJ.arrivalTime > carryOverJ.arrivalTimeEnd)      //not a valid class so nothing to pushed back.
        return;
    compareIntvlsMWFMinHeap intvlCompareObjForHeap;
    intervalInfo carryOverIntvl;
    carryOverIntvl.u = carryOverJ.prevNode;
    carryOverIntvl.v = v;
    carryOverIntvl.intvlStart = carryOverJ.arrivalTime-carryOverJ.lastTravelTime;
    carryOverIntvl.intvlEnd = carryOverJ.arrivalTimeEnd - carryOverJ.lastTravelTime;
    carryOverIntvl.lambda = carryOverJ.lastTravelTime;
    carryOverIntvl.nbrIndexFor_v = carryOverJ.prevNodeNbrIndex;
    carryOverIntvl.prevJourneyIndex = -1; carryOverIntvl.intvlId=-1;    //setting prev journeyIndex to -1 as there may have been some new prev journeys by now.
    listOfAdHocIntvls.push_back(carryOverIntvl);
    std::push_heap(listOfAdHocIntvls.begin(), listOfAdHocIntvls.end(), intvlCompareObjForHeap);
//    if (listJourneys[carryOverJ.prevNode][carryOverJ.prevJourneyIndex].arrivalTimeEnd == carryOverIntvl.intvlEnd)
//        get<1>(listJourneys[carryOverJ.prevNode][carryOverJ.prevJourneyIndex].
//           lastExpandedAt[carryOverJ.prevNodeNbrIndex]) = -1;   //Since the last portion of the journey was pushed back, it was not expanded.

}




void GraphDualCriteria::printmwfClassesResultsTest2(int source)
{
    int rv = 0;
//    ofstream earliestMWFOut(mwfResults);
    cout << V << "\n";
    for (int i = 0; i < finalMWFJourneyClass.size(); i++)
    {
        if (finalMWFJourneyClass[i].arrivalTimeStart >= infinity)
            continue;
        cout << i << " " << finalMWFJourneyClass[i].arrivalTimeStart << "  "  <<  finalMWFJourneyClass[i].wtTime <<"\n";
        if (finalMWFJourneyClass[i].arrivalTimeStart < infinity)
        {
            rv++;
            mwfJourneyClass currJourneyClass=finalMWFJourneyClass[i];
            cout << "("<< i <<"," <<finalMWFJourneyClass[i].arrivalTimeStart << "," <<finalMWFJourneyClass[i].wtTime<< ")<--";
            while (currJourneyClass.prevNode != -1)
            {
                int prevNode = currJourneyClass.prevNode, depTime = currJourneyClass.prevDepTime;
                int prevIndex=currJourneyClass.prevJourneyIndex;
                currJourneyClass = listJourneyClasses[prevNode].at(prevIndex);
                cout << depTime<< "("<< prevNode << "," << currJourneyClass.arrivalTimeStart <<"," << currJourneyClass.wtTime << ")<--";
            }
            cout << source << endl;
        }
    }
    cout << "Source: " << source << "\n";
    cout << "Num reachable vertices: " << rv << "\n";
    cout << "Ratio of reachable vertices " << (float)rv/V << "\n";
//    avgHops = sumHops/rv;
//    cout << "Max Hops = " << maxHopCount << "\n";
//    cout << "Avg # hops = " << avgHops << "\n";
}



void GraphDualCriteria::printmwfResultsTest2(int source)
{
    int rv = 0;
//    ofstream earliestMWFOut(mwfResults);
    cout << V << "\n";
    for (int i = 0; i < finalMWFJourneys.size(); i++)
    {
//        if (final >= infinity)
//            continue;
        cout << i << " " << finalMWFJourneys[i].arrivalTime << "  "  <<  finalMWFJourneys[i].wtTime <<"\n";
        if (finalMWFJourneys[i].arrivalTime < infinity)
        {
            rv++;
/*            mwfJourney currJourney=finalMWFJourneys[i];
            cout << "("<< i <<"," <<finalMWFJourneys[i].arrivalTime << "," <<finalMWFJourneys[i].wtTime<< ")<--";
            while (currJourney.prevNode != -1)
            {
                int prevNode = currJourney.prevNode, depTime = currJourney.prevDepTime;
                int prevIndex=currJourney.prevJourneyIndex;
                currJourney = listJourneys[prevNode].at(prevIndex);
                cout << depTime<< "("<< prevNode << "," << currJourney.arrivalTime <<"," << currJourney.wtTime << ")<--";
            }
            cout << source << endl;*/
        }
    }
    cout << "Source: " << source << "\n";
    cout << "Num reachable vertices: " << rv << "\n";
    cout << "Ratio of reachable vertices " << (float)rv/V << "\n";
//    avgHops = sumHops/rv;
//    cout << "Max Hops = " << maxHopCount << "\n";
//    cout << "Avg # hops = " << avgHops << "\n";
}


int GraphDualCriteria::searchPrevJourney(int node, int beforeTime)  //Returns index of journey at node before beforeTime.
{
    if (listJourneys[node].empty())
        return -1;
    int low = 0, high = (int)listJourneys[node].size()-1, mid = high/2, journeyFound = 0, retVal = -1;
    vector<mwfJourney> *nodeJourneys = &listJourneys[node];
    
    //Check last 2 journeys.
    if (nodeJourneys->at(high).arrivalTime <= beforeTime)
        return high;
    else if (high >=1)
    {
        if (nodeJourneys->at(high-1).arrivalTime <= beforeTime)
            return (high-1);
    }
    else if (high==0)
        return -1;
        
    
    while (!journeyFound)
    {
        //cout << "Went into !journeyFound loop" << endl;
        if (nodeJourneys->at(mid).arrivalTime <= beforeTime)
        {
            if (beforeTime == nodeJourneys->at(mid).arrivalTime)     //t is within the intvl.
            {
                journeyFound = 1;
                retVal = mid;
            }
            else
            {
                if ((high - low) > 1)
                {
                    low = mid;
                    mid = (high + low)/2;
                }
                else if (mid == high)
                {
                    journeyFound = 1;
                    retVal = mid;
                }
                else        //Difference in low and high is 1 and mid was low. This means, beforeTime occurs after arrival of low. Lets check if high satisfies beforeTime, if yes return high, else return mid/low
                {
                    if (nodeJourneys->at(high).arrivalTime <= beforeTime)
                    {
                        journeyFound = 1;
                        retVal = high;
                    }
                    else
                    {
                        journeyFound = 1;
                        retVal = mid;
                    }
                }
            }
        }
        else        //mid was after beforeTime. This means journey at mid is too late.
        {
            if ((high - low) > 1)
            {
                high = mid-1;        //mid can be safely excluded cuz it doesn't satisfy transmission after t.
                mid = (low + high)/2;
            }
            else        //High and low are adjacent or same. mid is too late, so only possiblity is low.
            {
                if (mid == high)
                {
                    if (nodeJourneys->at(low).arrivalTime <= beforeTime)
                    {
                        journeyFound = 1;
                        retVal = low;
                    }
                    else
                    {
                        journeyFound = 1;
                        retVal = -1;            //There is no valid interval.
                    }
                }
                else
                {
                    journeyFound = 1;
                    retVal = -1;            //There is no valid interval.
                }
            }
        }
    }
    return retVal;
}

int GraphDualCriteria::searchPrevJourneyClass(int node, int beforeTime)  //Returns index of journey at node before beforeTime.
{
    if (listJourneyClasses[node].empty())
        return -1;
    int low = 0, high = (int)listJourneyClasses[node].size()-1, mid = high/2, journeyFound = 0, retVal = -1;
    vector<mwfJourneyClass> *nodeJourneys = &listJourneyClasses[node];
    
    //Check last 2 journeys.
    if (nodeJourneys->at(high).arrivalTimeStart <= beforeTime)
        return high;
    else if (high >=1)
    {
        if (nodeJourneys->at(high-1).arrivalTimeStart <= beforeTime)
            return (high-1);
    }
    else if (high==0)
        return -1;
        
    
    while (!journeyFound)
    {
        //cout << "Went into !journeyFound loop" << endl;
        if (nodeJourneys->at(mid).arrivalTimeStart <= beforeTime)
        {
            if (beforeTime == nodeJourneys->at(mid).arrivalTimeStart)     //t is within the intvl.
            {
                journeyFound = 1;
                retVal = mid;
            }
            else
            {
                if ((high - low) > 1)
                {
                    low = mid;
                    mid = (high + low)/2;
                }
                else if (mid == high)
                {
                    journeyFound = 1;
                    retVal = mid;
                }
                else        //Difference in low and high is 1 and mid was low. This means, beforeTime occurs after arrival of low. Lets check if high satisfies beforeTime, if yes return high, else return mid/low
                {
                    if (nodeJourneys->at(high).arrivalTimeStart <= beforeTime)
                    {
                        journeyFound = 1;
                        retVal = high;
                    }
                    else
                    {
                        journeyFound = 1;
                        retVal = mid;
                    }
                }
            }
        }
        else        //mid was after beforeTime. This means journey at mid is too late.
        {
            if ((high - low) > 1)
            {
                high = mid-1;        //mid can be safely excluded cuz it doesn't satisfy transmission after t.
                mid = (low + high)/2;
            }
            else        //High and low are adjacent or same. mid is too late, so only possiblity is low.
            {
                if (mid == high)
                {
                    if (nodeJourneys->at(low).arrivalTimeStart <= beforeTime)
                    {
                        journeyFound = 1;
                        retVal = low;
                    }
                    else
                    {
                        journeyFound = 1;
                        retVal = -1;            //There is no valid interval.
                    }
                }
                else
                {
                    journeyFound = 1;
                    retVal = -1;            //There is no valid interval.
                }
            }
        }
    }
    return retVal;
}




int GraphDualCriteria::removeMinIntvl(intervalInfo &newIntvl, int& indexPreKnownIntvls)
{
    int newIntvlFrom = 0;
    compareIntvlsMWF intvlCompareObj;
    compareIntvlsMWFMinHeap intvlCompareObjForHeap;

    if ((indexPreKnownIntvls == listOfPreKnownIntvls.size()) && listOfAdHocIntvls.empty())
        return -1;
    
    if (listOfAdHocIntvls.empty())
    {
        newIntvl = listOfPreKnownIntvls[indexPreKnownIntvls++];
        newIntvlFrom = 0;
    }
    else if (intvlCompareObj(listOfPreKnownIntvls[indexPreKnownIntvls], listOfAdHocIntvls.front()))
    {
        newIntvl = listOfPreKnownIntvls[indexPreKnownIntvls++];
        newIntvlFrom = 0;
    }
    else
    {
        newIntvl = listOfAdHocIntvls.front();
        //cout << "Intvl is to be used from Ad Hoc list" << endl;
        std::pop_heap(listOfAdHocIntvls.begin(), listOfAdHocIntvls.end(), intvlCompareObjForHeap); listOfAdHocIntvls.pop_back();
        newIntvlFrom = 1;
    }
/*    if ((indexPreKnownIntvls % 10000000) == 0)
    {
        cout << "static intvls covered= " << indexPreKnownIntvls << ", Size of adHoc " << listOfAdHocIntvls.size() << endl;
    }*/
    return newIntvlFrom;
}
