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
        listJourneys.resize(V);
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

void GraphDualCriteria::run_mwf()
{
    time_sum=0;
    
    for(int i = 0 ;i < sources.size() ;i ++)
    {
        initial_ds_ea();
        //earliest_arrival(sources[i]);
        //initial_ds_ea();earliest_arrival_pair(sources[i]);
        //initial_ds_ea();
//        shortest_path(sources[i]);
        //minWaitForemost(sources[i]);
 //       minWaitForemostPrioritized(sources[i]);
        //minWaitForemostPrioritizedNoSet(sources[i]);
        mwfStreamingIntvls(sources[i]);
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


/*************************
 Represent each interval as {u,v,u_nbrIndx,(s,e),lambda,prevJourneyIndex} in a list Intvls[]        //u_nbrIndx represents index of v in nbrs list of u.
 journey is represented as (a,w,prevNode,lastExpAt) where lastExpAt is (s,lambda) indexed by nbr or NULL.
 Sort Intvls[] in order of s+lambda
 The whole graph is also stored in verts[]. For each node, list of nbrs. And for each nbr, list of intervals.
 PQ mergeIntvls = {}            //PQ with arrival time (s+lambda) as the key. When there is a tie, s is the key.
 listJ[u] is a set sorted list (by arrival time) of journeys arrived at node u
 for each u { openIntvl[u] = null; listJ[u] = {} }
 listJ[s][0]={t_start,0}
 newJourney = listJ[s].last();
 for each nbr w of s
     nextIntvl = nextFn(s,w,t_start)
     If (nextIntvl.s >= newJourney.a)
         vert[v][w][nextIntvl].prevJourneyIndex = listJ[s].last();     //Can be stored as Id in listJ[v]
     else If (listJ[v].lastJ.arr > nextIntvl.s and < nextIntvl.e)
         e = nextIntvl.e;
         newIntvlCreated=(s,w,listJ[u].end().arr,e,nextIntvl.lambda);
         newIntvlCreated.prevJourney=listJ[s].last();
         Insert(mergeIntvls, {newIntvlCreated})
 numNodesRchd = 1;      //source has been rchd.
 newIntvl = min (top(Intvls, mergeIntvls));     //Comparison of arrTime= s+lambda and secondary start time.
 while ( (Intvls[] || mergeIntvls[]) && (numNodesRchd < reachable) )
 {
    currArrivalTime = newIntvl.arrTime
    numNewNodesRchd=0;
    while (newIntvl != NULL && currArrialTime == newIntvl.arrTime)
    {
        newIntvl=removeMin(top(Intvls,mergeIntvls)
        If (newIntvl.Id != -1)      //this means this intvl came from static set of intvls.
            prevJourneyIndex=vert[u][v][newIntvl.Id].prevJourneyIndex;
        else        //This means this intvl came from mergeIntvls list and was created by breaking up existing interval.
            prevJourneyIndex=newIntvl.prevJourneyIndex
        if (prevJuorneyIndex==-1)
        {
            If (listJ[u].last().arr <= newIntvl.strt)
                prevJourneyIndex=listJ[u].last()
            else
                search prevJourney in listJ[u] (journey with arrTm <= newIntvl.strtTime)
        }
        If no previousJourneyFound
            continue;
        If (listJ[u][prevJourneyIndex].expandedIntvl[v].lambda <  newIntvl.lambda or listJ[u][prevJourneyIndex].expandedIntvl[v] == null)       //Based on intvl dominance criteria. (lambda > lastLambda )
        {
            newJourney.a=newIntvl.s+newIntvl.lambda;            //This will be extension of the journey at prevJourneyIndex
            newJourney.w=listJ[u][prevJourneyIndex].w+newIntvl.s-listJ[u][prevJourneyIndex].arr;      //This will be extension of the journey at prevJourneyIndex
            listJ[u][prevJourneyIndex].expandedIntvl[v] = {newIntvl.s,newIntvl.lambda};
            If (newJourney not dominated by listJ[v].last() || listJ[v] == null)     //Dominance based on mwf criteria.
                If (mwf[v] == NULL)
                    mwf[v]=newJourney;
                    numNewNodesRchd++;
                else if ((newJourney.arrTm,newJourney.wtTime) < (mwf[v].arrTm,mwf[v].wtTime))
                    mwf[v]=newJourney;
                listJ[v].append(newJourney);
                for each nbr w of v
                    nextIntvl = nextFn(v,w,newJourney.arr)        //This operates on static input intvl temporal graph
                    If ( newJourney.a <= nextIntvl.s)
                        vert[v][w][nextIntvl].prevJourneyIndex = listJ[v].last();     //Can be stored as Id in listJ[v]
                        newJourney.nexFnInt[w]=nextIntvl
                    else If (newJourney.arr > nextIntvl.s and < nextIntvl.e)
                        e = nextIntvl.e;
                        newIntvlCreated=(v,w,v_nbrIndexFor_w,newJourney.arr,e,nextIntvl.lambda);        //Sub-interval of the static intvl created as all journeys exit at start of an interval.
                        newIntvlCreated.prevJourneyIndex= listJ[v].last();
                        Insert(mergeIntvls, {newIntvlCreated})
        }
        else
            ignore intvl for expansion.
        newIntvl = min (top(Intvls, mergeIntvls));     //Comparison of arrTime= s+lambda and secondary start time.
    }
    numNodesRchd+= numNewNodesRchd;
 }
 
 ****************************/

void GraphDualCriteria::mwfStreamingIntvls(int source)
{
    mwfJourney newJourney;
    intervalInfo newIntvl;
    int indexPreKnownIntvls = 0;
    int prevJourneyIndex = -1;
    Timer t;

    int nextIntvl = 0, nextTravelTime = 0; int newIntvlFrom = 0;
    int numNodesRchd = 0;

    tuple<int,int> nodesRchable_maxFmstTime = earliest_arrival_pair(source,1); //4092652
    int nodesReachable= get<0>(nodesRchable_maxFmstTime), maxFmstTime=get<1>(nodesRchable_maxFmstTime);
    cout << "Num Nodes Reachable: " << nodesReachable << endl;
    finalMWFJourneys.resize(V);
    for (int i=0; i< V; i++)
    {
        finalMWFJourneys[i].arrivalTime = infinity;
        finalMWFJourneys[i].wtTime = infinity;
        finalMWFJourneys[i].prevNode=-1;
    }
    compareIntvlsMWFMinHeap intvlCompareObjForHeap;
    newJourney.arrivalTime = t_start;
    newJourney.prevNode = -1; newJourney.wtTime=0;
    newJourney.expandedAt.resize(vertices[source].numNbrs);
    listJourneys[source].push_back(newJourney);        //known journeys so far at source.
    finalMWFJourneys[source] = newJourney;
    
    bit_queue closedNodes((int)vertices.size());
    closedNodes.queue_add_bit(source);
    numNodesRchd++;
    cout << "Total num Intvls: " << listOfPreKnownIntvls.size() << endl;
    
    t.start();
    //setupNewJourney(source, newJourney.arrivalTime);      //Commented as may not be reqd for csg graphs. TBD
    newIntvlFrom = getMinIntvl(newIntvl, indexPreKnownIntvls);
    int u=0,v=0,nbrIndex=0,intvlId=0;
    while ( ((!listOfAdHocIntvls.empty()) || (indexPreKnownIntvls < listOfPreKnownIntvls.size()))
           && (numNodesRchd < nodesReachable) && (newIntvlFrom != -1) )
    {
//        if (indexPreKnownIntvls % 1000 == 0)
//            cout << "Intervals processed: " << indexPreKnownIntvls << "Num nodes reached: " << numNodesRchd << endl;
        int currArrivalTime = newIntvl.intvlStart+newIntvl.lambda;
        int numNewNodesRchd=0;
        while ((currArrivalTime == newIntvl.intvlStart+newIntvl.lambda) && (newIntvlFrom != -1))
        {
            if (newIntvlFrom == 0)      //Advance to next preKnownIntvl
                indexPreKnownIntvls++;
            else if (newIntvlFrom == 1)                      //remove min of AdHoc Interval.
            {
                cout << "Intvl is to be used from Ad Hoc list" << endl;
                std::pop_heap(listOfAdHocIntvls.begin(), listOfAdHocIntvls.end(), intvlCompareObjForHeap); listOfAdHocIntvls.pop_back();
            }
            else
                break;
//            if (indexPreKnownIntvls % 1000 == 0)
//                cout << "Intervals processed: " << indexPreKnownIntvls << "Num nodes reached: " << numNodesRchd << endl;
            u = newIntvl.u; v=newIntvl.v; nbrIndex=newIntvl.nbrIndexFor_v; intvlId=newIntvl.intvlId;
            if (newIntvl.intvlId != -1)
            {
                prevJourneyIndex = vertices[u].neighbors[nbrIndex].edgeSchedules[intvlId].prevJourneyIndex;
            }
            else
                prevJourneyIndex = newIntvl.prevJourneyIndex;
            if (prevJourneyIndex == -1)
            {
                if (listJourneys[u].empty())
                {
                    newIntvlFrom = getMinIntvl(newIntvl, indexPreKnownIntvls);
                    continue;
                }
                else if (listJourneys[u].at(listJourneys[u].size()-1).arrivalTime <= newIntvl.intvlStart)
                    prevJourneyIndex = (int)listJourneys[u].size()-1;
                else
                    prevJourneyIndex = searchPrevJourney(u,newIntvl.intvlStart);
            }
            if (prevJourneyIndex == -1)
            {
                newIntvlFrom = getMinIntvl(newIntvl, indexPreKnownIntvls);  //Move to the next interval.
                continue;
            }
            int prevJourneyLastExtLmbda = get<1>(listJourneys[u][prevJourneyIndex].expandedAt[nbrIndex]);
            int whehterExpanded = get<2>(listJourneys[u][prevJourneyIndex].expandedAt[nbrIndex]);
            if ((prevJourneyLastExtLmbda < newIntvl.lambda) || (whehterExpanded == 0))
            {
                newJourney.arrivalTime = newIntvl.intvlStart+newIntvl.lambda;
                if ((u== source) && (listJourneys[u][prevJourneyIndex].arrivalTime == t_start))
                    newJourney.wtTime = 0;
                else
                    newJourney.wtTime = listJourneys[u][prevJourneyIndex].wtTime + newIntvl.intvlStart-listJourneys[u][prevJourneyIndex].arrivalTime;
                newJourney.prevNode=u;
                listJourneys[u][prevJourneyIndex].expandedAt[nbrIndex] = make_tuple(newIntvl.intvlStart,newIntvl.lambda,1);

                newJourney.expandedAt.clear();
                newJourney.expandedAt.resize(vertices[v].numNbrs);
                if (!listJourneys[v].empty())
                {
                    tuple<int,int> j1; tuple<int,int> j2;
                    get<0>(j1) = listJourneys[v][listJourneys[v].size()-1].arrivalTime;
                    get<1>(j1) = listJourneys[v][listJourneys[v].size()-1].wtTime;
                    get<0>(j2) = newJourney.arrivalTime;
                    get<1>(j2) = newJourney.wtTime;
                    if (!checkDominance(j1, j2))
                    {
                        if (get<0>(j1) < get<0>(j2))
                            listJourneys[v].push_back(newJourney);
                        else
                            listJourneys[v][listJourneys[v].size()-1] = newJourney; //Replace last journey as arrTm is same and wt time more.
                    }
                    else        //The new journey from this interval is dominated so ignore it.
                    {
                        newIntvlFrom = getMinIntvl(newIntvl, indexPreKnownIntvls);
                        continue;
                    }
                }
                else
                {
                    numNewNodesRchd++;
                    listJourneys[v].push_back(newJourney);
                    finalMWFJourneys[v] = newJourney;
                    closedNodes.queue_add_bit(v);
                }
                if ((finalMWFJourneys[v].arrivalTime == newJourney.arrivalTime) &&
                    (finalMWFJourneys[v].wtTime > newJourney.wtTime))
                    finalMWFJourneys[v] = newJourney;
                
                //setupNewJourney(v, newJourney.arrivalTime);
            }
            newIntvlFrom = getMinIntvl(newIntvl, indexPreKnownIntvls);
        }
        numNodesRchd += numNewNodesRchd;
    }
    t.stop();
    time_sum += t.GetRuntime();
    printmwfResultsTest2(source);
}

void GraphDualCriteria::setupNewJourney(int v, int arrivalTime)
{
    int nextTravelTime=0,nextIntvl=0;
    compareIntvlsMWFMinHeap intvlCompareObjForHeap;
    for (int i=0; i< vertices[v].numNbrs;i++)
    {
        listJourneys[v][(int)listJourneys[v].size()-1].expandedAt[i] = make_tuple(-1,-1,0); //starting journey isnt expanded yet
        nextTravelTime = earliestUseEdgeAfterT(v, vertices[v].neighbors[i], arrivalTime, nextIntvl);
        if ((nextIntvl == -1) || (nextTravelTime >= infinity))
            continue;
        if (arrivalTime <= vertices[v].neighbors[i].edgeSchedules[nextIntvl].intvlStart)
        {
            vertices[v].neighbors[i].edgeSchedules[nextIntvl].prevJourneyIndex=(int)listJourneys[v].size()-1;
        }
        else if ( (vertices[v].neighbors[i].edgeSchedules[nextIntvl].intvlStart < arrivalTime)
                 && (arrivalTime <= vertices[v].neighbors[i].edgeSchedules[nextIntvl].intvlEnd) )
        {
            intervalInfo adHocIntvl;
            adHocIntvl.intvlEnd = vertices[v].neighbors[i].edgeSchedules[nextIntvl].intvlEnd;
            adHocIntvl.intvlStart = nextTravelTime; adHocIntvl.intvlId=-1;
            adHocIntvl.lambda= vertices[v].neighbors[i].edgeSchedules[nextIntvl].traveTime;
            adHocIntvl.u = v; adHocIntvl.nbrIndexFor_v=i; adHocIntvl.v = vertices[v].neighbors[i].nbrId;
            adHocIntvl.prevJourneyIndex = (int)listJourneys[v].size()-1;
            listOfAdHocIntvls.push_back(adHocIntvl); std::push_heap(listOfAdHocIntvls.begin(), listOfAdHocIntvls.end(), intvlCompareObjForHeap);
        }
    }
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
            rv++;
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
    
    if (high >=1)
    {
        if (nodeJourneys->at(high-1).arrivalTime <= beforeTime)
            return (high-1);
    }
    else if (high==0)
        return -1;
        
    
    while (!journeyFound)
    {
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

int GraphDualCriteria::getMinIntvl(intervalInfo &newIntvl, int indexPreKnownIntvls)
{
    int newIntvlFrom = 0;
    compareIntvlsMWF intvlCompareObj;
    if ((indexPreKnownIntvls == listOfPreKnownIntvls.size()) && listOfAdHocIntvls.empty())
        return -1;
    
    if (listOfAdHocIntvls.empty())
    {
        newIntvl = listOfPreKnownIntvls[indexPreKnownIntvls];
        newIntvlFrom = 0;
    }
    else if (intvlCompareObj(listOfPreKnownIntvls[indexPreKnownIntvls], listOfAdHocIntvls.front()))
    {
        newIntvl = listOfPreKnownIntvls[indexPreKnownIntvls];
        newIntvlFrom = 0;
    }
    else
    {
        newIntvl = listOfAdHocIntvls.front();
        newIntvlFrom = 1;
    }
    return newIntvlFrom;
}
