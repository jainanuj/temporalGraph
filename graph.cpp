//
//  graph.cpp
//  dummylib
//
//  Created by Anuj Jain on 11/20/20.
//

#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include "cppfibonacci/fibonacci.hpp"

Graph::Graph(const char* filePath)      //Reads the edges of the graph in interval format. Edges will be specified as u, v, intvlCount, [intvls]. Each intvl is specified as (start, end, travelTime). End is last time instant at which edge can be used.
{
    FILE* file;
    int x;
    int u, lambda, intvlStart, intvlEnd, adjustedEnd;
    Nbrs nbr;
    Node node;
    
    file = fopen(filePath,"r");
    
    x=fscanf(file, "%d %d",&V, &dynamic_E);
    vertices.resize(V);
    for (int i = 0; i < V; i++)
    {
        vertices[i].numNbrs = 0;
        vertices[i].nodeId = i;
    }
    for(int i = 0; i < dynamic_E; i ++)
    {
        x=fscanf(file, "%d %d %d", &u, &nbr.nbrId, &nbr.numIntvls); //(u, v, lamda, numIntvls)
        vertices[u].numNbrs++;
        nbr.edgeSchedules.resize(nbr.numIntvls);
        int intvls = 0; int intvlRead;
        for (intvlRead = 0; intvlRead <  nbr.numIntvls; intvlRead++)
        {
            x=fscanf(file, "%d %d %d", &intvlStart, &intvlEnd, &lambda);
            //adjustedEnd = intvlEnd - nbr.nbr_lambda;      This is for input model where end is time till edge is active.
            adjustedEnd = intvlEnd;     //This is for input model where end is last instant at which end can be used.
            //nbr.edgeSchedules[intvls].adjustedEnd = nbr.edgeSchedules[intvls].intvlEnd - nbr.nbr_lambda;
            if (adjustedEnd >= intvlStart)       //else ignore this intvl because it can't really be used.
            {
                nbr.edgeSchedules[intvls].intvlStart = intvlStart; nbr.edgeSchedules[intvls].intvlEnd = intvlEnd;
                nbr.edgeSchedules[intvls].adjustedEnd = adjustedEnd; nbr.edgeSchedules[intvls].traveTime = lambda;
                intvls++;
            }
        }
        nbr.adjustedNumIntvls = intvls;
        vertices[u].neighbors.push_back(nbr);
    }
    
    arr_time.resize(V);
    father.resize(V);
    f_time.resize(V);
    
    fclose(file);
}

void Graph::wuGraph(const char* filePath, int noL)
{
    FILE* file = fopen(filePath,"r");
    int x, wuEdges, vertices;
    x=fscanf(file, "%d %d",&vertices, &wuEdges);
    vector<tuple<int, int, int, int>> inputRows;
    compareTuple tupleComp;
    int u, v, t, lambda;
    inputRows.resize(wuEdges);
    const char* intermediateFile = "/Users/anujjain/research/temporalGraph/WuTemporalGraph/tempath/intermediate.txt";
    
    lambda = 1;
    for(int i = 0; i < wuEdges; i ++)
    {
        if (noL)
        {
            x=fscanf(file, "%d %d %d",&u, &v, &t);
        }
        else
            x=fscanf(file, "%d %d %d %d",&u, &v, &t, &lambda);
        inputRows[i] = std::make_tuple(u, v, t, lambda);
    }
    fclose(file);
    std::sort<vector<tuple<int, int, int, int>>::iterator, compareTuple>(inputRows.begin(), inputRows.end(), tupleComp);

    FILE* fileToWrite = fopen(intermediateFile, "w");
    fprintf(fileToWrite, "%d %d\n", vertices, wuEdges);
    int currentRunningU, currentRunningV, intvlCount = 0, numEdgesCollapsed = -1;
    vector<tuple<int, int>> currentIntvls;      //It is a vector of tuples. Each tuple consists of t, lambda.
    currentRunningU = std::get<0>(inputRows[0]); currentRunningV = std::get<1>(inputRows[0]);
    for (int i = 0; i < inputRows.size(); i++)
    {
        u = std::get<0>(inputRows[i]); v = std::get<1>(inputRows[i]);
        if ( (u == currentRunningU) && (v == currentRunningV) )
        {
            numEdgesCollapsed++;
            intvlCount++;
            currentIntvls.push_back(std::make_tuple(std::get<2>(inputRows[i]), std::get<3>(inputRows[i]) ) );       //Interval in intermediate file is t, lambda.
        }
        else
        {
            if (currentRunningU != currentRunningV)
            {
                fprintf(fileToWrite, "%d %d %d", currentRunningU, currentRunningV, intvlCount);     //Output format is u, v, intvlCount, (intervals). Intervals are of the form t, lambda.
                for (int intC = 0; intC < intvlCount; intC++)
                    fprintf(fileToWrite, " %d %d", std::get<0>(currentIntvls[intC]), std::get<1>(currentIntvls[intC]) );
                fprintf(fileToWrite, "\n");
            }
            else
                numEdgesCollapsed++;        //Completely ignore this edge as start and end node are same so it's a loop.
            currentIntvls.clear();      //Reset intvls for the new edge.
            intvlCount = 0;
            //Prepare to start looking at next edge.
            currentRunningU = u; currentRunningV = v;
            intvlCount++;
            currentIntvls.push_back(std::make_tuple(std::get<2>(inputRows[i]), std::get<3>(inputRows[i]) ) );
        }
    }
    //now print out the consolidated last edge.
    if (currentRunningU != currentRunningV)
    {
        fprintf(fileToWrite, "%d %d %d", currentRunningU, currentRunningV, intvlCount);
        for (int intC = 0; intC < intvlCount; intC++)
            fprintf(fileToWrite, " %d %d", std::get<0>(currentIntvls[intC]), std::get<1>(currentIntvls[intC]) );
        fprintf(fileToWrite, "\n");
    }
    else
        numEdgesCollapsed++;        //Completely ignore this edge as start and end node are same so it's a loop.
    fclose(fileToWrite);
    
    fileToWrite = fopen(intermediateFile, "r+b");
    string vertString = to_string(vertices);
    string edgesBefCollaps = to_string(wuEdges);
    string edgesAftCollaps = to_string(wuEdges - numEdgesCollapsed);
    int bufSpaceToFill = (int)(edgesBefCollaps.length() - edgesAftCollaps.length());
    for (int i = 0; i < bufSpaceToFill; i++ )
        edgesAftCollaps += " ";                     //Pad spaces if necessary to make sure that the older value of numEdges is completely removed from the file.
    int seek_pos = (int)vertString.length() + 1;        //+1 to account for the space between vertices and edges.
    fseek(fileToWrite, seek_pos, SEEK_SET);
    fputs(edgesAftCollaps.c_str(), fileToWrite);
    fclose(fileToWrite);
    
    collapseIntervalsWriteOuput(intermediateFile);      //Write the final output. Each Edge is in format ==> (u, v, intvlCount, [s, e, d])
}

void Graph::collapseIntervalsWriteOuput(const char *filePath)
{
    const char* opFile = "/Users/anujjain/research/temporalGraph/WuTemporalGraph/tempath/xuanOp.txt";
    FILE *intmdFile = fopen(filePath, "r");
    FILE *outputFile = fopen(opFile, "w");
    int x, xuanEdges, vertices, nextWuContact, intvlStart, intvlEnd, xuanIntvlCount, travelTime, nextTravelTime, adjustedIntvlCount;
    
    x=fscanf(intmdFile, "%d %d",&vertices, &xuanEdges);
    x=fprintf(outputFile, "%d %d\n", vertices, xuanEdges);
    int u, v, wuContactsCount;
    
    vector<tuple<int, int, int>> xuanIntvlsVector;      //Interval vector is (s, e, d). s=start time. e=end time. d=travel time or lambda.
    vector<tuple<int, int, int>> adjustedIntervals;
    
    //
    for(int i = 0; i < xuanEdges; i ++)
    {
        x=fscanf(intmdFile, "%d %d %d",&u, &v, &wuContactsCount);
        x=fscanf(intmdFile, "%d %d", &intvlStart, &travelTime);         //Starting intvl.
        intvlEnd = intvlStart; xuanIntvlCount = 1;
        xuanIntvlsVector.clear();
        for (int intC = 1; intC < wuContactsCount; intC++)
        {
            x=fscanf(intmdFile, "%d %d", &nextWuContact, &nextTravelTime);
            if ( ( nextWuContact == (intvlEnd+1)) && (nextTravelTime == travelTime) )       //Can merge the two contact instances into one interval.
                intvlEnd = nextWuContact;
            else
            {
                xuanIntvlCount++;
                xuanIntvlsVector.push_back(std::make_tuple(intvlStart, intvlEnd, travelTime));  //interval is a triplet - s, e, d. s=start, e=end, d=traveltime.
                intvlStart = nextWuContact; intvlEnd = intvlStart; travelTime = nextTravelTime;
            }
        }
        xuanIntvlsVector.push_back(std::make_tuple(intvlStart, intvlEnd, travelTime));      //Add the last intvl as well.
        adjustedIntervals.clear();
        adjustedIntervals = adjustSlowIntervals(xuanIntvlsVector);
        
        //Now write out the edge to the xuanOutputFile.
        adjustedIntvlCount = (int)adjustedIntervals.size();
        fprintf(outputFile, "%d %d %d", u, v, adjustedIntvlCount);
        for (int xuanIntvls = adjustedIntvlCount-1; xuanIntvls >= 0; xuanIntvls--)
        {
            fprintf(outputFile, " %d %d %d", std::get<0>(adjustedIntervals[xuanIntvls]), std::get<1>(adjustedIntervals[xuanIntvls]), std::get<2>(adjustedIntervals[xuanIntvls]));
        }
        fprintf(outputFile, "\n");
    }
    fclose(intmdFile);
    fclose(outputFile);
}

vector<tuple<int, int, int>> Graph::adjustSlowIntervals(std::vector<tuple<int, int, int> > &intervalVector)
{
    int len = (int) intervalVector.size();
    int j = 0, s, e, d, s2, e2, d2, newE;
    vector<tuple<int, int, int>> newIntvlsVector;      //Interval vector is (s, e, d). s=start time. e=end time. d=travel time or lambda.
    newIntvlsVector.push_back(intervalVector[len-1]); j++;

    for (int i = len - 2; i >= 0; i--)
    {
        s = std::get<0>(intervalVector[i]); e = std::get<1>(intervalVector[i]); d = std::get<2>(intervalVector[i]);
        s2 = std::get<0>(newIntvlsVector[j-1]); e2 = std::get<1>(newIntvlsVector[j-1]); d2 = std::get<2>(newIntvlsVector[j-1]);
        if ((e + d) > (s2 + d2))
        {
            newE = s2 + d2 - d;
            if (newE >= s)      //push the adjusted interval as s, newE, d. Else skip this interval.
            {
                newIntvlsVector.push_back(std::make_tuple(s, newE, d)); j++;
            }
        }
        else
        {
            newIntvlsVector.push_back(intervalVector[i]); j++;
        }
    }
    return newIntvlsVector;
}

int Graph::earliestUseEdgeAfterT(int u, Nbrs& v, int t, int &intvlID)
{
    int low = 0, high = v.adjustedNumIntvls-1, mid = high/2, intvlFound = 0, retVal = -1;
    intvlID = -1;
    
    while (!intvlFound)
    {
        if (v.edgeSchedules[mid].isValidFor(t))
        {
            if (t>= v.edgeSchedules[mid].intvlStart)     //t is within the intvl.
            {
                intvlFound = 1;
                retVal = t;
                intvlID = mid;
            }
            else if (t < v.edgeSchedules[mid].intvlStart)
            {
                if ((high - low) > 1)
                {
                    high = mid;
                    mid = (high + low)/2;
                }
                else if (mid == low)
                {
                    intvlFound = 1;
                    retVal = v.edgeSchedules[mid].intvlStart;   //t is lower than the low interval. safe to assume start of this intvl is earliest x-mission time.
                    intvlID = mid;
                }
                else        //Difference in low and high is 1 and mid was high. This means, t occurs before start of high intvl. Lets check if low satisfies t, if yes return low, else return mid.
                {
                    if (v.edgeSchedules[low].isValidFor(t))
                    {
                        intvlFound = 1;
                        retVal = (t > v.edgeSchedules[low].intvlStart) ? t : v.edgeSchedules[low].intvlStart;
                        intvlID = low;
                    }
                    else
                    {
                        intvlFound = 1;
                        retVal = v.edgeSchedules[mid].intvlStart;
                        intvlID = mid;
                    }
                }
            }
        }
        else        //mid was not a valid interval for transmission after t. This means valid interval is later than mid.
        {
            if ((high - low) > 1)
            {
                low = mid+1;        //mid can be safely excluded cuz it doesn't satisfy transmission after t.
                mid = (low + high)/2;
            }
            else        //High and low are adjacent or same. The mid didn't contain valid intvl. Only possiblity is to check the intvl at 'high' position.
            {
                if (mid == low)
                {
                    if (v.edgeSchedules[high].isValidFor(t))
                    {
                        intvlFound = 1;
                        retVal = (t > v.edgeSchedules[high].intvlStart) ? t : v.edgeSchedules[high].intvlStart;
                        intvlID = high;
                    }
                    else
                    {
                        intvlFound = 1;
                        retVal = -1;            //There is no valid interval.
                        intvlID = -1;
                    }
                }
                else
                {
                    intvlFound = 1;
                    retVal = -1;            //There is no valid interval.
                    intvlID = -1;
                }
            }
        }
    }
    return retVal;
}


void Graph::initial_query()
{
    t_start = 0;
    t_end = infinity;
    
    //int s;
    sources.push_back(1);
/*    for(int i = 0 ;i < 10 ;i ++)
    {
        s=rand()%V;
        sources.push_back(s);
    }*/

}


void Graph::initial_query(const char* filePath)
{
    t_start = 0;
    t_end = infinity;
    
    FILE* file = fopen(filePath,"r");
    int s, x;
    for(int i = 0 ;i < 100 ;i ++)
    {
        x=fscanf(file,"%d",&s);
    //    int y;
    //    x=fscanf(file,"%d%d",&s, &y);
        sources.push_back(s);
    }

}


void Graph::initial_ds_ea()
{
    for(int i=0; i<V; i++){
        arr_time[i]= infinity;
        father[i] = -1;
    }

}

void Graph::initial_ds_ld()
{
    for(int i=0; i<V; i++){
        arr_time[i]= -1;
    }

}

void Graph::initial_ds_f()
{
    for(int i=0; i<V; i++){
        f_time[i] = infinity;
    }

}

void Graph::initial_ds_s()
{
    for(int i=0; i<V; i++){
        f_time[i] = infinity;
    }

}


void Graph::run_earliest_arrival()
{
    time_sum=0;
    
    for(int i = 0 ;i < sources.size() ;i ++)
    {
        //initial_ds_ea();
        //earliest_arrival(sources[i]);
        //initial_ds_ea();earliest_arrival_pair(sources[i]);
        initial_ds_ea();earliest_arrival_fibo(sources[i]);
    }
    
    print_avg_time();
}

void Graph::earliest_arrival(int source)
{
    Timer t;
    t.start();
    Node *u;
    int nodeID;
    int tDepart, intvlID = -1;;
    vector<int> minHeap;    //Fix the pq to make it a min heap instead of max heap.
    nodeComparison heapCompFn(this);

    bit_queue nodesInHeap((int)vertices.size());
    bit_queue closedNodes((int)vertices.size());
    
        
    arr_time[source]=t_start;
    father[source] = source;
    minHeap.push_back(source);
    
    while (!minHeap.empty())
    {
        nodeID = minHeap.front();
        std::pop_heap<vector<int>::iterator, nodeComparison>(minHeap.begin(), minHeap.end(), heapCompFn); minHeap.pop_back();      //moves the min element to last and then removes it from heap.
        if (closedNodes.check_bit_obj_present(nodeID))
            continue;
        closedNodes.queue_add_bit(nodeID);         //Add it to the closedNodes.
        u = &vertices[nodeID];
        for (int i = 0; i < u->numNbrs; i++)
        {
            if (! closedNodes.check_bit_obj_present( u->neighbors[i].nbrId))
            {
                tDepart = earliestUseEdgeAfterT(nodeID, u->neighbors[i], arr_time[nodeID], intvlID);
                if ( (tDepart + u->neighbors[i].edgeSchedules[intvlID].traveTime) < arr_time[u->neighbors[i].nbrId])
                {
                    arr_time[u->neighbors[i].nbrId] = tDepart + u->neighbors[i].edgeSchedules[intvlID].traveTime;
                    father[u->neighbors[i].nbrId] = u->nodeId;
                    //if (! nodesInHeap.check_bit_obj_present(u->neighbors[i].nbrId) )
                    minHeap.push_back(u->neighbors[i].nbrId);
                    std::push_heap(minHeap.begin(), minHeap.end(), heapCompFn);
                }
            }
        }
    }
    t.stop();
    printResults(source);
    time_sum += t.GetRuntime();
}

void Graph::earliest_arrival_pair(int source)
{
    Timer t;
    t.start();
    Node *u;
    int nodeID;
    int tDepart, intvlID = -1;
    vector<pair<int, int>> minHeap;    //Fix the pq to make it a min heap instead of max heap.
    nodeComparison2 heapCompFn;

    bit_queue closedNodes((int)vertices.size());
        
    arr_time[source]=t_start;
    father[source] = source;
    minHeap.push_back(std::make_pair(arr_time[source], source));
    
    while (!minHeap.empty())
    {
        nodeID = minHeap.front().second;
        std::pop_heap<vector<pair<int, int>>::iterator, nodeComparison2>(minHeap.begin(), minHeap.end(), heapCompFn); minHeap.pop_back();      //moves the min element to last and then removes it from heap.
        if (closedNodes.check_bit_obj_present(nodeID))      //This was already closed.
            continue;
        closedNodes.queue_add_bit(nodeID);         //Add it to the closedNodes.
        u = &vertices[nodeID];
        for (int i = 0; i < u->numNbrs; i++)
        {
            if (! closedNodes.check_bit_obj_present( u->neighbors[i].nbrId))
            {
                tDepart = earliestUseEdgeAfterT(nodeID, u->neighbors[i], arr_time[nodeID], intvlID );
                if ( (tDepart + u->neighbors[i].edgeSchedules[intvlID].traveTime) < arr_time[u->neighbors[i].nbrId])
                {
                    arr_time[u->neighbors[i].nbrId] = tDepart + u->neighbors[i].edgeSchedules[intvlID].traveTime;
                    father[u->neighbors[i].nbrId] = u->nodeId;
                    minHeap.push_back(std::make_pair(arr_time[u->neighbors[i].nbrId], u->neighbors[i].nbrId));
                    std::push_heap(minHeap.begin(), minHeap.end(), heapCompFn);
                }
            }
        }
    }
    t.stop();
    printResults(source);
    time_sum += t.GetRuntime();
}

    

void Graph::earliest_arrival_fibo(int source)
{
    Timer t;
    t.start();
    Node *u;
    int nodeID, nodeArrTime;
    int tDepart, intvlID = -1;
    fibonacci_heap<int, int> minHeap;    //Fix the pq to make it a min heap instead of max heap.
    vector<fibonacci_heap<int, int>::node> fibNodes;
    fibNodes.resize(vertices.size());

    bit_queue nodesInHeap((int)vertices.size());
    bit_queue closedNodes((int)vertices.size());
            
    arr_time[source]=t_start;
    father[source] = source;
    fibNodes[source] = minHeap.insert(arr_time[source], source);
    
    while (minHeap.size())
    {
        nodeID = minHeap.top().data();
        nodeArrTime = minHeap.top().key();
        minHeap.remove();
        if (closedNodes.check_bit_obj_present(nodeID))      //This shouldn't happen.
            continue;
        closedNodes.queue_add_bit(nodeID);         //Add it to the closedNodes.
        u = &vertices[nodeID];
        for (int i = 0; i < u->numNbrs; i++)
        {
            if (! closedNodes.check_bit_obj_present( u->neighbors[i].nbrId))
            {
                tDepart = earliestUseEdgeAfterT(nodeID, u->neighbors[i], arr_time[nodeID], intvlID );
                if ( (tDepart + u->neighbors[i].edgeSchedules[intvlID].traveTime) < arr_time[u->neighbors[i].nbrId])
                {
                    arr_time[u->neighbors[i].nbrId] = tDepart + u->neighbors[i].edgeSchedules[intvlID].traveTime;
                    father[u->neighbors[i].nbrId] = u->nodeId;
                    if (nodesInHeap.check_bit_obj_present(u->neighbors[i].nbrId))
                        minHeap.decrease_key(fibNodes[u->neighbors[i].nbrId], arr_time[u->neighbors[i].nbrId]);
                    else
                        fibNodes[u->neighbors[i].nbrId] = minHeap.insert(arr_time[u->neighbors[i].nbrId], u->neighbors[i].nbrId);
                }
            }
        }
    }
    t.stop();
    printResults(source);
    time_sum += t.GetRuntime();
}


void Graph::printResults(int source)
{
    int x;
    cout << "Source is: " <<source << "\n";
    for (int i = 0; i < arr_time.size(); i++)
    {
        cout << "Arrival time at node " << i << " is: " << arr_time[i] << "\n";
        cout << "Path in rev is: ";
        x = i;
        while (father[x] != x)
        {
            cout << x <<", ";
            x = father[x];
        }
        cout << x << "\n";
    }
    
}

void Graph::run_earliest_arrival(const char* filePath)
{
    time_sum=0;
    
    FILE* file = fopen(filePath,"w");
    
    for(int i = 0 ;i < sources.size() ;i ++)
    {
        initial_ds_ea();
        earliest_arrival(sources[i], file);
    }
    
    fclose(file);
    
    print_avg_time();

}

void Graph::earliest_arrival(int source, FILE * file)
{
    Timer t;
    t.start();
    t.stop();
    time_sum += t.GetRuntime();
    
    print_result(source, arr_time, file);
}


void Graph::run_shortest()
{
    time_sum=0;
    
    for(int i = 0 ;i < sources.size() ;i ++)
    {
        //initial_ds_ea();
        //earliest_arrival(sources[i]);
        //initial_ds_ea();earliest_arrival_pair(sources[i]);
        initial_ds_ea(); shortest_path(sources[i]);
    }
    
//    print_avg_time();
}

void Graph::shortest_path(int source)
{
    vector<vector<incrementalJourney>> allHopJourneys; //At each hop there is a vector of foremost inremental journeys, discovered in that hop.
    allHopJourneys.resize(V);   //Max number of hops is the number of nodes.
    vector<int> newArrivalsInHop;
    newArrivalsInHop.resize(V);
    vector<int> currentTimeArrival;                //Arrival time at this node.
    int intvlID = -1, departTime, newArrivTime, hopCount = 0, nextNode;
    currentTimeArrival.resize(V);
    currentTimeArrival.assign(V, infinity);
    currentTimeArrival[source] = 0;
    incrementalJourney journeyIncrement(source, 0, -1, -1, -1, -1, -1);   //current node, arrTime, prevNode, prevEdgeId, prevIntvlId, pDep, indexInPrevHopV
    allHopJourneys[hopCount].push_back(journeyIncrement);

    shortestJourneys.resize(V);                     //Shortest Journeys for all nodes needs to be found.
    shortestJourneys[source].rPath.push_back(std::make_tuple(source, -1, -1));       //start node, nbr ID.
    shortestJourneys[source].sigmaSchedule.push_back(std::make_tuple(0, -1, -1));        //depart time, IntvlId
    int  numNodesSeen = 1, numNewNodesInCurrentHop = 1;     //number of new foremost nodes in current hop.
    while ( (hopCount < (V-1)) && (numNewNodesInCurrentHop > 0) && (numNodesSeen < V))
    {
        hopCount++;
        numNewNodesInCurrentHop = 0;    //number of foremost nodes
        newArrivalsInHop.assign(V, -1);         //this array is needed to avoid duplicates in the current hop.
        for (int i = 0; i < allHopJourneys[hopCount -1].size(); i++)        //extend all journeys from prev hop
        {
            incrementalJourney* incToExtend = &(allHopJourneys[hopCount-1][i]);
            int node = incToExtend->currentNodeID;
            int time = incToExtend->currentArrivalTime;
            
            for (int j = 0; j < vertices[node].numNbrs; j++)
            {
                departTime = earliestUseEdgeAfterT(node, vertices[node].neighbors[j], time, intvlID);
                if (departTime == -1)
                    continue;
                nextNode = vertices[node].neighbors[j].nbrId;
                newArrivTime = departTime + vertices[node].neighbors[j].edgeSchedules[intvlID].traveTime;
                if (newArrivTime < currentTimeArrival[nextNode])     //New foremost time to nextNode Discovered.
                {
                    if (newArrivalsInHop[nextNode] == -1)
                    {
                        incrementalJourney newIncrement(nextNode, newArrivTime, node, j, intvlID, departTime, i);
                        allHopJourneys[hopCount].push_back(newIncrement);
                        newArrivalsInHop[nextNode] = numNewNodesInCurrentHop++;
                    }
                    else
                    {
                        incrementalJourney *pNewIncrement = &(allHopJourneys[hopCount][newArrivalsInHop[nextNode]]);
                        pNewIncrement->currentArrivalTime = newArrivTime; pNewIncrement->prevNodeID = node; pNewIncrement->prevEdgeID = j;
                        pNewIncrement->prevIntvlID = intvlID; pNewIncrement->prevDepartureTime = departTime; pNewIncrement->indexPrevIncrement = i;
                    }
                    if (currentTimeArrival[nextNode] == infinity)      //This is shortest journey to nextNode as no path found yet.
                    {
                        numNodesSeen++;
                        int incrIndex = newArrivalsInHop[nextNode];
                        for (int shortestHopCount = hopCount; shortestHopCount > 0; shortestHopCount--)  //BuildShortestJourney for nextNode from Incremental Vectors.
                        {
                            incrementalJourney *pJourneyPath = &(allHopJourneys[shortestHopCount][incrIndex]);    //start node, it's nbr id.
                            int currentNode = pJourneyPath->currentNodeID;int startNode = pJourneyPath->prevNodeID; int edgeId = pJourneyPath->prevEdgeID;
                            int currentArrivalTime = pJourneyPath->currentArrivalTime;int depart = pJourneyPath->prevDepartureTime; int intvlUsed = pJourneyPath->prevIntvlID;
                            shortestJourneys[nextNode].rPath.push_back(std::make_tuple(currentNode, startNode, edgeId));
                            shortestJourneys[nextNode].sigmaSchedule.push_back(std::make_tuple(currentArrivalTime, depart, intvlUsed));    //dep. time on edge (node, nextnode), intvl used.
                            incrIndex = pJourneyPath->indexPrevIncrement;
                        }
                    }
                    currentTimeArrival[nextNode] = newArrivTime;
                }       //found a better arrival time.
            }       //Looking at all the nbrs of a node to extend in current hop count.
        }       //for loop extending journey from hopCount-1.
    }       //While loop to cover all hop counts.
    print_shortest_paths(source);
}

void Graph::print_shortest_paths(int source)
{
    cout << "Shortest paths from the source node " << source << " are as follows:\n";
    
    for (int i = 0; i < V; i++)
    {
        if ( (i == source) || (shortestJourneys[i].rPath.size() > 0) )
            cout << "\nShortest path to " << i << " has hop count=" << shortestJourneys[i].rPath.size() << ". Path is: \n";
        else
            continue;
        if (i == source)
        {
            cout << "(" << source << ")";
            continue;
        }
        for (int hopCount = (int)(shortestJourneys[i].rPath.size()) - 1; hopCount >= 0; hopCount--)
        {
            std::tuple<int, int, int> *pCurrentPrevNodesEdge = &(shortestJourneys[i].rPath[hopCount]);
            std::tuple<int, int, int> *pArriveDepartIntvlUsed = &(shortestJourneys[i].sigmaSchedule[hopCount]);
            int toNode = std::get<0>(*pCurrentPrevNodesEdge);
            int prevNode = std::get<1>(*pCurrentPrevNodesEdge);

            int arrivalTime = std::get<0>(*pArriveDepartIntvlUsed);
            int departTime = std::get<1>(*pArriveDepartIntvlUsed);

            cout << "(" << (prevNode) << "-->" << toNode << ")"
                << ", (dep at: " << departTime << ", arr at:" << arrivalTime << ");;";
        }
    }
    
}








void Graph::print_result(const int source, const vector<int>& t_time, FILE * file)
{
    
    for(int i = 0 ;i < V ;i ++)
    {
        if(t_time[i]!=infinity){
            fprintf(file, "%d %d %d %d %d\n", source, i, t_start, t_end, t_time[i]);
        }
    }

}


void Graph::print_result_ld(const int source, const vector<int>& t_time, FILE * file)
{
    
    for(int i = 0 ;i < V ;i ++)
    {
        if(t_time[i]!=-1){
            fprintf(file, "%d %d %d %d %d\n", source, i, t_start, t_end, t_time[i]);
        }
    }

}



void Graph::print_avg_time()
{
    cout<<"Average time: " << time_sum/100 <<endl;
}

void Graph::print_avg_time(const char* filePath1, const char* filePath2)
{
    FILE* file = fopen(filePath2,"a");
    
    fprintf(file, "%s\t%f\n", filePath1, time_sum/100);;
    fclose(file);
}
