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
//#include "cppfibonacci/fibonacci.hpp"
#include <fstream>
#include "fibheap.h"
#include "regHeap.hpp"
#include <string>

//#define MEASUREHEAP_DET 1
#define __TEST__
#define RUN_COUNT 100
#define HOPCHECK 1000
#define JOURCHECK 1000000


#define sourceId 0

string Graph::getOuptuFile(string inFile,string appendRes)
{
    string outFile = inFile;
    unsigned int len = (int)appendRes.length();
    outFile.replace(outFile.length()-4, len, appendRes);
    return outFile;
}

Graph::Graph(const char* filePath, int contactSeq = 0, const char * option="any")      //Reads the edges of the graph in interval format. Edges will be specified as u, v, intvlCount, [intvls]. Each intvl is specified as (start, end, travelTime). End is last time instant at which edge can be used.
{
    FILE* file;
    int x;
    int u, lambda, intvlStart, intvlEnd, adjustedEnd;
    Nbrs nbr;
    Node node;
    
    string inFile(filePath);
    string opFile("./");
    
    std::size_t lastSlash = inFile.find_last_of('/');
    
    opFile = opFile+inFile.substr(lastSlash+1);
    
    earliestResults= getOuptuFile(opFile, "_earliestResults.txt");
    mhfResults= getOuptuFile(opFile, "_mhfResults.txt");
    shortestResults= getOuptuFile(opFile, "_shortestResults.txt");
    hbhShrtstResults= getOuptuFile(opFile, "_hbhshrtstResuts.txt");
    minHeapMonitor= getOuptuFile(opFile, "_minHeapMon.txt");
    
    file = fopen(filePath,"r");
    
    x=fscanf(file, "%d %d",&V, &dynamic_E);
    vertices.resize(V);
    for (int i = 0; i < V; i++)
    {
        vertices[i].numNbrs = 0;
        vertices[i].nodeId = i;
        vertices[i].inDegree = 0;
    }
    for(int i = 0; i < dynamic_E; i ++)
    {
        x=fscanf(file, "%d %d %d", &u, &nbr.nbrId, &nbr.numIntvls); //(u, v, numIntvls)
        vertices[u].numNbrs++;
        vertices[nbr.nbrId].inDegree += nbr.numIntvls;
        nbr.edgeSchedules.resize(nbr.numIntvls);
        int intvls = 0; int intvlRead;
        for (intvlRead = 0; intvlRead <  nbr.numIntvls; intvlRead++)
        {
            if (contactSeq == 0)
                x=fscanf(file, "%d %d %d", &intvlStart, &intvlEnd, &lambda);
            else
            {
                x=fscanf(file, "%d %d", &intvlStart, &lambda);
                intvlEnd = intvlStart;
            }
            //adjustedEnd = intvlEnd - nbr.nbr_lambda;      This is for input model where end is time till edge is active.
            adjustedEnd = intvlEnd;     //This is for input model where end is last instant at which end can be used.
            //nbr.edgeSchedules[intvls].adjustedEnd = nbr.edgeSchedules[intvls].intvlEnd - nbr.nbr_lambda;
            if (adjustedEnd >= intvlStart)       //else ignore this intvl because it can't really be used.
            {
                nbr.edgeSchedules[intvls].intvlStart = intvlStart; nbr.edgeSchedules[intvls].intvlEnd = intvlEnd;
                nbr.edgeSchedules[intvls].adjustedEnd = adjustedEnd; nbr.edgeSchedules[intvls].traveTime = lambda;
                nbr.edgeSchedules[intvls].divTime = -1;   //This is used only when interval is divided for mwf.
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

void Graph::wuGraph(const char* filePath, int noL, int numDrop, int normalizeWrite)
{
//    FILE* file = fopen(filePath,"r");
    
    ifstream inputFile(filePath);
    string inputLine;
    int drop = 0;
    int minNodeId = infinity;
    int vertexIdAdj = 0;
    
    int x, wuEdges = 0, vertices = 0;
//    x=fscanf(file, "%d %d",&vertices, &wuEdges);
    vector<tuple<int, int, int, int>> inputRows;
    compareTupleWu tupleCompWu;
    int u, v, t, lambda;
//    inputRows.resize(wuEdges);
    
    string wuOpEnd = "_wu.txt";
    std::string wuFile(filePath);
    unsigned long len = wuOpEnd.length();
    wuFile.replace(wuFile.length()-4, len, wuOpEnd);

        
    lambda = 1;
    for (int dx=0; dx<numDrop; dx++)    //Drop the first numDrop lines as they are comments.
        getline(inputFile, inputLine);
    while (getline(inputFile, inputLine))
//    for(int i = 0; i < wuEdges; i ++)
    {
        if (noL == 1)
        {
            x = sscanf(inputLine.c_str(), "%d %d %d", &u, &v, &t);
//            x=fscanf(file, "%d %d %d",&u, &v, &t);
        }
        else if (noL == 2)      //drop the element after u, v. Input is: u, v, {some wt.}, t.
        {
            x = sscanf(inputLine.c_str(), "%d %d %d %d", &u, &v, &drop, &t);
        }
        else
        {
            x = sscanf(inputLine.c_str(), "%d %d %d %d", &u, &v, &t, &lambda);
//            x=fscanf(file, "%d %d %d %d",&u, &v, &t, &lambda);
        }
        inputRows.push_back(std::make_tuple(u, v, t, lambda));
        wuEdges++;
        if (u < minNodeId)
            minNodeId = u;
        if (v < minNodeId)
            minNodeId = v;

        if (u > vertices-1)
            vertices = u+1;
        if (v > vertices-1)
            vertices = v+1;
    }
    if (minNodeId > 0)
    {
        vertices -= minNodeId;
        vertexIdAdj=minNodeId;
    }
//    fclose(file);
    std::sort<vector<tuple<int, int, int, int>>::iterator, compareTupleWu>(inputRows.begin(), inputRows.end(), tupleCompWu);
    FILE* wuOutput = fopen(wuFile.c_str(), "w");
    fprintf(wuOutput, "%d %d\n", vertices, wuEdges);
    int normalizeSub = 0;
    if (normalizeWrite)
        normalizeSub = get<2>(inputRows[0]);
    else
        normalizeSub = 0;
    for (int i = 0; i < inputRows.size(); i++)
    {
        get<0>(inputRows[i]) = get<0>(inputRows[i]) - vertexIdAdj;
        get<1>(inputRows[i]) = get<1>(inputRows[i]) - vertexIdAdj;
        get<2>(inputRows[i]) = get<2>(inputRows[i]) - normalizeSub;
        fprintf(wuOutput, "%d %d %d %d\n", get<0>(inputRows[i]), get<1>(inputRows[i]), get<2>(inputRows[i]), get<3>(inputRows[i]) );
    }
    fclose(wuOutput);
    
    buildXuanGraph(filePath, inputRows, vertices, wuEdges);
}

void Graph::readWuFile(const char* filePath)
{
    ifstream inputFile(filePath);
    string inputLine;
    int x, wuEdges = 0, vertices = 0; int u, v, t, lambda;

    vector<tuple<int, int, int, int>> inputRows;

    getline(inputFile, inputLine);
    x = sscanf(inputLine.c_str(), "%d %d", &vertices, &wuEdges);
    while (getline(inputFile, inputLine))
    {
        {
            x = sscanf(inputLine.c_str(), "%d %d %d %d", &u, &v, &t, &lambda);
        }
        inputRows.push_back(std::make_tuple(u, v, t, lambda));
        //wuEdges++;
    }
    buildXuanGraph(filePath, inputRows, vertices, wuEdges);
}

void Graph::buildXuanGraph(const char* filePath, vector<tuple<int, int, int, int>>& inputRows, int vertices, int wuEdges)
{
    compareTuple tupleCompXuan;
    //const char* intermediateFile = "/Users/anujjain/research/temporalGraph/WuTemporalGraph/tempath/intermediate.txt";
    const char* intermediateFile = "./intermediate.txt";
    int u, v, t, lambda;

    string xuanOpEnd = "_xuanOp.txt";
    std::string opFile(filePath);
    unsigned long len = xuanOpEnd.length();
    opFile.replace(opFile.length()-4, len, xuanOpEnd);

    std::sort<vector<tuple<int, int, int, int>>::iterator, compareTuple>(inputRows.begin(), inputRows.end(), tupleCompXuan);
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
//            if (currentRunningU != currentRunningV)
//            {
                fprintf(fileToWrite, "%d %d %d", currentRunningU, currentRunningV, intvlCount);     //Output format is u, v, intvlCount, (intervals). Intervals are of the form t, lambda.
                for (int intC = 0; intC < intvlCount; intC++)
                    fprintf(fileToWrite, " %d %d", std::get<0>(currentIntvls[intC]), std::get<1>(currentIntvls[intC]) );
                fprintf(fileToWrite, "\n");
//            }
//            else
//                numEdgesCollapsed++;        //Completely ignore this edge as start and end node are same so it's a loop.
            currentIntvls.clear();      //Reset intvls for the new edge.
            intvlCount = 0;
            //Prepare to start looking at next edge.
            currentRunningU = u; currentRunningV = v;
            intvlCount++;
            currentIntvls.push_back(std::make_tuple(std::get<2>(inputRows[i]), std::get<3>(inputRows[i]) ) );
        }
    }
    //now print out the consolidated last edge.
//    if (currentRunningU != currentRunningV)
//    {
        fprintf(fileToWrite, "%d %d %d", currentRunningU, currentRunningV, intvlCount);
        for (int intC = 0; intC < intvlCount; intC++)
            fprintf(fileToWrite, " %d %d", std::get<0>(currentIntvls[intC]), std::get<1>(currentIntvls[intC]) );
        fprintf(fileToWrite, "\n");
//    }
//    else
//        numEdgesCollapsed++;        //Completely ignore this edge as start and end node are same so it's a loop.
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
    
    collapseIntervalsWriteOuput(intermediateFile, opFile.c_str());      //Write the final output. Each Edge is in format ==> (u, v, intvlCount, [s, e, d])
}


void Graph::collapseIntervalsWriteOuput(const char *filePath, const char *opFile)
{
//    const char* opFile = "/Users/anujjain/research/temporalGraph/WuTemporalGraph/tempath/xuanOp.txt";
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
    int low = 0, high = v.adjustedNumIntvls-1, mid = high/2, intvlFound = 0, retVal = infinity;
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
                        retVal = infinity;            //There is no valid interval.
                        intvlID = -1;
                    }
                }
                else
                {
                    intvlFound = 1;
                    retVal = infinity;            //There is no valid interval.
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
    
    int s;
#ifdef __TEST__
    sources.push_back(sourceId);
#else
    for(int i = 0 ;i < RUN_COUNT ;i ++)
    {
        s=rand()%V;
        sources.push_back(s);
    }
#endif
}


void Graph::initial_query(const char* filePath)
{
    t_start = 0;
    t_end = infinity;
    
    ifstream inputFile(filePath);
    string inputLine;

    int s;
    while (getline(inputFile,inputLine))
    {
        sscanf(inputLine.c_str(), "%d",&s);
        sources.push_back(s%V);
    }
}


void Graph::initial_ds_ea()
{
    for(int i=0; i<V; i++){
        arr_time[i]= infinity;
        father[i] = make_tuple(i,-1,0);
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
        initial_ds_ea();
//        earliest_arrival(sources[i]);
        earliest_arrival_pair(sources[i]);
//        earliest_arrival_fibo(sources[i]);
//        earliest_arrival_fibo_external(sources[i]);
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
    father[source] = make_tuple(source,-1,0);
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
                    father[u->neighbors[i].nbrId] = make_tuple(u->nodeId, intvlID, tDepart) ;
                    //if (! nodesInHeap.check_bit_obj_present(u->neighbors[i].nbrId) )
                    minHeap.push_back(u->neighbors[i].nbrId);
                    std::push_heap(minHeap.begin(), minHeap.end(), heapCompFn);
                }
            }
        }
    }
    t.stop();
//    printResults(source);
    time_sum += t.GetRuntime();
}

tuple<int,int> Graph::earliest_arrival_pair(int source,int retRchd)
{
    ofstream minHeapMonitorOutput(minHeapMonitor);
    unsigned long maxHeapSize = 0, avgHeapSize = 0, count = 0, currentAvg = 0;
    int numNodesReached=0;
    Timer t;
#ifdef MEASUREHEAP_DET
    clock_t ticks;
    clock_t insertTimer = 0; //clock() - clock();
    clock_t remMinTimer = 0; //clock() - clock();
#endif
    Node *u;
    int nodeID;
    int tDepart, intvlID = -1;
    vector<pair<int, int>> minHeap;    //Fix the pq to make it a min heap instead of max heap.
    int maxFmstTime = 0;
    nodeComparison2 heapCompFn;

    bit_queue closedNodes((int)vertices.size());
        
    arr_time[source]=t_start;
    father[source] = make_tuple(source,-1,0);
    minHeap.push_back(std::make_pair(arr_time[source], source));
    int i = 0, numRepeatedNodes = 0, numInserts = 0;
    
    t.start();

    while (!minHeap.empty())
    {
        i++;
        nodeID = minHeap.front().second;
#ifdef MEASUREHEAP_DET
        ticks = clock();
#endif
        std::pop_heap<vector<pair<int, int>>::iterator, nodeComparison2>(minHeap.begin(), minHeap.end(), heapCompFn); minHeap.pop_back();      //moves the min element to last and then removes it from heap.
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
        if (maxFmstTime < arr_time[nodeID])
            maxFmstTime = arr_time[nodeID];
        u = &vertices[nodeID];
        for (int i = 0; i < u->numNbrs; i++)
        {
            if (! closedNodes.check_bit_obj_present( u->neighbors[i].nbrId))
            {
                tDepart = earliestUseEdgeAfterT(nodeID, u->neighbors[i], arr_time[nodeID], intvlID );
                if (intvlID >= 0)
                {
                    if ( (tDepart + u->neighbors[i].edgeSchedules[intvlID].traveTime) < arr_time[u->neighbors[i].nbrId])
                    {
                        arr_time[u->neighbors[i].nbrId] = tDepart + u->neighbors[i].edgeSchedules[intvlID].traveTime;
                        father[u->neighbors[i].nbrId] = make_tuple(u->nodeId, intvlID, tDepart) ;
#ifdef MEASUREHEAP_DET
                        ticks = clock();
#endif
                        minHeap.push_back(std::make_pair(arr_time[u->neighbors[i].nbrId], u->neighbors[i].nbrId));
                        std::push_heap(minHeap.begin(), minHeap.end(), heapCompFn);
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
    if (retRchd == 1)
        return make_tuple(numNodesReached,maxFmstTime);
    time_sum += t.GetRuntime();
#ifdef __TEST__
    printEarliestResultsTest(source);
#endif
    return make_tuple(0,0);
}

    

void Graph::earliest_arrival_fibo(int source)
{
    ofstream minHeapMonitorOutput(minHeapMonitor);
    unsigned long maxHeapSize = 0, avgHeapSize = 0, count = 0, currentAvg = 0;
    Timer t;
#ifdef MEASUREHEAP_DET
    clock_t ticks;
    clock_t insertTimer = 0; //clock() - clock();
    clock_t decKeyTimer = 0; //clock() - clock();
    clock_t remMinTimer = 0; //clock() - clock();
#endif
    Node *u;
    int nodeID, nodeArrTime;
    int tDepart, intvlID = -1;
    fibHeap minHeap(V);
//    regHeap minHeap(V);

    bit_queue nodesInHeap((int)vertices.size());
    bit_queue closedNodes((int)vertices.size());
            
    arr_time[source]=t_start;
    father[source] = make_tuple(source,-1,0);
    minHeap.insert(source, arr_time[source]);
    int i = 0, numRepeatedNodes = 0, numDecreases = 0, numInserts = 0;
    
    t.start();
    while (!minHeap.empty())
    {
        i++;
#ifdef MEASUREHEAP_DET
        ticks = clock();
#endif
        pair<int, int> minNK = minHeap.removeMin();
#ifdef MEASUREHEAP_DET
        remMinTimer += (clock() - ticks);
#endif
        nodeID = get<0>(minNK);
        nodeArrTime = get<1>(minNK);

        if (nodeArrTime < 0)
        {
            cout << "Fib Heap returned, -ve key: " << nodeArrTime << " for node:" << nodeID << ". Exiting \n";
            exit(1);
        }

        if (closedNodes.check_bit_obj_present(nodeID))      //This shouldn't happen.
        {
            numRepeatedNodes++;
            continue;
        }
        closedNodes.queue_add_bit(nodeID);         //Add it to the closedNodes.
        u = &vertices[nodeID];
        for (int i = 0; i < u->numNbrs; i++)
        {
            if (! closedNodes.check_bit_obj_present( u->neighbors[i].nbrId))
            {
                tDepart = earliestUseEdgeAfterT(nodeID, u->neighbors[i], arr_time[nodeID], intvlID );
                if (intvlID != -1) {
                    if ( (tDepart + u->neighbors[i].edgeSchedules[intvlID].traveTime) < arr_time[u->neighbors[i].nbrId])
                    {
                        arr_time[u->neighbors[i].nbrId] = tDepart + u->neighbors[i].edgeSchedules[intvlID].traveTime;
                        father[u->neighbors[i].nbrId] = make_tuple(u->nodeId, intvlID, tDepart) ;
                        //if (nodesInHeap.check_bit_obj_present(u->neighbors[i].nbrId))
                        if (minHeap.partOfHeap(u->neighbors[i].nbrId))
                        {
#ifdef MEASUREHEAP_DET
                            ticks = clock();
#endif
                            minHeap.decrease_key(u->neighbors[i].nbrId, arr_time[u->neighbors[i].nbrId]);
#ifdef MEASUREHEAP_DET
                            decKeyTimer += (clock() - ticks);
#endif
                            numDecreases++;
                        }
                        else
                        {
#ifdef MEASUREHEAP_DET
                            ticks = clock();
#endif
                            minHeap.insert(u->neighbors[i].nbrId, arr_time[u->neighbors[i].nbrId]); numInserts++;
#ifdef MEASUREHEAP_DET
                            insertTimer += (clock() - ticks);
#endif
                        }
                    }
                }
            }
        }
/*        if (i % 50 == 0)
        {
            minHeapMonitorOutput << minHeap.numNodes << "; " ;
            if (minHeap.numNodes > maxHeapSize)
                maxHeapSize = (int)minHeap.numNodes;
            avgHeapSize = (currentAvg*(count) + (unsigned long)minHeap.numNodes )/(count + 1);
            currentAvg = avgHeapSize; count++;
            minHeapMonitorOutput << avgHeapSize << "\n";
        }*/

    }
    t.stop();
/*    minHeapMonitorOutput << "Max Heap size: " << maxHeapSize << "\n";
    minHeapMonitorOutput << "Avg Heap size: " << avgHeapSize << "\n";
    minHeapMonitorOutput << "Num Repeated Nodes: " << numRepeatedNodes << "\n";
    minHeapMonitorOutput << "Num Decreases: " << numDecreases << "\n";
    minHeapMonitorOutput << "Num Inserts: " << numInserts << "\n";
#ifdef MEASUREHEAP_DET
    minHeapMonitorOutput << "Insert Time: " << insertTimer << "\n";
    minHeapMonitorOutput << "Rem Min Time: " << remMinTimer << "\n";
    minHeapMonitorOutput << "DecKey Time: " << decKeyTimer << "\n";
#endif
*/
//    printResults(source);
    time_sum += t.GetRuntime();
    printEarliestResultsTest(source);
}

/*void Graph::earliest_arrival_fibo_external(int source)
{
    ofstream minHeapMonitorOutput(minHeapMonitor);
    unsigned long maxHeapSize = 0, avgHeapSize = 0, count = 0, currentAvg = 0;
    Timer t;
    clock_t ticks;
    clock_t insertTimer = 0; //clock() - clock();
    clock_t decKeyTimer = 0; //clock() - clock();
    clock_t remMinTimer = 0; //clock() - clock();
    t.start();
    Node *u;
    int nodeID, nodeArrTime;
    int tDepart, intvlID = -1;
    fibonacci_heap<int, int> minHeap;
    vector<fibonacci_heap<int, int>::node> fibNodes;
    fibNodes.resize(vertices.size());
    
    bit_queue nodesInHeap((int)vertices.size());
    bit_queue closedNodes((int)vertices.size());
            
    arr_time[source]=t_start;
    father[source] = make_tuple(source,-1,0);
    fibNodes[source] = minHeap.insert(arr_time[source], source);        //insert(key, value)
    int i = 0, numRepeatedNodes = 0, numDecreases = 0;
    
    while (minHeap.size() > 0)
    {
        i++;
        nodeID = minHeap.top().data();
        nodeArrTime = minHeap.top().key();
        ticks = clock();
        minHeap.remove();
        remMinTimer += (clock() - ticks);
        if (closedNodes.check_bit_obj_present(nodeID))      //This shouldn't happen.
        {
            numRepeatedNodes++;
            continue;
        }
        closedNodes.queue_add_bit(nodeID);         //Add it to the closedNodes.
        u = &vertices[nodeID];
        for (int i = 0; i < u->numNbrs; i++)
        {
            if (! closedNodes.check_bit_obj_present( u->neighbors[i].nbrId))
            {
                tDepart = earliestUseEdgeAfterT(nodeID, u->neighbors[i], arr_time[nodeID], intvlID );
                if (intvlID != -1) {
                    if ( (tDepart + u->neighbors[i].edgeSchedules[intvlID].traveTime) < arr_time[u->neighbors[i].nbrId])
                    {
                        arr_time[u->neighbors[i].nbrId] = tDepart + u->neighbors[i].edgeSchedules[intvlID].traveTime;
                        father[u->neighbors[i].nbrId] = make_tuple(u->nodeId, intvlID, tDepart) ;
                        if (nodesInHeap.check_bit_obj_present(u->neighbors[i].nbrId))
                        {
                            ticks = clock();
                            minHeap.decrease_key(fibNodes[u->neighbors[i].nbrId], arr_time[u->neighbors[i].nbrId]);
                            decKeyTimer += (clock() - ticks);
                            numDecreases++;
                        }
                        else
                        {
                            ticks = clock();
                            fibNodes[u->neighbors[i].nbrId] = minHeap.insert(arr_time[u->neighbors[i].nbrId], u->neighbors[i].nbrId);     //insert(key, value)
                            insertTimer += (clock() - ticks);
                            nodesInHeap.queue_add_bit(u->neighbors[i].nbrId);
                        }
                    }
                }
            }
        }
        if (i % 50 == 0)
        {
            minHeapMonitorOutput << minHeap.size() << "; " ;
            if (minHeap.size() > maxHeapSize)
                maxHeapSize = (int)minHeap.size();
            avgHeapSize = (currentAvg*(count) + (unsigned long)minHeap.size() )/(count + 1);
            currentAvg = avgHeapSize; count++;
            minHeapMonitorOutput << avgHeapSize << "\n";
        }

    }
    t.stop();
    minHeapMonitorOutput << "Max Heap size: " << maxHeapSize << "\n";
    minHeapMonitorOutput << "Avg Heap size: " << avgHeapSize << "\n";
    minHeapMonitorOutput << "Num Repeated Nodes: " << numRepeatedNodes << "\n";
    minHeapMonitorOutput << "Num Decreases: " << numDecreases << "\n";
    minHeapMonitorOutput << "Insert Time: " << insertTimer/CLOCKS_PER_SEC << "\n";
    minHeapMonitorOutput << "Rem Min Time: " << remMinTimer/CLOCKS_PER_SEC << "\n";
    minHeapMonitorOutput << "DecKey Time: " << decKeyTimer/CLOCKS_PER_SEC << "\n";

//    printResults(source);
    time_sum += t.GetRuntime();
    printEarliestResultsTest(source);
}*/


void Graph::printResults(int source)
{
    int x;
    cout << "Source is: " <<source << "\n";
    for (int i = 0; i < arr_time.size(); i++)
    {
        if (arr_time[i] >= infinity)
            continue;
        cout << "Arrival time at node " << i << " is: " << arr_time[i] << "\n";
        cout << "Path in rev is: ";
        x = i;
        while (get<0>(father[x]) != x)
        {
            cout << x <<"(" << get<1>(father[x]) << " " << get<2>(father[x]) << ") ";
            x = get<0>(father[x]);
        }
        cout << x << "\n";
    }
    
}

void Graph::printEarliestResultsTest(int source)
{
    int rv = 0;
    ofstream earliestOut(earliestResults);
    earliestOut << V << "\n";
    for (int i = 0; i < arr_time.size(); i++)
    {
//        if (arr_time[i] >= infinity)
//            continue;
        earliestOut << i << " " << arr_time[i] << "  "  << "\n";
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
        initial_ds_ea();
        shortest_path(sources[i]);
//        shortest_path_xuan(sources[i]);
    }
    
    print_avg_time();
}

void Graph::run_mwf()
{
    time_sum=0;
    
    for(int i = 0 ;i < sources.size() ;i ++)
    {
        cout << "prioritized mwf called"<<endl;
        initial_ds_ea();
        //earliest_arrival(sources[i]);
        //initial_ds_ea();earliest_arrival_pair(sources[i]);
        //initial_ds_ea();
//        shortest_path(sources[i]);
        //minWaitForemost(sources[i]);
        minWaitForemostPrioritized(sources[i]);
 //       minWaitForemostPrioritizedNoSet(sources[i]);
    }
    
    print_avg_time();
}


void Graph::shortest_path(int source)
{
    Timer t;

    vector<vector<incrementalJourney>> allHopJourneys; //At each hop there is a vector of foremost incremental journeys, discovered in that hop.
    vector<std::pair<int, int>> shortestJourneyPointer;     // <hop count where journey ended, index in allHopJourneys on that hop>
    vector<std::tuple<int, int, int>> earliestKnownTimeArrival;              // <earliestArrivalTime, hop in which this time achieved, index in allHopJourneys[thishop] vector.>

    allHopJourneys.resize(V);   //Max number of hops is the number of nodes.
    shortestJourneyPointer.resize(V);
    earliestKnownTimeArrival.resize(V);
    earliestKnownTimeArrival.assign(V, std::make_tuple(infinity, 0, -1));
    incrementalJourney journeyIncrement(source, 0, -1, -1, -1, -1, -1);   //current node, arrTime, prevNode, prevEdgeId, prevIntvlId, pDep, indexInPrevHopV

    t.start();
    int intvlID = -1, departTime, newArrivTime, hopCount = 0, nextNode;
    earliestKnownTimeArrival[source] = std::make_tuple(0, 0, 0);
    allHopJourneys[hopCount].push_back(journeyIncrement);

    int  numNodesSeen = 1, numNewNodesInCurrentHop = 1;     //number of new foremost nodes in current hop.
    while ( (hopCount < (V-1)) && (numNewNodesInCurrentHop > 0) && (numNodesSeen < V))
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
                    {
                        numNodesSeen++;
                        int incrIndex = std::get<2>(earliestKnownTimeArrival[nextNode]);
                        shortestJourneyPointer[nextNode] = std::make_pair(hopCount, incrIndex);
                    }
                    std::get<0>(earliestKnownTimeArrival[nextNode]) = newArrivTime;
                }       //found a better arrival time.
            }       //Looking at all the nbrs of a node to extend in current hop count.
        }       //for loop extending journey from hopCount-1.
    }       //While loop to cover all hop counts.

    t.stop();
    time_sum += t.GetRuntime();

#ifdef __TEST__
    build_shortestJourneys(source, shortestJourneyPointer, allHopJourneys);
    print_shortest_results_test(source);
#endif

//    print_shortest_paths(source);
}


int Graph::edgeAndScheduleSel(vector<std::tuple<int, int, int>>& e_min, vector<int>& t_min, vector<int>& t_LBD)
{
    Timer t;
    int numUpdated = 0;
    static int numTimesCalled = 0;
    int numNbrs = 0;
    numTimesCalled++;
    vector<int> t_arrival; t_arrival.resize(V);
    for (int i=0; i < V; i++)
    {
        get<0>(e_min[i]) = -1; get<1>(e_min[i]) = -1; get<2>(e_min[i]) = -1;
        t_min[i] = infinity;
        t_arrival[i] = t_LBD[i];
    }
    t.start();
    int departTime = 0, intvlID = -1, newArrivTime = -1, nbr=-1;
    for (int i=0;i<V;i++)
    {
        for (int j=0; j<vertices[i].numNbrs;j++)
        {
            nbr = vertices[i].neighbors[j].nbrId;
            departTime = earliestUseEdgeAfterT(i, vertices[i].neighbors[j], t_LBD[i], intvlID);
            if ( (departTime == -1) || (departTime >= infinity) || (intvlID == -1) )
                continue;
            newArrivTime = departTime + vertices[i].neighbors[j].edgeSchedules[intvlID].traveTime;
            if (newArrivTime < t_arrival[nbr])
            {
                get<0>(e_min[nbr]) = i;     //prevNode
                get<1>(e_min[nbr]) = j;
                get<2>(e_min[nbr]) = intvlID;
                t_min[nbr] = departTime;
                t_arrival[nbr] = newArrivTime;
                numUpdated++;
            }
            numNbrs++;
        }
    }
    t.stop();
    cout << "It took " << t.GetRuntime() << " secs to get through all edges. Num times called: " << numTimesCalled <<". Num Nbrs="<< numNbrs << endl;
    return numUpdated;
}

void Graph::shortest_path_xuan(int source)
{
    Timer t;
    vector<vector<incrementalJourney>> allHopJourneys; allHopJourneys.resize(V); //At each hop there is a vector of foremost incremental journeys, discovered in that hop.
    vector<std::pair<int, int>> shortestJourneyPointer; shortestJourneyPointer.resize(V);    // <hop count where journey ended, index in allHopJourneys on that hop>
    vector<std::pair<int, int>> journey1; journey1.resize(V);    // <hop count where journey ended, index in allHopJourneys on that hop>  This is the J[v].
    vector<std::pair<int, int>> journey2; journey2.resize(V);    // <hop count where journey ended, index in allHopJourneys on that hop>  This is the J[v].
    vector<std::pair<int, int>>& currentJourney = journey1;     // references to journey, so they can be swaped every hop.
    vector<std::pair<int, int>>& nextJourney = journey2;     // references to journey, so they can be swaped every hop.
    vector<std::pair<int, int>>& tempSwap = journey1;     // references to journey, so they can be swaped every hop.
    vector<int> t_LBD; t_LBD.resize(V);        // <earliestArrivalTime, hop in which this time achieved, index in allHopJourneys[thishop] vector.>
    vector<std::tuple<int, int, int>> e_min; e_min.resize(V);   //prevNode, nbrIndex in nbrList, intvlId on the edge(prevNode,current node). current node is the index in vector.
    vector<int> t_min; t_min.resize(V);
    
    t_LBD.assign(V, infinity);
    journey1.assign(V, std::make_pair(-1, -1));
    journey2.assign(V, std::make_pair(-1, -1));
    shortestJourneyPointer.assign(V, std::make_pair(-1, -1));

    incrementalJourney journeyIncrement(source, 0, -1, -1, -1, -1, -1);   //current node, arrTime, prevNode, prevEdgeId, prevIntvlId, pDep, indexInPrevHopV

    t.start();
    int hopCount = 0;
    t_LBD[source] = 0;
    allHopJourneys[hopCount].push_back(journeyIncrement);
    get<0>(currentJourney[source]) = hopCount; get<1>(currentJourney[source]) = (int)allHopJourneys[hopCount].size()-1;
    shortestJourneyPointer[source] = currentJourney[source];

    int  numNodesSeen = 1;      //Total number of nodes so far that don't have t_LBD as infinity.
    while ( (hopCount < (V-2)) && (numNodesSeen < V-1))
    {
        hopCount++;
        int numUpdated = edgeAndScheduleSel(e_min,t_min, t_LBD);
        if (numUpdated == 0)
            break;
        for (int i= 0; i < V; i++)
        {
            if (get<0>(e_min[i]) != -1)     //e_min[i] (0) is the prevNode to i.
            {
                int prevNode = get<0>(e_min[i]);     //prevNode is u.
                int nbrIndex = get<1>(e_min[i]); int intvlId = get<2>(e_min[i]);
                int arrTime = t_min[i] + vertices[prevNode].neighbors[nbrIndex].edgeSchedules[intvlId].traveTime;
                int hc_prevNode = get<0>(currentJourney[prevNode]);
                int index_prevNode = get<1>(currentJourney[prevNode]);
                incrementalJourney newIncrement(i, arrTime, prevNode, nbrIndex, intvlId, t_min[i], index_prevNode);
                allHopJourneys[hopCount].push_back(newIncrement);
                get<0>(nextJourney[i]) = hopCount; get<1>(nextJourney[i]) = (int)allHopJourneys[hopCount].size()-1;
                if (t_LBD[i] >= infinity)
                {
                    shortestJourneyPointer[i] = nextJourney[i];
                    numNodesSeen++;
                }
                t_LBD[i] = arrTime;
            }
            else
            {
                nextJourney[i] = currentJourney[i];
            }
        }
        tempSwap = currentJourney;
        currentJourney = nextJourney;
        nextJourney = tempSwap;
    }
    t.stop();
    time_sum += t.GetRuntime();

#ifdef __TEST__
    build_shortestJourneys(source, shortestJourneyPointer, allHopJourneys);
    print_shortest_results_test(source);
#endif

//    print_shortest_paths(source);
}

void Graph::minWaitForemost(int source)
{
    Timer t;
    vecFullList.resize(V);
    toExpandList.push_back(source);

    t.start();
    get<1>(vecFullList[source]).insert(make_tuple(0,0,source, 0, 0));
    get<0>(vecFullList[source]).insert(make_tuple(0,0,0));        //journeys terminating at source (arr_tm, wait time)
    mNumJourExtInst = 0;
    
    int exploreNode;
    while (!toExpandList.empty())
    {
        exploreNode = toExpandList.front();
        extendMWFJourneys(exploreNode);
        toExpandList.pop_front();
    }
    t.stop();
    time_sum += t.GetRuntime();
    printMWFWalks(source);
}

void Graph::printMWFWalks(int source)
{
    int prevNode, departTime, arrTime, waitTime, numRchableVerts = 0;
    set <tuple<int, int, int, int, int>, compareMFSetElements >::iterator itMWFJourney;
    cout << "output format: dest(arr,wait_time)<--(depTime)node..\n";
    for (int i=0; i < V; i++)
    {
        if (i==source)
        {
            numRchableVerts++;
            cout << source << "\n";
            continue;
        }
        prevNode = i;
        if ((get<1>(vecFullList[prevNode])).empty())
            continue;
        numRchableVerts++;
        itMWFJourney = (get<1>(vecFullList[prevNode])).begin();
        cout << "Walk: "<< i << " ";// << "("<< get<0>(*itMWFJourney) <<","<<get<1>(*itMWFJourney)<<")";//\n";
        arrTime = get<0>(*itMWFJourney); waitTime = get<1>(*itMWFJourney);
        //cout << arrTime*100+waitTime << "\n";
        prevNode = get<2>(*itMWFJourney); departTime = get<3>(*itMWFJourney);
        while (prevNode != source)
        {
            cout << "("<< arrTime <<","<<waitTime<<")";
            cout << "<--"  << "(" << departTime <<")"<< prevNode;
            itMWFJourney = (get<1>(vecFullList[prevNode])).lower_bound(make_tuple(departTime,0,0,0,0));
            if (itMWFJourney == (get<1>(vecFullList[prevNode])).end())
                --itMWFJourney;
            else if (departTime < get<0>(*itMWFJourney))
                --itMWFJourney;
            
            arrTime = get<0>(*itMWFJourney);
            waitTime = get<1>(*itMWFJourney);
            prevNode = get<2>(*itMWFJourney);
            departTime = get<3>(*itMWFJourney);
            if((get<1>(vecFullList[prevNode])).empty())
            {
                cout << "Something went wrong\n\n";
                break;
            }
        }
        cout << "("<< arrTime <<","<<waitTime<<")";
        cout << "<--"  << "(" << departTime << ")" << prevNode << "\n\n";
    }
    cout << "Num Reachable: " << numRchableVerts << endl;
}

void Graph::extendMWFJourneys(int exploreNode)
{
    set<tuple<int, int, int>>::iterator itNewJourneySet;
    tuple<int, int, int> journeyToExpand;
    while (!get<0>(vecFullList[exploreNode]).empty())
    {
        itNewJourneySet = (get<0>(vecFullList[exploreNode])).begin();
        journeyToExpand = *itNewJourneySet;
        extendJourney(journeyToExpand, exploreNode);
        (get<0>(vecFullList[exploreNode])).erase(itNewJourneySet);
    }
}

/************************************************************************************************************************************
 while TEL not empty do
   (u) = TEL.delete(); //  u has atleast one new journey.
   while FL[u]<1> not empty do. //FL[u]<1> represents second set in the tuple of wo sets stored at FL[u].
      (t,m) = FL[u]<1>.delete(); //t is arrival time at u, m is wait time upto u.
       for each (u,v) in E do
           generate the necessary one-edge expansions (v,t',m') of (u,t,m) one by one;
           if (t',m') is not dominated by a tuple in FL[v]<0>
                {FL[v]<0>.insert(t',m',u); // insert eliminates tuples in FL[v]<0> dominated by (t',m'). u is the previous node that is only used when tracing back paths.
                 FL[v]<1>.insert(t',m');  //set of new journeys. this should also eliminate dominated tuples
                  TEL.insert(v);} //This indicates v has some new journeys to be explored in the next iteration.
 **********************************************************************************************************************************/

/************************************************************************************************************************************
t_mwf[s] = (0,0); t_mwf[x] = (inf,inf) for all x neq s;
 Each node has 2 lists: FullList (FL) and ToExpandList (TEL). both lists are empty to begin with except FL[s] = (0,0) and TEL[s] = (0.0);
 PQ.insert(s,0,0)u=s;
 numNodesClosed++;
while ((numNodesClosed < N) && (PQ.notEmpty))
 {      (u,a,w) = PQ.removeMin();      //This PQ is on min a. If tie on min a, min w among all min a is picked.
    If ( u is not closed already)
    {       Mark u closed //(mwf walk has been found to u);
        numNodesClosed++;
    }
    else if (u,a,w) is marked dominated.  //This can only happen if u was already closed.
        continue;
    for each neighbor v of u
    {       generate necessary expansions (v,a',w') from u to v, using journeys in TEL[u];
        If (v,a',w') is not dominated by journeys in FL[v]
        {        Insert (v,a',w') in FL[v];
            Mark all the subsequent journeys FL[v] that are dominated by (a',w') as dominated;
            Insert (v,a',w') in the PQ;  //This can be done after all journeys to v have been expanded to avoid inserting dominated journeys in PQ.
        }
    }
 }
 **********************************************************************************************************************************/

/******************
 u=s; nodesReached=1; a(u)=0;w(u)=0
 mwf[s]=(tstart,0)
 An array lastJourney[] is mantained that stores only the last non dominated journey at every node.
while (nodesReached < N)
 {
    for each nbr of u
    {
        insert in min heap all necessary expansions as [arrivalTime,waitTime,nbr]
    }
    nonDominatedJFound=False;
    while (!nonDominatedJFound && !heap.empty() )
     {
        remove min as (a,w, v).        (based on arrivalTime; Tie-breaker with waitTime) //a=arrivalTime, w=waitTime, v=nodeId
        if (v,a,w) not dominated by the last journey at (v)
            lastJourney[v] = (a,w)                  //Replace last journey at v with (v,a,w)
            u=v,a(u)=a,w(u)=w
            nonDominatedJFound=true;
        else
            continue;
        If (v is reached first time.)
            mwf[v] = (a,w)
            nodesReached++
     }
 }
***********************/


/*************************
 Represent each interval as {u,v,u_nbrIndx,(s,e),lambda,prevJourneyIndex} in a list Intvls[]        //u_nbrIndx represents index of v in nbrs list of u.
 journey is represented as (a,w,prevNode,lastExpAt) where lastExpAt is (s,lambda) or NULL.
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
        removeMin(top(Intvls,mergeIntvls)
        If (newIntvl.Id != -1)      //this means this intvl came from static set of intvls.
            prevJourneyIndex=vert[u][v][newIntvl.Id].prevJourneyIndex;
        else        //This means this intvl came from mergeIntvls list and was created by breaking up existing interval.
            prevJourneyIndex=newIntvl.prevJourneyIndex
        if (prevJuorneyIndex==-1)
        {
            If (listJ[u].last().arr < newIntvl.strt)
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
            listJ[u][prevJourneyIndex].expandedIntvl[v] = {newIntvl};
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
                        newIntvlCreated=(v,w,newJourney.arr,e,nextIntvl.lambda);        //Sub-interval of the static intvl created as all journeys exit at start of an interval.
                        newIntvlCreated.prevJourneyIndex=newJourney;
                        Insert(mergeIntvls, {newIntvlCreated})
        }
        else
            ignore intvl for expansion.
        newIntvl = min (top(Intvls, mergeIntvls));     //Comparison of arrTime= s+lambda and secondary start time.
    }
    numNodesRchd+= numNewNodesRchd;
 }
 
 ****************************/

void Graph::minWaitForemostPrioritizedNoSet(int source)
{
    Timer t;
    int numNodesFinalized = 1;
    int numHeapAdded=0, numJourneysFinalized=0;
    vecTuple.resize(V);
    mwfJourneys.resize(V);
    for (int i=0;i<V;i++)
    {
        vecTuple[i]=make_tuple(-1,-1,-1,-1);
        mwfJourneys[i]=make_tuple(-1,-1,-1,-1);
    }
    nodeCompareArrWt heapCompFn;

    bit_queue closedNodes((int)vertices.size());
    closedNodes.queue_add_bit(source);
    tuple<int,int,int,int> currJourney = make_tuple(t_start,0,0,source);
    mwfJourneys[source] =currJourney;
    vecTuple[source] = currJourney;
    int nbr, intvlID, departTime, newArrivTime, newWaitTime;
    bool isDominated;
    bool bFoundNewViableJourney=false;
    t.start();
    tuple<int,int> nodesRchable_maxFmstTime = earliest_arrival_pair(source,1); //4092652
    int nodesReachable= get<0>(nodesRchable_maxFmstTime), maxFmstTime=get<1>(nodesRchable_maxFmstTime);
    cout << "Num Nodes Reachable: " << nodesReachable << endl;
    while (numNodesFinalized < nodesReachable)
    {
//        if (numHeapAdded % 1000000 == 0)
//            cout << "NumHeapAdded: " << numHeapAdded << " numFinal: " << numJourneysFinalized << " numNodesFinal: " << numNodesFinalized << endl;
        int node = get<3>(currJourney), arrTime = get<0>(currJourney), wtTime = get<1>(currJourney), numHops;
        for (int j=0; j<vertices[node].numNbrs;j++)
        {
            nbr = vertices[node].neighbors[j].nbrId;
            departTime = earliestUseEdgeAfterT(node, vertices[node].neighbors[j], arrTime, intvlID);
            if ( (departTime == -1) || (departTime >= infinity) || (intvlID == -1) ) //or departTime > max fmst arrival time. EARLYTMNT
                continue;
//            do {      // UNCOMMENT        For konect graphs all travel times are 1 so only first possible interval needs to be considered.
                newArrivTime = departTime + vertices[node].neighbors[j].edgeSchedules[intvlID].traveTime;
                if (newArrivTime > maxFmstTime)     //This can be simply discarded as arrTime is too high.
                    continue;
                if ((node == source) && (arrTime==t_start))  //this is nbr of source & start of the journey
                    newWaitTime = 0;
                else
                    newWaitTime = departTime-arrTime+wtTime;   //get<0> is arrival time at node and get<1> is wait time till node.
                numHops=get<2>(currJourney)+1;
                if (get<3>(vecTuple[nbr]) != -1)
                {
                    isDominated= checkDominance(make_tuple(get<0>(vecTuple[nbr]), get<1>(vecTuple[nbr])), make_tuple(newArrivTime,newWaitTime));
                    if (isDominated)
                        continue;
                }
                minHeap.push_back(make_tuple(newArrivTime,newWaitTime,numHops,nbr));
                std::push_heap(minHeap.begin(), minHeap.end(), heapCompFn);
                numHeapAdded++;
                intvlID++;
                if (intvlID < vertices[node].neighbors[j].numIntvls)
                    departTime = vertices[node].neighbors[j].edgeSchedules[intvlID].intvlStart;
                    
  //          } while (intvlID < vertices[node].neighbors[j].numIntvls);
        }
        //remove min
        bFoundNewViableJourney=false;
        while ((!bFoundNewViableJourney) && !(minHeap.empty()))
        {
            currJourney = minHeap.front();
            std::pop_heap(minHeap.begin(), minHeap.end(), heapCompFn); minHeap.pop_back();
            node = get<3>(currJourney); arrTime = get<0>(currJourney); wtTime = get<1>(currJourney);numHops=get<2>(currJourney);
            if (get<3>(vecTuple[node]) != -1)
            {
                isDominated= checkDominance(make_tuple(get<0>(vecTuple[node]), get<1>(vecTuple[node])), make_tuple(arrTime,wtTime));
                if (!isDominated)
                {
                    bFoundNewViableJourney = true;
                    vecTuple[node] = currJourney;
                    numJourneysFinalized++;
                }
                else
                    continue;
            }
            else if (!closedNodes.check_bit_obj_present(node))
            {
                closedNodes.queue_add_bit(node);
                mwfJourneys[node]=currJourney;
                numNodesFinalized++;
                bFoundNewViableJourney = true;
                vecTuple[node] = currJourney;
                numJourneysFinalized++;
            }
        }
    }
    t.stop();
    time_sum += t.GetRuntime();
    printMWFWalksPrioritizedNoSet(source);
    cout << "NumHpItems: " << numHeapAdded << " NumJourneysF: " << numJourneysFinalized << " NumNodesRchable: " << numNodesFinalized << endl;
}

void Graph::printMWFWalksPrioritizedNoSet(int source)
{
    int numRchd=0;
    for (int i=0;i<V;i++)
    {
//        if (get<3>(mwfJourneys[i]) != -1)
//        {
            cout << i << " " << get<0>(mwfJourneys[i]) << " " << get<1>(mwfJourneys[i]) << " " << get<2>(mwfJourneys[i]) << endl;
            numRchd++;
//        }
    }
//    cout << "Rchd: " << numRchd << endl;
}


void Graph::minWaitForemostPrioritized(int source)
{
    Timer t;
    int numNodesFinalized = 0;
    tuple<int,int> nodesRchable_maxFmstTime = earliest_arrival_pair(source,1); //4092652
    int nodesReachable= get<0>(nodesRchable_maxFmstTime), maxFmstTime=get<1>(nodesRchable_maxFmstTime);
    vecFullListPriority.resize(V);
    nodeCompareArrWt heapCompFn;

    bit_queue closedNodes((int)vertices.size());
    vecFullListPriority[source].insert(make_tuple(t_start,0,source,0,0,false));       //arr,wt,prevNd,prevDprt,numHops
    minHeap.push_back(std::make_tuple(t_start, 0,0, source));

    t.start();
    mNumJourExtInst = 0;
    tuple<int,int,int,int> exploreJourney;
    tuple<int,int,int,int,int,bool> journeyToFind;
    set <tuple<int, int, int, int, int, bool>, compareMFSetElementsPrioritized >::iterator itPrioritizedSet;
    while ((numNodesFinalized<nodesReachable) && !(minHeap.empty()))        //Make sure numReachable is set appropriately. TBD
    {
        exploreJourney = minHeap.front();
        std::pop_heap<vector<tuple<int,int,int,int>>::iterator, nodeCompareArrWt>(minHeap.begin(), minHeap.end(), heapCompFn);
        minHeap.pop_back();      //moves the min element to last and then removes it from heap.
        int nodeID=get<3>(exploreJourney);
        get<0>(journeyToFind)=get<0>(exploreJourney);
        get<1>(journeyToFind)=get<1>(exploreJourney);get<2>(journeyToFind)=0;
        get<3>(journeyToFind)=0;get<4>(journeyToFind)=0;get<5>(journeyToFind)=false;
        itPrioritizedSet = vecFullListPriority[nodeID].find(journeyToFind);
        if (itPrioritizedSet == vecFullListPriority[nodeID].end()) {
//            cout <<"Journey in PQ was not found on the node. This means journey was dominated"<<endl;
            if (!closedNodes.check_bit_obj_present(nodeID))
                cout <<"Node is not closed but journey in PQ was dominated. This should never happen. Node:"<<nodeID<<endl;
            continue;
        }
        if (!closedNodes.check_bit_obj_present(nodeID))
        {
            closedNodes.queue_add_bit(nodeID);
            numNodesFinalized++;
        }      //Num finalized journeys should be checked imm. after this, so unnecessary are not expanded. should use do-while TBD.
        extendPrioritizedJourney(exploreJourney, nodeID, source);
/**
 for each neighbor v of u
 {       generate necessary expansions (v,a',w') from u to v;
     If (v,a',w') is not dominated by journeys in FL[v]
     {        Insert (v,a',w') in FL[v];
         Mark all the subsequent journeys FL[v] that are dominated by (a',w') as dominated;
         Insert (v,a',w') in the PQ;
     }
 }

 */
    }
    t.stop();
    time_sum += t.GetRuntime();
    printMWFWalksPrioritized(source);
}

void Graph::printMWFWalksPrioritized(int source)
{
    int prevNode, departTime, arrTime, waitTime, currNode, numRchableVerts = 0;
    set <tuple<int, int, int, int, int, bool>, compareMFSetElementsPrioritized >::iterator itMWFJourney;
    cout << "output format: dest(arr,wait_time)<--(depTime)node..\n";
    for (int i=0; i < V; i++)
    {
        if (i==source)
        {
            numRchableVerts++;
            cout << source << " 0 0" << endl;
            continue;
        }
        currNode = i;
        if (vecFullListPriority[currNode].empty())
            continue;
        numRchableVerts++;
        itMWFJourney = vecFullListPriority[currNode].begin();
        cout << "Walk: "<< currNode << " ";// << "("<< get<0>(*itMWFJourney) <<","<<get<1>(*itMWFJourney)<<")";//\n";
        arrTime = get<0>(*itMWFJourney); waitTime = get<1>(*itMWFJourney);
        //cout << arrTime*100+waitTime << "\n";
        prevNode = get<2>(*itMWFJourney); departTime = get<3>(*itMWFJourney);
        //cout << "("<< arrTime <<","<<waitTime<<")";
        //cout << "<--"  << "(" << departTime <<")"<< prevNode;
        
        //cout <<  i << " " << arrTime << " " << waitTime << " " << endl;
        while (!((prevNode == source) && (currNode == source)))
        {
            cout << "("<< arrTime <<","<<waitTime<<")";
            cout << "<--"  << "(" << departTime <<")"<< prevNode;
            currNode=prevNode;
            itMWFJourney = vecFullListPriority[currNode].lower_bound(make_tuple(departTime,0,0,0,0,false));
            if (itMWFJourney == vecFullListPriority[prevNode].end())
                --itMWFJourney;
            else if (departTime < get<0>(*itMWFJourney))
                --itMWFJourney;
            
            arrTime = get<0>(*itMWFJourney);
            waitTime = get<1>(*itMWFJourney);
            prevNode = get<2>(*itMWFJourney);
            departTime = get<3>(*itMWFJourney);
            if(vecFullListPriority[prevNode].empty())
            {
                cout << "Something went wrong\n\n";
                break;
            }
        }
        cout << "("<< arrTime <<","<<waitTime<<")";
        cout << "<--"  << "(" << departTime << ")" << prevNode << "\n\n";
    }
    cout << "Num Reachable: " << numRchableVerts << endl;
}


void Graph::extendPrioritizedJourney(tuple<int,int, int,int> inJourney, int node, int source)
{
    int nbr, intvlID, departTime, newArrivTime, newWaitTime, numHops;
    bool isDominated;
    for (int j=0; j<vertices[node].numNbrs;j++)
    {
        nbr = vertices[node].neighbors[j].nbrId;
        departTime = earliestUseEdgeAfterT(node, vertices[node].neighbors[j], get<0>(inJourney), intvlID);
        if ( (departTime == -1) || (departTime >= infinity) || (intvlID == -1) ) //or departTime > max fmst arrival time. EARLYTMNT
            continue;
//        do {      // UNCOMMENT
            newArrivTime = departTime + vertices[node].neighbors[j].edgeSchedules[intvlID].traveTime;
            if ((node == source) && (get<0>(inJourney)==t_start))  //this is nbr of source & start of the journey
                newWaitTime = 0;
            else
                newWaitTime = departTime-get<0>(inJourney)+get<1>(inJourney);   //get<0> is arrival time at node and get<1> is wait time till node.
            numHops = get<2>(inJourney)+1;
            isDominated = checkDominatedAndInsertPrioritized(nbr,make_tuple(newArrivTime,newWaitTime, node, departTime, numHops,false));     //Check if this is dominated by the journeys already present at nbr.
//            if (!isDominated)
  //              toExpandList.insert(nbr);
            intvlID++;
            if (intvlID < vertices[node].neighbors[j].numIntvls)
                departTime = vertices[node].neighbors[j].edgeSchedules[intvlID].intvlStart;
                
  //      } while (intvlID < vertices[node].neighbors[j].numIntvls);
    }
}

bool Graph::checkDominatedAndInsertPrioritized(int destNode,tuple<int, int, int, int, int, bool> newJourney)
{
    if (vecFullListPriority[destNode].empty())
    {
        insertPrioritizedInJourneySets(destNode, newJourney, vecFullListPriority[destNode].begin());
        return false;
    }
    set <tuple<int, int, int, int, int, bool>, compareMFSetElementsPrioritized >::key_compare compUsed = vecFullListPriority[destNode].key_comp();
    set <tuple<int, int, int, int, int, bool>, compareMFSetElementsPrioritized >::iterator itJourneySet = vecFullListPriority[destNode].lower_bound(newJourney);
    
    if (itJourneySet == vecFullListPriority[destNode].end())
    {
        --itJourneySet;
        if (checkDominance(make_tuple(get<0>(*itJourneySet), get<1>(*itJourneySet)), make_tuple(get<0>(newJourney), get<1>(newJourney)) ) )
            return true;
        else
        {
            insertPrioritizedInJourneySets(destNode, newJourney, vecFullListPriority[destNode].end());
            return false;
        }
    }
    if (compUsed(newJourney,*itJourneySet)) //newJourney < *itjourneyset. This ensures if same journey present it is not inserted again.
    {
        if (itJourneySet==vecFullListPriority[destNode].begin())        //First journey is arriving later, so new journey can't be dominated.
        {
            insertPrioritizedInJourneySets(destNode, newJourney, vecFullListPriority[destNode].begin());
            return false;
        }
        else
        {
            --itJourneySet;     //Journey previous to the one that was first to be bigger than newJourney.
            if ( checkDominance(make_tuple(get<0>(*itJourneySet), get<1>(*itJourneySet)), make_tuple(get<0>(newJourney), get<1>(newJourney)) ) )
                return true;
            else
            {
                insertPrioritizedInJourneySets(destNode, newJourney, ++itJourneySet);          //Insert before the iterator itJourneySet.
                return false;
            }
        }
    }
    else        //Another journey with exact same arrival time and wait time already present. So consider this dominated as it is not a new journey.
        return true;
}

//Insert in both sets
//Erase all subsequent journeys that are now dominated by this.
void Graph::insertPrioritizedInJourneySets(int destNode, tuple<int,int,int,int, int,bool> newJourney, set <tuple<int, int, int,int, int, bool>, compareMFSetElementsPrioritized >::iterator itInsertPos)
{
    nodeCompareArrWt heapCompFn;
    set <tuple<int, int, int,int, int, bool>, compareMFSetElementsPrioritized >::iterator nextPos =  vecFullListPriority[destNode].insert(itInsertPos, newJourney);
    //Also insert it into the priority queue.//****************//
    minHeap.push_back(make_tuple(get<0>(newJourney), get<1>(newJourney), get<4>(newJourney), destNode));
    std::push_heap(minHeap.begin(), minHeap.end(), heapCompFn);
    
    tuple<int, int> newJourneySansPrev = make_tuple(get<0>(newJourney), get<1>(newJourney));
    bool bDominates = true;
    tuple<int, int, int,int, int, bool> nextJourney;
    ++nextPos;
    while ((nextPos != vecFullListPriority[destNode].end()) && bDominates)
    {
        nextJourney = *nextPos;
        if (checkDominance(newJourneySansPrev,
                           make_tuple(get<0>(nextJourney), get<1>(nextJourney)) ) )
            nextPos = vecFullListPriority[destNode].erase(nextPos); //If they are erased from the set. They won't be found when examined from heap. This means they were dominated.
        else
            bDominates = false;
    }
#ifdef __TEST__
    /*
    if (vecFullListPriority[destNode].size() > vertices[destNode].inDegree)
        cout << "Too many in journeys at: " << destNode << ". InDegree=" << vertices[destNode].inDegree << ". in journeys="  << get<1>(vecFullList[destNode]).size() << endl;
     */
    mNumJourExtInst++;
    int numHops=get<4>(newJourney);
/*    if ((mNumJourExtInst % JOURCHECK) == 0)
    {
        cout << "At jour: " << mNumJourExtInst << "At dest: " << destNode << " NumHops=" << numHops << " PrevNode:"<<get<2>(newJourney)<<" DepTm:"<<get<3>(newJourney)
        <<" ArrTm:"<<get<0>(newJourney)<<" WtTime:"<<get<1>(newJourney)<<endl;
    }
 */

//    if ((numHops % HOPCHECK) == 0)
//        cout << "At dest: " << destNode << " NumHops=" << numHops << " PrevNode:"<<get<2>(newJourney)<<" DepTm:"<<get<3>(newJourney)
//                <<" ArrTm:"<<get<0>(newJourney)<<" WtTime:"<<get<1>(newJourney)<<endl;
    
#endif

}

























//Every journey has number of hops associated with it.
//So append #of hops for every journey.
//When journey is picked up for extension, check #of hops and print every million multiple hops.
void Graph::extendJourney(tuple<int,int, int> inJourney, int node)
{
    int nbr, intvlID, departTime, newArrivTime, newWaitTime, numHops;
    bool isDominated;
    for (int j=0; j<vertices[node].numNbrs;j++)
    {
        nbr = vertices[node].neighbors[j].nbrId;
        departTime = earliestUseEdgeAfterT(node, vertices[node].neighbors[j], get<0>(inJourney), intvlID);
        if ( (departTime == -1) || (departTime >= infinity) || (intvlID == -1) ) //or departTime > max fmst arrival time. EARLYTMNT
            continue;
        do {      //TBD UNCOMMENT
            newArrivTime = departTime + vertices[node].neighbors[j].edgeSchedules[intvlID].traveTime;
            newWaitTime = departTime-get<0>(inJourney)+get<1>(inJourney);   //get<0> is arrival time at node and get<1> is wait time till node.
            numHops = get<2>(inJourney)+1;
            isDominated = checkDominatedAndInsert(nbr,make_tuple(newArrivTime,newWaitTime, node, departTime, numHops));     //Check if this is dominated by the journeys already present at nbr.
//            if (!isDominated)
  //              toExpandList.insert(nbr);
            intvlID++;
            if (intvlID < vertices[node].neighbors[j].numIntvls)
                departTime = vertices[node].neighbors[j].edgeSchedules[intvlID].intvlStart;
                
        } while (intvlID < vertices[node].neighbors[j].numIntvls);
    }
}

bool Graph::checkDominatedAndInsert(int destNode,tuple<int, int, int, int, int> newJourney)
{
    if (get<1>(vecFullList[destNode]).empty())
    {
        insertInJourneySets(destNode, newJourney, get<1>(vecFullList[destNode]).begin());
        return false;
    }
    set <tuple<int, int, int, int, int>, compareMFSetElements >::key_compare compUsed = get<1>(vecFullList[destNode]).key_comp();
    set <tuple<int, int, int, int, int>, compareMFSetElements >::iterator itJourneySet = get<1>(vecFullList[destNode]).lower_bound(newJourney);
    if (itJourneySet == get<1>(vecFullList[destNode]).end())
    {
        --itJourneySet;
        if (checkDominance(make_tuple(get<0>(*itJourneySet), get<1>(*itJourneySet)), make_tuple(get<0>(newJourney), get<1>(newJourney)) ) )
            return true;
        else
        {
            insertInJourneySets(destNode, newJourney, get<1>(vecFullList[destNode]).end());
            return false;
        }
    }
    if (compUsed(newJourney,*itJourneySet))     //newJourney < *itjourneyset.
    {
        if (itJourneySet==get<1>(vecFullList[destNode]).begin())        //First journey is arriving later, so new journey can't be dominated.
        {
            insertInJourneySets(destNode, newJourney, get<1>(vecFullList[destNode]).begin());
            return false;
        }
        else
        {
            --itJourneySet;     //Journey previous to the one that was first to be bigger than newJourney.
            if ( checkDominance(make_tuple(get<0>(*itJourneySet), get<1>(*itJourneySet)), make_tuple(get<0>(newJourney), get<1>(newJourney)) ) )
                return true;
            else
            {
                insertInJourneySets(destNode, newJourney, ++itJourneySet);          //Insert before the iterator itJourneySet.
                return false;
            }
        }
    }
    else        //Another journey with exact same arrival time and wait time already present. So consider this dominated as it is not a new journey.
        return true;
}

//Insert in both sets
//Erase all subsequent journeys that are now dominated by this.
void Graph::insertInJourneySets(int destNode, tuple<int,int,int,int, int> newJourney, set <tuple<int, int, int,int, int>, compareMFSetElements >::iterator itInsertPos)
{
    set <tuple<int, int, int,int, int>, compareMFSetElements >::iterator nextPos =  get<1>(vecFullList[destNode]).insert(itInsertPos, newJourney);
    tuple<int, int> newJourneySansPrev = make_tuple(get<0>(newJourney), get<1>(newJourney));
    tuple<int, int, int> newJourneySansPrevWHops = make_tuple(get<0>(newJourney), get<1>(newJourney), get<4>(newJourney));
    bool bDominates = true;
    tuple<int, int, int,int, int> nextJourney;
    ++nextPos;
    while ((nextPos != get<1>(vecFullList[destNode]).end()) && bDominates)
    {
        nextJourney = *nextPos;
        if (checkDominance(newJourneySansPrev,
                           make_tuple(get<0>(nextJourney), get<1>(nextJourney)) ) )
            nextPos = get<1>(vecFullList[destNode]).erase(nextPos);
        else
            bDominates = false;
    }
#ifdef __TEST__
    /*
    if (get<1>(vecFullList[destNode]).size() > vertices[destNode].inDegree)
        cout << "Too many in journeys at: " << destNode << ". InDegree=" << vertices[destNode].inDegree << ". in journeys="  << get<1>(vecFullList[destNode]).size() << endl;
     */
    mNumJourExtInst++;
    if ((mNumJourExtInst % JOURCHECK) == 0)
    {
        cout << "At dest: " << mNumJourExtInst <<endl;
    }

    int numHops=get<4>(newJourney);
    if ((numHops % HOPCHECK) == 0)
        cout << "At dest: " << destNode << " NumHops=" << numHops << " PrevNode:"<<get<2>(newJourney)<<" DepTm:"<<get<3>(newJourney)
                <<" ArrTm:"<<get<0>(newJourney)<<" WtTime:"<<get<1>(newJourney)<<endl;
    
#endif

    if (get<0>(vecFullList[destNode]).empty())
        toExpandList.push_back(destNode);
    pair<set <tuple<int, int, int>>::iterator, bool> newJourenysInsertRet = get<0>(vecFullList[destNode]).insert(newJourneySansPrevWHops);
    set<tuple<int, int, int>>::iterator itNextNewJourney = get<0>(newJourenysInsertRet);
    ++itNextNewJourney;
    tuple<int, int, int> nextNewJourney;
    bDominates = true;
    if (get<1>(newJourenysInsertRet))
    {
        while ((itNextNewJourney != get<0>(vecFullList[destNode]).end()) && bDominates)
        {
            nextNewJourney = *itNextNewJourney;
            if (checkDominance(newJourneySansPrev, make_tuple(get<0>(nextNewJourney), get<1>(nextNewJourney)) ))
                itNextNewJourney = get<0>(vecFullList[destNode]).erase(itNextNewJourney);
            else
                bDominates = false;
        }
    }
}

bool Graph::checkDominance(tuple<int, int> journey1,  tuple<int, int> journey2)
{
    int tA = get<0>(journey1), wA = get<1>(journey1), tB = get<0>(journey2), wB = get<1>(journey2);
    return (tB-tA+wA <= wB);
}














//----------------------------
void Graph::build_shortestJourneys(int source, vector<std::pair<int, int>>& shortestJourneyPointer, vector<vector<incrementalJourney>>& allHopJourneys)
{
    shortestJourneys.resize(V);                     //Shortest Journeys for all nodes needs to be found.
    shortestJourneys[source].rPath.push_back(std::make_tuple(source, -1, -1));       //start node, nbr ID.
    shortestJourneys[source].sigmaSchedule.push_back(std::make_tuple(0, -1, -1));        //depart time, IntvlId
    for (int node = 0; node < V; node++)
    {
        if (node == source)
            continue;
        int hopCount = std::get<0>(shortestJourneyPointer[node]);
        int incrIndex = std::get<1>(shortestJourneyPointer[node]);
        for (int shortestHopCount = hopCount; shortestHopCount > 0; shortestHopCount--)  //BuildShortestJourney for nextNode from Incremental Vectors.
        {
            incrementalJourney *pJourneyPath = &(allHopJourneys[shortestHopCount][incrIndex]);    //start node, it's nbr id.
            int currentNode = pJourneyPath->currentNodeID;int startNode = pJourneyPath->prevNodeID; int edgeId = pJourneyPath->prevEdgeID;
            int currentArrivalTime = pJourneyPath->currentArrivalTime;int depart = pJourneyPath->prevDepartureTime; int intvlUsed = pJourneyPath->prevIntvlID;
            shortestJourneys[node].rPath.push_back(std::make_tuple(currentNode, startNode, edgeId));
            shortestJourneys[node].sigmaSchedule.push_back(std::make_tuple(currentArrivalTime, depart, intvlUsed));    //dep. time on edge (node, nextnode), intvl used.
            incrIndex = pJourneyPath->indexPrevIncrement;
        }
    }
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

void Graph::print_shortest_results_test(int source)
{
    ofstream shortestOut(shortestResults);
    int maxHopCount = 0; int vertWHops = 1; int sumHops = 0; int avgHops; int arrivalTime; int numHops;
    shortestOut << V << "\n";
    for (int i = 0; i < V; i++)
    {
        if ( (i == source) || (shortestJourneys[i].rPath.size() > 0) )
        {
            numHops = (int)shortestJourneys[i].rPath.size()-1;
            arrivalTime = get<0>(shortestJourneys[i].sigmaSchedule[0]);
            shortestOut <<  i << "  " << arrivalTime << " " << shortestJourneys[i].rPath.size() <<  "\n";
            
        }
        else
            continue;
        sumHops += (int)shortestJourneys[i].rPath.size();
        vertWHops++;
        if (shortestJourneys[i].rPath.size() > maxHopCount)
            maxHopCount = (int)shortestJourneys[i].rPath.size();
/*#ifdef __TEST__
        for (int hopCount = (int)(shortestJourneys[i].rPath.size()) - 1; hopCount >= 0; hopCount--)
        {
            std::tuple<int, int, int> *pCurrentPrevNodesEdge = &(shortestJourneys[i].rPath[hopCount]);
            std::tuple<int, int, int> *pArriveDepartIntvlUsed = &(shortestJourneys[i].sigmaSchedule[hopCount]);
            int toNode = std::get<0>(*pCurrentPrevNodesEdge);
            int prevNode = std::get<1>(*pCurrentPrevNodesEdge);

            int arrivalTime = std::get<0>(*pArriveDepartIntvlUsed);
            int departTime = std::get<1>(*pArriveDepartIntvlUsed);

            shortestOut << " (" << (prevNode) << "-->" << toNode << ")"
                << ", (d: " << departTime << ", a:" << arrivalTime << ") ";
        }
        shortestOut << "\n";
#endif*/
    }
    avgHops = sumHops/vertWHops;
    cout << "Max Hops = " << maxHopCount << "\n";
    cout << "Avg # hops = " << avgHops << "\n";
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
    cout <<"Num of Sources tried on: " << sources.size() << endl;
    cout<<"Average time: " << time_sum/sources.size() <<endl;
}

void Graph::print_avg_time(const char* filePath1, const char* filePath2)
{
    FILE* file = fopen(filePath2,"a");
    
    fprintf(file, "%s\t%f\n", filePath1, time_sum/100);;
    fclose(file);
}
