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
    int divTime;
};

struct mwfJourney {
    int arrivalTime;
    int wtTime;
    int prevNode;
    int prevJourneyIndex;
    int prevDepTime;
    /** required only for classes*/
    bool isClass;
    int arrivalTimeEnd;
    int lastTravelTime;
    int prevNodeNbrIndex;
//    vector<tuple<int, int, int>> lastExpandedAt;     //(intvl.startTime,lambda) of intvl in which this journey was expanded on each nbr. Last element is flag whether it was expanded or not.
};

//journeyClass:
//start, end, prevLambda, prevNbrIndex,
struct mwfJourneyClass {
    int arrivalTimeStart;
    int arrivalTimeEnd;
    int prevNodeNbrIndex;
    int lastTravelTime;
    int wtTime;
    int prevNode;
    int prevJourneyIndex;
    int prevDepTime;
    vector<tuple<int, int, int>> lastExpandedAt;     //(intvl.startTime,lambda) of intvl in which this journey was expanded on each nbr. Last element is flag whether it was expanded or not.
};

//Carryover. - shouldn't say prevJourneyIndex. (cuz it may have changed).

class GraphDualCriteria  : public Graph
{
public:
    GraphDualCriteria(const char* filePath, int contactSeq, const char* option);
    bool lessCompArrHop(std::pair<int, int> newArrHop, std::pair<int,int> oldArrHop);
    void initial_ds_eha();
    void initial_ds_ewa();
    void initial_ds_ewa_classes();
    void initial_ds_hbhshrtst();

    int earliest_arrival_minHop_pair(int source);
    void run_mhf();
    void run_mwf();
    void run_hbh_shortest();
    void printmhfResultsTest(int source);
    void mhfHopByHop(int source);
    void shortestHopByHop(int source);
    void mwfStreamingIntvls(int source);
    void mwfStreamingIntvlsWJourneyClasses(int source);
    void mwfStreamingIntvlsV2(int source);
    
    void build_mhf_Journeys(int source, vector<std::tuple<int, int, int>>& mhfJourneyPointer, vector<vector<incrementalJourney>>& allHopJourneys);
    void printmhfResultsTest2(int source, vector<std::tuple<int, int, int>>& mhfJourney);
    void printmwfResultsTest2(int source);
    void printmwfClassesResultsTest2(int source);
    int removeMinIntvl(intervalInfo &newIntvl, int& indexPreKnownIntvls);
    int searchPrevJourney(int node, int beforeTime);
    int searchPrevJourneyClass(int node, int beforeTime);
    void setupNewJourney(int v, int arrivalTime);
    void setupNewJourneyClass(int v, int arrivalTime);
    int getPrevJourney(intervalInfo& intvl);
    int getPrevJourneyClass(intervalInfo& intvl);
    int checkNewJourneyAndInsert(mwfJourney& newJourney, int v);
    int checkNewJourneyAndInsertV2(mwfJourney& newJourney, int v);
    int checkNewJourneyClassAndInsert(mwfJourneyClass& newJourneyClass, int v);
    bool resolveOverlap(mwfJourneyClass& lastJClass, mwfJourneyClass& newJourneyClass, mwfJourneyClass& carryOverJClass);
    bool checkJourneyClassDominance(mwfJourneyClass& firstJClass, mwfJourneyClass& nextJClass);
    void buildAndPushIntvlInHeap(mwfJourneyClass& carryOverJClass, int v);
    void buildAndPushIntvlInHeapV2(mwfJourney &carryOverJ, int v);
    int createNewJourneyClass(mwfJourneyClass& prevJourneyClass, intervalInfo& intervalToExpand, mwfJourneyClass& newJourenyClass, int prevJourenyClassIndex);
    void recordFinalJourney(int v, mwfJourneyClass& newJourneyClass);
    
    bool checkShrtstDominanceAndPush(list<incrementalShortestJourney>& listShrtstJrnys, incrementalShortestJourney latestShrtstJrny,
                                     list<incrementalShortestJourney>::iterator posList);
    bool mergeJourneys(list<incrementalShortestJourney>& toList, list<incrementalShortestJourney>& fromList);
    void print_shortest_hbh_results(int source);

    
public:
    vector <pair<int,int>> arr_hop_time, f_time;

    vector<vector<mwfJourney>> listJourneys;    //List of journeys at each node.
    vector<intervalInfo> listOfPreKnownIntvls;
    vector<intervalInfo> listOfAdHocIntvls;
    vector<mwfJourney> finalMWFJourneys;
    vector<mwfJourneyClass> finalMWFJourneyClass;

    vector<vector<mwfJourneyClass>> listJourneyClasses;    //List of journeys at each node.
    int max_runTime, max_Src;
    
    vector<incrementalShortestJourney> shrtstPathAllVertices; //At each vertex there is a vector of shortest journeys discovered by the last hop.
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

/*
((u==10122) && (v==7239)) || ((u==2043) && (v==10122)) || ((u==12028) && (v==2043)) || ((u==3504) && (v==12028)) || ((u==73470) && (v==3504)) || ((u==6415) && (v==73470)) || ((u==67355) && (v==6415)) || ((u==0) && (v==67355)) || ((u==58668) && (v==0)) || ((u==113768) && (v==58668)) || ((u==7372) && (v==113768)) || ((u==18428) && (v==7372)) || ((u==200712) && (v==18428)) || ((u==9571) && (v==200712)) || ((u==30189) && (v==9571)) || ((u==196932) && (v==30189)) || ((u==3574) && (v==196932)) || ((u==11408) && (v==3574)) || ((u==6977) && (v==11408)) || ((u==1931) && (v==6977)) || ((u==1092) && (v==1931)) || ((u==11408) && (v==1092)) || ((u==7709) && (v==11408)) || ((u==276) && (v==7709)) || ((u==16100) && (v==276)) || ((u==14487) && (v==16100)) || ((u==6) && (v==14487)) || ((u==1926) && (v==6)) || ((u==2164) && (v==1926)) || ((u==7710) && (v==2164)) || ((u==7722) && (v==7710)) || ((u==8453) && (v==7722)) || ((u==956) && (v==8453)) || ((u==4590) && (v==956)) || ((u==10647) && (v==4590)) || ((u==1934) && (v==10647)) || ((u==20971) && (v==1934)) || ((u==1919) && (v==20971)) || ((u==7476) && (v==1919)) || ((u==9362) && (v==7476)) || ((u==24957) && (v==9362)) || ((u==11720) && (v==24957)) || ((u==29265) && (v==11720)) || ((u==3055) && (v==29265)) || ((u==4147) && (v==3055)) || ((u==16) && (v==4147)) || ((u==6963) && (v==16)) || ((u==145611) && (v==6963)) || ((u==xxxx) && (v==145611)) || ((u==24684) && (v==145611)) || ((u==51129) && (v==24684)) || ((u==5481) && (v==51129)) || ((u==1709) && (v==5481)) || ((u==29386) && (v==1709)) || ((u==8273) && (v==29386)) || ((u==0) && (v==8273)) || ((u==1) && (v==0)) || ((u==0) && (v==1)) )
*/

/*((u==10260) && (v==15260)) || ((u==14557) && (v==10260)) || ((u==10122) && (v==7239)) || ((u==2043) && (v==10122)) || ((u==12028) && (v==2043)) || ((u==3504) && (v==12028)) || ((u==73470) && (v==3504)) || ((u==6415) && (v==73470)) || ((u==67355) && (v==6415)) || ((u==0) && (v==67355)) || ((u==58668) && (v==0)) || ((u==113768) && (v==58668)) || ((u==7372) && (v==113768)) || ((u==18428) && (v==7372)) || ((u==200712) && (v==18428)) || ((u==9571) && (v==200712)) || ((u==30189) && (v==9571)) || ((u==196932) && (v==30189)) || ((u==3574) && (v==196932)) || ((u==11408) && (v==3574)) || ((u==6977) && (v==11408)) || ((u==1931) && (v==6977)) || ((u==1092) && (v==1931)) || ((u==11408) && (v==1092)) || ((u==7709) && (v==11408)) || ((u==276) && (v==7709)) || ((u==16100) && (v==276)) || ((u==14487) && (v==16100)) || ((u==6) && (v==14487)) || ((u==1926) && (v==6)) || ((u==2164) && (v==1926)) || ((u==7710) && (v==2164)) || ((u==7722) && (v==7710)) || ((u==8453) && (v==7722)) || ((u==956) && (v==8453)) || ((u==4590) && (v==956)) || ((u==10647) && (v==4590)) || ((u==1934) && (v==10647)) || ((u==20971) && (v==1934)) || ((u==1919) && (v==20971)) || ((u==7476) && (v==1919)) || ((u==9362) && (v==7476)) || ((u==24957) && (v==9362)) || ((u==11720) && (v==24957)) || ((u==29265) && (v==11720)) || ((u==3055) && (v==29265)) || ((u==4147) && (v==3055)) || ((u==16) && (v==4147)) || ((u==6963) && (v==16)) || ((u==145611) && (v==6963)) || ((u==24684) && (v==145611)) || ((u==51129) && (v==24684)) || ((u==5481) && (v==51129)) || ((u==1709) && (v==5481)) || ((u==29386) && (v==1709)) || ((u==8273) && (v==29386)) || ((u==0) && (v==8273)) || ((u==1) && (v==0)) || ((u==0) && (v==1))*/


/*
<--(51404400)8273(51231601,0)<--(51231600)0(50540401,50540399)<--(50540400)1(1,0)<--(0)0(0,0)<--(0)
*/


#endif /* graphDualCriteria_hpp */
