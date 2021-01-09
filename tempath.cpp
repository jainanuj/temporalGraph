#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include "graph.h"
using namespace std;

//./temppath earliest

int main(int argc, char* argv[]){
	
	char* option = argv[1];
	
	
	Graph g(argv[2]);
	g.initial_query();


	if(!strcmp(option,"earliest"))
    {
		g.run_earliest_arrival();
    }
/*    else if (!strcmp(option,"latest"))
    {
		g.run_latest_departure();
    }
    else if (!strcmp(option,"shortest"))
    {
		g.run_shortest();
    }
    else if(!strcmp(option,"fastest"))
    {
		g.run_fastest();
    }
	*/
	cout<<"hello world!!!" <<endl;
	
	return 0;
}



