For generating graphs with wu model and Xuan model from Koblenz networks. First line will b dropped as comment: <br>
./XuantemporalGraph wu \<filename\> <br>

For the koblenz networks where first 2 lines have to be dropped. <br>
./XuantemporalGraph wu \<filename\> \<1/2\>(drop num Lines to drop) \<0/1\>(normalize or not) <br>

If only first line has to be dropped than the first command is good enough. <br>
Normalize option means that the min timestamp starts at 0. So the min timestamp is subratcted from all timestamps. <br>

This runs algorithms for finding different types of paths in the temporal graph.
Graph is expressed in interval temporal format.
The different algorithms supported right now: <br>
earliest <br>
shortest <br>
mhf <br>
mwf <br>
hbhshrtst  <br>

shortest is actually min-hop. <br>

Graph can be in intvl format as:
u v numIntvls s1 e1 d1 s2 e2 d2 .....
or in contact sequence aggregate format as:
u v numIntvls c1 d1 c2 d2 ....

example, intvl format:
1 2 3 1 4 1 5 8 1 10 12 2

cs format:
1 2 3 1 1 5 1 10 2

command to generate CSG format and Intvl format from downloaded koblenz datasets
./XuantemporalGraph wu \<fileName\> \<1/2\>(drop num Lines to drop from top)  \<0/1\>(normalize or not) <br>

command to generat underlying static graph given an intvl temporal graph
./XuantemporalGraph static \<fileName_xuanOp.txt\> <br>

command to run with intvl format: <br>
./XuantemporalGraph earliest|shortest|mhf|mwf|hbhshrtst \<filename\> <br>

command to run with cs format: <br>
./XuantemporalGraph earliest|shortest \<filename\> cs <br>

To run with Test mode, for testign with just 1 source and comparing the outputs.
Uncomment #define __TEST__ on the top of both Wu and XuanTemporalGrpahs.
Compile and then just run with the same commands as above
