check push <br>
./XuantemporalGraph wu \<filename\> //generate graphs with wu model and Xuan model from Koblenz n/ws. First line will b dropped as comment:  <br> 

./XuantemporalGraph wu \<filename\> \<1/2\>(drop num Lines to drop) \<0/1\>(normalize or not) //same as above with num lines dropped and possible ts normalized output <br>

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

./XuantemporalGraph static \<fileName_xuanOp.txt\>  //generate underlying static graph given an intvl temporal graph<br>

./XuantemporalGraph earliest|shortest|mhf|mwf|hbhshrtst \<filename\> //run algos with intvl format <br>

./XuantemporalGraph earliest|shortest \<filename\> cs //run algos with cs intvl format <br>

To run with Test mode, for testign with just 1 source and comparing the outputs.
Uncomment #define __TEST__ on the top of both Wu and XuanTemporalGrpahs.
Compile and then just run with the same commands as above
