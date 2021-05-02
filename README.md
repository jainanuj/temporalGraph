For generating graphs with wu model and Xuan model from Koblenz networks:
./XuantemporalGraph wu <filename>
For the koblenz networks where first 2 lines have to be dropped.
./XuantemporalGraph wu <filename> 2

If only first line has to be dropped than the first command is good enough.


This runs algorithms for finding different types of paths in the temporal graph.
Graph is expressed in interval temporal format.
The different algorithms supported right now:
earliest
shortest

Graph can be in intvl format as:
u v numIntvls s1 e1 d1 s2 e2 d2 .....
or in contact sequence aggregate format as:
u v numIntvls c1 d1 c2 d2 ....

example, intvl format:
1 2 3 1 4 1 5 8 1 10 12 2

cs format:
1 2 3 1 1 5 1 10 2

command to run with intvl format:
./XuantemporalGraph earliest|shortest <filename>

command to run with cs format:
./XuantemporalGraph earliest|shortest <filename> cs

To run with Test mode, for testign with just 1 source and comparing the outputs.
Uncomment #define __TEST__ on the top of both Wu and XuanTemporalGrpahs.
Compile and then just run with the same commands as above
