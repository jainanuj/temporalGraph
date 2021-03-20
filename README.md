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

command to run with intvl format:
./XuantemporalGraph earliest|shortest <filename> cs
