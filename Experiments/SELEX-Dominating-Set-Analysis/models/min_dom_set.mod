# A model for 
# Minimum dominating set of a graph
# http://www.nada.kth.se/~viggo/wwwcompendium/node11.html
# To run solver using the test data below
# glpsol -m min_dom_set.mod -o test0.sol
# To run solver with separare data file
# glpsol -m min_dom_set.mod -d ../data/test.dat -o test0.sol

/* sets */ 
set NODES;

/* parameters */
param AdjTable {i in NODES, j in NODES};

/* decision variables: v[i]: whether to include a node in the dominating set */
var v {j in NODES} >= 0, integer;

/* objective function */
/* Minimize the size of the dominating set */
minimize z: sum{i in NODES} v[i];

/* constraints */
/* dominating: for each node, either itself of one of each neighbors has to be in the cover */
s.t. dom {i in NODES} : sum{j in NODES} AdjTable[i,j]*v[j] >= 1;


data;

set NODES      := "A" "B" "C" "D" "E";
param AdjTable:"A"	"B"	"C"	"D"	"E":= 
"A"	1	1	0	1	1
"B"	1	1	1	0	0
"C"	0	1	1	1	0
"D"	1	0	1	1	1
"E"	1	0	0	1	1;



 
