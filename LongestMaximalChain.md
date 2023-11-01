## Rough algo

1. For each point, find all spacelike points to it. A spacelike point is any such that neither are related by $\prec$.

This is just the problem of finding the longest path in an undirected graph.

Turns out that a 2 BFSs are sufficient:
https://www.geeksforgeeks.org/longest-path-undirected-tree/

Turns out it isn't - DFS from everywhere.