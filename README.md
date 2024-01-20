# ABCD: Algorithm For Balanced Component Discovery in Signed Networks

Authors: <em> Muhieddine Shebaro, Jelena Tešić </em>

![Highland Tribes Execution!](/images/updated_example1.png "Example")

**ABCD** (2023) is a heuristic searching algorithm in signed graphs that tackles the NP-hard problem of finding the largest balanced subgraph. It builds on the scalable discovery of fundamental cycles and utilizes the graph's node density distribution and near-optimal balanced states to minimize the number of vertices removed from the balanced sub-graph.

## Citation
Please cite the following publication: Shebaro, M., & Tesic, J. (2023). ABCD: Algorithm for Balanced Component Discovery in Signed Networks. ArXiv, abs/2311.00848.

**BibTeX entry:**
```
@misc{shebaro2023abcd,
      title={ABCD: Algorithm for Balanced Component Discovery in Signed Networks}, 
      author={Muhieddine Shebaro and Jelena Tešić},
      year={2023},
      eprint={2311.00848},
      archivePrefix={arXiv},
      primaryClass={cs.SI}
}
```

## How to Run the Code 

* Simply download ABCD.cpp compile the source file like the following:

```
user:~$ g++ -fopenmp ABCD.cpp
```
The signed graph must be in the following format to be compatible with graphC (src,dst,sign).
Preprocessing of the signed graph is embedded.

* To execute the compiled file, ABCD utilizes 2 parameters in this order (iteration_count K):
```
user:~$ ./a.out input.txt 1000 500
```
where K < iteration_count.
 
* **Output:** ABCD prints the vertex cardinality of the largest balanced subgraph found.

**Note:** if you run into errors related to the stack memory please run this command before executing the code:
```
user:~$ ulimit -s unlimited
```

[Data Lab @ TXST](DataLab12.github.io)


