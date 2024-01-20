# GraphB++: Scaling Frustration Index and Balanced State Discovery for Real Signed Graphs

Authors: <em> Muhieddine Shebaro, Jelena Tešić </em>

![Highland Tribes Execution!](/images/fig22.png "Example")

**GraphB++** (2023) extends GraphB+’s code to scale the approximation of the frustration index (NP-hard) for any real world signed network of any size or density based on efficiently finding a fundamental cycle basis. It also extends the frustration cloud to a (key,value) tuple collection F = B:(C, S). The states with their associated frequency and edge switches will be stored as a tuple in the memory-bound frustrated cloud in a scalable manner. We provide benchmark comparison of seven spanning tree sampling methods, including a breakthrough study of frustration index computation and timing for networks.

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

* Simply download one of the .cpp files and compile the source file like the following:

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


