# GraphB++: Scaling Frustration Index and Balanced State Discovery for Real Signed Graphs

Authors: <em> Muhieddine Shebaro, Jelena Tešić </em>

![Highland Tribes Execution!](/images/fig22.png "Example")

**GraphB++** (2023) extends GraphB+’s code to scale the approximation of the frustration index (NP-hard) for any real world signed network of any size or density based on efficiently finding a fundamental cycle basis. It also extends the frustration cloud to a (key,value) tuple collection F = B:(C, S). The states with their associated frequency and edge switches will be stored as a tuple in the memory-bound frustrated cloud in a scalable manner. We provide benchmark comparison of seven spanning tree sampling methods, including a breakthrough study of frustration index computation and timing for networks.

## Citation
Please cite the following publication: Shebaro, M., & Tesic, J. (2023). Scaling Frustration Index and Corresponding Balanced State Discovery for Real Signed Graphs. ArXiv, abs/2311.00869.

**BibTeX entry:**
```
@misc{shebaro2023scaling,
      title={Scaling Frustration Index and Corresponding Balanced State Discovery for Real Signed Graphs}, 
      author={Muhieddine Shebaro and Jelena Tešić},
      year={2023},
      eprint={2311.00869},
      archivePrefix={arXiv},
      primaryClass={cs.SI}
}
```

## How to Run the Code 

* Simply download one of the .cpp files (each source file has a different tree-sampling technique) and compile the source file like the following:

```
user:~$ g++ -fopenmp graphBplus_BFS.cpp
```
The signed graph must be in the following format to be compatible with graphC (src,dst,sign).
Preprocessing of the signed graph is embedded.

* To execute the compiled file, GraphB++ utilizes 1 parameter (iteration_count):
```
user:~$ ./a.out input.txt 1000
```
 
* **Output:** graphBpp is going to execute and eventually create two files *_balancedStates.csv and *_frustrationindex where they will encompass each balanced state with its associated frequency and frustration index respectively.

**Note:** if you run into errors related to the stack memory please run this command before executing the code:
```
user:~$ ulimit -s unlimited
```

[Data Lab @ TXST](DataLab12.github.io)


