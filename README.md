# Scaling Frustration Index and Balanced State Discovery in Large Signed Networks

Authors: <em> Muhieddine Shebaro, Jelena Tešić </em>

![Example!](balance.gif "Example")

**GraphB++** (2023) extends GraphB+’s code to scale the approximation of the frustration index (NP-hard) for any real world signed network of any size or density based on efficiently finding a fundamental cycle basis. It also extends the frustration cloud to a (key,value) tuple collection F = B:(C, S). The states with their associated frequency and edge switches will be stored as a tuple in the memory-bound frustrated cloud in a scalable manner. We provide benchmark comparison of seven spanning tree sampling methods, including a breakthrough study of frustration index computation and timing for networks.

**GraphL** (2024) is a scalable, powerful, and linear algorithm that discovers a “nearer” nearest balanced state of a signed network using gradient descent and fully signed networks.


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

* Simply download one of the algorithms (all the source code files correspond to GraphB++ with different tree-sampling techniques used except GraphL.cpp) and compile the former source file like the following:

```
user:~$ g++ -fopenmp graphBplus_BFS.cpp
```
The signed graph must be in the following format (src,dst,sign).
Preprocessing of the signed graph is embedded and neutral edges are treated as positive.

* To execute the compiled file if it's GraphB++ (any of the .cpp code files except GraphL), GraphB++ utilizes 2 parameters in this order (input_file_name iteration_count):
```
user:~$ ./a.out input.txt 5000
```
* To execute the compiled file if it's GraphBL, GraphBL utilizes 2 parameters in this order (input_file_name iteration_count):
```
user:~$ ./a.out input.txt 5000
```
By default, learning rate is 0.001, you can change that in the cpp code.

* GraphB++ simultaneously balances the signed graph and computes the frustration index. Eventually, it creates two files *_balancedStates.csv and *_frustrationindex where they will encompass each balanced state with its associated frequency and frustration index respectively. Whereas GraphL will print the frustration index in the terminal at the end of the execution.
* Enhanced_GraphBp.cpp file is another improved version of graphB+ that is linear and utilizes BFS sampling to compute the frustration index.
* GraphBplus_Consensus.cpp computes the consensus features of input signed graph.

**Note:** If you run into errors related to the stack memory please run this command before executing the code:
```
user:~$ ulimit -s unlimited
```

[Data Lab @ TXST](DataLab12.github.io)


