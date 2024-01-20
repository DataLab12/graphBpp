/*
graphB++ balancing algorithm for signed social network graphs

Copyright 2023, Texas State University. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
   * Neither the name of Texas State University nor the names of its
     contributors may be used to endorse or promote products derived from
     this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL TEXAS STATE UNIVERSITY BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Authors: Ghadeer Alabandi, Martin Burtscher, and Muhieddine Shebaro.
*/


#include <unistd.h>
#include <ios>
#include <fstream>
#include <cstdio>
#include <climits>
#include <algorithm>
#include <set>
#include <map>
#include <sys/time.h>
#include <iostream>
#include <random>
#include <list>
#include <vector>
#include <set>
#include<queue>
#include <omp.h>
#include <unordered_set>
#include <iterator>
static const bool verify = false;  // set to false for better performance
struct EdgeInfo {
  int beg;  // beginning of range (shifted by 1) | is range inverted or not
  int end;  // end of range (shifted by 1) | plus or minus (1 = minus, 0 = plus or zero)
};

struct Graph {
  int nodes;
  int edges;
  int* nindex;  // first CSR array
  int* nlist;  // second CSR array
  int* eweight;  // edge weights (-1, 0, 1)
  int* origID;  // original node IDs
};

static void freeGraph(Graph &g)
{
  g.nodes = 0;
  g.edges = 0;
  delete [] g.nindex;
  delete [] g.nlist;
  delete [] g.eweight;
  delete [] g.origID;
  g.nindex = NULL;
  g.nlist = NULL;
  g.eweight = NULL;
  g.origID = NULL;
}

struct CPUTimer
{
  timeval beg, end;
  CPUTimer() {}
  ~CPUTimer() {}
  void start() {gettimeofday(&beg, NULL);}
  double elapsed() {gettimeofday(&end, NULL); return end.tv_sec - beg.tv_sec + (end.tv_usec - beg.tv_usec) / 1000000.0;}
};

static inline int representative(const int idx, int* const label)
{
  int curr = label[idx];
  if (curr != idx) {
    int next, prev = idx;
    while (curr > (next = label[curr])) {
      label[prev] = next;
      prev = curr;
      curr = next;
    }
  }
  return curr;
}

static Graph readGraph(const char* const name)
{
  // read input from file
  FILE* fin = fopen(name, "rt");
  if (fin == NULL) {printf("ERROR: could not open input file %s\n", name); exit(-1);}
  size_t linesize = 256;
  char buf[linesize];
  char* ptr = buf;
  getline(&ptr, &linesize, fin);  // skip first line

  int selfedges = 0, wrongweights = 0, duplicates = 0, inconsistent = 0, line = 1, cnt = 0;
  int src, dst, wei;
  std::map<int, int> map;  // map node IDs to contiguous IDs
  std::set<std::pair<int, int>> set2;
  std::set<std::tuple<int, int, int>> set3;
  while (fscanf(fin, "%d,%d,%d", &src, &dst, &wei) == 3) {
    if (src == dst) {
      selfedges++;
    } else if ((wei < -1) || (wei > 1)) {
      wrongweights++;
    } else if (set2.find(std::make_pair(std::min(src, dst), std::max(src, dst))) != set2.end()) {
      if (set3.find(std::make_tuple(std::min(src, dst), std::max(src, dst), wei)) != set3.end()) {
        duplicates++;
      } else {
        inconsistent++;
      }
    } else {
      set2.insert(std::make_pair(std::min(src, dst), std::max(src, dst)));
      set3.insert(std::make_tuple(std::min(src, dst), std::max(src, dst), wei));
      if (map.find(src) == map.end()) {
        map[src] = cnt++;
      }
      if (map.find(dst) == map.end()) {
        map[dst] = cnt++;
      }
    }
    line++;
  }
  fclose(fin);

  // print stats
  //printf("  read %d lines\n", line);
  //if (selfedges > 0) printf("  skipped %d self-edges\n", selfedges);
  //if (wrongweights > 0) printf("  skipped %d edges with out-of-range weights\n", wrongweights);
  //if (duplicates > 0) printf("  skipped %d duplicate edges\n", duplicates);
  //if (inconsistent > 0) printf("  skipped %d inconsistent edges\n", inconsistent);
  /*if (verify) {
    if ((int)map.size() != cnt) {printf("ERROR: wrong node count\n"); exit(-1);}
    printf("  number of unique nodes: %d\n", (int)map.size());
    printf("  number of unique edges: %d\n", (int)set3.size());
  }*/

  // compute CCs with union find
  int* const label = new int [cnt];
  for (int v = 0; v < cnt; v++) {
    label[v] = v;
  }
  for (auto ele: set3) {
    const int src = map[std::get<0>(ele)];
    const int dst = map[std::get<1>(ele)];
    const int vstat = representative(src, label);
    const int ostat = representative(dst, label);
    if (vstat != ostat) {
      if (vstat < ostat) {
        label[ostat] = vstat;
      } else {
        label[vstat] = ostat;
      }
    }
  }
  for (int v = 0; v < cnt; v++) {
    int next, vstat = label[v];
    while (vstat > (next = label[vstat])) {
      vstat = next;
    }
    label[v] = vstat;
  }

  // determine CC sizes
  int* const size = new int [cnt];
  for (int v = 0; v < cnt; v++) {
    size[v] = 0;
  }
  for (int v = 0; v < cnt; v++) {
    size[label[v]]++;
  }

  // find largest CC
  int hi = 0;
  for (int v = 1; v < cnt; v++) {
    if (size[hi] < size[v]) hi = v;
  }

  // keep if in largest CC and convert graph into set format
  Graph g;
  g.origID = new int [cnt];  // upper bound on size
  int nodes = 0, edges = 0;
  std::map<int, int> newmap;  // map node IDs to contiguous IDs
  std::set<std::pair<int, int>>* const node = new std::set<std::pair<int, int>> [cnt];  // upper bound on size
  for (auto ele: set3) {
    const int src = std::get<0>(ele);
    const int dst = std::get<1>(ele);
    const int wei = std::get<2>(ele);
    if (label[map[src]] == hi) {  // in largest CC
      if (newmap.find(src) == newmap.end()) {
        g.origID[nodes] = src;
        newmap[src] = nodes++;
      }
      if (newmap.find(dst) == newmap.end()) {
        g.origID[nodes] = dst;
        newmap[dst] = nodes++;
      }
      node[newmap[src]].insert(std::make_pair(newmap[dst], wei));
      node[newmap[dst]].insert(std::make_pair(newmap[src], wei));
      edges += 2;
    }
  }
  if (verify) {
    if (nodes > cnt) {printf("ERROR: too many nodes\n"); exit(-1);}
    if (edges > (int)set3.size() * 2) {printf("ERROR: too many edges\n"); exit(-1);}
  }

  // create graph in CSR format
  g.nodes = nodes;
  g.edges = edges;
  g.nindex = new int [g.nodes + 1];
  g.nlist = new int [g.edges];
  g.eweight = new int [g.edges];
  int acc = 0;
  for (int v = 0; v < g.nodes; v++) {
    g.nindex[v] = acc;
    for (auto ele: node[v]) {
      const int dst = ele.first;
      const int wei = ele.second;
      g.nlist[acc] = dst;
      g.eweight[acc] = wei;
      acc++;
    }
  }
  g.nindex[g.nodes] = acc;
  if (verify) {
    if (acc != edges) {printf("ERROR: wrong edge count in final graph\n"); exit(-1);}
  }

  delete [] label;
  delete [] size;
  delete [] node;

  return g;
}

// source of hash function: https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
static inline unsigned int hash(unsigned int val)
{
  val = ((val >> 16) ^ val) * 0x45d9f3b;
  val = ((val >> 16) ^ val) * 0x45d9f3b;
  return (val >> 16) ^ val;
}

static void init(const Graph& g, double* const inCC, EdgeInfo* const einfo, int* const inTree, int* const negCnt,int* end)
{
  // shift nlist
  #pragma omp parallel for default(none) shared(g)
  for (int j = 0; j < g.edges; j++) {
    g.nlist[j] <<= 1;
  }


  // set minus if graph weight is -1
  #pragma omp parallel for default(none) shared(g, einfo,end)
  for (int j = 0; j < g.edges; j++) {
    einfo[j].end = (g.eweight[j] == -1) ? 1 : 0;
    end[j]=(g.eweight[j] == -1) ? 1 : 0;
  }

  // zero out inTree and negCnt
  #pragma omp parallel for default(none) shared(g, inTree, negCnt)
  for (int j = 0; j < g.edges; j++) {
    inTree[j] = 0;
    negCnt[j] = 0;
  }
}


static double initMinus(const Graph& g, const EdgeInfo* const einfo, bool* const minus)
{
  CPUTimer timer;
  timer.start();

  // set minus info to true
  #pragma omp parallel for default(none) shared(g, minus)
  for (int j = 0; j < g.edges; j++) {
    minus[j] = true;
  }

  // copy minus info of tree edges
  #pragma omp parallel for default(none) shared(g, einfo, minus) schedule(dynamic, 64)
  for (int i = 0; i < g.nodes; i++) {
    int j = g.nindex[i];
    while ((j < g.nindex[i + 1]) && (g.nlist[j] & 1)) {
      minus[j] = einfo[j].end & 1;
      j++;
    }
  }

  return timer.elapsed();
}

static double processCycles(const Graph& g, const int* const label, EdgeInfo* const einfo, bool* const minus)
{
  CPUTimer timer;
  timer.start();
  #pragma omp parallel for default(none) shared(g, label, einfo, minus) schedule(dynamic, 64)
  for (int i = 0; i < g.nodes; i++) {
    const int target0 = label[i];
    const int target1 = target0 | 1;
    int j = g.nindex[i + 1] - 1;
    while ((j >= g.nindex[i]) && !(g.nlist[j] & 1)) {
      int curr = g.nlist[j] >> 1;
      if (curr > i) {  // only process edges in one direction
        int sum = 0;
        while (label[curr] != target0) {
          int k = g.nindex[curr];
          while ((einfo[k].beg & 1) == ((einfo[k].beg <= target1) && (target0 <= einfo[k].end))){ k++;   }

          sum += einfo[k].end & 1;
          curr = g.nlist[k] >> 1;
        }
        minus[j] = sum & 1;  // make cycle have even number of minuses
      }
      j--;
    }
  }

  return timer.elapsed();
}
void mem_usage(double& vm_usage, double& resident_set) { //https://www.tutorialspoint.com/how-to-get-memory-usage-at-runtime-using-cplusplus
   vm_usage = 0.0;
   resident_set = 0.0;
   std::ifstream stat_stream("/proc/self/stat",std::ios_base::in); //get info from proc directory
   //create some variables to get info
   std::string pid, comm, state, ppid, pgrp, session, tty_nr;
   std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   std::string utime, stime, cutime, cstime, priority, nice;
   std::string O, itrealvalue, starttime;
   unsigned long vsize;
   long rss;
   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
   >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
   >> utime >> stime >> cutime >> cstime >> priority >> nice
   >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
   stat_stream.close();
   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // for x86-64 is configured to use 2MB pages
   vm_usage = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}
unsigned long get_mem_total() { //https://stackoverflow.com/questions/349889/how-do-you-determine-the-amount-of-linux-system-ram-in-c
    std::string token;
    std::ifstream file("/proc/meminfo");
    while(file >> token) {
        if(token == "MemFree:") {
            unsigned long mem;
            if(file >> mem) {
                return mem;
            } else {
                return 0;
            }
        }
        // Ignore the rest of the line
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return 0; // Nothing found
}

std::map<std::string,int> span_dfs;

int markk=0;
int min_edge_switch=INT_MAX;
class Graphdfs
{
	int V; 

public:
	Graphdfs(int V);
std::map<int, bool> visiteddfs;
    std::map<int, std::list<int> > adjdfs;
      std::vector< std::pair<int, std::pair<int,int>> > edges;
            std::map<std::string, int > edges_prim;


	void addEdgedfs(int v, int w);
  void CC(std::vector<int>* cliques);
  void DFS_FOR_CONNECTED_COMPONENTS(int v, bool visited[],std::vector<int>* cliques);
      int kruskalMST();
      void DFSmark(int v,int prev,int* labell);
void addEdgekrus(int u, int v, int w);
void addEdgeprim(int u, int v, int w);

int minKey(int key[], bool mstSet[])
{
	// Initialize min value
	int min = INT_MAX, min_index;
	for (int v = 0; v < V; v++)
		if (mstSet[v] == false && key[v] < min)
			min = key[v], min_index = v;

	return min_index;
}

void printMST(int parent[])
{
	for (int i = 1; i < V; i++){
			std::string total_string="";
            total_string.append(std::to_string(parent[i]));
            total_string.append(",");
            total_string.append(std::to_string(i));
            std::string total_string1="";
            total_string1.append(std::to_string(i));
            total_string1.append(",");
            total_string1.append(std::to_string(parent[i]));
    span_dfs[total_string]=1;
    span_dfs[total_string1]=1;
  }
}


void primMST()
{
	// Array to store constructed MST
	int parent[V];

	// Key values used to pick minimum weight edge in cut
	int key[V];

	// To represent set of vertices included in MST
	bool mstSet[V];

	// Initialize all keys as INFINITE
	for (int i = 0; i < V; i++)
		key[i] = INT_MAX, mstSet[i] = false;

	// Always include first 1st vertex in MST.
	// Make key 0 so that this vertex is picked as first
	// vertex.
	key[0] = 0;
	// First node is always root of MST
	parent[0] = -1;
	// The MST will have V vertices
	for (int count = 0; count < V - 1; count++) {

		// Pick the minimum key vertex from the
		// set of vertices not yet included in MST
		int u = minKey(key, mstSet);

		// Add the picked vertex to the MST Set
		mstSet[u] = true;


		for (int v = 0; v < V; v++){
       
      std::string total_string="";
            total_string.append(std::to_string(u));
            total_string.append(",");
            total_string.append(std::to_string(v));
            bool flagg=false;
            if(edges_prim.find(total_string)!=edges_prim.end()){
flagg=true;
            }
			if (edges_prim[total_string] && mstSet[v] == false
				&& edges_prim[total_string] < key[v]){
				parent[v] = u, key[v] = edges_prim[total_string];
   
        }
        if(flagg==false){
          edges_prim.erase(total_string);
        }

    }
	}

	// Print the constructed MST
	printMST(parent);
}


// This code is contributed by rathbhupendra



void DFS(int v,int prev,int* labell);
};

Graphdfs::Graphdfs(int V)
{
	this->V = V;
}

void Graphdfs::addEdgedfs(int v, int w)
{
    	adjdfs[v].push_back(w); 

}
void Graphdfs::addEdgekrus(int u, int v, int w)
    {
        edges.push_back({w, {u, v}});
    }
    void Graphdfs::addEdgeprim(int u, int v, int w)
    {
        std::string total_string="";
            total_string.append(std::to_string(u));
            total_string.append(",");
            total_string.append(std::to_string(v));
            edges_prim[total_string]=w;
    }
  
  struct DisjointSets
{
    int *parent, *rnk;
    int n;
  
    // Constructor.
    DisjointSets(int n)
    {
        // Allocate memory
        this->n = n;
        parent = new int[n+1];
        rnk = new int[n+1];
  
        // Initially, all vertices are in
        // different sets and have rank 0.
        for (int i = 0; i <= n; i++)
        {
            rnk[i] = 0;
  
            //every element is parent of itself
            parent[i] = i;
        }
    }
  
    // Find the parent of a node 'u'
    // Path Compression
    int find(int u)
    {
        /* Make the parent of the nodes in the path
        from u--> parent[u] point to parent[u] */
        if (u != parent[u])
            parent[u] = find(parent[u]);
        return parent[u];
    }
  
    // Union by rank
    void merge(int x, int y)
    {
        x = find(x), y = find(y);
  
        /* Make tree with smaller height
        a subtree of the other tree */
        if (rnk[x] > rnk[y])
            parent[y] = x;
        else // If rnk[x] <= rnk[y]
            parent[x] = y;
  
        if (rnk[x] == rnk[y])
            rnk[y]++;
    }
};
  
/* Functions returns weight of the MST*/
  
int Graphdfs::kruskalMST()
{
    int mst_wt = 0; // Initialize result
  
    // Sort edges in increasing order on basis of cost
    std::sort(edges.begin(), edges.end());
  
    // Create disjoint sets
    DisjointSets ds(V);
  
    // Iterate through all sorted edges
    std::vector< std::pair<int, std::pair<int,int>> >::iterator it;
    for (it=edges.begin(); it!=edges.end(); it++)
    {
        int u = it->second.first;
        int v = it->second.second;
  
        int set_u = ds.find(u);
        int set_v = ds.find(v);
  
        // Check if the selected edge is creating
        // a cycle or not (Cycle is created if u
        // and v belong to same set)
        if (set_u != set_v)
        {
             std::pair<int,int> p(u,v);
            std::string total_string="";
            total_string.append(std::to_string(p.first));
            total_string.append(",");
            total_string.append(std::to_string(p.second));
            std::pair<int,int> p1(v,u);
            std::string total_string1="";
            total_string1.append(std::to_string(p1.first));
            total_string1.append(",");
            total_string1.append(std::to_string(p1.second));
    span_dfs[total_string]=1;
    span_dfs[total_string1]=1;
  
            // Update MST weight
            mst_wt += it->first;
  
            // Merge two sets
            ds.merge(set_u, set_v);
        }
    }
  
    return mst_wt;
}

int mark=0;
void Graphdfs::DFS(int v,int prev,int* labell)
{

    visiteddfs[v] = true;
    labell[v]=mark;

    std::list<int>::iterator i;
    for (i = adjdfs[v].begin(); i != adjdfs[v].end(); ++i)
        if (!visiteddfs[*i]){
            std::pair<int,int> p(*i,v);
            std::string total_string="";
            total_string.append(std::to_string(p.first));
            total_string.append(",");
            total_string.append(std::to_string(p.second));
            std::pair<int,int> p1(v,*i);
            std::string total_string1="";
            total_string1.append(std::to_string(p1.first));
            total_string1.append(",");
            total_string1.append(std::to_string(p1.second));
    span_dfs[total_string]=1;
    span_dfs[total_string1]=1;
        mark++;

            DFS(*i,v,labell);
        }
}
void Graphdfs::CC(std::vector<int>* cliques)
{
    bool* visited = new bool[V];
    for (int v = 0; v < V; v++)
        visited[v] = false;
 
    for (int v = 0; v < V; v++) {
        if (visited[v] == false) {
            DFS_FOR_CONNECTED_COMPONENTS(v, visited,cliques);
            markk++;
         }
    }
    delete[] visited;
}
 
void Graphdfs::DFS_FOR_CONNECTED_COMPONENTS(int v, bool visited[],std::vector<int>* cliques)
{
    visited[v] = true;
    cliques[markk].push_back(v);
 
    std::list<int>::iterator i;
    for (i = adjdfs[v].begin(); i != adjdfs[v].end(); ++i)
        if (!visited[*i])
            DFS_FOR_CONNECTED_COMPONENTS(*i, visited,cliques);
}

void Graphdfs::DFSmark(int v,int prev,int* labell)
{

    visiteddfs[v] = true;
    labell[v]=mark;
    std::list<int>::iterator i;
    for (i = adjdfs[v].begin(); i != adjdfs[v].end(); ++i){
        if (!visiteddfs[*i]){
        mark++;
            DFSmark(*i,v,labell);
        }
    }

  
}
std::vector<std::string> split_string(std::string str)
{
  std::vector<std::string> splits;
    std::string ss = "";
    for (auto x : str){
        if (x == ','){
            splits.push_back(ss);
            ss = "";
        }
        else {
            ss = ss + x;
        }
    }
            splits.push_back(ss);
            return splits;
}
static void determineCCs(const Graph& g, int* const label, const bool* const minus, int* const count, double* const inCC, int* const negCnt,std::vector<int> *edges_balanced_src,std::vector<int> *edges_balanced_dst,std::vector<bool>*weight_balanced,std::map<std::string,bool> copy)
{
  // init CCs
  #pragma omp parallel for default(none) shared(g, label)
  for (int v = 0; v < g.nodes; v++) {
    label[v] = v;
  }

  // compute CCs with union find
  int x=0;
  for (int v = 0; v < g.nodes; v++) {
    const int beg = g.nindex[v];
    const int end = g.nindex[v + 1];
    int vstat = representative(v, label);
    for (int j = beg; j < end; j++) {
      const int nli = g.nlist[j] >> 1;
      if(v<nli){
        edges_balanced_src->push_back(v);
        edges_balanced_dst->push_back(nli);
        weight_balanced->push_back(minus[j]);
        
        std::string total_string="";
        total_string.append(std::to_string(v));
        total_string.append(",");
        total_string.append(std::to_string(nli));
        //#pragma omp parallel for default(none) shared(resolution,weight_balanced_copy,edges_balanced_src_copy,edges_balanced_dst_copy,minus,nli,v,j,g)
        //for(int q=0;q<g.edges/2;q++){
    

      }
      if (minus[j]) {
        negCnt[j]++;
      } else {
        int ostat = representative(nli, label);
        bool repeat;
        do {
          repeat = false;
          if (vstat != ostat) {
            int ret;
            if (vstat < ostat) {
              if ((ret = __sync_val_compare_and_swap(&label[ostat], ostat, vstat)) != ostat) {
                ostat = ret;
                repeat = true;
              }
            } else {
              if ((ret = __sync_val_compare_and_swap(&label[vstat], vstat, ostat)) != vstat) {
                vstat = ret;
                repeat = true;
              }
            }
          }
        } while (repeat);
      }
    }
  }

  // finalize CCs
  #pragma omp parallel for default(none) shared(g, label)
  for (int v = 0; v < g.nodes; v++) {
    int next, vstat = label[v];
    const int old = vstat;
    while (vstat > (next = label[vstat])) {
      vstat = next;
    }
    if (old != vstat) label[v] = vstat;
  }

  // determine CC sizes
  #pragma omp parallel for default(none) shared(g, count)
  for (int v = 0; v < g.nodes; v++) {
    count[v] = 0;
  }
  #pragma omp parallel for default(none) shared(g, count, label)
  for (int v = 0; v < g.nodes; v++) {
    #pragma omp atomic
    count[label[v]]++;
  }
  bool equal_partition=false;

  // find largest CC (source CC)
  int hi = 0;
  #pragma omp parallel for default(none) shared(g, hi, count)
  for (int v = 1; v < g.nodes; v++) {
    if (count[hi] < count[v]) {
      #pragma omp critical
      if (count[hi] < count[v]) hi = v;
    }
  }


  // init CC hop count (distance) from source CC, populate workset of edges that cross CCs
  std::set< std::pair<int, int> > ws;
  for (int v = 0; v < g.nodes; v++) {
    const int lblv = label[v];
    if (lblv == v) {
      count[lblv] = (lblv == hi) ? 0 : INT_MAX - 1;  // init count
    }
    for (int j = g.nindex[v]; j < g.nindex[v + 1]; j++) {
      const int nli = g.nlist[j] >> 1;
      const int lbln = label[nli];
      if (lblv < lbln) {  // only one direction
        ws.insert(std::make_pair(lblv, lbln));
      }
    }
  }

  // use Bellman Ford to compute distances
  bool changed;
  do {
    changed = false;
    for (auto p: ws) {
      const int lblv = p.first;
      const int lbln = p.second;
      const int distv = count[lblv];
      const int distn = count[lbln];
      if (distv + 1 < distn) {
        count[lbln] = distv + 1;
        changed = true;
      } else if (distn + 1 < distv) {
        count[lblv] = distn + 1;
        changed = true;
      }
    }
  } while (changed);

  // increment inCC if node is at even hop count from source CC
  #pragma omp parallel for default(none) shared(g, hi, label, count, inCC)
  for (int v = 0; v < g.nodes; v++) {
    inCC[v] += (count[label[v]] % 2) ^ 1;
  }

}



int main(int argc, char* argv[])
{
  printf("graphB++ balancing code for signed social network graphs (%s) -- All in One\n", __FILE__);
  printf("Copyright 2022 Texas State University\n");

  CPUTimer overall;
  overall.start();
  Graph g;
  Graph gspare;
  // process command line and read input
  if (argc != 3) {printf("USAGE: %s input_file_name iteration_count\n", argv[0]); exit(-1);}
  CPUTimer timer;
  timer.start();
  printf("verification: %s\n", verify ? "on" : "off");
  printf("input: %s\n", argv[1]);


  const int iterations = atoi(argv[2]);

  int* const switchesI = new int [iterations];
  std::vector<int> edges_balanced_src;
  std::vector<int> edges_balanced_dst;
  std::vector<bool> weight_balanced;

  timer.start();
 
  //printf("init time:  %.6f s\n", timer.elapsed());
  g = readGraph(argv[1]);
  gspare=readGraph(argv[1]);
  double graphBtime = 0;
  //int* const agreement = new int [g.edges/2];
  //int* const resolution = new int [g.edges/2];
   double* const inCC = new double [g.nodes];  // how often node was in largest CC or at an even distance from largest CC
   #pragma omp parallel for default(none) shared(g, inCC)
  for (int v = 0; v < g.nodes; v++) {
    inCC[v] = 0;
  }
   std::map<std::string,bool> copy;
 std::map<std::string, std::pair<int,int>> map;

  std::map<std::string,bool> mappy;
   for (int v = 0; v < g.nodes; v++) {
    int index=g.nindex[v];
    int index_max=g.nindex[v+1];
    while(index<index_max){
      if(v<g.nlist[index]){
     std::string total_string="";
          total_string.append(std::to_string(v));
          total_string.append(",");
          total_string.append(std::to_string(g.nlist[index]));
      if(g.eweight[index]==-1){
       mappy[total_string]=1;
          copy[total_string]=1;
      }else{
mappy[total_string]=0;
          copy[total_string]=0;
                }
      }
    index=index+1;
    }
  }


  for (int iter = 0; iter < iterations; iter++) {
    std::cout<<"Iteration "<<iter+1<<" Started!"<<std::endl;
    mark=0;

  g.edges=gspare.edges;
  g.nodes=gspare.nodes;
  #pragma parallel for default(none) shared(g,gspare)
 for(int i=0;i<g.edges;i++){
   g.nlist[i]=gspare.nlist[i];
   g.eweight[i]=gspare.eweight[i];
 }
  #pragma parallel for default(none) shared(g,gspare)
  for(int i=0;i<g.nodes;i++){
   g.origID[i]=gspare.origID[i];
 }


        bool* const minus = new bool [g.edges];
  int* const parent = new int [g.nodes];
  int* const queue = new int [g.nodes];  // first used as queue, then as CC size
  int*  label = new int [g.nodes];  // first used as count, then as label, and finally as CC label
    int*  labell = new int [g.nodes];  // first used as count, then as label, and finally as CC label

  int* const border = new int [g.nodes + 2];  // maybe make smaller
  int* const inTree = new int [g.edges];  // how often edge was in tree
  int* const negCnt = new int [g.edges];  // how often edge was negative
    int* const beg = new int [g.edges];
        int* const end = new int [g.edges];
        int* const copyedges = new int [g.edges];
                int* const finaledges = new int [g.edges];

  #pragma omp parallel for default(none) shared(g,copyedges)
  for(int i=0;i<g.edges;i++){
      copyedges[i]=g.nlist[i]<<1;
    }

  EdgeInfo* const einfo = new EdgeInfo [g.edges];
  int* const root = new int [g.nodes];  // tree roots
  #pragma omp parallel for default(none) shared(g,root)
 for (int i = 0; i < g.nodes; i++){
    root[i] = i;
 }

  std::partial_sort(root, root + std::min(iterations, g.nodes), root + g.nodes, [&](int a, int b) {
    return (g.nindex[a + 1] - g.nindex[a]) > (g.nindex[b + 1] - g.nindex[b]);
  });
  init(g, inCC, einfo, inTree, negCnt,end);
        edges_balanced_src.clear();
    weight_balanced.clear();
    edges_balanced_dst.clear();

    // generate tree
    timer.start();
  std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(1, g.edges);

        Graphdfs gdfs(g.nodes);
 
     for (int i = 0; i < g.nodes; i++) {
            int current = g.nindex[i];
            int stop = g.nindex[i + 1];
            for (int j = current; j < stop; j++) {
	                gdfs.addEdgeprim(i, g.nlist[j]>>1,distrib(gen));	 
            }

        }

        gdfs.primMST();
                std::map<std::string,int> spanningtree_edges=span_dfs;


          span_dfs.clear();
          std::cout<<"ggg"<<std::endl;
           Graphdfs glabel(g.nodes);
        for (auto const& x :spanningtree_edges)
{
  std::vector<std::string> splits=split_string(x.first);


    	                glabel.addEdgedfs(std::stoi(splits[0]), std::stoi(splits[1]));	 

}
glabel.DFSmark(iter%g.nodes,-1,labell);
    #pragma omp parallel for default(none) shared(label,labell,g)
      for(int i=0;i<g.nodes;i++){
      label[i]=labell[i]<<1;
    }
  for (int j = 0; j < g.edges; j++) copyedges[j] &= ~1;  // edge is not in tree
        #pragma omp parallel for default(none) shared(beg,spanningtree_edges,copyedges,finaledges,g,label)
        for (int i = 0; i < g.nodes; i++) {
            int current = g.nindex[i];
            int stop = g.nindex[i + 1];
            for (int j = current; j < stop; j++) {
              bool found=false;
                    std::string total_string="";
            total_string.append(std::to_string(i));
            total_string.append(",");
            total_string.append(std::to_string(copyedges[j]>>1));
                        if(spanningtree_edges.find(total_string)!=spanningtree_edges.end()){
                            finaledges[j]=copyedges[j]|1;

                            found=true;
                            if(label[i]<label[copyedges[j]>>1]){
                                                            beg[j]=label[copyedges[j]>>1];

                            }else{
                                                          beg[j]=label[i]|1;

                            }

                        }

              if(found==false){
                    finaledges[j]=copyedges[j];
                    beg[j]=beg[j]&0;
            }
            }

        }
int* represent=new int[g.nodes];
#pragma omp parallel for default(none) shared(represent,g,spanningtree_edges,label,copyedges)
for (int i = 0; i < g.nodes; i++) {
   int current = g.nindex[i];
            int stop = g.nindex[i + 1];
            int max_label=label[i];
            int thisnode=i;
            for (int j = current; j < stop; j++) {
               bool found=false;
            std::string total_string="";
            total_string.append(std::to_string(i));
            total_string.append(",");
            total_string.append(std::to_string(copyedges[j]>>1));
                        if(spanningtree_edges.find(total_string)!=spanningtree_edges.end()){
                            found=true;
                        }
            
              if(found==false){
                continue;
              }
              int current_label=label[copyedges[j]>>1];
              if(current_label>max_label){
                 thisnode=copyedges[j]>>1;
                 max_label=current_label;
              }
            }
            represent[i]=thisnode;
}

#pragma omp parallel for default(none) shared(end,g,spanningtree_edges,label,copyedges,represent)
for (int i = 0; i < g.nodes; i++) {
            int current = g.nindex[i];
            int stop = g.nindex[i + 1];
            for (int j = current; j < stop; j++) {
              bool found=false;
                      std::string total_string="";
            total_string.append(std::to_string(i));
            total_string.append(",");
            total_string.append(std::to_string(copyedges[j]>>1));
                        if(spanningtree_edges.find(total_string)!=spanningtree_edges.end()){
                            found=true;
                        }
            
              if(found==false){
                end[j]=(end[j] & 1);
                continue;
              }
              int source=label[i];
              int dest=label[copyedges[j]>>1];
              int label_end=-1;
              if(source<dest){
                int q=copyedges[j]>>1;
                int next=represent[q];
                while(label[next]>label[q]){
                    q=next;
                    next=represent[q];
                }
              label_end=label[q];
              end[j]=(end[j] & 1) | label_end;
              continue;

              }
                else{
                int q=i;
                int next=represent[q];
                while(label[next]>label[q]){
                    q=next;
                    next=represent[q];
                }
              label_end=label[q];
                    

              end[j]=(end[j] & 1) | label_end;
              continue;

              }


              
            }
}
//swap
int* parentt=new int[g.edges];
#pragma omp parallel for default(none) shared(parentt,label,copyedges,g)
for (int i = 0; i < g.nodes; i++) {
            int current = g.nindex[i];
            int stop = g.nindex[i + 1];
            for (int j = current; j < stop; j++) {
                if(label[i]>label[copyedges[j]>>1]){
                  parentt[j]=1;
                }
            }
}
#pragma omp parallel for default(none) shared(beg,finaledges,copyedges,end,inTree,negCnt,parentt,g,spanningtree_edges)
for (int i = 0; i < g.nodes; i++) {
            int current = g.nindex[i];
            int stop = g.nindex[i + 1];
            int temp1=g.nindex[i];
            for (int j = current; j < stop; j++) {
              bool found=false;
                                 std::string total_string="";
            total_string.append(std::to_string(i));
            total_string.append(",");
            total_string.append(std::to_string(copyedges[j]>>1));
                        if(spanningtree_edges.find(total_string)!=spanningtree_edges.end()){     
                                                        found=true;
                                                      }
                            if(found==false){
                              if(j==(stop-1)){
                                continue;
                              }
                                                                                                           std::swap(copyedges[j], copyedges[stop-1]);
                                          std::swap(finaledges[j], finaledges[stop-1]);
                                                                                    std::swap(beg[j], beg[stop-1]);
                                                                                   std::swap(end[j], end[stop-1]);
                                                                                   std::swap(inTree[j], inTree[stop-1]);
                                                                                   std::swap(negCnt[j], negCnt[stop-1]);
                                                           std::swap(parentt[j], parentt[stop-1]);
                                                stop=stop-1;
                                                j=g.nindex[i]-1;
                                          continue;
                            }else{
                                 if(parentt[j]==1){
                                   if(j<=temp1){
                                     continue;
                                   }
                                        std::swap(finaledges[j], finaledges[temp1]);
                                                                                std::swap(copyedges[j], copyedges[temp1]);

                                                                                    std::swap(beg[j], beg[temp1]);
                                                                                   std::swap(end[j], end[temp1]);
                                                                                   std::swap(negCnt[j], negCnt[temp1]);
                                                                                   std::swap(inTree[j], inTree[temp1]);
                                                                                                                                              std::swap(parentt[j], parentt[temp1]);    
                                                                                   temp1++;
                                                                                                                 }
                            }
               
            }
}




#pragma omp parallel for default(none) shared(beg,end,einfo,g,finaledges,inTree)
for(int i=0;i<g.edges;i++){
  einfo[i].beg=beg[i]; 
  einfo[i].end=end[i];
  g.nlist[i]=finaledges[i];
  inTree[i]=(g.nlist[i]&1);
}

    // initialize plus/minus
    timer.start();
    graphBtime += initMinus(g, einfo, minus);
    // find cycles
    timer.start();
    graphBtime += processCycles(g, label, einfo, minus);
    // determine connected components
    timer.start();

    determineCCs(g, label, minus, queue, inCC, negCnt,&edges_balanced_src,&edges_balanced_dst,&weight_balanced,copy);
            std::cout<<"Iteration "<<iter+1<<" Completed!"<<std::endl;
              
  //#pragma omp parallel for default(none) shared(g,mapped,edges_balanced_src,edges_balanced_dst,weight_balanced)
              int switches=0;

           std::set<std::tuple<int,int,int>> tempset;
            for(int temp=0;temp<g.edges/2;temp++){
      int src=edges_balanced_src[temp];
      int dst=edges_balanced_dst[temp];
      bool weight=weight_balanced[temp];
      std::string total_string="";
          total_string.append(std::to_string(src));
          total_string.append(",");
          total_string.append(std::to_string(dst));
          if(weight!=mappy[total_string]){
            //std::cout<<total_string.c_str()<<std::endl;
            //std::cout<<mapped[total_string]<<"x"<<mappy[total_string]<<std::endl;
            switches++;
          }
                      std::tuple<int,int,int>p(g.origID[src],g.origID[dst],weight);
                      tempset.insert(p);
  }
  switchesI[iter]=switches;

            if(switches<min_edge_switch){
              min_edge_switch=switches;
            }
       std::string total_string = "";
  std::set<std::tuple<int,int,int>>::iterator itr;
  // Displaying set elements
  for (itr = tempset.begin();
       itr != tempset.end(); itr++)
  {
    total_string.append(std::to_string(std::get<0>(*itr)));
    total_string.append("-->");
        total_string.append(std::to_string(std::get<1>(*itr)));
        total_string.append(":");
            total_string.append(std::to_string(std::get<2>(*itr)));
            total_string.append(" | ");

    }
    if (map.find(total_string) == map.end()) {
      std::pair<int,int> p(switchesI[iter],1);
        map[total_string] = p;
} else {
      std::pair<int,int> p(switchesI[iter],map[total_string].second+1);
        map[total_string] = p;
}


  spanningtree_edges.clear();
  delete [] minus;
  delete [] einfo;
  delete [] parent;
  delete [] queue;
  delete [] label;
  delete [] border;
  delete [] inTree;
  delete [] negCnt;
  delete [] root;
  delete [] copyedges;
  delete [] finaledges;
  delete [] beg;
  delete [] end;
  delete [] labell;
  delete [] represent;
  }
  printf("graph bal time: %.6f s\n", graphBtime);
   std::string namei(argv[1]);
    namei.pop_back();
    namei.pop_back();
    namei.pop_back();
    namei.pop_back();
    namei.append("_FrustrationIndex.txt");
  FILE *findex = fopen(namei.c_str(), "wt");
    fprintf(findex,"-------------------------GraphB++-------------------------\n");
  fprintf(findex,"Frustration Index (%d iterations): %d",iterations, min_edge_switch);
  
      int counterI=0;
        

  freeGraph(g);
  freeGraph(gspare);
   g=readGraph(argv[1]);
   
 

  for (auto const& x : map)
{
    if(x.second.first==min_edge_switch){
      counterI++;
    }
}

    fprintf(findex,"\n# of Unique Balanced States With Minimum Number of Switches: %d", counterI);
    fprintf(findex,"\n");
    fprintf(findex,"\nDetailed Analysis:");

    for (int iter = 0; iter < iterations; iter++) {
                    fprintf(findex,"\n# of Switches of iteration %d: %d",iter+1, switchesI[iter]);
    }
      fprintf(findex,"\nMinimum # of Switches after %d iteration(s): %d",iterations, min_edge_switch);
  fclose(findex);
  std::cout << "Saved to \""<<namei<<"\"."<< std::endl;

    std::string name(argv[1]);
    name.pop_back();
    name.pop_back();
    name.pop_back();
    name.pop_back();
    name.append("_balancedStates.csv");

  FILE *balanced = fopen(name.c_str(), "wt");
                  fprintf(balanced, "Balanced State ID,Frequency,Frustration Index\n");
      int tempC=1;
     for (auto const& x : map) {
      fprintf(balanced, "%d,%d",tempC,x.second.second);
            fprintf(balanced, ",%d\n",x.second.first);
            tempC++;
    }
      fclose(balanced);
        std::cout << "Saved to \""<<name<<"\"."<< std::endl;


  delete [] inCC;

     freeGraph(g);

  printf("overall runtime with I/O: %.6f s\n", overall.elapsed());
}
