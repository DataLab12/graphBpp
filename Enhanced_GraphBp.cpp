/* GraphB++, 2024
 Author: Muhieddine Shebaro @ DataLab, Texas State University*/


#include <cstdio>
#include <climits>
#include <algorithm>
#include <set>
#include <map>
#include <sys/time.h>
#include <iostream>
#include <list>
#include <vector>
#include <set>
#include <queue>
#include <omp.h>
#include <random>
#include <chrono>

struct Graph {
  int nodes;
  int edges;
  int* nindex;  // first CSR array
  int* nlist;  // second CSR array
  int* eweight;  // edge weights (-1, 0, 1)
  bool* switched;
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
  delete [] g.switched;
  g.nindex = NULL;
  g.nlist = NULL;
  g.eweight = NULL;
  g.origID = NULL;
  g.switched = NULL;
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
  printf("  read %d lines\n", line);
  if (selfedges > 0) printf("  skipped %d self-edges\n", selfedges);
  if (wrongweights > 0) printf("  skipped %d edges with out-of-range weights\n", wrongweights);
  if (duplicates > 0) printf("  skipped %d duplicate edges\n", duplicates);
  if (inconsistent > 0) printf("  skipped %d inconsistent edges\n", inconsistent);
 

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

  // create graph in CSR format
  g.nodes = nodes;
  g.edges = edges;
  g.nindex = new int [g.nodes + 1];
  g.nlist = new int [g.edges];
  g.eweight = new int [g.edges];
  g.switched = new bool [g.edges];
  int acc = 0;
  for (int v = 0; v < g.nodes; v++) {
    g.nindex[v] = acc;
    for (auto ele: node[v]) {
      const int dst = ele.first;
      const int wei = ele.second;
      g.nlist[acc] = dst;
      g.switched[acc]=0;
      if(wei==0){
        g.eweight[acc] = 1;

      }else{
        g.eweight[acc] = wei;
      }
      acc++;
    }
  }
    g.nindex[g.nodes] = acc;

  delete [] label;
  delete [] size;
  delete [] node;

  return g;
}





int main(int argc, char* argv[])
{
  printf("graphB++ balancing code for signed social network graphs (%s)\n", __FILE__);
  printf("Copyright 2024 Texas State University\n");

  CPUTimer overall;
  overall.start();
  Graph g;
  // process command line and read input
  if (argc != 3) {printf("USAGE: %s input_file_name iteration_count\n", argv[0]); exit(-1);}
  CPUTimer timer;
  timer.start();
  printf("input: %s\n", argv[1]);
  g=readGraph(argv[1]);
  printf("input time: %.6f s\n", timer.elapsed());
  int iterations=atoi(argv[2]);
  std::vector<int> history;
  int frustration_max=INT_MAX;
int iter=0;

while(iter<iterations){
  int frustration=0.0;
    int* S=new int[g.nodes]; 
    int s=iter%g.nodes;
    std::vector<bool> visited;
    visited.resize(g.nodes,false);
 
    std::list<int> queue;
 
    visited[s] = true;
    queue.push_back(s);
    S[s]=1;
    while(!queue.empty())
    {
        s = queue.front();        
        queue.pop_front();
        int index=g.nindex[s];
        int index_max=g.nindex[s+1];
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(g.nlist + index, g.nlist + index_max,std::default_random_engine(seed));
        std::shuffle(g.eweight + index, g.eweight + index_max, std::default_random_engine(seed));

        while(index<index_max){
        {     
            if (!visited[g.nlist[index]])
            {   
                S[g.nlist[index]]=S[s]*g.eweight[index];
                visited[g.nlist[index]] = true;
                queue.push_back(g.nlist[index]);
            }else{

        if(S[s]*S[g.nlist[index]]*g.eweight[index]==-1){
          g.switched[index]=1;
          frustration++;
        }else{
          g.switched[index]=0;
        }


            }
            index=index+1;
        }
    }
    }
  frustration=frustration/2;
  history.push_back(frustration);
  if(frustration_max>frustration){
    frustration_max=frustration;
  }
   iter++;
  delete [] S;
std::cout<<"Iteration "<<iter<<" completed!"<<std::endl;
}
  std::cout<<frustration_max<<std::endl;

FILE *f = fopen("history.txt", "wt");

for(const int& i : history){
 fprintf(f, "%d\n", i);
} 
fclose(f);
  freeGraph(g);

  printf("overall runtime with I/O: %.6f s\n", overall.elapsed());
}
