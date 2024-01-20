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

Authors: Ghadeer Alabandi and Martin Burtscher and Muhieddine Shebaro
*/

#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

#include <cstdio>
#include <climits>
#include <algorithm>
#include <set>
#include <map>
#include <vector>
#include <map>
#include <sys/time.h>
#include <iostream>

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
  printf("  read %d lines\n", line);
  if (selfedges > 0) printf("  skipped %d self-edges\n", selfedges);
  if (wrongweights > 0) printf("  skipped %d edges with out-of-range weights\n", wrongweights);
  if (duplicates > 0) printf("  skipped %d duplicate edges\n", duplicates);
  if (inconsistent > 0) printf("  skipped %d inconsistent edges\n", inconsistent);
  if (verify) {
    if ((int)map.size() != cnt) {printf("ERROR: wrong node count\n"); exit(-1);}
    printf("  number of unique nodes: %d\n", (int)map.size());
    printf("  number of unique edges: %d\n", (int)set3.size());
  }

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

static void init(const Graph& g, int* const inCC, EdgeInfo* const einfo, int* const inTree, int* const negCnt)
{
  // shift nlist
  #pragma omp parallel for default(none) shared(g)
  for (int j = 0; j < g.edges; j++) {
    g.nlist[j] <<= 1;
  }

  // zero out inCC
  #pragma omp parallel for default(none) shared(g, inCC)
  for (int v = 0; v < g.nodes; v++) {
    inCC[v] = 0;
  }

  // set minus if graph weight is -1
  #pragma omp parallel for default(none) shared(g, einfo)
  for (int j = 0; j < g.edges; j++) {
    einfo[j].end = (g.eweight[j] == -1) ? 1 : 0;
  }

  // zero out inTree and negCnt
  #pragma omp parallel for default(none) shared(g, inTree, negCnt)
  for (int j = 0; j < g.edges; j++) {
    inTree[j] = 0;
    negCnt[j] = 0;
  }
}

static double generateSpanningTree(const Graph& g, const int root, const int seed, EdgeInfo* const einfo, int* const parent, int* const queue, int* const border, int* const label, int* const inTree, int* const negCnt)
{
  const int seed2 = seed * seed + seed;

  // initialize
  #pragma omp parallel for default(none) shared(g)
  for (int j = 0; j < g.edges; j++) g.nlist[j] &= ~1;  // edge is not in tree
  #pragma omp parallel for default(none) shared(g, parent)
  for (int i = 0; i < g.nodes; i++) parent[i] = -1;
  int tail = 1;
  parent[root] = INT_MAX & ~3;
  queue[0] = root;

  // BFS traversal
  int level = 0;
  border[0] = 0;
  border[1] = tail;
  while (border[level + 1] < g.nodes) {  // skipping last iteration
    const int bit = (level & 1) | 2;
    #pragma omp parallel for default(none) shared(g, level, tail, queue, parent, seed2, border, bit)
    for (int i = border[level]; i < border[level + 1]; i++) {
      const int node = queue[i];
      const int me = (node << 2) | bit;
      #pragma omp atomic write
      parent[node] = parent[node] & ~3;
      for (int j = g.nindex[node]; j < g.nindex[node + 1]; j++) {
        const int neighbor = g.nlist[j] >> 1;
        const int seed3 = neighbor ^ seed2;
        const int hash_me = hash(me ^ seed3);
        int val, hash_val;
        do {  // pick parent deterministically
          #pragma omp atomic read
          val = parent[neighbor];
          hash_val = hash(val ^ seed3);
        } while (((val < 0) || (((val & 3) == bit) && ((hash_val < hash_me) || ((hash_val == hash_me) && (val < me))))) && (__sync_val_compare_and_swap(&parent[neighbor], val, me) != val));
        if (val < 0) {
          #pragma omp atomic capture
          val = tail++;
          queue[val] = neighbor;
        }
      }
    }
    level++;
    if (border[level] == tail) {printf("ERROR: input appears to have multiple connected components; terminating program\n"); exit(-1);}
    border[level + 1] = tail;
  }
  const int levels = level + 1;
  if (verify) {
    if (border[levels] != tail) {printf("ERROR: head mismatch\n"); exit(-1);}
    if (tail != g.nodes) {printf("ERROR: tail mismatch\n"); exit(-1);}
    for (int i = 0; i < g.nodes; i++) {
      if (parent[i] < 0) {printf("ERROR: found unvisited node %d\n", i); exit(-1);}
    }
  }

  CPUTimer timer;
  timer.start();

  // bottom up: push counts
  #pragma omp parallel for default(none) shared(g, label, border)
  for (int i = 0; i < g.nodes; i++) label[i] = 1;
  for (int level = levels - 1; level > 0; level--) {  // skip level 0
    #pragma omp parallel for default(none) shared(level, queue, label, parent, border)
    for (int i = border[level]; i < border[level + 1]; i++) {
      const int node = queue[i];
      #pragma omp atomic
      label[parent[node] >> 2] += label[node];
    }
  }
  if (verify) {
    if (label[root] != g.nodes) {printf("ERROR: root count mismatch\n"); exit(-1);}
  }
  // top down: label tree + set nlist flag + set edge info + move tree nodes to front + make parent edge first in list
  label[root] = 0;
  for (int level = 0; level < levels; level++) {
    #pragma omp parallel for default(none) shared(g, level, queue, parent, label, einfo, inTree, negCnt, border)
    for (int i = border[level]; i < border[level + 1]; i++) {
      const int node = queue[i];
      const int par = parent[node] >> 2;
      const int nodelabel = label[node];
      const int beg = g.nindex[node];
      int paredge = -1;
      int lbl = (nodelabel >> 1) + 1;
      int pos = beg;
      for (int j = beg; j < g.nindex[node + 1]; j++) {
        const int neighbor = g.nlist[j] >> 1;
        if (neighbor == par) {
          paredge = j;
        } else if ((parent[neighbor] >> 2) == node) {
          const int count = label[neighbor];
          label[neighbor] = lbl << 1;
          lbl += count;
          // set child edge info
          einfo[j].beg = label[neighbor];
          einfo[j].end = (einfo[j].end & 1) | ((lbl - 1) << 1);
          g.nlist[j] |= 1;  // child edge is in tree
          // swap
          if (pos < j) {
            std::swap(g.nlist[pos], g.nlist[j]);
            std::swap(einfo[pos], einfo[j]);
            std::swap(inTree[pos], inTree[j]);
            std::swap(negCnt[pos], negCnt[j]);
            if (paredge == pos) paredge = j;
          }
          pos++;
        }
      }
      if (paredge >= 0) {
        // set parent edge info
        einfo[paredge].beg = nodelabel | 1;
        einfo[paredge].end = (einfo[paredge].end & 1) | ((lbl - 1) << 1);
        g.nlist[paredge] |= 1;  // parent edge is in tree
        // move parent edge to front of list
        if (paredge != beg) {
          if (paredge != pos) {
            std::swap(g.nlist[pos], g.nlist[paredge]);
            std::swap(einfo[pos], einfo[paredge]);
            std::swap(inTree[pos], inTree[paredge]);
            std::swap(negCnt[pos], negCnt[paredge]);
            paredge = pos;
          }
          if (paredge != beg) {
            std::swap(g.nlist[beg], g.nlist[paredge]);
            std::swap(einfo[beg], einfo[paredge]);
            std::swap(inTree[beg], inTree[paredge]);
            std::swap(negCnt[beg], negCnt[paredge]);
          }
        }
      }

      if (verify) {
        if (i == 0) {
          if (lbl != g.nodes) {printf("ERROR: lbl mismatch\n"); exit(-1);}
        }
        int j = beg;
        while ((j < g.nindex[node + 1]) && (g.nlist[j] & 1)) j++;
        while ((j < g.nindex[node + 1]) && !(g.nlist[j] & 1)) j++;
        if (j != g.nindex[node + 1]) {printf("ERROR: not moved %d %d %d\n", beg, j, g.nindex[node + 1]); exit(-1);}
      }
    }
  }

  const double rt = timer.elapsed();

  // update inTree
  #pragma omp parallel for default(none) shared(g, inTree)
  for (int j = 0; j < g.edges; j++) {
    inTree[j] += g.nlist[j] & 1;
  }

  return rt;
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
          while ((einfo[k].beg & 1) == ((einfo[k].beg <= target1) && (target0 <= einfo[k].end))){ k++;};
          if (verify) {
            if ((k >= g.nindex[curr + 1]) || !(g.nlist[k] & 1)) {printf("ERROR: couldn't find path\n"); exit(-1);}
          }
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

static void determineCCs(const Graph& g, int* const label, const bool* const minus, int* const count, int* const inCC, int* const negCnt,std::vector<int> *edges_balanced_src,std::vector<int> *edges_balanced_dst,std::vector<bool>*weight_balanced)
{
  // init CCs
  #pragma omp parallel for default(none) shared(g, label)
  for (int v = 0; v < g.nodes; v++) {
    label[v] = v;
  }

  // compute CCs with union find
  //#pragma omp parallel for default(none) shared(g, label, minus, negCnt,edges_balanced_src,edges_balanced_dst,weight_balanced) schedule(dynamic, 64)
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


int min_edge_switch=INT_MAX;

int main(int argc, char* argv[])
{
  printf("graphB++ balancing code for signed social network graphs (%s)\n", __FILE__);
  printf("Copyright 2021 Texas State University\n");

  CPUTimer overall;
  overall.start();

  // process command line and read input
  if (argc != 3) {printf("USAGE: %s input_file_name iteration_count\n", argv[0]); exit(-1);}
  CPUTimer timer;
  timer.start();
  printf("verification: %s\n", verify ? "on" : "off");
  printf("input: %s\n", argv[1]);
  Graph g = readGraph(argv[1]);
  printf("nodes: %d\n", g.nodes);
  printf("edges: %d\n", g.edges);
  const int iterations = atoi(argv[2]);
  printf("input time: %.6f s\n", timer.elapsed());
  std::vector<int> edges_balanced_src;
  std::vector<int> edges_balanced_dst;
  std::vector<bool> weight_balanced;

      std::map<std::string, bool> mappy;

   std::map<std::string, std::pair<int,int>> map;

  // allocate all memory
  bool* const minus = new bool [g.edges];
  int* const parent = new int [g.nodes];
  int* const queue = new int [g.nodes];  // first used as queue, then as CC size
  int* const label = new int [g.nodes];  // first used as count, then as label, and finally as CC label
  int* const border = new int [g.nodes + 2];  // maybe make smaller
  int* const inCC = new int [g.nodes];  // how often node was in largest CC or at an even distance from largest CC
  int* const inTree = new int [g.edges];  // how often edge was in tree
  int* const negCnt = new int [g.edges];  // how often edge was negative
  EdgeInfo* const einfo = new EdgeInfo [g.edges];
  int* const root = new int [g.nodes];  // tree roots

  timer.start();
   for (int v = 0; v < g.nodes; v++) {
    int index=g.nindex[v];
    int index_max=g.nindex[v+1];
    while(index<index_max){
      if(v<g.nlist[index]){
          std::string total_string = "";
          total_string.append(std::to_string(v));
          total_string.append(",");
          total_string.append(std::to_string(g.nlist[index]));
      if(g.eweight[index]==-1){
          mappy[total_string] = 1;
      }else{
          mappy[total_string] = 0;
      }
      }
    index=index+1;
    }
  }
  for (int i = 0; i < g.nodes; i++) root[i] = i;
  std::partial_sort(root, root + std::min(iterations, g.nodes), root + g.nodes, [&](int a, int b) {
    return (g.nindex[a + 1] - g.nindex[a]) > (g.nindex[b + 1] - g.nindex[b]);
  });
  init(g, inCC, einfo, inTree, negCnt);
  //printf("init time:  %.6f s\n", timer.elapsed());
double avgg=0;
  double graphBtime = 0;
  for (int iter = 0; iter < iterations; iter++) {
    edges_balanced_src.clear();
    weight_balanced.clear();
    edges_balanced_dst.clear();
    //printf("tree %d\n", iter);

    // generate tree
    timer.start();
    graphBtime += generateSpanningTree(g, root[iter % g.nodes], iter + 17, einfo, parent, queue, border, label, inTree, negCnt);
    //if (iter == 0) printf("tree time:  %.6f s\n", timer.elapsed());

    // initialize plus/minus
    timer.start();
    graphBtime += initMinus(g, einfo, minus);
    //if (iter == 0) printf("pm time:    %.6f s\n", timer.elapsed());

    // find cycles
    timer.start();
            int switches=0;
    double temp= processCycles(g, label, einfo, minus);
    //if (iter == 0) printf("cycle time: %.6f s\n", timer.elapsed());


    // determine connected components
    timer.start();

    determineCCs(g, label, minus, queue, inCC, negCnt,&edges_balanced_src,&edges_balanced_dst,&weight_balanced);
            std::cout<<"Iteration "<<iter+1<<" Completed!"<<std::endl;
              
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
            switches++;
          }
            std::tuple<int,int,int>p(g.origID[src],g.origID[dst],weight);
              tempset.insert(p);
  }
            if(switches<min_edge_switch){
              min_edge_switch=switches;
            }
      std::string total_string = "";
  std::set<std::tuple<int,int,int>>::iterator itr;
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
   double vm, rss;
   mem_usage(vm, rss);
   if((get_mem_total()-get_mem_total()*0.25) <= vm){
     std::map<std::string, std::pair<int,int>>::iterator largest = std::max_element(map.begin(),map.end(),[] (const std::pair<std::string, std::pair<int,int>>& a, const std::pair<std::string, std::pair<int,int>>& b)->bool{ return a.second.first <= b.second.first; } );
     if(largest->second.first<=switches){
       continue;
     }
        map.erase(largest->first);
   }
     if (map.find(total_string) == map.end()) {
      std::pair<int,int> p(switches,1);
        map[total_string] = p;
} else {
      std::pair<int,int> p(switches,map[total_string].second+1);
        map[total_string] = p;
}
    //if (iter == 0) printf("CC time:    %.6f s\n", timer.elapsed());
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
  for (auto const& x : map)
{
    if(x.second.first==min_edge_switch){
      counterI++;
    }
}

    fprintf(findex,"\n# of Unique Balanced States With Minimum Number of Switches: %d", counterI);
    fprintf(findex,"\n");
    fprintf(findex,"\nDetailed Analysis:");
    int iter=1;
for (auto const& x : map)
{
            fprintf(findex,"\n# of Switches of Balanced State #%d: %d",iter, x.second.first);
            iter++;

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

  // finalize
  freeGraph(g);
  delete [] minus;
  delete [] einfo;
  delete [] parent;
  delete [] queue;
  delete [] label;
  delete [] border;
  delete [] inCC;
  delete [] inTree;
  delete [] negCnt;
  delete [] root;

  printf("overall runtime with I/O: %.6f s\n", overall.elapsed());
}
