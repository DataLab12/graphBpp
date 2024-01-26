/*
Harary Clustering based on Balanced Laplacian and GraphB++ for signed social network graphs

Copyright 2024, Texas State University. All rights reserved.

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
#include <list>


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

static void init(const Graph& g, int* const inCC, EdgeInfo* const einfo, int* const inTree, int* const negCnt,float* const influence,float* const strength)
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
 #pragma omp parallel for default(none) shared(g, influence)
  for (int v = 0; v < g.nodes; v++) {
    influence[v] = 0;
  }
  #pragma omp parallel for default(none) shared(g, strength)
  for (int v = 0; v < g.nodes; v++) {
    strength[v] = 0;
  }
  // zero out inTree and negCnt
  #pragma omp parallel for default(none) shared(g, inTree, negCnt)
  for (int j = 0; j < g.edges; j++) {
    inTree[j] = 0;
    negCnt[j] = 0;
  }
}

static void generateSpanningTree(const Graph& g, const int root, const int seed, EdgeInfo* const einfo, int* const parent, int* const queue, int* const border, int* const label, int* const inTree, int* const negCnt)
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

  // update inTree
  #pragma omp parallel for default(none) shared(g, inTree)
  for (int j = 0; j < g.edges; j++) {
    inTree[j] += g.nlist[j] & 1;
  }

}

static void initMinus(const Graph& g, const EdgeInfo* const einfo, bool* const minus)
{

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

}

static void processCycles(const Graph& g, const int* const label, EdgeInfo* const einfo, bool* const minus)
{
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
          sum += einfo[k].end & 1;
          curr = g.nlist[k] >> 1;
        }
        minus[j] = sum & 1;  // make cycle have even number of minuses
      }
      j--;
    }
  }
}

static void determineCCs(const Graph& g, int* const label, const bool* const minus, int* const count, int* const inCC, int* const negCnt,std::vector<int> *edges_balanced_src,std::vector<int> *edges_balanced_dst,std::vector<bool>*weight_balanced,std::map<std::string,double>*  agreement,std::map<std::string,int>* resolution,std::map<std::string,bool> copy)
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
        std::string total_string="";
        total_string.append(std::to_string(v));
        total_string.append(",");
        total_string.append(std::to_string(nli));
        if((copy)[total_string]==minus[j]){
            (*resolution)[total_string]=(*resolution)[total_string]+1;


        }
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

int src,dst;
bool src_flag,dst_flag;
 for (int v = 0; v < g.nodes; v++) {
    const int beg = g.nindex[v];
    const int end = g.nindex[v + 1];
    for (int j = beg; j < end; j++) {
      const int nli = g.nlist[j] >> 1;
      if(v<nli){
        std::string total_string="";
        total_string.append(std::to_string(v));
        total_string.append(",");
        total_string.append(std::to_string(nli));

        src=v;
        dst=nli;
        src_flag=false;
        dst_flag=false;
        if(((count[label[src]] % 2) ^ 1) == 1){
          src_flag=true;
        }
        if(((count[label[dst]] % 2) ^ 1) == 1){
          dst_flag=true;
        }
    if(src_flag==true && dst_flag==true){
      (*agreement)[total_string]=(*agreement)[total_string]+1;
    }

      }
    }

}

}



std::map<int,std::vector<float>> extract_features(const char* name, int i)
{

  Graph g = readGraph(name);
  const int iterations = i;
  std::vector<int> edges_balanced_src;
  std::vector<int> edges_balanced_dst;
  std::vector<bool> weight_balanced;

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
   
  float* const influence = new float [g.nodes];
  float* const strength = new float [g.nodes];
  std::map<std::string,double> agreement;
  std::map<std::string,int> resolution;
    std::map<std::string,bool> copy;
    int* const degree = new int [g.nodes];
    for (int v = 0; v < g.nodes; v++) {
    degree[v] = 0;
  }

for (int v = 0; v < g.nodes; v++) {
    int index=g.nindex[v];
    int index_max=g.nindex[v+1];
    while(index<index_max){
      if(v<g.nlist[index]){
        degree[v]=degree[v]+1;
        degree[g.nlist[index]]=degree[g.nlist[index]]+1;

      if(g.eweight[index]==-1){
       std::string total_string="";
    total_string.append(std::to_string(v));
    total_string.append(",");
    total_string.append(std::to_string(g.nlist[index]));
        copy[total_string]=1;
      }else{
        std::string total_string="";
    total_string.append(std::to_string(v));
    total_string.append(",");
    total_string.append(std::to_string(g.nlist[index]));
        copy[total_string]=0;
      }
      }
    index=index+1;
    }
  }
   

  


  for (int i = 0; i < g.nodes; i++) root[i] = i;
  std::partial_sort(root, root + std::min(iterations, g.nodes), root + g.nodes, [&](int a, int b) {
    return (g.nindex[a + 1] - g.nindex[a]) > (g.nindex[b + 1] - g.nindex[b]);
  });
  init(g, inCC, einfo, inTree, negCnt,influence,strength);
  for (int iter = 0; iter < iterations; iter++) {
    edges_balanced_src.clear();
    weight_balanced.clear();
    edges_balanced_dst.clear();
    //printf("tree %d\n", iter);

    // generate tree
    generateSpanningTree(g, root[iter % g.nodes], iter + 17, einfo, parent, queue, border, label, inTree, negCnt);

    // initialize plus/minus
    initMinus(g, einfo, minus);
    // find cycles
    processCycles(g, label, einfo, minus);


    // determine connected components

    determineCCs(g, label, minus, queue, inCC, negCnt,&edges_balanced_src,&edges_balanced_dst,&weight_balanced,&agreement,&resolution,copy);
   

    }



for (int v = 0; v < g.nodes; v++) {
    int index=g.nindex[v];
    int index_max=g.nindex[v+1];
    while(index<index_max){
      const int nli = g.nlist[index] >> 1;
      if(v<nli){
          int src=v;
      int dst=nli;
       std::string total_string="";
      total_string.append(std::to_string(src));
            total_string.append(",");
            total_string.append(std::to_string(dst));
      influence[src]= 1.0*influence[src]+ (1.0*agreement[total_string]/iterations);
      strength[src]=1.0*strength[src]+ (1.0*resolution[total_string]/iterations);
      influence[dst]= 1.0*influence[dst]+ (1.0*agreement[total_string]/iterations);
      strength[dst]=1.0*strength[dst]+ (1.0*resolution[total_string]/iterations);
      }
    index=index+1;
    }
  }

std::map<int,std::vector<float>> node_features;
  
  for (int i = 0; i < g.nodes; i++) {
        influence[i]= (1.0*influence[i])/(1.0*degree[i]);
        strength[i]=(1.0*strength[i])/(1.0*degree[i]);
        float value=1.0* inCC[i] / iterations;
        node_features[g.origID[i]].push_back(value);
        float value1=1.0* influence[i];
        node_features[g.origID[i]].push_back(value1);
        float value2=1.0*strength[i];
        node_features[g.origID[i]].push_back(value2);
  }




  // finalize
  freeGraph(g);
  delete [] minus;
  delete [] einfo;
  delete [] parent;
  delete [] queue;
  delete [] border;
  delete [] inCC;
  delete [] inTree;
  delete [] negCnt;
  delete [] root;
  delete [] label;
  delete [] influence;
  delete [] strength;
  delete [] degree;


  return node_features;
}
