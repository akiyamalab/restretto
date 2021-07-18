#ifndef CALC_MCFP_H_
#define CALC_MCFP_H_

#include "common.hpp"
#include <vector>
#include <queue>

namespace MCFP {
  using namespace std;
  struct node {
    int to, cost;
    node(int to, int cost): to(to), cost(cost) {}
  };

  class LeftBackSSP {
  public:
    const int INFTY = 1000000007;
    static const int CASE_BACKWARD = -1;
    static const int CASE_STEPBACKWARD = -2;
    static const int CASE_STEPFORWARD = -3;
    const int CASE_UNKNOWN = -4;
    struct fnode {
      int to;      //a node id where this node point to
      int from;    //a node id where this node is pointed from
      int cap;     //capacity, this number contains 0 or 1
      int rev_cap; //reverse-capacity
      int rev_one_cap; //reverse-capacity to a node which has previous id number
    };

    vector<node> graph;
    vector<fnode> flow_nodes;
  private:
    vector<int> prev;
    void initializeFlowNodes();
    int iterativeFlow();

  public:
    LeftBackSSP(const vector<node> &graph): graph(graph), flow_nodes(vector<fnode>(graph.size())) {}
    int run(int supply);
  };

  int runLeftBackSSP(const vector<node> &graph, int supply, vector<int>& ret);
}

#endif
