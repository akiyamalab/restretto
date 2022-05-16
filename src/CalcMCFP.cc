#include "CalcMCFP.hpp"

namespace MCFP {
  using namespace std;
  void LeftBackSSP::initializeFlowNodes() {
    for (int i = 0; i < (int)graph.size(); i++) {
      fnode &fnd = flow_nodes[i];
      fnd.to = graph[i].to; // i"から"伸びている辺
      fnd.rev_cap = 0;
      fnd.rev_one_cap = 0;
      if (fnd.to >= 0) {
        fnd.cap = 1;
        flow_nodes[fnd.to].from = i;
      } else {
        fnd.cap = 0;
      }
    }
  }

  int LeftBackSSP::iterativeFlow() {
    if (prev.size() == 0) {
      prev = vector<int>(graph.size(), CASE_UNKNOWN);
    }
    vector<int> dist(graph.size(), INFTY);
    dist[0] = 0;
    for (int now = 0; now < (int)graph.size();) {
      int nxt = now + 1;
      int cost = dist[now];
      if (flow_nodes[now].cap) {
        int to = flow_nodes[now].to;
        int to_cost = cost - graph[now].cost;
        if (dist[to] > to_cost) {
          dist[to] = to_cost;
          prev[to] = now;
        }
      }
      if (now + 1 < (int)graph.size()) {
        int to = now + 1;
        int to_cost = cost;
        if (dist[to] > to_cost) {
          dist[to] = to_cost;
          prev[to] = CASE_STEPFORWARD;
        }
      }
      if (prev[now] != CASE_STEPFORWARD and flow_nodes[now].rev_one_cap) {
        int to = now - 1;
        int to_cost = cost;
        if (dist[to] > to_cost) {
          dist[to] = to_cost;
          prev[to] = CASE_STEPBACKWARD;
          nxt = to;
        }
      }
      if (flow_nodes[now].rev_cap) {
        int to = flow_nodes[now].from;
        int to_cost = cost + graph[to].cost;
        if (dist[to] > to_cost) {
          dist[to] = to_cost;
          prev[to] = CASE_BACKWARD;
          nxt = to;
        }
      }
      now = nxt;
    }

    // flow
    int now = graph.size() - 1;
    while (now != 0) {
      switch (prev[now]) {
      case LeftBackSSP::CASE_BACKWARD: // back
        flow_nodes[now].cap++;
        now = flow_nodes[now].to;
        flow_nodes[now].rev_cap--;
        break;
      case LeftBackSSP::CASE_STEPBACKWARD: // step back
        now = now + 1;
        flow_nodes[now].rev_one_cap--;
        break;
      case LeftBackSSP::CASE_STEPFORWARD: // step forward
        flow_nodes[now].rev_one_cap++;
        now = now - 1;
        break;
      default: // forward
        flow_nodes[now].rev_cap++;
        now = prev[now];
        flow_nodes[now].cap--;
      }
    }

    return dist[graph.size()-1];
  }

  int LeftBackSSP::run(int supply) {
    int ret = 0;
    initializeFlowNodes();
    for (int loop = 0; loop < supply; loop++) {
      ret += iterativeFlow();
    }
    return ret;
  }

  int runLeftBackSSP(const vector<node> &graph, int supply, vector<int>& ret) { //supply: メモリに一度に収まるフラグメントグリッドの個数
    LeftBackSSP solver(graph);
    int c = -solver.run(supply - 1); // IPSJ2018「キャッシュメモリの容量が M のオフラインキャッシュ問題のメモリ戦略の最適化は，このグラフに流量 M−1 のフローを流すことに対応する」なぜこうなるのかまでは引用論文をよく見ないと分からない．
    int sz = graph.size(); //データセットにおけるフラグメントの総個数(=重複を許す)
    ret.clear();
    ret.resize(sz, -1);
    queue<int> empty;
    for (int i = 0; i < supply; ++i) {
      empty.push(i);
    }
    for (int i = 0; i < sz; ++i) {
      if (solver.flow_nodes[i].to == i) {
        ret[i] = ret[i - 1];
      }
      else if (solver.flow_nodes[i].rev_cap) {
        ret[i] = ret[solver.flow_nodes[i].from - 1];
      }
      else {
        assert(!empty.empty());
        ret[i] = empty.front();
        empty.pop();
      }
      // not flow or no out edge
      // be careful about self-loop
      if (i == sz - 1 || (solver.flow_nodes[i + 1].to != i + 1 && (solver.flow_nodes[i + 1].cap == 1 || solver.flow_nodes[i + 1].to < 0))) {
        empty.push(ret[i]);
      }
    }
    return c;
  }
}
