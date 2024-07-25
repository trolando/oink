/*
 * File prepared by Paweł Parys
 */

#include <stack>
#include <queue>
#include <cassert>
#include <climits>
#include <numeric> // for iota

#include "zlkpp.hpp"

namespace pg {

ZLKPPSolver::ZLKPPSolver(Oink& oink, Game& game, int variant) : Solver(oink, game), variant(variant) {}

bool ZLKPPSolver::get_attractor(int player, std::vector<int> &nodes) {
    // initially is_in_attractor[v] = 0 and num_successors[v] = -1 for all v

    std::queue<int> Q;
    std::vector<int> to_be_cleaned;
    bool changed = false;

    for (int v : nodes) {
        is_in_attractor[v] = true;
        Q.push(v);
    }
    
    while (!Q.empty()) {
        int v = Q.front();
        Q.pop();
        for (const int *ptr2 = ins(v); *ptr2 >= 0; ++ptr2) {
            int v2 = *ptr2; // predecessor of v
            if (!cur_nodes_bm[v2] || is_in_attractor[v2])
                continue;
            if (owner(v2) != player) {
                if (num_successors[v2] < 0) {
                    // initially num_successors is -1, so at the end v is automatically subtracted
                    for (const int *ptr3 = outs(v2); *ptr3 >= 0; ++ptr3)
                        if (cur_nodes_bm[*ptr3])
                            ++num_successors[v2];
                    assert(num_successors[v2] >= 0);
                    to_be_cleaned.push_back(v2);
                }
                else
                    --num_successors[v2];
                if (num_successors[v2])
                    continue;
            }
            else
                strategy[v2] = v;
            is_in_attractor[v2] = true;
            nodes.push_back(v2);
            changed = true;
            Q.push(v2);
        }
    }

    // cleanup
    for (int v : nodes)
        is_in_attractor[v] = false;
    for (int v : to_be_cleaned)
        num_successors[v] = -1;

    // returned value = have we added any nodes to the attractor
    return changed;
}

void ZLKPPSolver::remove_nodes(const std::vector<int> nodes) { // "nodes" need not to be sorted
    for (int v : nodes) {
        cur_nodes_bm[v] = false;
        cur_nodes_next[cur_nodes_prev[v]] = cur_nodes_next[v];
        cur_nodes_prev[cur_nodes_next[v]] = cur_nodes_prev[v];
        if (v == cur_first_node)
            cur_first_node = cur_nodes_next[v];
    }
    cur_num_nodes -= nodes.size();
}

void ZLKPPSolver::restore_nodes(const std::vector<int> nodes) {
    // we assume that "nodes" were previously removed by "remove_nodes" (and that they were not used in between, so links in these nodes are correct)
    for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
        int v = *it;
        cur_nodes_bm[v] = true;
        cur_nodes_next[cur_nodes_prev[v]] = v;
        cur_nodes_prev[cur_nodes_next[v]] = v;
        if (cur_nodes_prev[v] >= v)
            cur_first_node = v;
    }
    cur_num_nodes += nodes.size();
}

std::vector<int> ZLKPPSolver::get_cur_nodes() {
    std::vector<int> ret;
    if (cur_num_nodes) {
        int v = cur_first_node;
        for (;;) {
            ret.push_back(v);
            v = cur_nodes_next[v];
            if (v == cur_first_node)
                break;
        }
    }
    return ret;
}

void ZLKPPSolver::set_cur_nodes(const std::vector<int> nodes) {  // it assumes that the "current game" is a subset of "nodes", and that "nodes" are sorted
    assert(!nodes.empty()); // we assume that "nodes" is nonempty
    for (unsigned int a = 0; a < nodes.size(); ++a) {
        int v = nodes[a], v2 = nodes[(a + 1) % nodes.size()];
        cur_nodes_bm[v] = true;
        cur_nodes_next[v] = v2;
        cur_nodes_prev[v2] = v;
    }
    cur_first_node = nodes.front();
    cur_num_nodes = nodes.size();
}

std::vector<int> ZLKPPSolver::get_nodes_of_max_priority(int max_priority) {
    std::vector<int> top;
    if (!cur_num_nodes)
        return top;
    int v = cur_nodes_prev[cur_first_node];
    while (priority(v) == max_priority) {
        top.push_back(v);
        const int *ptr2 = outs(v);
        for (;; ++ptr2) {
            assert(*ptr2 >= 0);
            if (cur_nodes_bm[*ptr2])
                break;
        }
        strategy[v] = *ptr2; // we can go to an arbitrary node that is currently in the game
        if (v == cur_first_node)
            break;
        v = cur_nodes_prev[v];
    }
    return top; // the returned vector is sorted
}

bool ZLKPPSolver::do_step(int max_priority, int prec_cur, int prec_opo, int &num_nodes_in_h, bool &reached_bottom_cur, bool &reached_bottom_opo) {
    num_nodes_in_h = INT_MAX;
    if (!cur_num_nodes)
        return false;
    if (prec_opo < min_dominion) {
        reached_bottom_opo = true;
        return false;
    }
    ++iterations;
    assert(cur_num_nodes); // we assume that the current game is nonempty
    auto top = get_nodes_of_max_priority(max_priority);
    get_attractor(max_priority % 2, top);
    remove_nodes(top);
    auto lose = solve(max_priority - 1, prec_opo, prec_cur, reached_bottom_opo, reached_bottom_cur);
    num_nodes_in_h = cur_num_nodes;
    restore_nodes(top);
    bool changed = get_attractor((max_priority - 1) % 2, lose);
    remove_nodes(lose);
    return changed;
}

void ZLKPPSolver::solve_liverpool(int max_priority, int prec_cur, int prec_opo, bool &reached_bottom_cur, bool &reached_bottom_opo) {
    int initial_num_nodes = cur_num_nodes, dummy;
    if (prec_opo / 2 >= min_dominion)
        solve_liverpool(max_priority, prec_cur, prec_opo / 2, reached_bottom_cur, reached_bottom_opo);
    if (initial_num_nodes <= prec_opo / 2)
        return;
    if (do_step(max_priority, prec_cur, prec_opo, dummy, reached_bottom_cur, reached_bottom_opo) && prec_opo / 2 >= min_dominion)
        solve_liverpool(max_priority, prec_cur, prec_opo / 2, reached_bottom_cur, reached_bottom_opo);
}

std::vector<int> ZLKPPSolver::solve(int max_priority, int prec_cur, int prec_opo, bool &reached_bottom_cur, bool &reached_bottom_opo) {
    if (!cur_num_nodes)
        return std::vector<int>();
    auto saved_nodes = get_cur_nodes();
    int num_nodes_in_h;
    switch (variant) {

    case ZLK_STANDARD:
        while (do_step(max_priority, prec_cur, prec_opo, num_nodes_in_h, reached_bottom_cur, reached_bottom_opo));
        break;

    case ZLK_WARSAW:
        bool reached_bottom_opo_last, atr_greater;
        atr_greater = true;
        while (atr_greater) {
            reached_bottom_opo_last = false;
            atr_greater = do_step(max_priority, prec_cur, prec_opo / 2, num_nodes_in_h, reached_bottom_cur, reached_bottom_opo_last);
            reached_bottom_opo |= reached_bottom_opo_last;
        }
        if (num_nodes_in_h <= prec_opo / 2 || !reached_bottom_opo_last)
            break;
        if (do_step(max_priority, prec_cur, prec_opo, num_nodes_in_h, reached_bottom_cur, reached_bottom_opo))
            while (do_step(max_priority, prec_cur, prec_opo / 2, num_nodes_in_h, reached_bottom_cur, reached_bottom_opo));
        break;
    
    case ZLK_LIVERPOOL:
        solve_liverpool(max_priority, prec_cur, prec_opo, reached_bottom_cur, reached_bottom_opo);
        break;
    }
    
    // we should set strategy in nodes of maximal priority to an arbitrary node that remains in the subgame now:
    get_nodes_of_max_priority(max_priority);

    auto winning_region = get_cur_nodes();
    set_cur_nodes(saved_nodes); // the current game remains unchanged at the end
    return winning_region;
}

void ZLKPPSolver::run() {
    min_dominion = 2;
    for (int v = 0; v < nodecount(); ++v)
        for (const int *ptr2 = outs(v); *ptr2 >= 0; ++ptr2)
            if (*ptr2 == v) {
                min_dominion = 1;
                goto finish;
            }
    finish:
    
    iterations = 0;

    cur_nodes_bm = new bool[nodecount()];
    std::fill(cur_nodes_bm, cur_nodes_bm + nodecount(), true);
    
    cur_nodes_next = new int[nodecount()];
    std::iota(cur_nodes_next, cur_nodes_next + nodecount() - 1, 1);
    cur_nodes_next[nodecount() - 1] = 0;

    cur_nodes_prev = new int[nodecount()];
    std::iota(cur_nodes_prev + 1, cur_nodes_prev + nodecount(), 0);
    cur_nodes_prev[0] = nodecount() - 1;

    cur_first_node = 0;

    cur_num_nodes = nodecount();

    num_successors = new int[nodecount()];
    std::fill(num_successors, num_successors + nodecount(), -1);

    is_in_attractor = new bool[nodecount()];
    std::fill(is_in_attractor, is_in_attractor + nodecount(), 0);
    
    strategy = new int[nodecount()];
    
    // remove disabled nodes (they could be disabled by preprocessing)
    for (int v = 0; v < nodecount(); ++v)
        if (disabled[v])
            remove_nodes(std::vector<int>(1, v));

    assert(cur_num_nodes);
    int max_priority = priority(cur_nodes_prev[cur_first_node]);
    int pow2 = 0;
    while ((1 << pow2) - 1 < nodecount())
        ++pow2;
    bool reached_bottom_cur, reached_bottom_opo;
    auto win = solve(max_priority, (1 << pow2) - 1, (1 << pow2) - 1, reached_bottom_cur, reached_bottom_opo);

    int player = max_priority % 2;
    for (int v = cur_nodes_prev[cur_first_node];;) {
        int who;
        if (!win.empty() && win.back() == v) {
            who = player;
            win.pop_back();
        } else
            who = !player;
        Solver::solve(v, who, strategy[v]);
        if (v == cur_first_node)
            break;
        v = cur_nodes_prev[v];
    }

    logger << "solved with " << iterations << " iterations." << std::endl;

    delete[] cur_nodes_bm;
    delete[] cur_nodes_next;
    delete[] cur_nodes_prev;
    delete[] num_successors;
    delete[] is_in_attractor;
    delete[] strategy;
}

}
