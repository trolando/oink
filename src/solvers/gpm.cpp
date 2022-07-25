#include <algorithm>

#include "gpm.hpp"
#include "oink/uintqueue.hpp"

namespace pg {

Small::Small(const Game &game, const int player) : game(game)
{
    this->player = player;
    this->n = game.nodecount();                  // n := number of vertices
    this->d = game.priority(game.nodecount()-1); // d := highest priority
    this->l = (d+2)/2;                           // l := length of the tuples

    // count how many vertices of each priority there are
    this->counts = new int[d+1];
    for (int p=0; p<=d; p++) this->counts[p] = 0;
    for (int v=0; v<n; v++) this->counts[game.priority(v)]++;

    // allocate and initialize the vertices
    data = new int[l * (n+4)];
    std::fill(data, data + (l*(n+4)), '\0');
    top.resize(n+4);
}

Small::Small(const Small &src) : game(src.game)
{
    this->player = src.player;
    this->n = src.n;
    this->d = src.d;
    this->l = src.l;
    this->counts = new int[d+1];
    this->data = new int[l * (n+4)];
    std::copy(src.counts, src.counts+(d+1), this->counts);
    std::copy(src.data, src.data+(l*(n+4)), this->data);
    this->top = src.top;
}

Small::~Small()
{
    delete[] data;
    delete[] counts;
}

void
Small::set(const Measures &other)
{
    auto src = dynamic_cast<const Small*>(&other);

    delete[] data;
    delete[] counts;

    this->player = src->player;
    this->n = src->n;
    this->d = src->d;
    this->l = src->l;
    this->counts = new int[d+1];
    this->data = new int[l * (n+4)];
    std::copy(src->counts, src->counts+(d+1), this->counts);
    std::copy(src->data, src->data+(l*(n+4)), this->data);
    this->top = src->top;
}

void Small::copy(int from, int to)
{
    // decrease counts?
    if (to >= 0 && top[to+4] == 0 && top[from+4] == 1) {
        int pr = game.priority(to);
        if ((pr&1)==player) {
            counts[pr]--;
        }
    }

    int *fp = data + (from+4)*l;
    int *tp = data + (to+4)*l;
    for (int i=0; i<l; i++) *tp++ = *fp++;
    top[to+4] = top[from+4];
}

void Small::copy_from(const Measures& mother, int from, int to)
{
    const Small *other = dynamic_cast<const Small*>(&mother);

    // decrease counts?
    if (to >= 0 && top[to+4] == 0 && other->top[from+4] == 1) {
        int pr = game.priority(to);
        if ((pr&1)==player) {
            counts[pr]--;
        }
    }

    int *fp = other->data + (from+4)*l;
    int *tp = data + (to+4)*l;
    for (int i=0; i<l; i++) *tp++ = *fp++;
    top[to+4] = other->top[from+4];
}

void Small::copy_bot(int to)
{
    int *tp = data + (to+4)*l;
    for (int i=0; i<l; i++) *tp++ = 0;
    top[to+4] = 0;
}

void Small::copy_top(int to)
{
    // decrease counts?
    if (to >= 0 && top[to+4] == 0) {
        int pr = game.priority(to);
        if ((pr&1)==player) {
            counts[pr]--;
        }
    }
            
    top[to+4] = 1;
}

bool Small::is_top(int tgt) const
{
    return top[tgt+4];
}

void Small::see(int tgt, int priority)
{
    if (top[tgt+4]) return; // don't care if already at top

    // for player 0 tuple of 0, 2, 4, 6, 8...
    // for player 1 tuple of 1, 3, 5, 7, 9...

    int *tp = data + (tgt+4)*l;
    for (int i=0; i<l; i++) {
        int tuple_p = this->player + 2*i;
        if (tuple_p < priority) {
            tp[i] = 0;
        } else if (tuple_p == priority) {
            // check if reached the max! if so, carry on
            if (tp[i] >= counts[tuple_p]) {
                tp[i] = 0;
                // on carry overflow, this is a \Top
                if ((priority+2) > d) {
                    top[tgt + 4] = 1;
                    return; // done
                } else {
                    priority+=2;
                }
            } else {
                tp[i]++;
                return; // done
            }
        } else {
            break;
        }
    }
}

void Small::inc(int tgt)
{
    if (top[tgt+4]) return; // nothing to do

    int *tp = data + (tgt+4)*l;
    for (int i=0; i<l; i++) {
        int tuple_p = this->player + 2*i;
        if (tp[i] == counts[tuple_p]) {
            tp[i] = 0;
        } else {
            tp[i]++;
            return;
        }
    }

    top[tgt+4] = true;
}

/**
 * Return -1 if value[i] < value[j]
 * Return  0 if value[i] = value[j]
 * Return  1 if value[i] > value[j]
 */
int Small::compare(int i, int j) const
{
    // Handle \Top cases
    if (top[i+4]) {
        if (top[j+4]) return 0;
        return 1;
    }
    if (top[j+4]) return -1;
    // Compare tuples, starting at the highest (last) value
    int *ip = data + (i+4)*l;
    int *jp = data + (j+4)*l;
    for (int x=l-1; x>=0; x--) {
        if (ip[x] < jp[x]) return -1;
        if (ip[x] > jp[x]) return 1;
    }
    // No difference? then equal
    return 0;
}

bool Small::eq(int i, int j) const
{
    if (top[i+4]) {
        if (top[j+4]) return true;
        return false;
    }
    if (top[j+4]) return false;


    int *ip = data + (i+4)*l;
    int *jp = data + (j+4)*l;
    for (int x = 0; x < l; x++) {
        if (ip[x] != jp[x]) return false;
    }
    return true;
}

void Small::stream(std::ostream &out, int i) const
{
    int *ip = data + (i+4)*l;
    if (top[i+4]) {
        out << "Top";
    } else {
        out << "[";
        for (int x=l-1; x >= 0; x--) out << " " << ip[x];
        out << " ]";
    }
}

Ordered::Ordered(const Game &game, const int player)
{
    this->player = player;
    this->n = game.nodecount();                  // n := number of vertices
    this->d = game.priority(game.nodecount()-1); // d := highest priority
    priorities.resize(d + 1);

    // find highest own priority
    this->max_pr = -1;
    this->max_opp_pr = -1;
    for (int i=game.nodecount()-1; i>=0; i--) {
        const int pr = game.priority(i);
        if ((pr&1) == player ) {
            if (pr > this->max_pr) {
                this->max_pr = pr;
            }
        } else if(pr > max_opp_pr) {
            this->max_opp_pr = pr;
        }
        priorities[pr] = true;
    }
    if (max_opp_pr == -1) {
        max_opp_pr = player;
    }

    // compute number of vertices with player's priority
    this->N = 0;
    for (auto v=0; v<game.nodecount(); v++) {
        if ((game.priority(v)&1) == player) N++;
    }

    // see the paper: use l = ceil(log2(N+1)) + 1
    l = 1;
    while ((1L<<(l-1)) <= N) l++;

    //std::cout << "n=" << N << std::endl;
    //std::cout << "l=" << l << std::endl;

    // count how many vertices of each priority there are
    this->counts = new int[d+1];
    for (int p=0; p<=d; p++) this->counts[p] = 0;
    for (int v=0; v<n; v++) this->counts[game.priority(v)]++;

    // allocate and initialize the vertices
    data = new int[l * (n+4)];
    std::fill(data, data + (l*(n+4)), -1);
    top.resize(n+4);
}

Ordered::Ordered(const Ordered &src)
{
    this->player = src.player;
    this->n = src.n;
    this->N = src.N;
    this->d = src.d;
    this->l = src.l;
    this->max_pr = src.max_pr;
    this->max_opp_pr = src.max_opp_pr;
    this->counts = new int[d+1];
    this->data = new int[l * (n+4)];
    std::copy(src.counts, src.counts+(d+1), this->counts);
    std::copy(src.data, src.data+(l*(n+4)), this->data);
    this->top = src.top;
    this->priorities = src.priorities;
}

Ordered::~Ordered()
{
    delete[] data;
    delete[] counts;
}

void
Ordered::set(const Measures &other)
{
    auto src = dynamic_cast<const Ordered*>(&other);

    delete[] data;
    delete[] counts;

    this->player = src->player;
    this->n = src->n;
    this->N = src->N;
    this->d = src->d;
    this->l = src->l;
    this->max_pr = src->max_pr;
    this->max_opp_pr = src->max_opp_pr;
    this->counts = new int[d+1];
    this->data = new int[l * (n+4)];
    std::copy(src->counts, src->counts+(d+1), this->counts);
    std::copy(src->data, src->data+(l*(n+4)), this->data);
    this->top = src->top;
    this->priorities = src->priorities;
}

void Ordered::copy(int from, int to)
{
    int *fp = data + (from+4)*l;
    int *tp = data + (to+4)*l;
    for (int i=0; i<l; i++) *tp++ = *fp++;
    top[to+4] = top[from+4];
}

void Ordered::copy_from(const Measures& mother, int from, int to)
{
    const Ordered *other = dynamic_cast<const Ordered*>(&mother);
    int *fp = other->data + (from+4)*l;
    int *tp = data + (to+4)*l;
    for (int i=0; i<l; i++) *tp++ = *fp++;
    top[to+4] = other->top[from+4];
}

void Ordered::copy_bot(int to)
{
    int *tp = data + (to+4)*l;
    for (int i=0; i<l; i++) *tp++ = -1;
    top[to+4] = 0;
}

void Ordered::copy_top(int to)
{
    top[to+4] = 1;
}

bool Ordered::is_top(int tgt) const
{
    return top[tgt+4];
}

int Ordered::val(int tgt) const
{
    int res = 0;
    int *tp = data + (tgt+4)*l;
    bool found = false;
    for (int x=l-1; x>=0; x--) {
        res *= 2;
        if (!found) {
            res += (tp[x] != -1) ? 1 : 0;
            if ((tp[x] & 1) == player) {
                found = true;
            }
        }
    }
    return res;
}

void Ordered::see(int tgt, int priority)
{
    if (top[tgt+4]) return; // nothing to do
    copy(tgt, -2);
    copy(tgt, -3);
    up(-2, priority);
    if (bump(-3)) {
        up(-3, priority);
        if (compare(-3, -2) < 0) {
            copy(-3, -2);
        }
    }
    copy(-2, tgt);
    int last = *(data + (tgt+4)*l + l - 1);
    if (last != -1 && (last & 1) == player) {
        top[tgt+4] = true;
    }
}

bool Ordered::up(int tgt, int p)
{
    int *tp = data + (tgt+4)*l;

    // return if already top
    if (top[tgt+4]) return false;

    // first set the first _ or non-<player> value to "priority"
    for (int i=0; i<l; i++) {
        if (tp[i] == -1 or ((tp[i] & 1) != player)) {
            tp[i] = p; // set it
            break;
        } else {
            tp[i] = -1;
        }

    }

    // next set the highest non-_ < p value to "p" and the rest to -1
    for (int i=l-1; i>=0; i--) {
        if (tp[i] == -1) continue;
        if ((p&1) == player) {
            if (tp[i] >= p) continue;
        } else {
            if (tp[i] > p) continue;
        }
        tp[i] = p;
        for (int j=i-1; j>=0; j--) tp[j] = -1; // set rest to -1
        return true;
    }

    return true;
}

bool Ordered::bump(int tgt)
{
    if (top[tgt+4]) return false; // nothing to do

    int *tp = data + (tgt+4)*l;

    for (int i = 1; i<l; i++) {
        if (tp[i] == -1) {
            // if the current is _, then the smallest increase is max_opp_pr
            tp[i] = max_opp_pr;
            for (int z=i-1; z>=0; z--) tp[z] = -1;
            return true;
        } else if ((tp[i] & 1) != player) {
            // if the current is of opponent parity, then subtract 2 OR change parity
            if (tp[i] > (1-player)) {
                // if odd != -1 / even != 0 then we can just subtract 2 and we are done
                tp[i] -= 2;
                for (int z=i-1; z>=0; z--) tp[z] = -1;
                return true;
            } else if (i < l-1 and tp[i+1] >= player) {
                // if is 1/0, then can set to 0/1 if the previous is set and we're not the highest...
                tp[i] = player;
                for (int z=i-1; z>=0; z--) tp[z] = -1;
                return true;
            }
        } else {
            // if it's of parity <player>, and the next one is higher, then
            // increase and done probably the only one that is relevant, since
            // this algorithm is used when an "even" sequence must be bumped
            if (tp[i] + 2 <= max_pr and i < l-1 and (tp[i+1] == -1 || tp[i+1] >= tp[i] + 2)) {
                tp[i] += 2;
                for (int z=i-1; z>=0; z--) tp[z] = -1;
                return true;
            }
        }
    }

    return false;
}

void Ordered::inc(int tgt)
{
    if (top[tgt+4]) return;

    int *tp = data + (tgt+4)*l;

    for (int i=0; i<l; i++) {
        int min = d;
        for (int j=i+1; j<l; j++) {
            if (tp[j] == -1) continue;
            if (tp[j] < min) min = tp[j];
        }
        if (tp[i] == -1) {
            tp[i] = (min&1) == player && min != 0 ? min-1 : min;
            return;
        } else if ((tp[i]&1) == player) {
            // only bump if the next higher allows it
            if ((tp[i]+2) <= min) {
                tp[i] += 2;
                return;
            }
        } else {
            // 5->3, 4->2, 3->1, 2->0, 1->0, 0->1
            if (tp[i] >= 2) {
                tp[i] -= 2;
                return;
            } else if (tp[i] == 1) {
                tp[i] = 0;
                return;
            } else if (min >= 1) {
                tp[i] = 1;
                return;
            }
        }
        tp[i] = -1;
    }

    top[tgt+4] = true; // cannot go higher, therefore top
}

/**
 * Return -1 if value[i] < value[j]
 * Return  0 if value[i] = value[j]
 * Return  1 if value[i] > value[j]
 */
int Ordered::compare(int i, int j) const
{
    // Handle \Top cases
    if (top[i+4]) {
        if (top[j+4]) return 0;
        return 1;
    }
    if (top[j+4]) return -1;

    // Compare tuples, starting at the highest (last) value
    int *ip = data + (i+4)*l;
    int *jp = data + (j+4)*l;
    for (int x=l-1; x>=0; x--) {
        if (ip[x] == jp[x]) continue;
        // cases where i.x is _ or j.x is _
        if (ip[x] == -1) return -1;
        if (jp[x] == -1) return 1;
        // case distinction based on parity
        int pr_i = ip[x];
        if ((pr_i&1) != player) pr_i = -pr_i;
        int pr_j = jp[x];
        if ((pr_j&1) != player) pr_j = -pr_j;
        if (pr_i < pr_j) return -1;
        if (pr_i > pr_j) return 1;
    }
    // No difference? then equal
    return 0;
}

bool Ordered::eq(int i, int j) const
{
    if (top[i+4]) {
        if (top[j+4]) return true;
        return false;
    }
    if (top[j+4]) return false;


    int *ip = data + (i+4)*l;
    int *jp = data + (j+4)*l;
    for (int x = 0; x < l; x++) {
        if (ip[x] != jp[x]) return false;
    }
    return true;
}

void Ordered::stream(std::ostream &out, int i) const
{
    int *ip = data + (i+4)*l;
    if (top[i+4]) {
        out << "Top";
    } else {
        out << "[";
        for (int x = l - 1; x >= 0; x--) {
            if (ip[x] == -1) out << " _";
            else out << " " << ip[x];
        }
        out << " ]";
    }
}

GPMSolver::GPMSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

GPMSolver::~GPMSolver()
{
}

bool
GPMSolver::lift(Measures &pm, const int player, int v) {
    bool updated = false;
    if (owner(v) == player) {
        // Take the highest successor
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            int to = *curedge;
            if (!G[to]) continue;
            pm.copy(to, -1);
            pm.see(-1, priority(v));
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "Lift option (takes highest) " << label_vertex(v) << " to " << label_vertex(to) << " ";
                pm.stream(logger, -1);
                logger << std::endl;
            }
#endif
            if ((!updated || pm.compare(-1, -4) > 0)) {
                pm.copy(-1, -4);
                updated = true;
            }
        }
    } else {
        // Take the lowest successor
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            int to = *curedge;
            if (!G[to]) continue;
            pm.copy(to, -1);
            pm.see(-1, priority(v));
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "Lift option (takes lowest) " << label_vertex(v) << " to " << label_vertex(to) << " ";
                pm.stream(logger, -1);
                logger << std::endl;
            }
#endif
            if ((!updated || pm.compare(-1, -4) <= 0)) {
                pm.copy(-1, -4);
                updated = true;
            }
        }
    }
    return updated;
}

bool
GPMSolver::update(Measures &pm, const int player) {
    auto &Q = player == 0 ? Q0 : Q1;
    //auto &R = player == 0 ? R0 : R1;
    auto &T = player == 0 ? T0 : T1;
    //auto &B = player == 0 ? B0 : B1;
    auto &L = player == 0 ? L0 : L1;

    QQ.reset();
    L.reset();
    for (auto v = Q.find_first(); v != bitset::npos; v = Q.find_next(v)) {
        lift_attempts++;
        if (lift(pm, player, v)) {
            auto res = pm.compare(-4, v);
            if (res < 0) {
#ifndef NDEBUG
                if (trace) {
                    logger << "\033[1m" << player << "-lifting\033[m \033[36m" << label_vertex(v) << "\033[m from \033[1m";
                    pm.stream(logger, v);
                    logger << " to ";
                    pm.stream(logger, -4);
                    logger << std::endl;
                }
#endif
                //LOGIC_ERROR;
            }
            if (res > 0) {
                lifts++;
#ifndef NDEBUG
                if (trace) {
                    logger << "\033[1m" << player << "-lifting\033[m \033[36m" << label_vertex(v) << "\033[m from \033[1m";
                    pm.stream(logger, v);
                    logger << "\033[m to \033[1m";
                    pm.stream(logger, -4);
                    logger << "\033[m" << std::endl;
                }
#endif
                L[v] = 1;
                // actually lift
                pm.copy(-4, v);
                // update B and R
                if (pm.is_top(-4)) {
                    T[v] = true;
                    //R[v] = false;
                }
                //B[v] = false;
                // update QQ
                for (auto curedge = ins(v); *curedge != -1; curedge++) {
                    int from = *curedge;
                    if (!T[from]) QQ[from] = true;
                }
            }
        }
    }

    Q.swap(QQ);
    return Q.any();
}


/**
 * Attract as player <pl> via <v> to <Z>, vertices in <R>. Only escape to <R>.
 * Add attracted vertices to <Z> and to queue <Q>.
 */
void
GPMSolver::attractVertices(int pl, int v, bitset &R, bitset &Z)
{
    // attract vertices with an edge to <v>
    for (auto curedge = ins(v); *curedge != -1; curedge++) {
        int from = *curedge;
        if (!Z[from] && R[from]) {
            if (owner(from) != pl) {
                // check if opponent can escape
                bool escapes = false;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (R[to] && !Z[to]) {
                        escapes = true;
                        break;
                    }
                }
                if (escapes) continue;
            }
            // attract
            Z[from] = true;
            Queue.push(from);
        }
    }
}


void
GPMSolver::shortcuts(const int player, Measures &pm0, Measures &pm1)
{
    auto& Q = player == 0 ? Q0 : Q1;
    auto& Qo = player == 0 ? Q1 : Q0;
    auto& T = player == 0 ? T0 : T1;
    auto& To = player == 0 ? T1 : T0;
    auto& pm = player == 0 ? pm0 : pm1;
    auto& pm_o = player == 0 ? pm1 : pm0;
    auto& L = player == 0 ? L0 : L1;

/*
    // Attract towards Bottom!
    QQ = B;
    for (auto v = QQ.find_first(); v != bitset::npos; v = QQ.find_next(v)) {
        logger << "bottom: " << label_vertex(v) << std::endl;
        Queue.push(v);
    }
    while (Queue.nonempty()) {
        const int v = Queue.pop();
        for (auto curedge = ins(v); *curedge != -1; curedge++) {
            int from = *curedge;
            if (QQ[from] || !R[from]) continue; // don't attract
            if (owner(from) == player) {
                // see we can escape (not have to play to bottom)
                bool escapes = false;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (!R[to] || QQ[to]) continue; // not a potential escape
                    // it is only a valid escape if the measure obtained by going from
                    // <from> to <to> actually is >= the current measure...
                    pm.copy(to, -1);
                    pm.see(-1, priority(from));
                    if (pm.compare(-1, from) >= 0) {
                        escapes = true;
                        break;
                    }
                }
                if (escapes) continue;
            }
            // attract
            QQ[from] = true;
            Queue.push(from);
        }
    }
    QQ ^= G;
    QQ -= T;
    for (auto v = QQ.find_first(); v != bitset::npos; v = QQ.find_next(v)) {
        logger << "player " << player << " should win " << label_vertex(v) << std::endl;
    }
*/




    // logger<<"running shortcuts for player " << player << std::endl; 

    // now find vertices that are "stable"
    // any vertex in L is lifted previously, so it's increasing for L
    // also add all in T
    L |= T;
    for (auto v = L.find_first(); v != bitset::npos; v = L.find_next(v)) {
        Queue.push(v);
    }
    while (Queue.nonempty()) {
        const int v = Queue.pop();
        for (auto curedge = ins(v); *curedge != -1; curedge++) {
            int from = *curedge;
            if (L[from] || !G[from]) continue; // don't attract
            if (owner(from) != player) {
                // see if opponent can escape
                bool escapes = false;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (!G[to] || L[to]) continue; // not a potential escape
                    // it is only a valid escape if the measure obtained by going from
                    // <from> to <to> actually is == the current measure...
                    pm.copy(to, -1);
                    pm.see(-1, priority(from));
                    if (pm.compare(-1, from) <= 0) {
                        escapes = true;
                        break;
                    }
                }
                if (escapes) continue;
            }
            // attract
            L[from] = true;
            Queue.push(from);
        }

    }
    L ^= G;
    for (auto v = L.find_first(); v != bitset::npos; v = L.find_next(v)) {
        // logger << "player " << (1-player) << " should win " << label_vertex(v) << std::endl;
        pm_o.copy_top(v);
        To[v] = true;
        //Ro[v] = false;
        Qo[v] = false;
        //Bo[v] = false;
        for (auto curedge = ins(v); *curedge != -1; curedge++) {
            int from = *curedge;
            if (G[from] && !To[from]) Qo[from] = true;
        }
    }
}


void
GPMSolver::solve(Measures &pm, const int player)
{
    // extract strategy
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        if (!pm.is_top(v)) {
            // won by ~X
            str[v] = -1;
            if (owner(v) == player) {
                // find lowest prog...
                for (auto curedge = outs(v); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (G[to] && !pm.is_top(to)) {
                        if (str[v] == -1) {
                            pm.copy(to, -1);
                            str[v] = to;
                        } else if (pm.compare(to, -1) < 0) {
                            pm.copy(to, -1);
                            str[v] = to;
                        }
                    }
                }

                if (str[v] == -1) LOGIC_ERROR;
            }
        }
    }
    // mark as solved
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        if (!pm.is_top(v)) {
#ifndef NDEBUG
            if(trace) {
                logger << "player " << player << " wins " << label_vertex(v) << " with strategy " << label_vertex(str[v]) << " and ";
                pm.stream(logger, v);
                logger << std::endl;
            }
#endif
            Solver::solve(v, player, str[v]);
            G[v] = false;
        }
    }
    Solver::flush();
}

bool
GPMSolver::parseOptions(std::string& opts) {
    std::optional<MeasureKind> parsed = parse_measure_kind(opts);
    if (parsed.has_value()) {
        measure_kind = *parsed;
        if (trace) {
            logger << "selected ";
            stream_measure_kind(logger, measure_kind);
            logger << " for GPM" << std::endl;
        }
        return true;
    } else {
        logger << "unrecognised option for GPM solver, possible options: ";
        stream_measure_kinds(logger);
        logger << std::endl;
        return false;
    }
}

void
GPMSolver::run()
{
    str = new int[nodecount()];
    G = disabled;
    G.flip();

    //Q = disabled;
    // B0 = B1 = G;
    // R0 = R1 = G;

    Q0 = Q1 = G;
    QQ.resize(nodecount());
    Queue.resize(nodecount());
    T0.resize(nodecount());
    T1.resize(nodecount());
    L0.resize(nodecount());
    L1.resize(nodecount());

    Measures *pm0 = new_measure(measure_kind, game, 0);
    Measures *pm1 = new_measure(measure_kind, game, 1);

    bool c0 = true;
    bool c1 = true;

    long iterations = 0;

    while (true) {
        iterations++;
        if (c0) {
            c0 = update(*pm0, 0);
            if (!c0 || (iterations%10) == 0) shortcuts(0, *pm0, *pm1);
        }

        if (c1) {
            c1 = update(*pm1, 1);
            if (!c1 || (iterations%10) == 0) shortcuts(1, *pm0, *pm1);
        }

        if (!c0 && !c1) break;
    }

    solve(*pm1, 0);
    solve(*pm0, 1);

    logger << "solved with " << lifts << " lifts, " << lift_attempts << " lift attempts, " << iterations << " iterations." << std::endl;

    delete[] str;
    delete pm0;
    delete pm1;
}

std::optional<MeasureKind> parse_measure_kind(std::string& opts) {
    if (opts.empty()) {
        return MeasureKind::Ordered;
    } else if (opts == "opm" || opts == "qpt" || opts == "ordered" || opts == "OPM" || opts == "QPT" || opts == "ORDERED") {
        return { MeasureKind::Ordered };
    } else if (opts == "spm" || opts == "small" || opts == "SPM" || opts == "SMALL") {
        return { MeasureKind::Small };
    }

    return std::nullopt;
}

Measures *new_measure(const MeasureKind kind, const Game &game, const int player) {
    switch (kind) {
        case MeasureKind::Small:
            return new Small(game, player);
        case MeasureKind::Ordered:
            return new Ordered(game, player);
        default:
            THROW_ERROR("Measure kind exists but is not implemented");
    }
}

void stream_measure_kind(std::ostream& logger, MeasureKind kind) {
    switch (kind) {
        case MeasureKind::Small:
            logger << "Small Progress Measures (SPM)";
            break;
        case MeasureKind::Ordered:
            logger << "Ordered Progress Measure (OPM)";
            break;
        default:
            THROW_ERROR("Measure kind exists but is not implemented");
            break;
    }
}

void stream_measure_kinds(std::ostream& logger) {
    logger << "spm, opm";
}

}
