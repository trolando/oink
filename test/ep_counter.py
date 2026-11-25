#!/usr/bin/env python3
"""Generate counterexamples for uncapped EP updates.

Use --energy with an even chain length for the endpoint-sweep family. The
default family tests the historical stopping-game policy-improvement rule.

The arena and initial Even moves come from Friedmann (2011), Figure 4 and
Sections 6.4 to 6.7, https://lmcs.episciences.org/1026/pdf. The three forced
vertices at the exit make the stopping-game scores order paths as in the proof.

Until the top bit closes, Friedmann's best Odd reply leads every counter
vertex to x. The two successor paths in such a profile share a suffix. Let
t be the value where they merge. The stopping value of x is 2047/33791, and
the sum of all stopping probabilities outside the exit is at most
1/(31*32**3), so t lies between 0.060 and 0.061. A prefix vertex of rank j
contributes (reward-t)*32**(-j) times the probability of reaching it. The
prefixes are disjoint, so all ranks below the highest differing vertex have
total weight at most 1/31 of its weight. That vertex decides the comparison,
in the direction of its parity, exactly as in Friedmann's path order.

For an Odd vertex u, an alternative whose path avoids u compares two sink
paths with u prepended to both. An alternative whose path revisits u would
close a cycle. Its maximum is even, since in the original counter Odd would
prefer any higher odd cycle to the priority-1 sink. The alternative's path
therefore adds a prefix with even maximum to the selected successor's path.
In both cases the original reply minimizes the stopping score, and the
subsequent Even switches agree with Friedmann's locally optimal rule. Before
the highest bit closes, that reply defeats the current Even strategy at the
genuinely Even-winning d_bits. Strategy verification rejects these policies
throughout the exponentially long prefix of that historical improvement rule.
"""

import argparse


def write_counter(bits):
    game = {}

    def add(name, priority, owner, *successors):
        game[name] = priority, owner, successors

    # EP takes the first edge as its initial strategy at each vertex.
    add("c", 8 * bits + 4, 0, "r", "s")
    add("s", 8 * bits + 6, 0, "x", *(f"f{i}" for i in range(1, bits + 1)))
    add("r", 8 * bits + 8, 0, "x", *(f"g{i}" for i in range(1, bits + 1)))
    for i in range(1, 2 * bits + 1):
        successors = ("c", "s", "r") if i == 1 else ("r", "s", f"t{i - 1}")
        add(f"t{i}", 4 * bits + 2 * i + 1, 0, *successors)
        add(f"a{i}", 4 * bits + 2 * i + 2, 1, f"t{i}")

    # The gate priorities follow Figure 4, which gives distinct priorities.
    # Odd has no exit from e_bits (Friedmann, Section 6.7). Even can win at
    # d_bits by choosing e_bits, so an earlier Even policy cannot certify.
    for i in range(1, bits + 1):
        lane = (f"a{j}" for j in range(1, 2 * i + 1))
        higher = (f"g{j}" for j in range(i + 1, bits + 1))
        exits = (f"d{i}",) if i == bits else (f"d{i}", f"h{i}")
        add(f"d{i}", 4 * i - 1, 0, "r", f"e{i}", "s", *lane)
        add(f"e{i}", 4 * i, 1, *exits)
        add(f"g{i}", 4 * i + 2, 0, f"k{i}", f"f{i}")
        add(f"k{i}", 8 * bits + 4 * i + 5, 0, "x", *higher)
        add(f"f{i}", 8 * bits + 4 * i + 7, 1, f"e{i}")
        add(f"h{i}", 8 * bits + 4 * i + 8, 1, f"k{i}")

    # At x, play enters an Odd-winning cycle. At ranks 1, 2, 3 respectively
    # the stopping rewards are 1, 0, 1, giving x the score 2047/33791.
    add("x", 12 * bits + 12, 1, "odd")
    add("odd", 12 * bits + 11, 0, "even")
    add("even", 12 * bits + 10, 1, "odd")

    vertex_id = {name: i for i, name in enumerate(game)}
    print(f"parity {len(game) - 1};")
    for name, (priority, owner, successors) in game.items():
        edges = ",".join(str(vertex_id[target]) for target in successors)
        print(f'{vertex_id[name]} {priority} {owner} {edges} "{name}";')


def write_energy_counter(length):
    """The first supported dominion appears after F_(length+1)+1 uncapped sweeps."""
    n = length + 6
    z, b, c, t, d, a = 0, 1, length + 2, length + 3, length + 4, length + 5
    rows = {
        z: (length + 5, 0, [c]),
        b: (length + 4, 1, [z, b, 2, t]),
        c: (2, 0, [b, a]),
        t: (1, 0, [t, z]),
        d: (1, 0, [d, c]),
        a: (0, 1, [a, d]),
    }
    for i in range(1, length + 1):
        v = length + 2 - i
        successors = [v, c, v - 1]
        if i > 1:
            successors.append(v + 1)
        rows[v] = (i + 2, 1 - i % 2, successors)
    print(f"parity {n - 1};")
    for v in range(n):
        priority, owner, successors = rows[v]
        print(f"{v} {priority} {owner} {','.join(map(str, successors))};")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("size", type=int, help="number of bits or energy-chain vertices")
    parser.add_argument("--energy", action="store_true", help="generate the slow energy-lift family")
    args = parser.parse_args()
    if args.energy:
        if args.size < 8 or args.size % 2:
            parser.error("the energy counter needs an even chain length of at least 8")
        write_energy_counter(args.size)
    else:
        if args.size < 1:
            parser.error("the counter needs at least one bit")
        write_counter(args.size)
