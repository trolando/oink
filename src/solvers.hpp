#include <iostream>
#include <functional>
#include <vector>

#ifndef SOLVERS_HPP
#define SOLVERS_HPP

namespace pg {

class Oink;
class Game;
class Solver;

class Solvers
{
public:
    Solvers();

    /**
     * Get number of solvers
     */
    unsigned count() { return labels.size(); }

    /**
     * Get label of solver <id>
     */
    std::string label(int id) { return labels[id]; }

    /**
     * Get description of solver <id>
     */
    std::string desc(int id) { return descriptions[id]; }

    /**
     * Get whether solver <id> can be run in parallel
     */
    bool isParallel(int id) { return ispar[id]; }

    /**
     * Construct solver <id> with the given parameters
     */
    Solver* construct(int id, Oink* oink, Game* game, std::ostream& lgr) { return constructors[id](oink, game, lgr); }

    /**
     * Obtain the id matching the given solver label
     */
    int id(std::string label);

    /**
     * Write a formatted list of all solvers to the given ostream
     */
    void list(std::ostream &out);

protected:
    std::vector<std::string> labels;
    std::vector<std::string> descriptions;
    std::vector<bool> ispar;
    std::vector<std::function<Solver*(Oink*, Game*, std::ostream&)>> constructors;
    void add(std::string, std::string, int, std::function<Solver*(Oink*, Game*, std::ostream&)>);
};

}

#endif 
