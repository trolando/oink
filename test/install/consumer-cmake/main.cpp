// Exercises the installed Oink public API and version header.
#include "oink/game.hpp"
#include "oink/io.hpp"
#include "oink/version.hpp"

#include <iostream>

int main()
{
    pg::Game g(3, 3);
    std::cout << "oink " << pg::version() << " — " << g.vertexcount()
              << " vertices, compressed input: " << pg::supported_input_formats()
              << "\n";
    return g.vertexcount() == 3 ? 0 : 1;
}
