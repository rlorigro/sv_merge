#include "Graph.hpp"


namespace sv_merge {
Node::Node(const string &name) :
        name(name) {}


Node::Node(string &name) :
        name(name) {}
}