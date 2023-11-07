#include "HeteroGraph.hpp"


namespace sv_merge{

HeteroNode::HeteroNode(const string& name):
        name(name),
        type('*')
{}


HeteroNode::HeteroNode(string& name):
        name(name),
        type('*')
{}


HeteroNode::HeteroNode(const string& name, char type):
        name(name),
        type(type)
{}


HeteroNode::HeteroNode(string& name, char type):
        name(name),
        type(type)
{}

}