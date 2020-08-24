/**
 * Simple Hypergraph - auxiliary data structure. 
 *
 * Jose Lugo-Martinez, jlugomar@andrew.cmu.edu
 * Computational Biology Department
 * School of Computer Science
 * Carnegie Mellon University
 *
 * Sep-23-2013
 *
 */


#ifndef __SIMPLE_HYPERGRAPH_H__
#define __SIMPLE_HYPERGRAPH_H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <iostream>
#include <vector>
#include <list>
#include <set>
using namespace std;


class SimpleHypergraph  {
public:
    SimpleHypergraph()  {}
    ~SimpleHypergraph()  {}

    /** Return all hyperedges that contain specified vertex. */
    vector<unsigned> get_incident_edges(unsigned);

    /** Read a target file and create a hypergraph. */
    static SimpleHypergraph read_hypergraph(const char*, const char*, const char*);

    /** Read a target file and create a dual hypergraph from input hypergraph. */
    static SimpleHypergraph read_dual_hypergraph(const char*, const char*, const char*);

    /** Prints a hypergraph into corresponding user-speified files. */
    void print_hypergraph(const char*, const char*, const char*);

    /** Transform an input hypergraph into a standard graph using clique expansion of a hypergraph. */
    void transform_clique_expansion(const char*, const char*);

    /** Transform an input hypergraph into a standard graph using star expansion of a hypergraph. */
    void transform_star_expansion(const char*, const char*);

    /** Breadth-first assignment of distances from a source node. */
    vector<unsigned> breadth_first_sort(unsigned) const;

    /** Gets neighboring nodes for each node in a hypergraph. */
    vector<set<unsigned> > get_neighbors() const;


    string node_labels;                     // Vertex (or node) labels.
	string edge_labels;					    // Hyperedge labels.
    vector<vector<unsigned> > vertex_set;   // Hyperedges with vertex list.
    vector<vector<unsigned> > edge_set;     // Vertices with hyperedge list.
};

#endif
