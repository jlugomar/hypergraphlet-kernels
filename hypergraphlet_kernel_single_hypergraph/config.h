/**
 *
 * Jose Lugo-Martinez, jlugomar@andrew.cmu.edu
 * Computational Biology Department
 * School of Computer Science
 * Carnegie Mellon University
 *
 * Sep-23-2013
 *
 * Last modified: Nov-01-2018
 *
 */

#ifndef __CONFIG_H__
#define __CONFIG_H__

#define HYPERGRAPHLET_SIZE 15
#define NODES_MAX_SIZE 4
#define HYPEREDGES_MAX_SIZE 11
#define HYPERGRAPHLETS_TYPES 472        // Maximum number of n-hypergraphlets up to N=4
#define HYPERGRAPHLETS_TYPES_SIZE 511
#define LOG_HYPERGRAPHLETS_TYPES_SIZE 9
#define NODES_ALPHABET_SIZE 31          // Note: User-defined vertices alphabet size.
#define LOG_NODES_ALPHABET_SIZE 5       // Note: Depends on NODES_ALPHABET_SIZE.
#define HYPEREDGES_ALPHABET_SIZE 31     // Note: User-defined hyperedges alphabet size.
#define LOG_HYPEREDGES_ALPHABET_SIZE 5  // Note: Depends on HYPEREDGES_ALPHABET_SIZE.
#define ZERO_CHAR '@'
#define NODES_SIMILARITY_THRESHOLD 1.0
#define HYPEREDGES_SIMILARITY_THRESHOLD 1.0
#define HYPERGRAPHLET_SIMILARITY_THRESHOLD 1.0

// Flags for setting best hypergraphlets combination.
#define HYPERGRAPHLETS_1 1  // If 1 count 1-hypergraphlets, otherwise ignore 1-hypergraphlets.
#define HYPERGRAPHLETS_2 1  // If 1 count 2-hypergraphlets, otherwise ignore 2-hypergraphlets.
#define HYPERGRAPHLETS_3 1  // If 1 count 3-hypergraphlets, otherwise ignore 3-hypergraphlets.
#define HYPERGRAPHLETS_4 1  // If 1 count 4-hypergraphlets, otherwise ignore 4-hypergraphlets.

#define ENABLE_INDUCED_HYPERGRAPHLETS 1	// If 1 enable counting of partial decompositions of large
                                        // hyperedges via induced subhypergraph definition (C. Berge, 1980),
                                        // otherwise ignore partial decomposition of hyperedges.

#define DISABLE_MULTIPLICITY 0	// If 1 restrict induced hypergraphlets from having equivalent hyperedges with same label.

#define OUTPUT_FORMAT 2 // If 0 print triangular kernel matrix in binary format (for efficient SVM^light);
                        // If 1 print triangular kernel matrix to standard output;
                        // Otherwise print square kernel matrix to standard output.

#include <stdint.h>
#include <fstream>
#include <utility>
#include <list>
#include <map>
#include <vector>
using namespace std;


typedef uint64_t Element;
typedef pair <Element, Element> Key;


// Enhanced data structure to accommodate for label substitutions and hyperedge insertions and deletions.
struct MismatchInfo
{
    float matches;                          // Contains counts for exact non-isomorphic vertex- and hyperedge-labeled hypergraphlets found on a given hypergraph.
    float mismatches;                       // Contains counts for inexact non-isomorphic labeled hypergraphlets obtained by hyperedge indels algorithm.
    map<Key, float> mismatchesHypergraph;   // Contains counts for each inexact non-isomorphic labeled hypergraphlet generated by label
                                            // substitutions algorithm and satisfying HYPERGRAPHLET_SIMILARITY_THRESHOLD (defined above).
};

#endif
