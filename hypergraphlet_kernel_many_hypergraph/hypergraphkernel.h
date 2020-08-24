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

#ifndef __HYPERGRAPHKERNELS_H__
#define __HYPERGRAPHKERNELS_H__

#include "utils.h"
#include "simhypergraph.h"
#include <fstream>
#include <utility>
#include <list>
#include <map>
#include <vector>
#include <set>
#include <sys/time.h>
using namespace std;


class HypergraphKernels  {
public:
    HypergraphKernels() : NORMALIZE(false), DUALITY(false), VERBOSE(false), SF(0.0), ELM(0), EM(0) {}
    ~HypergraphKernels()  {}
     
    /** Read vertex- and edge-labeled hypergraphs (or dual hypergraphs) and a list of vertex of interest. */ 
    void read_hypergraph(string nl_path, string el_path, string hg_path, const vector<string> &hg_filenames, const vector<unsigned> &vertices_of_interest);

    /** Read a probability similarity matrix for each vertex label so that one can weight each vertex label substitution. */
    void read_sim_vlm_matrix(string filename);

    /** Read a probability similarity matrix for each edge label so that one can weight each edge label substitution. */
    void read_sim_elm_matrix(string filename);

    /** Stores the class labels for each example, if provided. */
    inline void set_labels(const vector<int> &l)  { labels = l; }

    /** Read and compute cumulative random walk hypergraph kernel matrix. 
     *  This is a memory-efficient but slower version so that (dual) hypergraphs are not stored in main memory. */
    void read_and_compute_random_walk_cumulative_matrix(string nl_path, string el_path, string hg_path, const vector<string> &hg_filenames, const vector<unsigned> &vertices_of_interest, int steps, double restart);

    /** Compute cumulative random walk hypergraph kernel matrix. */
    void compute_random_walk_cumulative_matrix(int steps, double restart);

    /** Read and compute random walk hypergraph kernel matrix.
      *  This is a memory-efficient but slower version so that (dual) hypergraphs are not stored in main memory. */
    void read_and_compute_random_walk_matrix(string nl_path, string el_path, string hg_path, const vector<string> &hg_filenames, const vector<unsigned> &vertices_of_interest, int steps, double restart);

    /** Compute random walk hypergraph kernel matrix. */
    void compute_random_walk_matrix(int steps, double restart);

    /** Read and compute label substitutions hypergraph kernel matrix.
      *  This is a memory-efficient but slower version so that (dual) hypergraphs are not stored in main memory. */
    void read_and_compute_label_mismatch_matrix(string nl_path, string el_path, string hg_path, const vector<string> &hg_filenames, const vector<unsigned> &vertices_of_interest);

    /** Compute label substitutions hypergraph kernel matrix. */
    void compute_label_mismatch_matrix();

    /** Read and compute edge indels hypergraph kernel matrix.
      *  This is a memory-efficient but slower version so that (dual) hypergraphs are not stored in main memory. */
    void read_and_compute_edge_mismatch_matrix(string nl_path, string el_path, string hg_path, const vector<string> &hg_filenames, const vector<unsigned> &vertices_of_interest);

    /** Compute edge indels hypergraph kernel matrix. */
    void compute_edge_mismatch_matrix();

    /** Read and compute edit distance (with 1 operation) hypergraph kernel matrix.
      *  This is a memory-efficient but slower version so that (dual) hypergraphs are not stored in main memory. */
    void read_and_compute_edit_distance_matrix(string nl_path, string el_path, string hg_path, const vector<string> &hg_filenames, const vector<unsigned> &vertices_of_interest);

    /** Compute edit distance (with 1 operation) hypergraph kernel matrix. */
    void compute_edit_distance_matrix();

    /** Read and compute edit distance (with 2 operations) hypergraph kernel matrix.
      *  This is a memory-efficient but slower version so that (dual) hypergraphs are not stored in main memory. */
    void read_and_compute_edit_distance2_matrix(string nl_path, string el_path, string hg_path, const vector<string> &hg_filenames, const vector<unsigned> &vertices_of_interest);

    /** Compute edit distance (with 2 operations) hypergraph kernel matrix. */
    void compute_edit_distance2_matrix();

    /** Writes distance kernel matrix in either: 1:triangular binary form, 2:triangular standard output form, 3: squared-matrix standard output form. */
    void write_matrix(const char *file);

    /**  */
    void write_sparse_svml_lm(string nl_path, string el_path, string hg_path, const vector<string> &hg_filenames, const vector<unsigned> &vertices_of_interest, const char *file);

    /** */
    void write_sparse_svml_em(string nl_path, string el_path, string hg_path, const vector<string> &hg_filenames, const vector<unsigned> &vertices_of_interest, const char *file);

    /** */
    void write_sparse_svml_ed(string nl_path, string el_path, string hg_path, const vector<string> &hg_filenames, const vector<unsigned> &vertices_of_interest, const char *file);

    /** */
    void write_sparse_svml_ed2(string nl_path, string el_path, string hg_path, const vector<string> &hg_filenames, const vector<unsigned> &vertices_of_interest, const char *file);

    /** */
    void write_labels(const char *file);

    inline void set_normalize()  { NORMALIZE = true; } 

    inline void set_duality()  { DUALITY = true; } 

    inline void set_verbose()  { VERBOSE = true; }

    inline void set_number_vertex_label_mismatches(float fraction)  { SF = fraction; }
    
    inline void set_vertex_label_mismatches_alphabet(string alphabet)  { ALPHABET = alphabet; }

    inline void set_vertex_label_mismatches_root_alphabet(string alphabet)  { ALPHABET_ROOT = alphabet; }

    inline void set_edge_label_mismatches_alphabet(string alphabet)  { EDGES_ALPHABET = alphabet; }

    inline void set_number_edge_label_mismatches(int elm)  { ELM = elm; }

	inline void set_number_edges_mismatches(unsigned edges_mismatches)  { EM = edges_mismatches; }
   
private:
	/** Returns the cumulative random walk hypergraph kernel between two rooted neighborhoods. */
    float random_walk_cumulative(SimpleHypergraph &hg1, SimpleHypergraph &hg2, unsigned hg1_root, unsigned hg2_root, int steps, double restart);

	/** Returns the random walk hypergraph kernel between two rooted neighborhoods. */
    float random_walk(SimpleHypergraph &hg1, SimpleHypergraph &hg2, unsigned hg1_root, unsigned hg2_root, int steps, double restart);

	/** Insert the counts of induced labeled 2-hypergraphlets. */
	void insert_induced_hypergraphlets(vector<map<Key,MismatchInfo> > &hash, char root, char a, vector<char> e1, bool found_ri) ;

	/** Insert the counts of induced labeled 3-hypergraphlets. */
	void insert_induced_hypergraphlets(vector<map<Key,MismatchInfo> > &hash, char root, char a, char b, vector<char> e1, vector<char> e2, vector<char> e3, vector<char> e4, bool found_rij, bool found_ri, bool found_rj, bool found_ij);

	/** Insert the counts of induced labeled 4-hypergraphlets. */
	void insert_induced_hypergraphlets(vector<map<Key,MismatchInfo> > &hash, char root, char a, char b, char c, vector<char> e1, vector<char> e2, vector<char> e3, vector<char> e4, vector<char> e5, vector<char> e6, vector<char> e7, vector<char> e8, vector<char> e9, vector<char> e10, vector<char> e11, bool found_rijk, bool found_rij, bool found_rik, bool found_rjk, bool found_ijk, bool found_ri, bool found_rj, bool found_rk, bool found_ij, bool found_ik, bool found_jk);

	/** Returns the counts of induced labeled hypergraphlets on a rooted neighborhood. */
	vector<map<Key,MismatchInfo> > get_induced_hypergraphlets_counts(SimpleHypergraph &hg, unsigned hg_root);
    
    /** Returns the counts of labeled hypergraphlets on a rooted neighborhood. */
	vector<map<Key,MismatchInfo> > get_hypergraphlets_counts(SimpleHypergraph &hg, unsigned hg_root);

    /** Adds the counts for inexact hypergraphlets based on vertex label substitutions. */
    void add_vertex_label_mismatch_counts(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, unsigned long hg_type, int VLM, bool option);

    /** Adds the counts for inexact hypergraphlets based on edge label substitutions. */
    void add_edge_label_mismatch_counts(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, unsigned long hg_type, int ELM);

    /** Updates the vector of counts with inexact hypergraphlets based on vertex and edge label substitutions. **/
    void update_label_mismatch_counts(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, unsigned long hg_type, bool option, int VLM, int ELM, bool eq);

    /** Adds the counts for inexact hypergraphlets based on edge insertions and deletions. */
    void add_edge_mismatch_counts(vector<map<Key,MismatchInfo> > &hash);

    /** Adds the counts for inexact hypergraphlets based on 1-edge insertion and deletion. */
    void add_1_edge_mismatch_counts(vector<map<Key,MismatchInfo> > &hash);

    /** Adds the counts for inexact hypergraphlets based on 2-edge insertions and deletions. */
    void add_2_edge_mismatch_counts(vector<map<Key,MismatchInfo> > &hash);

    /** Normalizes the kernel matrix using the method for normalizing spectral kernel matrices. */
    void normalize_spectral(map<Key,MismatchInfo> &, unsigned long hg_type);

    /** Computes hypergraphlet distance between two vector of hypergraphlet counts. */
    float distance_hash_join(map<Key,MismatchInfo>, map<Key,MismatchInfo>, unsigned long hg_type);

    /** Normalizes the kernel matrix using the method for normalizing the spectral kernel matrix. */
    void normalize_spectral(vector<map<Key,MismatchInfo> >&);
     
    /** Computes hypergraphlet distance between two vector of hypergraphlet counts. */
    float distance_hash_join(vector<map<Key,MismatchInfo> >, vector<map<Key,MismatchInfo> >);

    // Data members.
    bool NORMALIZE, DUALITY, VERBOSE;
    float SF;
	int ELM;
    unsigned EM;
    string ALPHABET;
    string ALPHABET_ROOT;
    string EDGES_ALPHABET;

    vector<int>                 labels;
    vector<SimpleHypergraph>    hypergraphs;
    vector<unsigned>            roots;  // Vertices of interest.
    map<string,float>           sim_vlm_matrix;    
	map<string,float>           sim_elm_matrix;    
    vector<vector<map<Key,MismatchInfo> > > hashes;
    vector<vector<float> >  hyperkernel;
    map<Key, list<Key> >    vl_mismatch_neighborhood;
    map<Key, list<Key> >    el_mismatch_neighborhood;
};

#endif

