/**
 *
 * Jose Lugo-Martinez, jlugomar@andrew.cmu.edu
 * Computational Biology Department
 * School of Computer Science
 * Carnegie Mellon University
 *
 * Sep-23-2013
 *
 */

#ifndef UTILS_H
#define	UTILS_H

#include "config.h"
#include <vector>
#include <list>
using namespace std;


float compare_labels(char label1, char label2);

int randint(int max);

double randdouble();

list<unsigned> get_incident_edges(vector<vector<unsigned> > vertex_set, unsigned vertex);

unsigned long get_hypergraphlet_length(unsigned long hg_type);

unsigned long get_hyperedges_number(unsigned long hg_type);

void compare_two(char &a, char &b, char &a1, char &b1);

void compare_three(char &a, char &b, char &c, char &a1, char &b1, char &c1);

unsigned int set_k(unsigned long hg_type, float sf);

void insert_permutation(const Key &target, vector<Key> &mismatches_list);

void insert_hyperedge_label(const char &hyperedge_label, vector<char> &hyperedge_labels);

bool search_hyperedge_label(const char &hyperedge_label, vector<char> &hyperedge_labels);

string get_nkey(Element k);

string get_ekey(Element k);

string print_nkey(Element k);

string print_ekey(Element k);

void make_key(Key &k, const Element nk, const Element ek);

Element make_nodes_key(char root, char a, char b, char c, unsigned long hg_type);

Element make_hyperedges_key(char e1, char e2, char e3, char e4, char e5, char e6, char e7, char e8, char e9, char e10, char e11, unsigned long hg_type);

void initialize_vertices_labels(Element nkey, char &root, char &a, char &b, char &c);

void initialize_edges_labels(Element ekey, char &e1, char &e2, char &e3, char &e4, char &e5, char &e6, char &e7, char &e8, char &e9, char &e10, char &e11);

Element get_feature_id_nodes(Element nkey, unsigned long hg_type);

Element get_feature_id_hyperedges(Element ekey, unsigned long hg_type);

Key get_feature_id(Element nkey, Element ekey, unsigned long hg_type);

void increment_match_hash(map<Key,MismatchInfo> &hash, const Key &k, vector<Key> &mismatches_list);

void increment_match_hash(map<Key,MismatchInfo> &hash, const vector<Key> &keys, vector<vector<Key> > &mismatches_list);

float retrieve_exact_matches_count(map<Key,MismatchInfo> &hash, const Key &k);

float retrieve_edge_mismatch_count(map<Key,MismatchInfo> &hash, const Key &k);

float retrieve_label_mismatch_count(unsigned long hg_type, map<Key,MismatchInfo> &hash, const Key key);

void insert_mismatch_counts(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, const Key &k);

void insert_mismatches_hash(map<Key,MismatchInfo> &mismatch_hash, const Key &k, vector<Key> mismatches_list);

Key create_permutations_subset(vector<Key> &mismatches, char root, char a, char b, char c, unsigned long hg_type, char e1, char e2, char e3, char e4, char e5, char e6, char e7, char e8, char e9, char e10, char e11);

vector<Key> generate_labels(vector<vector<Key> > &mismatches, char root, char a, char b, char c, unsigned long hg_type, vector<char> e1s, vector<char> e2s, vector<char> e3s, vector<char> e4s, vector<char> e5s, vector<char> e6s, vector<char> e7s, vector<char> e8s, vector<char> e9s, vector<char> e10s, vector<char> e11s);

void insert_hypergraphlet_mismatch_neighborhood(list<Key> &neighborhood, Key k);

void generate_hypergraphlet_mismatch_neighborhood_m1(list<Key> &neighborhood, string ALPHABET_ROOT, string ALPHABET, map<string, float> &sim_matrix, unsigned long hg_type, Element nkey, Element ekey);

void generate_hypergraphlet_mismatch_neighborhood_m2(list<Key> &neighborhood, string ALPHABET_ROOT, string ALPHABET, map<string, float> &sim_matrix, unsigned long hg_type, Element nkey, Element ekey);

void generate_graphlet_mismatch_neighborhood_for_edges_m1(list<Key> &neighborhood, string EDGES_ALPHABET, map<string, float> &sim_elm_matrix, unsigned long hg_type, Element nkey, Element ekey);

void generate_graphlet_mismatch_neighborhood_for_edges_m2(list<Key> &neighborhood, string EDGES_ALPHABET, map<string, float> &sim_elm_matrix, unsigned long hg_type, Element nkey, Element ekey);

void generate_vertex_label_mismatch_hypergraphlets(map<Key, list<Key> > &vl_mismatch_neighborhood, map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, Key key, unsigned long hg_type, string ALPHABET_ROOT, string ALPHABET, map<string, float> &sim_matrix, int VLM);

void generate_edge_label_mismatch_hypergraphlets(map<Key, list<Key> > &el_mismatch_neighborhood, map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, Key key, unsigned long hg_type, string EDGES_ALPHABET, map<string, float> &sim_elm_matrix, int ELM);

void increment_mismatch_count(map<Key,MismatchInfo> &hash, const Key &key, const Key &mismatch_key, float sim_score, float mult_factor);

float compare_hypergraphlets_vertex_labels(Element nkey1, Element nkey2, unsigned long hg_type, map<string, float> &sim_vlm_matrix, float &sim_vlm_score);

float compare_hypergraphlets_edge_labels(Element ekey1, Element ekey2, unsigned long hg_type, map<string, float> &sim_elm_matrix, float &sim_elm_score);

void update_mismatch_count(map<Key,MismatchInfo> &hash, Key k, float mult_factor, unsigned long hg_type, map<string, float> sim_vlm_matrix, map<string, float> sim_elm_matrix, int VLM, int ELM, bool eq);

void increment_edge_mismatch_hash(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, const Key &k, float mult_factor, vector<Key> &mismatches_list);

void insert_edge_mismatch_hypergraphlet(list<pair <unsigned long, Key> > &EM_set, pair <unsigned long, Key> graphlet, unsigned index);

void update_edge_mismatch_count(vector<list<pair <unsigned long, Key> > > &EM_set, Key k, unsigned long hg_type, unsigned EDGE_MISMATCHES_ALLOWED, unsigned vindex, string EDGES_ALPHABET);

#endif
