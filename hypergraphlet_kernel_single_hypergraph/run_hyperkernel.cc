/**
 * Generates a kernel matrix according to a
 * user-specified hypergraph-based kernel 
 * method over a set of vertices of interest
 * within an input undirected vertex- and 
 * edge-labeled hypergraph.
 *
 * Modified from Vladimir Vacic graphlet kernel
 *
 * Jose Lugo-Martinez, jlugomar@andrew.cmu.edu
 * Computational Biology Department
 * School of Computer Science
 * Carnegie Mellon University
 *
 * Sep-23-2013
 *
 * Copyright (c) 2014 Jose Lugo-Martinez,
 * Vladimir Vacic, and Predrag Radivojac.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "hyperkernel.h"
#include "string.h"
#include <iostream>
#include <fstream>
using namespace std;


void print_help()  {
    cout << "Usage: run_hyperkernel -p FILE -n FILE -g HG_FILE -l L_FILE -e E_FILE -t TYPE -z TASK -[k|s] OUTPUT [...]\n";
    cout << "Options:\n\n";

    cout << "  -h         Displays this message.\n\n";

    cout << "  -t TYPE    Kernel type (0-Cumulative Random Walk, 1-Standard Random Walk, 2-Standard Hypergraphlet, 3-Label Substitutions Hypergraphlet, 4-Edge Indels Hypergraphlet, 5-Edit Distance Hypergraphlet).\n";
    cout << "             Defaults to standard hypergraphlet.\n\n";

    cout << "  -z TASK    Task type (0-Vertex classification, 1-Edge classification)\n";
    cout << "             Defaults to vertex classification.\n\n";

    cout << "  -p FILE    List of positive (vertices) examples.\n";
    cout << "  -n FILE    List of negative (vertices) examples.\n";
    cout << "  -g HG_FILE Hypergraph file without file extension.\n";
    cout << "  -l L_FILE  Vertex labels file for input hypergraph without file extension.\n\n";
    cout << "  -e E_FILE  Edge labels file for input hypergraph without file extension.\n\n";

    cout << "  -N         Normalize the kernel matrix.\n";
    cout << "             Defaults to false.\n\n";

    cout << "  -k KERNEL  Output file for the kernel matrix in standard output.\n";
    cout << "   or\n";
    cout << "  -s SPARSE  Output file for the sparse attribute matrix (SVML).\n";
    cout << "             Defaults to KERNEL.\n\n";

    cout << "  -I STEPS   Number of steps. (Needed for Random Walk Hyperkernels)\n";
    cout << "             Defaults to 100,000 steps.\n\n";

    cout << "  -R RESTART Restart probability. (Needed for Random Walk Hyperkernels)\n\n";
    cout << "             Defaults to 0.3\n\n";

    cout << "  -S SIMVMAT Similarity matrix for weighting each possible vertex label substitution. (Needed for Label Substitutions and Edit Distance Hyperkernels)\n";
    cout << "             Defaults to uniform weights (i.e. all label mismatches are equally weighted to 1).\n\n";

    cout << "  -P SIMEMAT Similarity matrix for weighting each possible edge label substitution. (Needed for Label Substitutions, Edge Indels and Edit Distance Hyperkernels)\n";
    cout << "             Defaults to uniform weights (i.e. all label mismatches are equally weighted to 1).\n\n";

    cout << "  -K VLABMIS Fraction of nodes in the n-hypergraphlet allowed to have vertex label mismatches. (Needed for Label Substitutions and Edit Distance Hyperkernels)\n";
    cout << "             For example, K=0.34 allows 1-vertex label mismatch for 3- and 4-hypergraphlets, whereas K=0.26 allows 1-vertex label mismatch for 4-hypergraphlets only.\n";
    cout << "             Defaults to 0.0 (i.e. no vertex label mismatches allowed).\n\n";

    cout << "  -M ELABMIS Total number of edge label mismacthes allowed between graphlets. (Needed for Label Substitutions, Edge Indels and Edit Distance Hyperkernels)\n";
    cout << "             Defaults to 0.\n\n";

    cout << "  -m EDGMIS  Total number of edge insertions or deletions allowed between hypergraphlets. (Needed for Edge Indels and Edit Distance Hyperkernels)\n";
    cout << "             Defaults to 0.\n\n";

    cout << "  -V ALPHA   Vertex labels alphabet over problem statement is defined (i.e. all possible labels for a vertex). (Needed for Label Substitutions and Edit Distance Hyperkernels)\n\n";

    cout << "  -E ALPHA   Edge labels alphabet over problem statement is defined (i.e. all possible labels for an edge). (Needed for Label Substitutions, Edge Indels and Edit Distance Hyperkernels)\n\n";

    cout << "  -c LABELS  Output file for each example class label.\n\n";   
   
    cout << "  -v         Verbose (prints progress messages).\n\n"; 
}

int main(int argc, char* argv[])  {
    typedef enum kerneltype  {
        RANDOM_WALK_CUMULATIVE,
        RANDOM_WALK,
        STANDARD_GRAPHLET,
        LABEL_MISMATCH,
        EDGE_MISMATCH,
        EDIT_DISTANCE
     } KernelType;

    typedef enum tasktype  {
        VERTEX_CLASSIFICATION,
        EDGE_CLASSIFICATION
    } TaskType;

    typedef enum outformat  {
        KERNEL,
        SPARSE_SVML
     } OutputFormat;

    string pos_file;
    string neg_file;
    string hg_path, l_path, e_path;
    string output_file;
    OutputFormat format(KERNEL);
    KernelType hk_type(STANDARD_GRAPHLET);
    TaskType task_type(VERTEX_CLASSIFICATION);
    string labels_file;
    bool normalize(false);
    bool verbose(false);

    // Random Walk Hyperkernel Parameters
    int steps(100000); 
    double restart(0.3);

    // Label Substitutions and Edit Distance Hyperkernels Parameters
    float vlabel_mismatches(0.0);
    int elabel_mismatches(1);
    string vertices_alphabet;
    string root_alphabet;
    string edges_alphabet;
    string sim_vlm_matrix_file, sim_elm_matrix_file;

    // Edge Indels and and Edit Distance Hyperkernels Parameter
    unsigned edgmis(1);


    // Parse command line arguments.
    for (int i=1; i<argc && (argv[i])[0] == '-'; i++)  {
        switch ((argv[i])[1])  {
            case 'h': print_help(); exit(0);
            case 't': 
                i++; 
                switch (to_i(argv[i]))  {
                    case 0:
                        hk_type=RANDOM_WALK_CUMULATIVE;
                        break;
                    case 1:
                        hk_type=RANDOM_WALK;
                        break;
                    case 2:
                        hk_type=STANDARD_GRAPHLET;
                        break;
                    case 3:
                        hk_type=LABEL_MISMATCH;
                        break;
                    case 4:
                        hk_type=EDGE_MISMATCH;
                        break;
                    case 5:
                        hk_type=EDIT_DISTANCE;
                        break;
                    default:
                        hk_type=STANDARD_GRAPHLET;
                }
                break;
            case 'z': 
                i++;
                switch (to_i(argv[i]))  {
                    case 0:
                        task_type=VERTEX_CLASSIFICATION;
                        break;
                    case 1:
                        task_type=EDGE_CLASSIFICATION;
                        break;
                    default:
                        task_type=VERTEX_CLASSIFICATION;
                }
                break;
            case 'p': i++; pos_file=argv[i]; break;
            case 'n': i++; neg_file=argv[i]; break;
            case 'g': i++; hg_path=argv[i]; break;
            case 'l': i++; l_path=argv[i]; break;
			case 'e': i++; e_path=argv[i]; break;
            case 'N': normalize=true; break;            
            case 'k': i++; format=KERNEL; output_file=argv[i]; break;
            case 's': i++; format=SPARSE_SVML; output_file=argv[i]; break;
            // Hyperkernel-specific parameters                    
            case 'I': i++; steps=to_i(argv[i]); break;
            case 'R': i++; restart=to_f(argv[i]); break;

            case 'S': i++; sim_vlm_matrix_file=argv[i]; break;
            case 'P': i++; sim_elm_matrix_file=argv[i]; break;
            case 'K': 
                i++; 
                vlabel_mismatches=to_f(argv[i]);
                if(vlabel_mismatches < 0.0 || vlabel_mismatches > 1.0)  { 
                    cerr << "ERROR: Fraction of nodes allowed to have vertex label mismatches must be 0<=K<=1, but you entered " << argv[i] << endl;
                    print_help();  exit(1);
                }
                break;
            case 'M':
                i++;
                elabel_mismatches=to_i(argv[i]);
                if(elabel_mismatches < 0 || elabel_mismatches > 2)  {
                    cerr << "ERROR: Total number of edge label mismatches must be either 1 or 2, but you entered " << argv[i] << endl;
                    print_help();  exit(1);
                }
                break;
            case 'm':
                i++;
                edgmis=to_i(argv[i]);
                if(edgmis > 2)  {
                    cerr << "ERROR: Total number of edge insertions and deletions must be either 1 or 2, but you entered " << argv[i] << endl;
                    print_help();  exit(1);
                }
                break;
            case 'V': i++; vertices_alphabet=argv[i]; root_alphabet=argv[i]; break;
            case 'E': i++; edges_alphabet=argv[i]; break;
            case 'c': i++; labels_file=argv[i]; break;
            case 'v': verbose=true; break;
            default: 
                cerr << "ERROR: Unknown option " << argv[i] << endl;
                print_help();  exit(1);
        }
    }

    if (0 == output_file.size())  {
        cerr << "ERROR: Output file name not specified." << endl;  print_help();  exit(1);
    }

    if (0 == sim_vlm_matrix_file.size() && (hk_type == LABEL_MISMATCH || hk_type == EDIT_DISTANCE))  {
        //Use default node labels similarity matrix
        sim_vlm_matrix_file = "user_defined.matrix";
    }

    if (0 == sim_elm_matrix_file.size() && (hk_type == LABEL_MISMATCH || hk_type == EDGE_MISMATCH || hk_type == EDIT_DISTANCE))  {
        //Use default hyperedge labels similarity matrix
        sim_elm_matrix_file = "user_defined.matrix";
    }

    if ((hk_type == LABEL_MISMATCH || hk_type == EDIT_DISTANCE) && vlabel_mismatches > 0.0 && (0 == vertices_alphabet.size() || 0 == root_alphabet.size()))  {
        cerr << "ERROR: Alphabet for the vertex labels not specified. It is required for selected hyperkernel type." << endl;  print_help();  exit(1);
    }

    HyperKernels hk;

    string line;
    vector<unsigned> roots;
    vector<int> labels;
    
    // Read lists of hypergraphs.
    ifstream p(pos_file.c_str(), ios::in);
    if (p.fail()) {
        cerr << "WARNING: Positive file " << pos_file << " cannot be opened." << endl;
    }
    else  {
        while(getline(p, line))  {
            vector<string> tokens = split(line, '\t');            
            roots.push_back(to_i(strip(tokens[0])));
            labels.push_back(1);
        }
    }
    p.close();

    ifstream n(neg_file.c_str(), ios::in);
    if (n.fail())  {
        cerr << "WARNING: Negative file " << neg_file << " cannot be opened." << endl;
    }
    else  {
        while(getline(n, line))  {
            vector<string> tokens = split(line, '\t');
            roots.push_back(to_i(strip(tokens[0])));
            labels.push_back(-1);
        }
    }
    n.close();

    if (roots.size() < 2)  {
        cerr << "ERROR: Too few vertices of interest." << endl << endl; print_help(); exit(1);
    }

    if (normalize)  hk.set_normalize();
    if (verbose)  hk.set_verbose();

    switch (task_type)  {
        case VERTEX_CLASSIFICATION:
            hk.read_hypergraph(l_path, e_path, hg_path, roots);
            break;
        case EDGE_CLASSIFICATION:
            hk.read_dual_hypergraph(l_path, e_path, hg_path, roots);
            break;
    }
    hk.set_labels(labels);
    
    switch (format)  {
        case KERNEL:
			switch (hk_type)  {
                case RANDOM_WALK_CUMULATIVE:
                    hk.compute_random_walk_cumulative_matrix(steps, restart);
                    break;
				case RANDOM_WALK:
					hk.compute_random_walk_matrix(steps, restart);
					break;
				case STANDARD_GRAPHLET:
					hk.set_number_vertex_label_mismatches(0.0);                    
					hk.compute_label_mismatch_matrix();
					break;
				case LABEL_MISMATCH:
					hk.set_number_vertex_label_mismatches(vlabel_mismatches);
                    hk.set_number_edge_label_mismatches(elabel_mismatches);
                    switch (task_type)  {
                        case VERTEX_CLASSIFICATION:
                            hk.set_vertex_label_mismatches_alphabet(vertices_alphabet);
                            hk.set_vertex_label_mismatches_root_alphabet(root_alphabet);
                            hk.set_edge_label_mismatches_alphabet(edges_alphabet);
                            break;
                        case EDGE_CLASSIFICATION:
                            hk.set_vertex_label_mismatches_alphabet(edges_alphabet);
                            hk.set_vertex_label_mismatches_root_alphabet(edges_alphabet);
                            hk.set_edge_label_mismatches_alphabet(vertices_alphabet);
                            break;
                    }
                    hk.read_sim_vlm_matrix(sim_vlm_matrix_file);
                    hk.read_sim_elm_matrix(sim_elm_matrix_file);
					hk.compute_label_mismatch_matrix();
                    break;
                case EDGE_MISMATCH:
					hk.set_number_edges_mismatches(edgmis);
                    switch (task_type)  {
                        case VERTEX_CLASSIFICATION:
                            hk.set_edge_label_mismatches_alphabet(edges_alphabet);
                            break;
                        case EDGE_CLASSIFICATION:
                            hk.set_edge_label_mismatches_alphabet(vertices_alphabet);
                            break;
                    }
                    hk.compute_edge_mismatch_matrix();
					break;
                case EDIT_DISTANCE:
                    hk.set_number_vertex_label_mismatches(vlabel_mismatches);
                    hk.set_number_edge_label_mismatches(elabel_mismatches);
                    switch (task_type)  {
                        case VERTEX_CLASSIFICATION:
                            hk.set_vertex_label_mismatches_alphabet(vertices_alphabet);
                            hk.set_vertex_label_mismatches_root_alphabet(root_alphabet);
                            hk.set_edge_label_mismatches_alphabet(edges_alphabet);
                            break;
                        case EDGE_CLASSIFICATION:
                            hk.set_vertex_label_mismatches_alphabet(edges_alphabet);
                            hk.set_vertex_label_mismatches_root_alphabet(edges_alphabet);
                            hk.set_edge_label_mismatches_alphabet(vertices_alphabet);
                            break;
                    }
                    hk.read_sim_vlm_matrix(sim_vlm_matrix_file);
                    hk.read_sim_elm_matrix(sim_elm_matrix_file);
                    hk.set_number_edges_mismatches(edgmis);
                    if (edgmis == 2)
                        hk.compute_edit_distance2_matrix();
                    else
                        hk.compute_edit_distance_matrix();                    
                    break;
			}
			hk.write_matrix(output_file.c_str());
            break;
        case SPARSE_SVML:
            switch (hk_type)  {
                case RANDOM_WALK_CUMULATIVE:
                    hk.compute_random_walk_cumulative_matrix(steps, restart);
                    hk.write_matrix(output_file.c_str());
                    break;
				case RANDOM_WALK:
					hk.compute_random_walk_matrix(steps, restart);
                    hk.write_matrix(output_file.c_str());
					break;
				case STANDARD_GRAPHLET:
					hk.set_number_vertex_label_mismatches(0.0);
					hk.write_sparse_svml_lm(output_file.c_str());
					break;
				case LABEL_MISMATCH:
					hk.set_number_vertex_label_mismatches(vlabel_mismatches);
                    hk.set_number_edge_label_mismatches(elabel_mismatches);
                    switch (task_type)  {
                        case VERTEX_CLASSIFICATION:
                            hk.set_vertex_label_mismatches_alphabet(vertices_alphabet);
                            hk.set_vertex_label_mismatches_root_alphabet(root_alphabet);
                            hk.set_edge_label_mismatches_alphabet(edges_alphabet);
                            break;
                        case EDGE_CLASSIFICATION:
                            hk.set_vertex_label_mismatches_alphabet(edges_alphabet);
                            hk.set_vertex_label_mismatches_root_alphabet(edges_alphabet);
                            hk.set_edge_label_mismatches_alphabet(vertices_alphabet);
                            break;
                    }
					hk.read_sim_vlm_matrix(sim_vlm_matrix_file);
                    hk.read_sim_elm_matrix(sim_elm_matrix_file);
					hk.write_sparse_svml_lm(output_file.c_str());
                    break;
                case EDGE_MISMATCH:
					hk.set_number_edges_mismatches(edgmis);
                    switch (task_type)  {
                        case VERTEX_CLASSIFICATION:
                            hk.set_edge_label_mismatches_alphabet(edges_alphabet);
                            break;
                        case EDGE_CLASSIFICATION:
                            hk.set_edge_label_mismatches_alphabet(vertices_alphabet);
                            break;
                    }
					hk.write_sparse_svml_em(output_file.c_str());
                    break;
                case EDIT_DISTANCE:
                    hk.set_number_vertex_label_mismatches(vlabel_mismatches);
                    hk.set_number_edge_label_mismatches(elabel_mismatches);
                    switch (task_type)  {
                        case VERTEX_CLASSIFICATION:
                            hk.set_vertex_label_mismatches_alphabet(vertices_alphabet);
                            hk.set_vertex_label_mismatches_root_alphabet(root_alphabet);
                            hk.set_edge_label_mismatches_alphabet(edges_alphabet);
                            break;
                        case EDGE_CLASSIFICATION:
                            hk.set_vertex_label_mismatches_alphabet(edges_alphabet);
                            hk.set_vertex_label_mismatches_root_alphabet(edges_alphabet);
                            hk.set_edge_label_mismatches_alphabet(vertices_alphabet);
                            break;
                    }
                    hk.read_sim_vlm_matrix(sim_vlm_matrix_file);
                    hk.read_sim_elm_matrix(sim_elm_matrix_file);
                    hk.set_number_edges_mismatches(edgmis);
                    if (edgmis == 2)
                        hk.write_sparse_svml_ed2(output_file.c_str());
                    else
                        hk.write_sparse_svml_ed(output_file.c_str()); 
            }
    }

    if (labels_file.size() > 0)
        hk.write_labels(labels_file.c_str());

    exit(0);
}

