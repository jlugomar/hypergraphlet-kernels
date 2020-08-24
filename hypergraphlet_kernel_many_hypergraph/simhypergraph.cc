#include "simhypergraph.h"
#include "string.h"
#include <fstream>
#include <set>
#include <queue>


vector<unsigned> SimpleHypergraph::get_incident_edges(unsigned vertex)  {
    if (vertex >= this->node_labels.size())  {
        cerr << "ERROR: Vertex index " << vertex << " >= hypergraph nodes size " << this->node_labels.size() << "." << endl; exit(1);
    }
    return this->edge_set[vertex];
}


SimpleHypergraph SimpleHypergraph::read_hypergraph(const char *nlabels_file, const char *elabels_file, const char *hypergraph_file)  {
    ifstream nlin(nlabels_file, ios::in);
	ifstream elin(elabels_file, ios::in);
    ifstream hgin(hypergraph_file, ios::in);

    if (nlin.fail())  {
        cerr << "ERROR: Hypergraph vertex labels  file " << nlabels_file << " could not be opened." << endl; exit(1);
    }

    if (elin.fail())  {
        cerr << "ERROR: Hypergraph edge labels  file " << elabels_file << " could not be opened." << endl; exit(1);
    }

    if (hgin.fail())  {
        cerr << "ERROR: Hypergraph file " << hypergraph_file << " could not be opened." << endl; exit(1);
    }

    SimpleHypergraph hg;
    
    // Read hypergraph node labels
    string line;
    if (getline(nlin,line))
        hg.node_labels = strip(line);

    // Read hypergraph edge labels
    if (getline(elin,line))
        hg.edge_labels = strip(line);

    // Set size of input hypergraph
    hg.vertex_set.resize(hg.edge_labels.size());
    hg.edge_set.resize(hg.node_labels.size());

    // Read hypergraph
	while(getline(hgin, line))  {
        vector<string> tokens = split(line, '\t');
        unsigned edge = to_i(tokens[0]); // edge id
        
        // Parse vertex set list
        for (unsigned i=1; i<tokens.size(); i++)  {
            unsigned vertex = to_i(tokens[i]); // vertex id
            
            if (vertex >= hg.node_labels.size())  {
                cerr << "ERROR: Node index " << vertex << " >= hypergraph nodes size " << hg.node_labels.size() << " in target hypergraph file " << hypergraph_file << "." << endl; exit(1);
            }

            if (edge >= hg.vertex_set.size()) {
                cerr << "ERROR: Edge index " << edge << " >= hypergraph edges size " << hg.vertex_set.size() << " in target hypergraph file " << hypergraph_file << "." << endl; exit(1);
            }

            hg.vertex_set[edge].push_back(vertex);
            hg.edge_set[vertex].push_back(edge);
        }
    }
    return hg;
}


SimpleHypergraph SimpleHypergraph::read_dual_hypergraph(const char *nlabels_file, const char *elabels_file, const char *hypergraph_file)  {
    ifstream nlin(nlabels_file, ios::in);
    ifstream elin(elabels_file, ios::in);
    ifstream hgin(hypergraph_file, ios::in);

    if (nlin.fail())  {
        cerr << "ERROR: Hypergraph vertex labels  file " << nlabels_file << " could not be opened." << endl; exit(1);
    }

    if (elin.fail())  {
        cerr << "ERROR: Hypergraph edge labels  file " << elabels_file << " could not be opened." << endl; exit(1);
    }

    if (hgin.fail())  {
        cerr << "ERROR: Hypergraph file " << hypergraph_file << " could not be opened." << endl; exit(1);
    }

    SimpleHypergraph hg, dual_hg;

    // Read hypergraph node labels and set dual hypergraph edge labels
    string line;
    if (getline(nlin,line))  { 
        hg.node_labels = strip(line);
        dual_hg.edge_labels = strip(line);
    }
        
    // Read hypergraph edge labels and set dual hypergraph node labels
    if (getline(elin,line))  {
        hg.edge_labels = strip(line);
        dual_hg.node_labels = strip(line);
    }

    // Set size of input hypergraph
    hg.vertex_set.resize(hg.edge_labels.size());
    hg.edge_set.resize(hg.node_labels.size());
    
    // Set size of dual hypergraph
    dual_hg.vertex_set.resize(dual_hg.edge_labels.size());
    dual_hg.edge_set.resize(dual_hg.node_labels.size());

    // Read hypergraph
	while(getline(hgin, line))  {
        vector<string> tokens = split(line, '\t');
        unsigned edge = to_i(tokens[0]); // edge id
		
		// Parse vertex set list
        for (unsigned i=1; i<tokens.size(); i++)  {
            unsigned vertex = to_i(tokens[i]); // vertex id

            if (vertex >= hg.node_labels.size())  {
                cerr << "ERROR: Node index " << vertex << " >= hypergraph nodes size " << hg.node_labels.size() << " in target hypergraph file " << hypergraph_file << "." << endl; exit(1);
            }

            if (edge >= hg.vertex_set.size()) {
                cerr << "ERROR: Edge index " << edge << " >= hypergraph edges size " << hg.vertex_set.size() << " in target hypergraph file " << hypergraph_file << "." << endl; exit(1);
            }

            hg.vertex_set[edge].push_back(vertex);
            hg.edge_set[vertex].push_back(edge);
        }
    }

    // Generate dual hypergraph
    for (unsigned epsilon=0; epsilon<hg.edge_set.size(); epsilon++)  {  // Iterate over input hypergraph nodes which are dual hypergraph hyperedges
        for (unsigned e_i=0; e_i<hg.edge_set[epsilon].size(); e_i++)  {
            unsigned edge = hg.edge_set[epsilon][e_i];
            dual_hg.vertex_set[epsilon].push_back(edge);
            dual_hg.edge_set[edge].push_back(epsilon);
            if (hg.edge_set[epsilon].size() == 1)  {
                dual_hg.node_labels.push_back(dual_hg.node_labels[edge]);
                dual_hg.edge_set.resize(dual_hg.node_labels.size());
                dual_hg.vertex_set[epsilon].push_back(dual_hg.node_labels.size() - 1);
                dual_hg.edge_set[(dual_hg.node_labels.size() - 1)].push_back(epsilon);
            }
        }
    }
    return dual_hg;
}


void SimpleHypergraph::print_hypergraph(const char* nlabels_filename, const char* elabels_filename, const char* hypergraph_filename)  {
    ofstream nlabels_file(nlabels_filename, ios::out);
    // Output vertex labels for hypergraph
    nlabels_file << node_labels << endl;
    nlabels_file.close();

    ofstream elabels_file(elabels_filename, ios::out);
    // Output hyperedge labels for hypergraph
    elabels_file << edge_labels << endl;
    elabels_file.close();

    ofstream hypergraph_file(hypergraph_filename, ios::out);
    // Output hypergraph
    for (unsigned hyperedge=0; hyperedge<vertex_set.size(); hyperedge++)  {  // Iterate over each hyperedge on hypergraph
        hypergraph_file << hyperedge;
        for (unsigned v_i=0; v_i<vertex_set[hyperedge].size(); v_i++)  {  // Iterate over each vertex in current hyperedge
            unsigned u = vertex_set[hyperedge][v_i];
            hypergraph_file << "\t" << u;
        }
        hypergraph_file << endl;
    }
    hypergraph_file.close();
}
    

inline void insert_edge(vector<pair <unsigned, char> > &adj, pair <unsigned, char> target)  {
    vector<pair <unsigned, char> > :: iterator it;
    for (it = adj.begin(); it != adj.end(); it++)  {
        if (it->first == target.first && it->second == target.second)
            return;
    }
    adj.push_back(target);
}


void SimpleHypergraph::transform_clique_expansion(const char* labels_filename, const char* graph_filename)  {
    vector<vector<pair<unsigned, char> > > adj; // Adjacency lists
    adj.resize(node_labels.size());
    for (unsigned hyperedge=0; hyperedge<vertex_set.size(); hyperedge++)  {  // Iterate over each hyperedge on hypergraph
        for (unsigned u_i=0; u_i<vertex_set[hyperedge].size(); u_i++)  {  // Iterate over each vertex in current hyperedge
            unsigned u = vertex_set[hyperedge][u_i];
            for (unsigned v_j=0; v_j<u_i; v_j++)  {
                unsigned v = vertex_set[hyperedge][v_j];
                pair <unsigned, char> p1 (v, edge_labels[hyperedge]);
                pair <unsigned, char> p2 (u, edge_labels[hyperedge]);
                insert_edge(adj[u], p1);
                insert_edge(adj[v], p2);
            }
        }
    }
    ofstream labels_file(labels_filename, ios::out);
    // Output node labels for transformed graph
    labels_file << node_labels << endl;
    labels_file.close();

    ofstream graph_file(graph_filename, ios::out);
    // Output transformed graph
    for (unsigned i=0; i<adj.size(); i++)  {
        graph_file << i;
        for (unsigned j=0; j<adj[i].size(); j++)  {
            graph_file << "\t" << adj[i][j].first << "," << adj[i][j].second;
        }
        graph_file << endl;
    }
    graph_file.close();
}


void SimpleHypergraph::transform_star_expansion(const char* labels_filename, const char* graph_filename)  {
    vector<vector<pair<unsigned, char> > > adj; // Adjacency lists
    adj.resize((node_labels.size() + edge_labels.size()));
    for (unsigned hyperedge=0; hyperedge<vertex_set.size(); hyperedge++)  {  // Iterate over each hyperedge on hypergraph
        unsigned edgeID = node_labels.size() + hyperedge;
        for (unsigned u_i=0; u_i<vertex_set[hyperedge].size(); u_i++)  {  // Iterate over each vertex in current hyperedge
            unsigned u = vertex_set[hyperedge][u_i];
            pair <unsigned, char> p1 (u, edge_labels[hyperedge]);
            pair <unsigned, char> p2 (edgeID, edge_labels[hyperedge]);
            insert_edge(adj[edgeID], p1);
            insert_edge(adj[u], p2);
        }
    }
    ofstream labels_file(labels_filename, ios::out);
    // Output node labels for transformed graph
    labels_file << node_labels << edge_labels << endl;
    labels_file.close();

    ofstream graph_file(graph_filename, ios::out);
    // Output transformed graph
    for (unsigned i=0; i<adj.size(); i++)  {
        graph_file << i;
        for (unsigned j=0; j<adj[i].size(); j++)  {
            graph_file << "\t" << adj[i][j].first << "," << adj[i][j].second;
        }
        graph_file << endl;
    }
    graph_file.close();
}


vector<unsigned> SimpleHypergraph::breadth_first_sort(unsigned hg_root) const  {
    vector<unsigned> dist(node_labels.size(), UINT_MAX);
    dist[hg_root] = 0;
   
    queue<unsigned> Q;
    Q.push(hg_root);
    
    while (!Q.empty())  {
        unsigned i = Q.front();

        for (unsigned j=0; j<edge_set[i].size(); j++) {
            unsigned e = edge_set[i][j];
            for (unsigned k=0; k<vertex_set[e].size(); k++)  {
                unsigned l = vertex_set[e][k];
            
                if (UINT_MAX==dist[l])  {
                    dist[l] = dist[i] + 1;
                    // Ignore vertices than cannot be reached by hypergraphlet enumeration algorithm
                    if (dist[l] <= 2)
                        Q.push(l);
                }
            }
        }
        Q.pop();
    }
    return dist;
}


vector<set<unsigned> > SimpleHypergraph::get_neighbors() const  {
    vector<set<unsigned> > neighbors;
    
    neighbors.resize(node_labels.size());
    
    for (unsigned hyperedge=0; hyperedge<vertex_set.size(); hyperedge++)  {  // Iterate over each hyperedge on hypergraph
        for (unsigned u_i=0; u_i<vertex_set[hyperedge].size(); u_i++)  {  // Iterate over each vertex in current hyperedge
            unsigned u = vertex_set[hyperedge][u_i];
            for (unsigned v_i=0; v_i<u_i; v_i++)  {  // Iterate over other vertices in current hyperedge
            	unsigned v = vertex_set[hyperedge][v_i];
            	neighbors[u].insert(v);
            	neighbors[v].insert(u);
            }
        }
    }
    return neighbors;
}
