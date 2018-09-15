#include <iostream>
#include <set>
#include <algorithm>
#include <limits>

#include "classification.h"
#include "globals.h"
#include "matrix.h"

using namespace std;
using namespace globals;

typedef matrix mat;
typedef vector<number> vec;

int flood_fill(vector<int>& res, const vector<pair<int, int>>& tree, int node, int number) {
    if (res[node] != 0)
        return 0;
    res[node] = number;

    int tot = 1;

    for (int i = 0; i < tree.size(); i++) {
        if (tree[i].first == node)
            tot += flood_fill(res, tree, tree[i].second, number);
        else if (tree[i].second == node)
            tot += flood_fill(res, tree, tree[i].first, number);
    }

    return tot;
}

void walk_tree(vector<int>& centrality, const vector<pair<int, int>>& tree, pair<int, int> edge, int current) {
    int other = edge.first == current ? edge.second : edge.first;

    for (int i = 0; i < tree.size(); i++)
        if ((tree[i].first == other && tree[i].second != current)
                || (tree[i].second == other && tree[i].first != current)) {
            walk_tree(centrality, tree, tree[i], other);
            centrality[current]++;
        }
}

vector<int> classification::group_models(const vector<Lambda>& hmms, int number_of_groups, bool optimal_guess) {
    //assert(number_of_groups > 0);
    //assert(number_of_groups > hmms.size());
    if (number_of_groups > hmms.size()) {
        cerr << "Too few hmms (" << hmms.size() << " vs " << number_of_groups << ")" << endl;
        vector<int> res;
        return res;
    }
    matrix distances(hmms.size(), hmms.size());

    number distance = 0;

    cerr << "Calculating distances" << endl;

    //calculate distances between models
    for (int i = 0; i < distances.getHeight(); i++) {
        for (int j = 0; j < i; j++) {
            distance = hmms[i].A.distance_squared(hmms[j].A) + hmms[i].B.distance_squared(hmms[j].B);
            distances.set(i, j, distance);
            distances.set(j, i, distance);
        }
    }

    cerr << distances << endl;

    cerr << "Calculating MST" << endl;

    //prim's algorithm to build a MST
    vector<pair<int, int>> tree;
    std::set<int> nodes_in_tree;
    std::set<int> nodes_not_in_tree;

    nodes_in_tree.insert(0);

    for (int i = 1; i < hmms.size(); i++)
        nodes_not_in_tree.insert(i);

    pair<int, int> minedge;
    number mindist;

    while (nodes_in_tree.size() < hmms.size()) {
        minedge = {-1, -1};
        mindist = numeric_limits<number>::max();

        for (auto from : nodes_in_tree) {
            for (auto to : nodes_not_in_tree) {
                if (distances.get(from, to) < mindist) {
                    minedge = {from, to};
                    mindist = distances.get(from, to);
                }
            }
        }

        tree.push_back(minedge);
        nodes_in_tree.insert(minedge.second);
        nodes_not_in_tree.erase(minedge.second);
    }

    cerr << "Sorting MST" << endl;
    
    //sort MST edges
    sort(tree.begin(), tree.end(),
        [distances](pair<int, int> a, pair<int, int> b) {
            return distances.get(a.first, a.second) < distances.get(b.first, b.second);
        }
    );

    for (auto p : tree)
        cerr << "(" << p.first << ", " << p.second << ")" << "\t";
    cerr << endl;

    cerr << "Pruning " << (number_of_groups - 1) << " edges from MST" << endl;

    //remove number_of_groups longest edges
    std::set<int> nodes_to_check;
    for (int i = 0; i < number_of_groups - 1 && tree.size() - i - 1 >= 0; i++) {
        pair<int, int> p = tree[tree.size() - i - 1];
        nodes_to_check.insert(p.first);
        nodes_to_check.insert(p.second);
    }

    for (int i = 0; i < number_of_groups - 1 && tree.size() - 1 >= 0; i++)
        tree.erase(tree.begin() + tree.size() - 1);

    cerr << "Assigning numbers to clusters" << endl;

    //flood fill the clusters with their numbers
    vector<int> assignment(hmms.size());
    int groupnumber = 1;

    for (auto node : nodes_to_check) {
        int groupsize = flood_fill(assignment, tree, node, groupnumber);
        if (groupsize > 0)
            groupnumber++;
    }

    cerr << "Group assignment: " << endl;
    for (auto number : assignment)
        cerr << number << " ";
    cerr << endl;
    
    if (optimal_guess) {
        vector<int> optimal_guesses(number_of_groups + 1);
        vector<int> centrality(hmms.size());
        
        for (int i = 0; i < hmms.size(); i++)
            for (auto edge : tree)
                if (edge.first == i || edge.second == i)
                    walk_tree(centrality, tree, edge, i);
        
        cerr << "Centrality: ";
        for (auto node : centrality)
            cerr << node << " ";
        cerr << endl;

        for (int group = 1; group < optimal_guesses.size(); group++) {
            int maxhmm = -1;
            int maxcentrality = -1;
            for (int i = 0; i < hmms.size(); i++)
                if (assignment[i] == group && centrality[i] > maxcentrality) {
                    maxcentrality = centrality[i];
                    maxhmm = i;
                }
            optimal_guesses[group] = maxhmm;
        }

        return optimal_guesses;
    } else
        return assignment;
}