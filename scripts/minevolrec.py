""" Inferring Transcript Phylogenies from Transcript Ortholog Clusters. 

Usage:
======
    python3 minclevol.py    -l argument0 -clus argument1 -nhx argument2 -map argument3 -matx argument4 -join argument5 
                            -outf argument6 -outp argument7 -c argument8

    example: (see execution_inferring_transcript_phylogenies.sh)
"""

#!/usr/bin/env python3
#!interpreter [optional-arg]
__author__ = "Wend Yam D D Ouedraogo"
__copyright__ = "Copyright 2023, CoBIUS Lab"
__credits__ = ["Wend Yam Donald Davy Ouedraogo"]
__license__ = "CoBIUS License"
__version__ = "1.0.0"
__maintainer__ = "Wend Yam Donald Davy Ouedraogo"
__email__ = "wend.yam.donald.davy.ouedraogo@usherbrooke.ca"
__status__ = "production"
__date___ = "2023-08-04"

########################
###### librairies ######
########################
import argparse
import pandas as pd
import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import copy as cp
import phylotreelib as pt
import os
import time
import random

def build_arguments_parser():
    ''' parsor program parameter function '''
    parser = argparse.ArgumentParser(description="parsor program parameter")
    parser.add_argument('-l', '--label', required=True, help='0 | 1, Boolean variable controlling the utilization of the labeled version of the algorithm (0) or the non-labeled version (1).')
    parser.add_argument('-clus', '--clusters', required=True, help='FASTA file containing cluster IDs in front of each id_transcript separated by semicolons')
    parser.add_argument('-nhx', '--nhx', required=True, help='txt file containg gene tree \t\tformat => NHX format')
    parser.add_argument('-map', '--mappings', required=True, help='FASTA file containg transcripts and their corresponding genes separated by semicolons.')
    parser.add_argument('-matx', '--matrix', required=True, help='CSV file containing matrix separated by comma \';\'(header: True, index: True)Note: Index must correspond to the header.')
    parser.add_argument('-forest', '--forest', required=False, default=1, help='0 | 1, boolean variable controls the reconstruction of a rooted binary transcript tree(1) or transcript forest(0).')
    parser.add_argument('-forest_threshold', '--forest_threshold', required=False, default=0.3, help='float variable represents threshold to cut the dendogram(minimum evolution tree)')
    parser.add_argument('-join', '--join', required=False, default=1, help='min(1) | mean(0) | max(2). By default(1)')
    parser.add_argument('-outf', '--output', required=False, default='.', help='output folder')
    parser.add_argument('-outp', '--prefix', required=False, default='iteration', help='prefix of output files')
    parser.add_argument('-c', '--compute', required=False, default=1, help='0 | 1, Compute all solutions(By default False(1))')
    return parser

def get_dataframe_structures_from_inputs(clusters_input, mappings_input):
    ''' get dataframe from user's inputs '''
    clusters_file = open(clusters_input, 'r')
    mappings_file = open(mappings_input, 'r')
    clusters_df = pd.DataFrame(columns=['id_transcript', 'cluster_id'])
    mappings_df = pd.DataFrame(columns=['id_transcript', 'id_gene'])

    clusters_lines = [str(_).split('\n')[0] for _ in clusters_file.readlines()]
    for clusters_line in clusters_lines:
        if clusters_line.startswith('>'):
            id_transcript = clusters_line.split(':')[0].split('>')[-1]
            id_cluster = clusters_line.split(':')[-1]
            clusters_df.loc[id_transcript] = [id_transcript, id_cluster]
        else:
            raise TypeError('Inappropriate argument -clus format')
    
    mappings_lines = [str(_).split('\n')[0] for _ in  mappings_file.readlines()]
    for mappings_line in mappings_lines:
        if mappings_line.startswith('>'):
            id_transcript = mappings_line.split('>')[-1].split(':')[0]
            id_gene = mappings_line.split(':')[-1]
            mappings_df.loc[id_transcript] = [id_transcript, id_gene]
        else:
            raise TypeError('Inappropriate argument -map format')

    mappings_file.close()
    clusters_file.close()
    mappings_df = mappings_df.set_index('id_transcript').join(clusters_df.set_index('id_transcript'))

    return mappings_df

def get_ortholog_trees(mappings_df, gene_tree, matrix):
    """
    Computes ortholog trees from orthologous groups of transcripts
    Parameters:
        mappings_df: a pandas dataframe(3 columns[id_transcript, id_gene, cluster_id])
        gene_tree: a phylogenetic Tree ete3 python object representing a gene tree, feature: gene_name
        matrix : Pairwise similarity matrix
    Returns:
        Dictionary [keys:clusters, values: transcripts subtree ete3 python object representing transcripts clusters]
    """
    orthologs_trees_dict = {}
    clusters_indexes = list(set(mappings_df.cluster_id.values))
    for cluster_index in clusters_indexes:
        # get the df that contain the transcripts in the same group
        sub_df = mappings_df[mappings_df['cluster_id']==cluster_index]
        # copy the gene tree and prune it
        copy_gene_tree = cp.deepcopy(gene_tree)
        genes_in_the_same_group = list(set(sub_df.id_gene.values))
        if len(genes_in_the_same_group) == 1:
            copy_gene_tree = Tree('{};'.format(genes_in_the_same_group[0]))
            observed_node = copy_gene_tree&genes_in_the_same_group[0]
            observed_node.add_features(gene_name=genes_in_the_same_group[0])
        else:
            copy_gene_tree.prune(genes_in_the_same_group)
        leaves = [leaf for leaf in copy_gene_tree.get_leaves()]
        for node in leaves:
            node_name = node.gene_name
            transcripts = sub_df[sub_df['id_gene']==node_name].index.values
            if len(transcripts) == 1:
                node.name = transcripts[0]
            else:
                distance_partial_matrix = get_partial_distance_matrix(transcripts, matrix)
                partial_dmat = pt.Distmatrix.from_distdict(distance_partial_matrix)
                partial_tree = partial_dmat.upgma()
                newick_partial_tree = partial_tree.newick()
                t = Tree(newick_partial_tree, format=1)
                t.resolve_polytomy(recursive=True)
                if node.is_root():
                    copy_gene_tree = t
                else:
                    node.add_child(t)
                    node.delete()

        ortholog_tree = copy_gene_tree
        labelled_ortholog_tree = label_ortholog_tree(ortholog_tree, gene_tree)
        orthologs_trees_dict['orthologs_group_{}'.format(cluster_index)] = labelled_ortholog_tree
    return orthologs_trees_dict

def label_ortholog_tree(t, gene_tree):
    """
    Labels the orthologous tree based on specified criteria.

    This function performs the labeling of an orthologous tree, considering the LCA
    reconciliation with respect to the Gene Tree given at the input.
   
    Parameters:
    - t: the ortholog tree to label.
    - gene_tree: a given NHX gene tree consistent with the transcripts in t.

    Returns:
    - a labelled ortholog tree(features=['label', 'lca_gene', 'enum']).
    """
    for enum, node in enumerate(t.traverse('postorder')):
        if node.is_leaf():
            name_gene = str(node.name).split('_')[-1]
            node.add_features(lca_gene=name_gene, label='speciation')
        else:
            children = node.children
            left_child = children[0]
            right_child = children[1]
            lca_left = left_child.lca_gene
            lca_right = right_child.lca_gene
            if lca_left == lca_right:
                node.add_features(lca_gene=left_child.lca_gene, label='creation', name=enum)
            else:
                genes = [_.lca_gene for _ in left_child.get_leaves()+right_child.get_leaves()]
                if len(list(set(genes))) == 1:
                    common_ancestor = gene_tree&genes[0]
                else:
                    common_ancestor = gene_tree.get_common_ancestor(list(set(genes)))
                label_name = None
                if hasattr(common_ancestor, 'D'):
                    evol_type = common_ancestor.D
                    if evol_type == 'N':
                        label_name = 'speciation'
                    elif evol_type == 'Y':
                        label_name = 'duplication'
                    else:
                        raise ValueError('Inappropriate NHX format! -D')
                elif hasattr(common_ancestor, 'DD'):
                    evol_type = common_ancestor.DD
                    if evol_type == 'N':
                        label_name = 'speciation'
                    elif evol_type == 'Y':
                        label_name = 'duplication'
                    else:
                        raise ValueError('Inappropriate NHX format! -DD')
                else:
                    ValueError('Inappropriate NHX format!')

                node.add_features(lca_gene=common_ancestor.gene_name, label=label_name, name=enum)
    return t

def get_me_tree(ortholog_trees_set, matrix, join_input, forest, threshold):
    """
    Retrieves the minimum evolution tree based on the provided orthologous tree set,
    distance matrix, and specified join input.

    This function takes an orthologous tree set, a distance matrix, and a join input as
    parameters to calculate and return the minimum evolution tree. The minimum evolution
    tree is constructed by iteratively joining the most closely related clusters until
    the entire tree is formed. The process follows a minimum evolution criterion.

    Parameters:
    - ortholog_trees_set (dict): A set of orthologous trees used in the tree construction.
    - matrix (pandas.DataFrame): A similarity matrix representing the evolutionary distances
      between clusters.
    - join_input (integer 0|1|2): The method or criteria used for joining clusters in the tree
      construction process.
    - forest (integer 0|1):
    - threshold (float):

    Returns:
    - min_evolution_tree (ete3 Tree): The minimum evolution tree constructed based on the
      provided parameters.
    """
    ortholog_trees_indexes = list(ortholog_trees_set.keys())
    me_trees = []
    me_tree = None
    if len(ortholog_trees_indexes) == 1:
        me_tree = Tree('{};'.format(ortholog_trees_indexes[0]), format=1)
        me_trees.append(me_tree)
    elif len(ortholog_trees_indexes) == 2:
        print(forest)
        if forest == 0:
            me_trees.extend([ortholog_trees_indexes[0], ortholog_trees_indexes[1]])
        elif forest == 1:
            me_tree = Tree('({}, {});'.format(ortholog_trees_indexes[0], ortholog_trees_indexes[1]), format=1)
            me_trees.append(me_tree)
        else:
            raise ValueError('Inappropriate forest argument. forest argument must beinteger 0||1')
    elif len(ortholog_trees_indexes) > 2:
        distance_matrix = get_dictdistances_of_ortholog_trees(matrix, ortholog_trees_set, join_input)
        dmat = pt.Distmatrix.from_distdict(distance_matrix)
        nj_t = dmat.nj()
        newick_t = nj_t.newick()
        me_tree = Tree(newick_t, format=1)
        me_tree.resolve_polytomy(recursive=True)
        if forest == 1:
            me_trees.append(me_tree)
        elif forest == 0:
            print(me_tree.write(format=6))
            cut_trees = cutting_dendogram(me_tree, threshold, [])
            update_cut_trees = get_all_cut_trees(cut_trees, me_tree)                           
            me_trees.extend(update_cut_trees)
            print(me_trees)
        else:
            raise ValueError('Inappropriate forest argument. forest argument must beinteger 0||1')
    else:
        raise ValueError('Inappropriate join argument. Join argument must be 0|1|2')
    
    return me_trees

def cutting_dendogram(node, threshold, M):
    
    if not node.is_leaf():
        children = node.children
        left_child = children[0]
        right_child = children[1]
        d_left = node.get_distance(left_child)
        #print(d_left)
        cut_left = cp.deepcopy(left_child).detach()
        if d_left >= threshold and d_left != 1.0:
            cut_left = cp.deepcopy(left_child).detach()
            M.append(cut_left)
            
        d_right = node.get_distance(right_child)
        cut_right = cp.deepcopy(right_child).detach()
        #print(d_right)
        if d_right >= threshold and d_right != 1.0:
            
            M.append(cut_right)
            
        cutting_dendogram(cut_left, threshold, M)
        cutting_dendogram(cut_right, threshold, M)
    return M

def get_all_cut_trees(M, dendogram):
    print(M)
    #'''
    copy_t = cp.deepcopy(dendogram)
    if len(M) != 0:
        leaves = []
        for cut_tree in M:
            leaves.extend([str(_.name) for _ in cut_tree.get_leaves()])
        #print(leaves)
        leaves_left = [leaf.name for leaf in dendogram.get_leaves() if str(leaf.name) not in leaves]
        #print(leaves_left)
        if len(leaves_left) != 0:
            copy_t.prune(leaves_left)
            M.append(copy_t)
    else:
        M.append(copy_t)
    #'''
    #print(M)
    return M
  
def get_dictdistances_of_ortholog_trees(matrix, ortholog_trees_set, method):
    ''' convert similarity matrix --matx to distance matrix for neighbor joining algorithm computation '''
    dict_distances = {}
    ids_clusters = list(ortholog_trees_set.keys())
    for id_cluster_x in ids_clusters:
        leaves_x = [str(_.name) for _ in ortholog_trees_set[id_cluster_x].get_leaves()]
        tmp_dict_distances = {}
        for id_cluster_y in ids_clusters:
            if id_cluster_x == id_cluster_y:
                tmp_dict_distances[id_cluster_y] = None
            else:
                leaves_y = [str(_.name) for _ in ortholog_trees_set[id_cluster_y].get_leaves()]
                score = 0.0
                if method == 1 :
                    for leaf_x in leaves_x:
                        for leaf_y in leaves_y:
                            current_score = matrix.loc[leaf_x][leaf_y]
                            if current_score > score:
                                score = current_score
                    tmp_dict_distances[id_cluster_y] = round(1.0 - float(score), 3)
                elif method == 2:
                    for leaf_x in leaves_x:
                        for leaf_y in leaves_y:
                            current_score = matrix.loc[leaf_x][leaf_y]
                            if current_score < score:
                                score = current_score
                    tmp_dict_distances[id_cluster_y] = round(1.0 - float(score), 3)
                elif method == 0:
                    for leaf_x in leaves_x:
                        for leaf_y in leaves_y:
                            score += matrix.loc[leaf_x][leaf_y]
                    tmp_dict_distances[id_cluster_y] = round(1.0 - float(score/(len(leaves_x)*len(leaves_y))), 3)
                else:
                    raise ValueError('Inappropriate join argument. Join argument must be 0|1|2')
        dict_distances[id_cluster_x] = tmp_dict_distances
    return dict_distances

def get_partial_distance_matrix(list_transcripts, matrix):
    partial_dict_distances = {}
    for transcript_x in list_transcripts:
        tmp_partial_dict_distances = {}
        for transcript_y in list_transcripts:
            if transcript_x == transcript_y:
                tmp_partial_dict_distances[transcript_y] = None
            else:
                tmp_partial_dict_distances[transcript_y] = round(1-float(matrix.loc[transcript_x][transcript_y]),3)
        partial_dict_distances[transcript_x] = tmp_partial_dict_distances
    return partial_dict_distances

def fix_locations(number, tree):
    """
    Adjusts the locations of elements in a tree based on a specified number.

    This function takes a numerical parameter and a tree structure as input. It performs
    adjustments to the locations of elements within the tree based on the specified number.

    Parameters:
    - number (int): A numerical parameter used to determine the adjustments.
    - tree (Tree): The tree structure that needs location adjustments.

    Returns:
    - adjusted_tree (Tree): The tree structure with updated element locations.
    """
    for node in tree.traverse('preorder'):
        number += 1
        node.add_features(location=number)
    return tree

def init_graph_of_recurrence(t1):
    """
    Initializes a graph of recurrence based on a given input.

    This function takes a parameter t1 and constructs an initial graph of recurrence.
    The graph represents relationships or dependencies among elements, and its construction
    is determined by the characteristics of the input parameter t1.

    Parameters:
    - t1 (ete3 Tree): Description of the input parameter t1.

    Returns:
    - recurrence_graph (Graph): The initialized graph of recurrence based on the provided input.
    """
    root_t1 = t1.get_tree_root().location
    g = Tree('{};'.format(root_t1))
    return g

def get_bipartitions(t):
    """
    Extracts bipartitions from a given tree.

    This function takes a tree structure as input and analyzes its topology to identify
    and extract bipartitions. Bipartitions represent the division of the tree into two
    disjoint sets of leaves. The extracted bipartitions provide insights into the hierarchical
    relationships among elements in the tree.

    Parameters:
    - tree (Tree): The input tree structure for bipartition analysis.

    Returns:
    - bipartition_list (list): A list containing the identified bipartitions in the tree.
      Each bipartition is represented as a tuple of two sets, indicating the partitioning
      of leaves in the tree.
    """
    if t.is_leaf():
        return [None, None]
    else:
        t_copy = cp.deepcopy(t)
        children = t_copy.children
        left_child = children[0].detach()
        right_child = children[0].detach()
        return [cp.deepcopy(left_child), cp.deepcopy(right_child)]

def minCLevol(t1, t2, graph_of_recurrence, mappings_df, gene_tree, inters_dict):
    """
    Computes the best reconciliation cost of a Minimum Evolution-Reconciliation 
    SuperTree(MinCLevol) between two trees.

    This function calculates the Minimum Evolution-Reconciliation SuperTree between two given trees,
    tree1 and tree2, utilizing a graph of recurrence, a DataFrame of mappings, a gene tree,
    and an pre-computed inters() stored in a dictionary.

    Parameters:
    - t1 (Tree): The first input tree for comparison.
    - t2 (Tree): The second input tree for comparison.
    - graph_of_recurrence (Graph): A graph representing recurrence relationships among elements.
    - mappings_df (DataFrame): A DataFrame containing mappings between elements in the trees.
    - gene_tree (Tree): The gene tree representing evolutionary relationships.
    - inters_dict (dict): A dictionary representing inter(x,y) between elements.

    Returns:
    - the best reconciliation cost after merge the two trees t1 and t2
    - the graph of recurrence.
    """
    if t1 == None and t2 == None:
        return [0, graph_of_recurrence]
    elif t1 == None and t2 != None:
        return [t2.get_tree_root().global_cost, graph_of_recurrence]
    elif t1 != None and t2 == None:
        return [t1.get_tree_root().global_cost, graph_of_recurrence]
    else:
        leaves_t2 = [str(_.name) for _ in t2.get_leaves()]
        leaves_t1 = [str(_.name) for _ in t1.get_leaves()]
        is_t1_a_leaf = len(leaves_t1) == 1

        # get bipartitions
        t1l = None
        t1r = None
        if not is_t1_a_leaf:
            bipartitions = get_bipartitions(t1)
            t1l = cp.deepcopy(bipartitions[0])
            t1r = cp.deepcopy(bipartitions[1])
        
        # update the graph
        name_node = t1.get_tree_root().location
        node_of_interest = graph_of_recurrence&name_node
        name_leftNode = None
        name_rightNode = None
        if not is_t1_a_leaf:
            name_leftNode = t1l.get_tree_root().location
            name_rightNode = t1r.get_tree_root().location
            for name in [name_leftNode, name_rightNode]:
                tmp = Tree('{};'.format(name))
                node_of_interest.add_child(tmp)
        node_of_interest.name = str(name_node) + 'fusion'
        tmp_node = Tree('{};'.format(name_node))
        node_of_interest.add_child(tmp_node)

        # first case <====> (T1, None, None, T2)
        local_merge_first_case = get_local_lca_reconciliation_cost(leaves_t1, leaves_t2, mappings_df, gene_tree, inters_dict)
        cost_firstCase = local_merge_first_case[0] + minCLevol(t1, None, graph_of_recurrence, mappings_df, gene_tree, inters_dict)[0] + minCLevol(None, t2, graph_of_recurrence, mappings_df, gene_tree, inters_dict)[0]
        evolution_first_case = local_merge_first_case[1]

        if not is_t1_a_leaf:
            leaves_t1l = [str(_.name) for _ in t1l.get_leaves()]
            leaves_t1r = [str(_.name) for _ in t1r.get_leaves()]

            # second case <====> (T1L, T2, T1R, None)
            local_merge_second_case = get_local_lca_reconciliation_cost(leaves_t1l+leaves_t2, leaves_t1r, mappings_df, gene_tree, inters_dict)
            cost_secondCase = local_merge_second_case[0] + minCLevol(t1l, t2, graph_of_recurrence, mappings_df, gene_tree, inters_dict)[0] + minCLevol(t1r, None, graph_of_recurrence, mappings_df, gene_tree, inters_dict)[0]
            evolution_second_case = local_merge_second_case[1]
            
            # third case <====> (T1L, None, T1R, T2)
            local_merge_third_case = get_local_lca_reconciliation_cost(leaves_t1l, leaves_t1r+leaves_t2, mappings_df, gene_tree, inters_dict)
            cost_thirdCase = local_merge_third_case[0] + minCLevol(t1l, None, graph_of_recurrence, mappings_df, gene_tree, inters_dict)[0] + minCLevol(t1r, t2, graph_of_recurrence, mappings_df, gene_tree, inters_dict)[0]
            evolution_third_case = local_merge_third_case[1]
        else:
            cost_secondCase = np.inf
            cost_thirdCase = np.inf
            evolution_third_case = 'nil'
            evolution_second_case = 'nil'

        # update the graph
        node_1stCase = graph_of_recurrence&name_node
        node_1stCase.add_features(cost=cost_firstCase, solution=name_node, all_costs={name_node: cost_firstCase}, evol_type=evolution_first_case)
        node_parent = node_1stCase.up
        solution = min(cost_firstCase, cost_secondCase, cost_thirdCase)
        evolution = {}
        for sol in [(cost_firstCase, name_node, evolution_first_case), (cost_secondCase, name_leftNode, evolution_second_case), (cost_thirdCase, name_rightNode, evolution_third_case)]:
            evolution[sol[1]]=sol[2]

        if name_leftNode != None and name_rightNode != None:
            node_parent.add_features(cost=solution, all_costs={name_node: cost_firstCase, name_leftNode: cost_secondCase, name_rightNode: cost_thirdCase}, evol_type=evolution)
        else:
            node_parent.add_features(cost=solution, all_costs={name_node: cost_firstCase}, evol_type=evolution)
        return [solution, graph_of_recurrence]

def backtracking_minCLevol(best_rec_cost, recursion_tree, t1, t2, matrix):
    """
    Performs backtracking to find the Minimum Evolution-Reconciliation 
    SuperTree.

    This function utilizes backtracking to explore possible solutions for Minimum Evolution-Reconciliation 
    SuperTree problem. It takes the best_rec_cost, recursion_tree,
    t1, t2, and a similarity matrix as input, and iteratively explores the solution space to
    determine the Minimum Evolution-Reconciliation SuperTree.

    Parameters:
    - best_rec_cost (float): The best known recurrence cost in the search process.
    - recursion_tree (Tree): The recursion tree guiding the backtracking process.
    - t1 (Tree): The first input tree for comparison.
    - t2 (Tree): The second input tree for comparison.
    - matrix (pandas DataFrame): A similarity matrix representing the evolutionary distances
      between clusters.

    Returns:
    - All Minimum Evolution-Reconciliation SuperTrees found through the backtracking process.
    """
    leaves = recursion_tree.get_leaves()
    solutions_index = []
    for leaf in leaves:
        path = [leaf] + leaf.get_ancestors() 
        cur_best = None
        for i in range(len(path), 0, -1):
            if i==len(path):
                ancestor = path[i-1]
                dict_costs = ancestor.all_costs
                min_indexes = [str(index) for index in dict_costs.keys() if dict_costs[index]==best_rec_cost]
                cur_best = min_indexes
            elif i==1:
                ancestor = path[i-1]
                ancestor_name = ancestor.name
                if ancestor_name in cur_best:
                    evol_type = ancestor.evol_type
                    min_score = 0.0
                    max_score = np.inf
                    mean_score = None
                    observed_node = t1.search_nodes(location=int(ancestor_name))[0]
                    transcripts_leaves_t1 = [node.name for node in observed_node.get_leaves()]
                    transcripts_leaves_t2 = [node.name for node in t2.get_leaves()]
                    scores = []
                    for transcript_t1 in transcripts_leaves_t1:
                        for transcript_t2 in transcripts_leaves_t2:
                            score = matrix.loc[transcript_t1][transcript_t2]
                            scores.append(score)
                            if score > min_score:
                                min_score = score
                            if score < max_score:
                                max_score = score
                    mean_score = round(np.mean(scores), 3)
                    solutions_index.append((best_rec_cost, ancestor_name, (min_score, max_score, mean_score), evol_type))
                else:
                    break
            else:
                ancestor = path[i-1]
                ancestor_name = ancestor.name
                ancestor_index = ancestor_name.split('fusion')[0]
                if ancestor_index in cur_best:
                    dict_costs = ancestor.all_costs
                    min_cost = ancestor.cost
                    min_indexes = [str(index) for index in dict_costs.keys() if dict_costs[index]==min_cost]
                    cur_best = min_indexes
                else:
                    break
    return [solutions_index, [t1, t2]]

def compute_solution(bool_label, name_node, t1_original, t2_original, i, evolution, mappings_df, gene_tree, inters_dict):
    t1 = cp.deepcopy(t1_original)
    t2 = cp.deepcopy(t2_original)
    node_of_interest_in_t1 = t1.search_nodes(location=int(name_node))
    node_of_interest_in_t2 = t2.search_nodes(location=int(name_node))

    tree = None
    tree_grafted = None
    node_of_interest = None
    if len(node_of_interest_in_t1) == 1:
        node_of_interest = node_of_interest_in_t1[0]
        tree = t1
        tree_grafted = t2
    elif len(node_of_interest_in_t2) == 1:
        node_of_interest = node_of_interest_in_t2[0]
        tree = t2
        tree_grafted = t1
    else:
        raise ValueError('Node not found !')
    
    root_t = tree.get_tree_root()
    induction_node = Tree('{};'.format(str(i)+'_induced'))
    if node_of_interest == root_t:
        induction_node.add_child(t1)
        induction_node.add_child(t2)
        # update
        children = induction_node.children
        leaves_t1 = [node.name for node in children[0].get_leaves()]
        leaves_t2 = [node.name for node in children[1].get_leaves()]
        local_score_1 = get_local_lca_reconciliation_cost(leaves_t1, leaves_t2, mappings_df, gene_tree, inters_dict)[0]

        global_score_1 = local_score_1 + children[0].global_cost + children[1].global_cost
        if bool_label == 0:
            induction_node.add_features(local_cost=local_score_1, global_cost=global_score_1, label=evolution)
        else:
            induction_node.add_features(local_cost=local_score_1, global_cost=global_score_1)
        return induction_node
    else:
        node_parent = node_of_interest.up
        induction_node.add_child(tree_grafted)
        tmp = node_of_interest.detach()
        induction_node.add_child(tmp)
        node_parent.add_child(induction_node)
        # update
        leaves_t1 = [node.name for node in tree_grafted.get_leaves()]
        leaves_t2 = [node.name for node in tmp.get_leaves()]
        local_score_2 = get_local_lca_reconciliation_cost(leaves_t1, leaves_t2, mappings_df, gene_tree, inters_dict)[0]
        global_score_2 = local_score_2 + tmp.global_cost + tree_grafted.global_cost
        if bool_label == 1:
            induction_node.add_features(local_cost=local_score_2, global_cost=global_score_2)
        else:
            induction_node.add_features(local_cost=local_score_2, global_cost=global_score_2, label=evolution)
        ancestors_nodes = [node for node in induction_node.get_ancestors()]
        for anc_node in ancestors_nodes:
            children = anc_node.children
            left_child = children[0]
            right_child = children[1]
            leaves_left = [node.name for node in left_child.get_leaves()]
            leaves_right = [node.name for node in right_child.get_leaves()]
            local_score_3 = get_local_lca_reconciliation_cost(leaves_left, leaves_right, mappings_df, gene_tree, inters_dict)[0]
            global_score_3 = local_score_3 + left_child.global_cost + right_child.global_cost
            anc_node.local_cost = local_score_3
            anc_node.global_cost = global_score_3
        return tree

def get_inter(tree):
    inters_dict = {}
    already_checked = []
    nodes = [node for node in tree.traverse('postorder')]
    for node_x in nodes:
        for node_y in nodes:
            pair = (node_x.gene_name, node_y.gene_name)
            pair_sym = (node_y.gene_name, node_x.gene_name)
            if pair not in already_checked and pair_sym not in already_checked:
                already_checked.append(pair_sym)
                already_checked.append(pair)
                if node_x.gene_name == node_y.gene_name:
                    inters_dict[pair] = 0
                    inters_dict[pair_sym] = 0
                else:
                    get_ancestors_x = [node.gene_name for node in node_x.get_ancestors()]
                    get_ancestors_y = [node_y.gene_name] + [node.gene_name for node in node_y.get_ancestors()]
                    if node_y.gene_name in get_ancestors_x:
                        number = 0
                        node_observed = get_ancestors_x[0]
                        while node_observed != node_y.gene_name:
                            node_observed = get_ancestors_x[number]
                            number += 1
                        if number >= 1:
                            number -= 1
                        inters_dict[pair_sym] = number
                        inters_dict[pair] = number
                    elif node_x.gene_name in get_ancestors_y:
                        number = 0
                        node_observed = get_ancestors_y[0]
                        while node_observed != node_x.gene_name:
                            node_observed = get_ancestors_y[number]
                            number += 1
                        if number >= 1:
                            number -= 1
                        inters_dict[pair_sym] = number
                        inters_dict[pair] = number
                    else:
                        # intersection is nil
                        pass
    return inters_dict

def get_local_lca_reconciliation_cost(leaves_t1, leaves_t2, mappings_df, gene_tree, inters_dict):
    all_lca_leaves_t1 = list(set([mappings_df.loc[_].id_gene for _ in  leaves_t1]))
    all_lca_leaves_t2 = list(set([mappings_df.loc[_].id_gene for _ in  leaves_t2]))
    all_lca_leaves_t = list(set([mappings_df.loc[_].id_gene for _ in  leaves_t1+leaves_t2]))

    s_t1 = None 
    if len(all_lca_leaves_t1) == 1:
        s_t1 = all_lca_leaves_t1[0]
    else:
        s_t1 = gene_tree.get_common_ancestor(all_lca_leaves_t1).gene_name

    s_t2 = None 
    if len(all_lca_leaves_t2) == 1:
        s_t2 = all_lca_leaves_t2[0]
    else:
        s_t2 = gene_tree.get_common_ancestor(all_lca_leaves_t2).gene_name

    s_t = None
    if len(all_lca_leaves_t) == 1:
        s_t = all_lca_leaves_t[0]
    else:
        s_t = gene_tree.get_common_ancestor(all_lca_leaves_t).gene_name

    if all([True if _==s_t else False for _ in [s_t, s_t1, s_t2]]):
        return (1, 'creation')
    elif s_t != s_t1 and s_t == s_t2:
        pair_1 = (s_t, s_t1)
        return (2 + inters_dict[pair_1], 'creation')
    elif s_t == s_t1 and s_t != s_t2:
        pair_2 = (s_t, s_t2)
        return (2 + inters_dict[pair_2], 'creation')
    elif s_t != s_t1 and s_t != s_t2:
        pair_1 = (s_t, s_t1)
        pair_2 = (s_t, s_t2)
        observed_node = gene_tree.search_nodes(gene_name=s_t)[0]
        evol_type = None

        if hasattr(observed_node, 'D'):
            if observed_node.D == 'Y':
                evol_type = 'duplication'
            elif observed_node.D == 'N':
                evol_type = 'speciation'
            else:
                raise ValueError('Inappropriate NHX format for the gene tree given at the input !')
        
        elif hasattr(observed_node, 'DD'):
            if observed_node.DD == 'Y':
                evol_type = 'duplication'
            elif observed_node.DD == 'N':
                evol_type = 'speciation'
            else:
                raise ValueError('Inappropriate NHX format for the gene tree given at the input !')
        else:
            raise ValueError('Inappropriate NHX format for the gene tree given at the input !')
        
        return (inters_dict[pair_1] + inters_dict[pair_2], evol_type)

#'''
def lca_reconciliation_cost(t, mappings_df, gene_tree, inters_dict):
    costs = []
    for node in t.traverse('postorder'):
        if node.is_leaf():
            costs.append(0)
        else:
            children = node.get_children()
            left_leaves = [str(_.name) for _ in children[0].get_leaves()]
            right_leaves = [str(_.name) for _ in children[1].get_leaves()]
            node_cost = get_local_lca_reconciliation_cost(left_leaves, right_leaves, mappings_df, gene_tree, inters_dict)[0]
            costs.append(node_cost)
    return sum(costs)
#'''

def init_costs(t, gene_tree, mappings_df, inters_dict):
    for node in t.traverse('postorder'):
        if node.is_leaf():
            node.add_features(local_cost=0, global_cost=0)
        else:
            children = node.children
            left_leaves = [_.name for _ in  children[0].get_leaves()]
            right_leaves = [_.name for _ in  children[1].get_leaves()]
            local_cost = get_local_lca_reconciliation_cost(left_leaves, right_leaves, mappings_df, gene_tree, inters_dict)[0]
            g_cost = local_cost + children[0].global_cost + children[1].global_cost
            node.add_features(local_cost=local_cost, global_cost=g_cost)
    return t

def map_transcripts_lca(solutions):
    solutions_mapped = []
    for solution in solutions:
        tree = solution
        for node in tree.traverse('preorder'):
            if not node.is_leaf():
                leaves = [str(_.name) for _ in node.get_leaves()]
                leaves.sort()
                node.add_features(t_lca='{}'.format('&'.join(leaves)))
        solutions_mapped.append(solution)
    return solutions_mapped

def compare_transcripts_lca(t1, t2):
    is_chiral = True
    leaves_t1 = [str(_.name) for _ in t1.get_leaves()]
    leaves_t2 = [str(_.name) for _ in t2.get_leaves()]

    leaves_t1.sort()
    leaves_t2.sort()
    if len(leaves_t1) != len(leaves_t2) or leaves_t1 != leaves_t2:
        raise ValueError('Inappropriate leaves!')
    
    for leave in leaves_t1:
        node_t1 = t1&leave
        node_t2 = t2&leave
        anc_t1 = [node.t_lca for node in node_t1.get_ancestors()]
        anc_t2 = [node.t_lca for node in node_t2.get_ancestors()]

        if len(anc_t1) != len(anc_t2) or anc_t1 != anc_t2:
            is_chiral = False
            break
    return is_chiral

def avoid_chirality(solutions):
    if len(solutions) == 1:
        return solutions
    else:
        all_solutions = []
        solutions_mapped = map_transcripts_lca(solutions)
        already_done = []
        for i in range(0, len(solutions_mapped)-1):
            if i not in already_done:
                solution_i = solutions_mapped[i]
                t_i = solution_i
                for j in range(i+1, len(solutions_mapped)):
                    if j not in already_done:
                        solution_j = solutions_mapped[j]
                        t_j = solution_j
                        is_chiral = compare_transcripts_lca(t_i, t_j)
                        if is_chiral:
                            already_done.append(j)
        for i, solution in enumerate(solutions_mapped):
            if i not in already_done:
                all_solutions.append(solution)
    return all_solutions
    
def transcript_tree_construction(ortholog_trees_set, me_tree, matrix, gene_tree, mappings_df, join_input, bool_label, bool_compute_all):

    leaves_me_tree = [node.name for node in me_tree.get_leaves()]
    INTERS_DICT = get_inter(gene_tree)
    if len(leaves_me_tree) == 1:
        transcript_tree = ortholog_trees_set[leaves_me_tree[0]]
        return [init_costs(transcript_tree, gene_tree, mappings_df, INTERS_DICT)]
    else:
        fusion_results = {}
        fusion_ordered_steps = [node for node in me_tree.traverse('postorder') if not node.is_leaf()]

        
        for step, fusion_ordered_step in enumerate(fusion_ordered_steps):
            ortholog_tree_a = None
            ortholog_tree_b = None

            children = fusion_ordered_step.children
            left_child = children[0]
            right_child = children[1]

            if left_child.is_leaf():
                tmp_tree = ortholog_trees_set[left_child.name]
                ortholog_trees_a = [init_costs(tmp_tree, gene_tree, mappings_df, INTERS_DICT)]
            else:
                if bool_compute_all == 0:
                    ortholog_trees_a = fusion_results[left_child]
                elif bool_compute_all == 1:
                    ortholog_trees_a = [fusion_results[left_child][0]]
                else:
                    raise ValueError('Inapproptiate value of -c')

            if right_child.is_leaf():
                tmp_tree = ortholog_trees_set[right_child.name]
                ortholog_trees_b = [init_costs(tmp_tree, gene_tree, mappings_df, INTERS_DICT)]
            else:
                if bool_compute_all == 0:
                    ortholog_trees_b = fusion_results[right_child]
                elif bool_compute_all == 1:
                    ortholog_trees_b = [fusion_results[right_child][0]]
                else:
                    raise ValueError('Inapproptiate value of -c')

            solutions_tmp = []
            for ortholog_tree_a in ortholog_trees_a:
                t1 = fix_locations(0, ortholog_tree_a)
                n_leaves_t1 = len(t1.get_leaves())
                for ortholog_tree_b in ortholog_trees_b:
                    t2 = fix_locations(n_leaves_t1, ortholog_tree_b)

                    graph_of_recurrence_1 = init_graph_of_recurrence(t1)
                    minimum_reconciliation_cost_principle_1 = minCLevol(t1, t2, graph_of_recurrence_1, mappings_df, gene_tree, INTERS_DICT)
                    graph_of_recurrence_1 = minimum_reconciliation_cost_principle_1[1]
                    min_lca_rec_cost_1 = minimum_reconciliation_cost_principle_1[0]

                    graph_of_recurrence_2 = init_graph_of_recurrence(t2)
                    minimum_reconciliation_cost_principle_2 = minCLevol(t2, t1, graph_of_recurrence_2, mappings_df, gene_tree, INTERS_DICT)
                    graph_of_recurrence_2 = minimum_reconciliation_cost_principle_2[1]
                    min_lca_rec_cost_2 = minimum_reconciliation_cost_principle_2[0]

                    plausible_solutions = None
                    if min_lca_rec_cost_1 < min_lca_rec_cost_2:
                        # we only consider min_lca_rec_cost_1
                        plausible_solutions = [backtracking_minCLevol(min_lca_rec_cost_1, graph_of_recurrence_1, t1, t2, matrix)]
                    elif min_lca_rec_cost_1 > min_lca_rec_cost_2:
                        # we only consider min_lca_rec_cost_2
                        plausible_solutions = [backtracking_minCLevol(min_lca_rec_cost_2, graph_of_recurrence_2, t2, t1, matrix)]
                    elif min_lca_rec_cost_1==min_lca_rec_cost_2:
                        # we consider both
                        plausible_solutions = [backtracking_minCLevol(min_lca_rec_cost_1, graph_of_recurrence_1, t1, t2, matrix)] + [backtracking_minCLevol(min_lca_rec_cost_2, graph_of_recurrence_2, t2, t1, matrix)]
                    else:
                        raise ValueError('no solution found !')
                    solutions_tmp.extend(plausible_solutions)
  
            min_rec = min([solution[0][0][0] for solution in solutions_tmp])
            solutions = [solution for solution in solutions_tmp if solution[0][0][0]==min_rec]

            solutions_trees = None
            if join_input == 1:
                minimum_score = max([solution[0][0][2][0] for solution in solutions])
                solutions_trees = [solution for solution in solutions if solution[0][0][2][0]==minimum_score]
            elif join_input == 2:
                maximum_score = min([solution[0][0][2][1] for solution in solutions])
                solutions_trees = [solution for solution in solutions if solution[0][0][2][1]==maximum_score]
            elif join_input == 0:
                mean_score = max([solution[0][0][2][2] for solution in solutions])
                solutions_trees = [solution for solution in solutions if solution[0][0][2][2]==mean_score]
            else:
                raise ValueError('Inappropriate value --join !')

            # compute solutions and update local_lca / global_lca
            exhaustive_solutions = []
            for solution_tree in solutions_trees:
                #print(solution_tree)
                #print('_-----------------------------_')
                exhaustive_solutions.append(compute_solution(bool_label, solution_tree[0][0][1], solution_tree[1][0], solution_tree[1][1], step, solution_tree[0][0][3], mappings_df, gene_tree, INTERS_DICT))

            # avoid chirality
            exhaustive_solutions_without_chirality = avoid_chirality(exhaustive_solutions)
            
            fusion_results[fusion_ordered_step] = exhaustive_solutions_without_chirality 

    trees = fusion_results[me_tree.get_tree_root()]
    solutions_to_return = []
    if bool_label == 1:
        for tree in trees:
            solutions_to_return.append(label_ortholog_tree(tree, gene_tree))
    else:
        solutions_to_return = trees
    return solutions_to_return

def get_orthologs_set_colored(ortholog_trees_set):
    colors_maps = {}
    if len(ortholog_trees_set.keys()) > 100:
        colors = [
        "red", "orange", "yellow", "green", "blue", "indigo", "violet",
        "white", "gray", "brown", "pink", "purple", "cyan", "magenta",
        "amber", "apricot", "aquamarine", "azure", "beige", "chartreuse", 
        "coral", "crimson", "cyan", "fuchsia", "gold", "goldenrod", "green-yellow", 
        "lavender", "lemon", "lime", "magenta-rose", "maroon", "mint", "olive", 
        "peach", "periwinkle", "plum", "rose", "russet", "salmon", "sapphire", 
        "scarlet", "silver", "tan", "teal", "turquoise", "vermilion", 
        "viridian", "wheat", "zaffre", "amaranth", "auburn", "cerulean", 
        "champagne", "cobalt", "coffee", "copper", "cream", "emerald", "garnet", 
        "ginger", "ivory", "jade", "jasmine", "lavender-blush", "lilac", "mahogany", 
        "melon", "midnight-blue", "navy-blue", "ochre", "pear", "rosewood", "ruby", 
        "sepia", "saffron", "smoky-black", "steel-blue", "tangerine", "topaz", 
        "ultramarine", "vanilla", "violet-red", "wisteria", "xanadu", "yellow-green", 
        "zinnwaldite", "alabaster", "almond", "amethyst", "asparagus", "bittersweet",
        "buff", "burgundy", "byzantium", "capri", "carmine", "celadon", "cerise"
        # Add more colors as needed
        ]
        for i, key_o in enumerate(ortholog_trees_set.keys()):
            color = colors[i]
            for transcript in [leaves.name for leaves in ortholog_trees_set[key_o].get_leaves()]:
                colors_maps[transcript] = color
    else:
        already_checked = []
        for i, key_o in enumerate(ortholog_trees_set.keys()):
            color = '#247e68'
            while color in already_checked:
                color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
            already_checked.append(color)
            for transcript in [leaves.name for leaves in ortholog_trees_set[key_o].get_leaves()]:
                colors_maps[transcript] = color

    return colors_maps

def f_export_data(enum_solution, solution_tree, ortholog_trees_set, colors_maps, minimum_evolution_tree, computing_time, enum_tree, output_prefix, output_folder):
    # Create directories
    try:
        #os.mkdir('{}/{}_solution_{}'.format(output_folder, output_prefix, enum_solution))
        #os.mkdir('{}/{}_solution_{}/ortholog_trees'.format(output_folder, output_prefix, enum_solution))
        #os.mkdir('{}/{}_solution_{}/ortholog_trees/nhx'.format(output_folder, output_prefix, enum_solution))
        
        #os.mkdir('{}/{}_solution_{}/dendogram'.format(output_folder, output_prefix, enum_solution))
        
        main_out__folder_path = os.path.join(output_folder, output_prefix)
        os.makedirs(main_out__folder_path, exist_ok=True)
        
        main_folder_path = os.path.join(output_folder , output_prefix, '{}_subtree'.format(enum_tree))
        os.makedirs(main_folder_path, exist_ok=True)
                
        sub_main_folder_path = os.path.join(output_folder , output_prefix, '{}_subtree'.format(enum_tree), 'solution_{}'.format(enum_solution))
        os.makedirs(sub_main_folder_path, exist_ok=True)
        
        main_ortholog_trees_path = os.path.join(output_folder , output_prefix, '{}_subtree'.format(enum_tree), 'solution_{}'.format(enum_solution), 'ortholog_trees')
        os.makedirs(main_ortholog_trees_path, exist_ok=True)
        
        main_nhx_path = os.path.join(output_folder , output_prefix, '{}_subtree'.format(enum_tree), 'solution_{}'.format(enum_solution), 'ortholog_trees', 'nhx')
        os.makedirs(main_nhx_path, exist_ok=True)
        
        main_dendogram_path = os.path.join(output_folder , output_prefix, '{}_subtree'.format(enum_tree), 'solution_{}'.format(enum_solution), 'dendogram')
        os.makedirs(main_dendogram_path, exist_ok=True)
        
        #os.makedirs('{}/{}_solution_{}'.format(output_folder, output_prefix, enum_solution), exist_ok=True)
        #os.makedirs('{}/{}_solution_{}/ortholog_trees'.format(output_folder, output_prefix, enum_solution), exist_ok=True)
        #os.makedirs('{}/{}_solution_{}/ortholog_trees/nhx'.format(output_folder, output_prefix, enum_solution), exist_ok=True)
        #os.makedirs('{}/{}_solution_{}/dendogram'.format(output_folder, output_prefix, enum_solution), exist_ok=True)
        # print(f"Folder '{folder_path}' created or overwritten successfully.")
    except Exception as e:
        print(f"Error creating or overwriting folder '{output_folder}': {str(e)}")


    # get steps
    for ortholog_tree in ortholog_trees_set.keys():
        tree = ortholog_trees_set[ortholog_tree]
        nhx_tree = tree.write(features=['label'])
        root_label = tree.label
        nhx_tree = nhx_tree.replace(';', '1:0[&&NHX:label={}];'.format(root_label))
        path_ortho = os.path.join(output_folder , output_prefix, '{}_subtree'.format(enum_tree), 'solution_{}'.format(enum_solution), 'ortholog_trees', 'nhx', '{}.nhx'.format(ortholog_tree))
        file_ortholog = open(path_ortho, 'w')
        file_ortholog.write(nhx_tree)
        file_ortholog.close()
        ts = viz_transcripts(tree, 'NA', False, colors_maps)
        tree_save_path = os.path.join(output_folder , output_prefix, '{}_subtree'.format(enum_tree), 'solution_{}'.format(enum_solution), 'ortholog_trees', 'nhx', '{}.svg'.format(ortholog_tree))
        tree.render(tree_save_path, w=183, units='mm', tree_style=ts)

    # get time computing
    path_computing_time = os.path.join(output_folder , output_prefix, '{}_subtree'.format(enum_tree), 'solution_{}'.format(enum_solution), 'computing_time.txt')
    file_computing_time = open(path_computing_time, 'w')
    file_computing_time.write(str(computing_time))
    file_computing_time.close()

    # get the ME tree
    path_dendogram_os = os.path.join(output_folder , output_prefix, '{}_subtree'.format(enum_tree), 'solution_{}'.format(enum_solution), 'dendogram', 'me_tree.nwk')
    minimum_evolution_tree.write(outfile=path_dendogram_os, format=5)

    # get the solution
    path_solution_os = os.path.join(output_folder , output_prefix, '{}_subtree'.format(enum_tree), 'solution_{}'.format(enum_solution), 'solution.nhx')
    file_solution = open(path_solution_os, 'w')
    nhx_solution_tree = solution_tree.write(features=['label'])
    solution_root_label = solution_tree.label
    nhx_solution_tree = nhx_solution_tree.replace(';', '1:0[&&NHX:label={}];'.format(solution_root_label))
    file_solution.write(nhx_solution_tree)
    file_solution.close()
    ts_solution = viz_transcripts(solution_tree, computing_time, True, colors_maps)
    path_viz_solution_os = os.path.join(output_folder , output_prefix, '{}_subtree'.format(enum_tree), 'solution_{}'.format(enum_solution), 'solution.svg')
    solution_tree.render(path_viz_solution_os, w=183, units='mm', tree_style=ts_solution)

    return True

def viz_transcripts(tree, computing_time, is_solution, colors):
    if is_solution:
        cost = tree.get_tree_root().global_cost
    ts = TreeStyle()
    ts.show_leaf_name = True
    number_creation = 0
    for n in tree.traverse():
        nstyle = NodeStyle()
        n.dist = 1.0
        type_evol = n.label
        if n.is_leaf():
            nstyle['shape'] = 'sphere'
            nstyle['fgcolor'] = 'gold'
            n.set_style(NodeStyle())
            n.img_style['bgcolor'] = colors[n.name]
        else:
            if type_evol == 'creation':
                nstyle['shape'] = 'sphere'
                nstyle['fgcolor'] = 'green'
                number_creation += 1
            elif type_evol == 'duplication':
                nstyle['shape'] = 'square'
                nstyle['fgcolor'] = 'red'
            elif type_evol == 'speciation':
                nstyle['shape'] = 'circle'
                nstyle['fgcolor'] = 'blue'
            else:
                raise ValueError('Error occured plotting the reconciled transcript tree!')
            n.set_style(nstyle)
    if is_solution:
        ts.title.add_face(TextFace("transcript phylo: cost = {} & number creation events = {} & computing time = {}".format(cost, number_creation, computing_time), fsize=5), column=0)
    return ts

def read_tree(nhx_input):
    try:
        GENE_TREE = Tree(nhx_input, format=1)
        for enum, node in enumerate(GENE_TREE.traverse('preorder')):
            if node.is_root():
                node.add_features(gene_name='root')
            elif not node.is_leaf():
                node.add_features(gene_name='ancestor_{}'.format(enum))
            else:
                # observed node is a leaf
                node.add_features(gene_name=node.name)
        return GENE_TREE
    except:
        raise TypeError("An error occurred while reading the GENE TREE. Please review the NHX TREE format.")

def successful_message(label, time):
    if label == 0:
        message = f"""
    🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟
    🌟                                                                                             🌟
    🌟   L-MinEvolRec was successfully executed! Check the results in the output directory         🌟
    🌟                                                                                             🌟
    🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟
    \n computing time : {time}
    """
    else:
        message = f"""
    🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟
    🌟                                                                                           🌟
    🌟   MinEvolRec was successfully executed! Check the results in the output directory         🌟
    🌟                                                                                           🌟
    🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟
    \n computing time : {time}
    """
    print(message)
    return True

def main_minevolrec(label, clusters, input_tree, input_matrix, input_mappings, join, forest, threshold, compute_solutions, prefix, output_path):
    # Format inputs to dataframes
    MAPPINGS_DF = get_dataframe_structures_from_inputs(clusters, input_mappings)
    #print('bonjour')
    # Read NHX tree and add a label on each internal node
    GENE_TREE = read_tree(input_tree)
    #print('bonsoir')
    # get the guide tree(MINIMUM EVOLUTION TREE) and ortholog trees 
    try:
        MATRIX = pd.read_csv(input_matrix, sep=';').set_index('transcripts')
        ORTHOLOG_TREES_SET = get_ortholog_trees(MAPPINGS_DF, GENE_TREE, MATRIX)
        MINIMUM_EVOLUTION_TREES = get_me_tree(ORTHOLOG_TREES_SET, MATRIX, join, forest, threshold)    
    except:
        raise TypeError('Inappropriate matrix format arg --matx')
    
    # compute transcript phylogenies 
    #print(ORTHOLOG_TREES_SET)
    #for ortholog in ORTHOLOG_TREES_SET:
    #    print(ORTHOLOG_TREES_SET[ortholog])
    #print(MINIMUM_EVOLUTION_TREES)
    #print('bienvenue')
    cpu_times = []
    for enum_tree, me in enumerate(MINIMUM_EVOLUTION_TREES):
        start_computing = time.time()
        all_transcript_phylogenies = transcript_tree_construction(ORTHOLOG_TREES_SET, me, MATRIX, GENE_TREE, MAPPINGS_DF, join, label, compute_solutions)
        stop_computing = time.time()
        CPU_TIME_COMPUTING = stop_computing - start_computing
        cpu_times.append(CPU_TIME_COMPUTING)
        # export data
        try:
            COLORS_MAPS = get_orthologs_set_colored(ORTHOLOG_TREES_SET)
            for enum_solution, solution_tree in enumerate(all_transcript_phylogenies):
                f_export_data(enum_solution, solution_tree, ORTHOLOG_TREES_SET, COLORS_MAPS, me, CPU_TIME_COMPUTING, enum_tree, prefix, output_path)
        except:
            raise('An Error Occurred During Data Export !')
        #print('aurevoir')
    successful_message(label, sum(cpu_times))


#############################
########### MAIN ############
#############################
    

if __name__ == '__main__':
    
    #####################################
    ############ User inputs ############
    #####################################
    parser = build_arguments_parser()
    args = parser.parse_args()
    BOOL_LABEL = int(args.label)
    CLUSTERS_INPUT = args.clusters
    NHX_INPUT = args.nhx
    MAPPINGS_INPUT = args.mappings
    MATRIX_INPUT = args.matrix
    BOOL_FOREST = int(args.forest)
    THRESHOLD_INPUT = float(args.forest_threshold)
    JOIN_INPUT = int(args.join)
    OUTPUT_FOLDER_INPUT = str(args.output)
    PREFIX_OUTPUT = str(args.prefix)
    BOOL_COMPUTE_SOLUTIONS = int(args.compute)
    
    
    ######################################
    ########### MAIN PROGRAM #############
    ######################################
    main_minevolrec(label=BOOL_LABEL, clusters=CLUSTERS_INPUT, input_tree=NHX_INPUT, input_matrix=MATRIX_INPUT, input_mappings=MAPPINGS_INPUT, join=JOIN_INPUT, forest=BOOL_FOREST, threshold=THRESHOLD_INPUT, compute_solutions=BOOL_COMPUTE_SOLUTIONS, prefix=PREFIX_OUTPUT, output_path=OUTPUT_FOLDER_INPUT)
  
    