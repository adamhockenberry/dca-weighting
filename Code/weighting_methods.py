from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import random

"""
Possible weighting schemes to add to this:
    1. Krogh-Mitchison max-ent method (it's slow/scales terribly but perhaps useful)
    2. A re-scaled ACL method to correct the problem (quirk?) of short branches close to
        the root getting an outsized influence

Finally, it should be noted to any users that a major problem with all of these methods is the
effective number of sequences. Which is to say, it's not clear at all what any of the weights
should sum to in any of these methods. Ideally, a method would return weights on a 0-1 scale and
the sum of all weights would be less than the number of terminal clades in a tree (unless the 
tree is a star phylogeny in which case they would be equal.
"""



###################
#GSC weights
###################
def calc_GSC_weights(rooted_tree):
    """
    This is a fairly straightforward and fast implementation of the Gerstein-Sonnhammer-Chothia
    weighting algorithm (GSC). I tried a "leaf first" method at first but found this "root first" method 
    (inspired from code I saw in PyCogent) to be so much cleaner/easier to work with for my purposes.

    Input(s):
    rooted_tree - a Bio.Phylo tree object that should be properly rooted for the numbers to make any sense
    
    Output(s):
    weights_dict - a dictionary object with key:val pairs that are Bio.Phylo clade objects (specifically 
            terminals)(keys) and floats that represent the weight (vals)
    """
    #Two simple assertions that will ensure that everything goes smoothly
    assert rooted_tree.rooted
    assert rooted_tree.is_bifurcating()
    
    #Instantiate the empty weights dict
    weights_dict = {}
    #Every clade starts with their own branch length as a weight
    for node in rooted_tree.get_terminals() + rooted_tree.get_nonterminals():
        weights_dict[node] = node.branch_length
    
    ########################
    #Call recursive function
    ########################
    weights_dict, finished = recursive_GSC_fxn(rooted_tree.root, weights_dict, [])
    
    #Remove non-terminal nodes from the weights_dict
    for i in rooted_tree.get_nonterminals():
        del weights_dict[i]
    #Ensure that the sum of the calculated weights are approximately equal to the total branch length
    assert np.isclose(np.sum(list(weights_dict.values())), rooted_tree.total_branch_length())
    return weights_dict

def recursive_GSC_fxn(node, weights_dict, finished):
    """
    This is the update function for GSC that will get recursively called (using
    a depth first search strategy) to compute GSC weights for all possible terminals
    in the tree. 
    
    Input(s):
    node - the current Bio.Phylo clade object
    weights_dict - a dictionary keeping track of all weights where key:val pairs are
            Bio.Phylo clade objects and the current weight for each given clade (float)
    finished - not entirely necessary, but I like to keep track of all the nodes that I've seen
            so this is a list containing all clades that have been rooted/tested

    Output(s):
    weights_dict - the updated weights_dicty object described above
    finished - the updated finished list described above
    """
    if len(node.clades) == 2:#Everything should be bifurcating
        lclade = node.clades[0]
        rclade = node.clades[1]
        if node.branch_length is not None:
            #Calculate the downstream branch length of each daughter clade
            l_tot = lclade.total_branch_length()
            r_tot = rclade.total_branch_length()
            if l_tot + r_tot > 0.:
                #Partition the current branch length in a weighted fashion to the daughters
                weights_dict[lclade] += weights_dict[node]*(l_tot/(l_tot + r_tot))
                weights_dict[rclade] += weights_dict[node]*(r_tot/(l_tot + r_tot))
            else:
                #In the special case of two zero branch length daughters, divide current weight
                #evenly amongst them
                weights_dict[lclade] += weights_dict[node]*0.5
                weights_dict[rclade] += weights_dict[node]*0.5
        #Call recursive function. This is depth first because we're always going down the lclades first
        weights_dict, finished = recursive_GSC_fxn(lclade, weights_dict, finished)
        weights_dict, finished = recursive_GSC_fxn(rclade, weights_dict, finished)
        finished.append(node)
    elif len(node.clades) == 0:#Am at a terminal node
        finished.append(node)
    return weights_dict, finished


def normalize_GSC_weights(weights_dict, rooted_tree):
    """
    This is a novel (I believe) normalization scheme that I came up with to re-scale 
    GSC weights. It simply divides the weight for each terminal by the overall depth of that 
    terminal and therefore can/should(?) be thought of as a % independence.

    Input(s):
    weights_dict - weights dictionary, output from "calc_GSC_weights"
    rooted_tree - a Bio.Phylo tree object that was used as input for "calc_GSC_weights"
    
    Output(s):
    normed_weights_dict - key:val pairs are Bio.Phylo clade objects (same as in the 
            weights_dict) (keys) and a single number that represents the normalized weight.
    """
    normed_weights_dict = {}
    depths_dict = rooted_tree.depths()
    for term in rooted_tree.get_terminals():
        if weights_dict[term] > 0:
            normed_weights_dict[term] = weights_dict[term]/depths_dict[term]
        else:
            normed_weights_dict[term] = 0. 

    return normed_weights_dict


###################
#ACL weights
###################
def calc_ACL_weights(rooted_tree):
    """
    This is an implemntation of  Altschul, Carroll and Lipman weighting (ACL) that represents the "current"
    flowing out from each terminal leaf in a tree assuming that a power source was plugged into the root. 
    Tree branches provide resistance, in the analogy. This implementation works well, and is pretty fast, 
    but I couldn't figure out how to avoid a matrix inversion which should supposedly be possible using 
    Felsenstein's pruning algorithm. Thus, np.linalg.inv makes a clunky showing with all it's associated problems.

    NOTE: I would really like to think of a way to normalize these values as a % independence-like metric
    
    Input(s):
    rooted_tree - a Bio.Phylo tree object that must be (intelligently) rooted 

    Output(s):
    weights_dict - the weights dictionary with key:val pairs of Bio.Phylo clade objects for each terminal (keys)
            and floats representing the weight (vals)
    rooted_tree - the rooted_tree object (returned because it's possible I have pruned zero branch length clades 
            and/or altered the length of the root clade
    """
    assert rooted_tree.is_bifurcating()
    assert rooted_tree.rooted
    
    bl_list = [term.branch_length for term in rooted_tree.get_terminals()]+\
            [internal.branch_length for internal in rooted_tree.get_nonterminals()]
    bl_array = np.array(bl_list)
    min_bl = np.min(bl_array[np.nonzero(bl_array)])
    if min_bl > 10e-7:
        min_bl = 10e-9
    else:
        min_bl = min_bl/100.
    #Remove zero branch length terminal clades entirely from the tree
    for term in rooted_tree.get_terminals():
        if term.branch_length == 0.0:
            term.branch_length = min_bl
    #zero_bl_clades = [term for term in rooted_tree.get_terminals() if term.branch_length==0.]
    #if len(zero_bl_clades) != 0:
    #    rooted_tree = trim_zero_bls(rooted_tree)
    #    print('Found some zero branch length terminals and have removed them. Note that '
    #            'returned tree will contain fewer taxa than original')
    
    if rooted_tree.root.branch_length:
        print('The passed tree contains a non-zero branch length root. This has been removed and'
                'at the end of this algorithm will be set to "None"')
    #Set root branch length to zero
    rooted_tree.root.branch_length = 0.
    #Get initial terminal order
    initial_order = rooted_tree.get_terminals()
    #Set up an empty matrix
    initial_matrix = np.zeros((len(initial_order),len(initial_order)))
    
    ##############################
    #Calling the recursive function
    ##############################
    vcv_matrix, finished_list = vcv_recursive(rooted_tree.root, initial_matrix, [])
    
    #And cleaning things up
    inv_vcv_matrix = np.linalg.inv(vcv_matrix)
    inv_weights = inv_vcv_matrix.sum(axis=1)/inv_vcv_matrix.sum()
    #Reset root branch length to None (as it should have beeen)
    rooted_tree.root.branch_length = None
    weights_dict = dict(zip(initial_order, inv_weights))
    return weights_dict, rooted_tree

def vcv_recursive(putative_root, vcv_matrix, finished):
    """
    This is a fast recursive function to calculate the variance co-variance matrix. 
    
    Input(s):
    putative_root - the root node (initially) and subsequently, other clades from a Bio.Phylo tree object
    vcv_matrix - the matrix which will store everything. Needs to be n x n where n is the 
            total number of terminals in the tree object where putative_root belongs. Should be zeroes to start
    finished - a list of all called nodes to keep track of everything for record keeping

    Output(s):
    vcv_matrix - the updated vcv_matrix described above
    finished - the updated finished list described above
    """
    terminals = putative_root.get_terminals()
    if not set(terminals).issubset(set(finished)):
        vcv_matrix[len(finished):len(finished)+len(terminals), len(finished):len(finished)+len(terminals)] += putative_root.branch_length
    if len(putative_root.clades) == 2:
            vcv_matrix, finished = vcv_recursive(putative_root.clades[0], vcv_matrix, finished)
            vcv_matrix, finished = vcv_recursive(putative_root.clades[1], vcv_matrix, finished)
    elif len(putative_root.clades) == 0:
        finished.append(putative_root)
    return vcv_matrix, finished

def trim_zero_bls(my_tree):
    """
    Certain weighting schemes (specifically ACL) really do not like zero branch length terminals so this code
    prunes them kind of randomly. Am still not entirely sure why this must be the case but c'est la vie.

    Input(s):
    my_tree - a Bio.Phylo tree object

    Output(s):
    my_tree - the same tree object with zero branch length terminals removed
    """
    zero_bls = [term for term in my_tree.get_terminals() if term.branch_length==0.]
    while len(zero_bls)>0:
        my_tree.prune(random.choice(zero_bls))
        zero_bls = [term for term in my_tree.get_terminals() if term.branch_length==0.]
    return my_tree

###################
#HH weights
###################
def calc_HH_weights(records_list):
    """
    This is a modified Henikoff and Henikoff algorithm that I developed to correct
    for gapped sequences. The entire thing proceeds as usual but instead of either
    ignoring gaps or treating them as a 21st character, I give all gaps a weight of
    zero and therefore downweight each column of the alignment according to whatever
    it's weight *would* be times the number of ungapped positions/all positions.

    Essentially, if a column has a ton of gaps those gaps all get zero and also the
    remaining characters get heavily downweighted as well because, in my thinking, the
    column is largely uninformative.

    Depending on the application, the normalization at the bottom is pretty critical and
    should be considered (divide all values by the max?, the mean?, etc.)
    
    Input(s):
    records_list - a list of Bio.Seq objects that should be aligned sequences

    Output(s):
    weights_dict - a dictionary of key:val pairs where key is the id (from record.id) 
            of sequences and the value is a float of the weight assigned to said
            sequence
    """
    seqs = np.array([list(record.seq) for record in records_list])
    ids = [record.id for record in records_list]
    seqs_T = seqs.T

    weights_T = []
    all_weights = []
    for i in seqs_T[:]:
        counter_dict = Counter(i)
        #Comment/delete the below line to treat gaps as a 21st character
        #However I leave them in there and then use them to downweight
        #the entire column as described in the docs
        del counter_dict['-']
        
        r = len(counter_dict.keys())
        weights_dict = {}
        for key, val in counter_dict.items():
            weights_dict[key] = 1./(r*val)
        
        #Note that implicitly here gapped sequences will get zero weight if gaps
        #were deleted from the counter_dict above
        temp_array = np.zeros(i.shape)
        for key, val in weights_dict.items():
            np.place(temp_array, i==key, [val])
        #Rescale to punish gapped columns    
        positions = np.sum(list(counter_dict.values())) #number of non-gapped positions
        temp_array = temp_array * (positions/seqs_T.shape[1])
        #Append
        weights_T.append(temp_array)
    weights_T = np.array(weights_T)
    all_weights = weights_T.T
    all_weights = np.sum(all_weights, axis=1)
    weights_dict = dict(zip(ids, all_weights))
    return weights_dict


#def GSC_adhock(rooted_tree, output_all=False):
#    """
#    DEPRECATED TIP FIRST ALGORITHM THAT IS TOO COMPLICATED

#    This is a fairly straightforward and fast implementation of the Gerstein-Sonnhammer-Chothia
#    weighting algorithm (GSC). Am sure it could be simplified and probably sped up but it scales
#    well and works.
#    
#    Input(s):
#    rooted_tree - a Bio.Phylo tree object that should be properly rooted for the numbers to make any sense
#    output_all - a boolean that when toggled to True (or really any non-false value) will return two
#            dictionaries. The basic weights_dict described below and an optional full_weights_dict
#            that keeps track of the weight of each node at each depth and is useful for certain applications
#            that require repeated calculation of weights for various root points on a tree
#    
#    Output(s):
#    weights_dict - a dictionary object with key:val pairs that are Bio.Phylo clade objects (specifically 
#            terminals)(keys) and floats that represent the weight (vals)
#    full_weights_dict(optional) - a dictionary object with key:val pairs that are Bio.Phylo clade objects
#            (specifically terminals)(keys) and lists of the weight assigned at each depth (vals). The final
#            entry in the list is the correct weight assigned to each terminal. 
#    """
#    assert rooted_tree.rooted
#    assert rooted_tree.is_bifurcating()
#    
#    #Initialize the weights dictionary with zeros in the starting list for each terminal
#    full_weights_dict = {}
#    for terminal in rooted_tree.get_terminals():
#        full_weights_dict[terminal] = [0.0]
#
#    #This can probably be cleaned up to just specify that the root shouldn't have any branch length
#    #because I can't figure out why on earth it ever would
#    #This code just gets the depths of each terminal
#    if rooted_tree.root.branch_length:
#        rooted_tree.root.branch_length = None
#        print("The root appeared to have branch length which was removed during the"
#                "calculation of GSC weights")
#    #Get the depths of all nodes in unit branch lengths
#    all_depths_dict = rooted_tree.depths(unit_branch_lengths=True)
#    #Convert this dictionary into a "reverse" dictionary where depth is the key
#    #and the value is a list of all nodes at that particular depth
#    all_depths_reverse_dict = defaultdict(list)
#    for key, val in all_depths_dict.items():
#        all_depths_reverse_dict[val].append(key)
#    #####################
#    #On to the algorithm.
#    #####################
#    #Start with the deepest nodes and work backwards (i.e. leaf to root)
#    for i in range(max(all_depths_reverse_dict.keys()), 0, -1):
#        #First copy over the previous value for each node
#        for key,val in full_weights_dict.items():
#            full_weights_dict[key].append(val[-1])
#        #Now iterate through the clades at that particular depth
#        for clade in all_depths_reverse_dict[i]:
#            #Get all the downstream terminals
#            try:
#                aterms = clade.clades[0].get_terminals()
#                bterms = clade.clades[1].get_terminals()
#            #This error will be raised for terminal nodes who subsequently get all
#            #their branch length without the need for any further calculations
#            except IndexError:
#                full_weights_dict[clade][-1] += clade.branch_length
#                continue
#            #But for NON-terminal nodes, the branch length must be divided between 
#            #downstream terminals in a weghted manner
#            a = [full_weights_dict[term][-1] for term in aterms]
#            b = [full_weights_dict[term][-1] for term in bterms]
#            tot_weights = np.sum(a+b)
#            
#            if tot_weights != 0:
#                #Divide the branch length according to the existing weight of each downstream terminal
#                for term in aterms + bterms:
#                    full_weights_dict[term][-1] += ((full_weights_dict[term][-1]/tot_weights)*clade.branch_length)
#            #If the weight of the downstream terminals is all zero, divide the branch length evenly
#            #between them
#            else:
#                for term in aterms + bterms:
#                    full_weights_dict[term][-1] += clade.branch_length/len(aterms+bterms)
#    #Make a new dictionary of just the final weights
#    weights_dict = {}
#    for key, val in full_weights_dict.items():
#        weights_dict[key] = val[-1]
#    #Output both dictionaries if requested
#    if output_all:
#        return weights_dict, full_weights_dict
#    #Otherwise just return the basic weights_dict with the final values
#    else:
#        return weights_dict
#
