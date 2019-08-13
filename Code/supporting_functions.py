
###########################################################################################
#Mid-point rooting
###########################################################################################
def MP_root(tree):
    """
    This function implements a mid-point rooting algorithm that is linear with tree
    size as opposed to the behavior of many mid-point rooting algorithms that are polynomial
    by requiring all pair-wise comparisons of terminal-terminal distances (this includes
    biopython's midpoint rooting method, strangely). 

    The only really strange behavior that I'm aware of right now are what to do when the 
    algorithm encounters ties. For now I'm just choosing one at random, where random is 
    actually whoever pops up first in the list (so maybe not fully random).

    Input(s):
    tree - a Bio.Phylo tree object

    Output(s):
    rooted-tree - the initial tree object re-rooted
    """

    initial_bl = tree.total_branch_length()
    initial_term_names = [i.name for i in tree.get_terminals()]
    #Root randomly with an outgroup at a terminal
    tree.root_with_outgroup(tree.get_terminals()[0], outgroup_branch_length=0.0)
    
    #Through some type of bug in Bio.Phylo, new terminals can pop into existence and
    #this is a hack to remove them. I think they come from non-zero root branch lengths
    #but this should be investigated and the code can be cleaned up to flag these cases
    pruned_bls = []
    for terminal in tree.get_terminals():
        if terminal.name not in initial_term_names:
            pruned_bls.append(terminal.branch_length)
            tree.prune(terminal)
            print('Pruned a strange terminal that popped into existence during midpoint rooting')
    
    #I'm not entirely sure how this algorithm works on non-bifurcating trees. Thus, this
    #initial assertion statement. It might be fine, but I'd have to think about it.
    assert tree.is_bifurcating()
    
    #Find which terminal/s is farthest away from my randomly selected root terminal
    depths = tree.depths()
    max_depth = max(list(depths.values()))
    clade1 = [i for i,j in depths.items() if j == max_depth]
    #Idealy this would be a list of length 1. If not, there are multiple terminals
    #that are maximally distant (could be due to polytomies or I suppose random chance?)
    if len(clade1) > 1:
        print('Potential for multiple midpoints. Choosing the first that I encounter')

    #Re-root at that farthest terminal (NOTE: just choosing the first clade in the list
    #that may be longer than length 1)
    tree.root_with_outgroup(clade1[0], outgroup_branch_length=0.0)
    
    #And finally find which terminal is farthest from THAT one to identify the farthest pair
    depths = tree.depths()
    max_depth = max(list(depths.values()))
    clade2 = [i for i,j in depths.items() if j == max_depth]
    #Same constraint/caveat applies as above with regard to "ties"
    if len(clade2) > 1:
        print('Potential for multiple midpoints. Choosing the first that I encounter')
    
    #Given the clade pairs, re-root the tree
    rooted_tree = final_root(clade1[0], clade2[0], depths[clade2[0]], tree)
    
    #Ensuring that I've fully conserved branch length after all these manipulations
    #because I've had problems with gaining / losing owing to what I think are issues
    #in Bio.Phylo that I think I've fully figured out.
    assert np.isclose(initial_bl-np.sum(pruned_bls), rooted_tree.total_branch_length())
    
    return rooted_tree



def rel_time_AJH(tree):
    """
    This closely (exactly?) follows the original implementation but note that the explanation
    of the original algorithm fails to mention what happens with zero length branches which of course
    give zero division errors
    """
    depth_dict = tree.depths(unit_branch_lengths=True)
    for key in tree.get_terminals():
        del depth_dict[key]

    inv_depth_dict = {}
    for key,val in depth_dict.items():
        try:
            inv_depth_dict[val].append(key)
        except KeyError:
            inv_depth_dict[val] = [key]


    for depth in range(max(list(inv_depth_dict.keys())), -1, -1):
        for clade in inv_depth_dict[depth]:
            temp = clade.depths()
            lens = [temp[term]-clade.branch_length for term in clade.get_terminals()]
            for ds_clade in clade.clades:
                ds_lens = [temp[term]-clade.branch_length for term in ds_clade.get_terminals()]
                if np.mean(lens) > 0:
                    ds_clade.rate = np.mean(ds_lens) / np.mean(lens)
                else:
                    ds_clade.rate = 0.
                for all_ds in ds_clade.get_terminals() + ds_clade.get_nonterminals():
                    if all_ds == ds_clade:
                        pass
                    else:
                        all_ds.rate = all_ds.rate*ds_clade.rate
    return tree
