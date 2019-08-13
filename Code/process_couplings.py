import pandas as pd
import numpy as np

def process_ccmpredpy(infile_loc):
    """
    This function will process the output files from CCMPredPy into a nicer looking
    and easier to work with pandas dataframe. 

    Note that the CCMPredPy files must first be processed by the bash script included
    in this code to remove their metadata
    """
    df_ccm_pivot = pd.read_csv(infile_loc, ' ', header=None, skipfooter=0)
    df_ccm_pivot = df_ccm_pivot.apply(pd.to_numeric)
    df_ccm_pivot.values[tuple([np.arange(df_ccm_pivot.shape[0])]*2)] = np.nan
    df_ccm = df_ccm_pivot.where(np.triu(np.ones(df_ccm_pivot.shape)).astype(np.bool))
    df_ccm = df_ccm.stack().reset_index()
    df_ccm.reset_index(drop=True, inplace=True)
    #Renaming the columns to something sensible
    df_ccm.columns = ['aa1_loc', 'aa2_loc', 'couplings']
    #And indexing the amino acids so that their numbering starts at 1
    #There was once a reason why I did this but it eludes me now
    df_ccm['aa1_loc'] = df_ccm['aa1_loc'] + 1
    df_ccm['aa2_loc'] = df_ccm['aa2_loc'] + 1
    assert all(df_ccm['aa1_loc'] < df_ccm['aa2_loc'])
    return df_ccm

def remove_close(df, primary_distance_cutoff):
    """
    This code will take what should ideally be a merged dataframe and cut any rows
    out according to the difference in their aa1_loc and aa2_loc columns. Which is
    to say, it will get rid of comparisons between amino acids that are too close
    to each other in primary sequence space since these are less interesting.

    Common parameters are probably 6 ish to isolate medium to long range interactions
    and 12 to isolate long range ones in particular. 
    """
    df['abs_diff_in_loc'] = df['aa1_loc'] - df['aa2_loc']
    df['abs_diff_in_loc'] = df['abs_diff_in_loc'].abs()
    ##Only consider amino acids separated by greater than some chain distance
    df = df[df['abs_diff_in_loc'] >=primary_distance_cutoff]
    df.reset_index(drop=True, inplace=True)
    return df

def process_contacts_df(df_contacts):
    """
    This is for processing the contacts dataframe and making the column names match up to the
    coupling dataframe for subsequent merging.
    """
    df_contacts.columns = df_contacts.columns.astype(int)
    df_contacts_stack = df_contacts.where(np.triu(np.ones(df_contacts.shape)).astype(np.bool))
    df_contacts_stack = df_contacts_stack.stack().reset_index()
    df_contacts_stack['abs_diff'] = df_contacts_stack['level_0'] - df_contacts_stack['level_1']
    df_contacts_stack['abs_diff'] = df_contacts_stack['abs_diff'].abs()
    df_contacts_stack = df_contacts_stack[df_contacts_stack['abs_diff'] >= 1]
    df_contacts_stack.reset_index(drop=True, inplace=True)
    df_contacts_stack.columns = ['aa1_loc', 'aa2_loc', 'distance', 'primary_chain_distance']
    return df_contacts, df_contacts_stack

def merge_contacts_couplings(df_contacts_stack, df_couplings_stack, seq):
    """
    This will merge the contact and coupling dataframes. The sequence makes things
    a bit easier so I request it as well.
    """
    assert all(df_couplings_stack['aa1_loc'] == df_contacts_stack['aa1_loc'])
    assert all(df_couplings_stack['aa2_loc'] == df_contacts_stack['aa2_loc'])
    merged_df = pd.concat([df_contacts_stack, df_couplings_stack[['couplings']]],\
                          axis=1, join_axes=[df_contacts_stack.index])
    #Note that this is critical for other code! The rows are sorted according
    #to their coupling scores!
    merged_df.sort_values('couplings', ascending=False, inplace=True)

    seq_dict = dict([(i+1,j) for i,j in list(enumerate(list(seq)))])
    merged_df['aa1_aa'] = merged_df['aa1_loc'].map(seq_dict)
    merged_df['aa2_aa'] = merged_df['aa2_loc'].map(seq_dict)
    return merged_df

def ppv_from_df(merged_df, number_to_test, length_cutoff=8):
    """
    Pretty straightforward calculation of positive predictive value given
    the merged_df input, n to test, and what to define as a contact (in
    units of angstroms)
    """
    temp_df = merged_df[:number_to_test]
    tps = temp_df[temp_df['distance']<=length_cutoff]['distance'].count()
    totals = temp_df['distance'].count()
    return  tps/totals, totals
