import pandas as pd
import numpy as np

def process_ccmpredpy(infile_loc):
    df_ccm_pivot = pd.read_csv(infile_loc, ' ', header=None, skipfooter=0)
    df_ccm_pivot = df_ccm_pivot.apply(pd.to_numeric)
    df_ccm_pivot.values[[np.arange(df_ccm_pivot.shape[0])]*2] = np.nan
    df_ccm = df_ccm_pivot.where(np.triu(np.ones(df_ccm_pivot.shape)).astype(np.bool))
    df_ccm = df_ccm.stack().reset_index()
    df_ccm.reset_index(drop=True, inplace=True)
    df_ccm.columns = ['aa1_loc', 'aa2_loc', 'couplings']
    df_ccm['aa1_loc'] = df_ccm['aa1_loc'] + 1
    df_ccm['aa2_loc'] = df_ccm['aa2_loc'] + 1
    assert all(df_ccm['aa1_loc'] < df_ccm['aa2_loc'])
    return df_ccm

def remove_close(df, primary_distance_cutoff):
    df['abs_diff_in_loc'] = df['aa1_loc'] - df['aa2_loc']
    df['abs_diff_in_loc'] = df['abs_diff_in_loc'].abs()
    ##Only consider amino acids separated by greater than some chain distance
    df = df[df['abs_diff_in_loc'] >=primary_distance_cutoff]
    df.reset_index(drop=True, inplace=True)
    return df

def process_contacts_df(df_contacts, primary_distance_cutoff):
    df_contacts.columns = df_contacts.columns.astype(int)
    df_contacts_stack = df_contacts.where(np.triu(np.ones(df_contacts.shape)).astype(np.bool))
    df_contacts_stack = df_contacts_stack.stack().reset_index()
    df_contacts_stack['abs_diff'] = df_contacts_stack['level_0'] - df_contacts_stack['level_1']
    df_contacts_stack['abs_diff'] = df_contacts_stack['abs_diff'].abs()
    df_contacts_stack = df_contacts_stack[df_contacts_stack['abs_diff'] >= primary_distance_cutoff]
    df_contacts_stack.reset_index(drop=True, inplace=True)
    df_contacts_stack.columns = ['aa1_loc', 'aa2_loc', 'distance', 'primary_chain_distance']
    return df_contacts, df_contacts_stack

def merge_contacts_couplings(df_contacts_stack, df_couplings_stack, seq):
    assert all(df_couplings_stack['aa1_loc'] == df_contacts_stack['aa1_loc'])
    assert all(df_couplings_stack['aa2_loc'] == df_contacts_stack['aa2_loc'])
    merged_df = pd.concat([df_contacts_stack, df_couplings_stack[['couplings']]],\
                          axis=1, join_axes=[df_contacts_stack.index])
    merged_df.sort_values('couplings', ascending=False, inplace=True)

    seq_dict = dict([(i+1,j) for i,j in list(enumerate(list(seq)))])
    merged_df['aa1_aa'] = merged_df['aa1_loc'].map(seq_dict)
    merged_df['aa2_aa'] = merged_df['aa2_loc'].map(seq_dict)
    return merged_df

def ppv_from_df(merged_df, number_to_test, length_cutoff=8):
    temp_df = merged_df[:number_to_test]
    tps = temp_df[temp_df['distance']<=length_cutoff]['distance'].count()
    totals = temp_df['distance'].count()
    return  tps/totals, totals
