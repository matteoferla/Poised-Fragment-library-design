import pandas as pd

def load_n_clean(filename):
    df = pd.read_pickle(filename)
    df['has_ring'] = df.mol.apply(AllChem.CalcNumRings).astype(bool)
    return df1.loc[df1.has_ring].copy()

def combine_safely(dfs, key, dtype=int):
    return pd.DataFrame({f'{key}{i}': df.groupby('inchi')[key].sum() for i, df in enumerate(dfs)})\
               .fillna(0).astype(dtype).sum(axis=1)

dfs = [load_n_clean(fn) for fn in ('Enamine_REAL_HAC11_1Msub_synthons.pkl.gz',
                                   'Enamine_REAL_1M_random_subsample2_synthons.pkl.gz',
                                   'LC_Stock_HTS.pkl.gz',
                                   'Mcule_ultimate_1Mchunk0.pkl.gz',
                                  )]

df = pd.DataFrame(dict(tally=combine_safely(dfs, 'counts'),
                  USRCAT=pd.concat(dfs).drop_duplicates('inchi').set_index('inchi').USRCAT
                 )
            )
df.to_pickle('common_synthons.pkl.gz')