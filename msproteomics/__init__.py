# coding=utf-8
# Copyright 2025 Thang V Pham
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
With contributions from Truong Xuan Nam and others at Thuy Loi University

Do Duy Luc, Tran Thanh Kien: MaxLFQ in Python
Tuong Dang Vuong Quoc, Nguyen Gia Bao, Truong Tuan Dung: C++ integration

"""

__author__ = 'Thang V Pham'

import numpy as np
import pandas as pd
from pandas.core.dtypes.common import is_numeric_dtype
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm
import warnings
import re


def read(file_path,
         primary_id_col = "PG.ProteinGroups",
         secondary_id_cols = np.array(["EG.ModifiedSequence", "FG.Charge", "F.FrgIon", "F.Charge"]),
         sample_id_col = "R.Condition",
         intensity_col = "F.PeakArea",
         other_cols = np.array(["F.ExcludedFromQuantification", "PG.Qvalue", "EG.Qvalue"])):
    
    cols = np.concatenate(([primary_id_col],
                            secondary_id_cols,
                            [sample_id_col],
                            [intensity_col],
                            other_cols),
                            axis=0)
    return pd.read_csv(file_path,
                       delimiter="\t",
                       usecols=cols)


def preprocess(quant_table,
               primary_id_col = "PG.ProteinGroups",
               secondary_id_cols = np.array(["EG.ModifiedSequence", "FG.Charge", "F.FrgIon", "F.Charge"]),
               sample_id_col = "R.Condition",
               intensity_col = "F.PeakArea",
               median_normalization = True,
               log2_intensity_cutoff = 0,
               pdf_out = "qc-plots.pdf",
               pdf_width = 12,
               pdf_height = 8):

    if isinstance(quant_table, pd.DataFrame):
        if not is_numeric_dtype(quant_table[intensity_col]):
            raise TypeError("Intensity column must be numeric")

        print("Concatenating secondary ids...")
        second_id = quant_table[secondary_id_cols[0]].astype(str)
        for col in range(1, len(secondary_id_cols)):
            second_id += quant_table[secondary_id_cols[col]].astype(str)
        df = pd.DataFrame({'protein_list': quant_table[primary_id_col],
                           'sample_list': quant_table[sample_id_col],
                           #'quant': np.log2(quant_table[intensity_col]),
                           'quant': np.log2(quant_table[intensity_col].values, out = np.full_like(quant_table[intensity_col], np.nan), where = quant_table[intensity_col] > 0),
                           'id': second_id})
        df.dropna(axis=0)
        # intensity cut off
        figs = []
        if log2_intensity_cutoff is not None:
            print('Removing low intensities...')
            if pdf_out is not None:
                # histogram
                fig1 = plt.figure(figsize=(pdf_width, pdf_height))
                n, bins, patch = plt.hist(df['quant'], bins=50, density=True)
                plt.ylabel('Density')
                plt.xlabel('log2 intensity')
                plt.title('Histogram of log2 intensities')
                #plt.arrow(log2_intensity_cutoff, 0, log2_intensity_cutoff, max(n)/2, color='r')
                plt.annotate('Cutoff',
                              xy =(log2_intensity_cutoff, 0),
                              xytext =(log2_intensity_cutoff, max(n)/2), 
                              arrowprops = dict(width = 0.001,
                              color = 'r',
                              shrink = 0.0),)
                figs.append(fig1)
            df = df[df["quant"] > log2_intensity_cutoff]
        samples = df['sample_list'].unique()

        if pdf_out is not None:
            dl = []
            m = []
            for index, sample in enumerate(samples):
                dl.append(df.loc[df['sample_list'] == sample, 'quant'])
                m.append(np.median(dl[index]))
            print("Barplotting raw data ...")
            fig2 = plt.figure(figsize=(pdf_width, pdf_height))
            y_pos = np.arange(1, len(samples) + 1)
            plt.boxplot(dl,
                        flierprops=dict(marker='o', markerfacecolor='blue', markersize=1, markeredgecolor='none'))
            plt.xticks(y_pos, samples, rotation=90)
            figs.append(fig2)

        if median_normalization is True:
            print("Median normalization ...")
            for sample, f in zip(samples, np.mean(m) - m):
                df.loc[df['sample_list'] == sample, 'quant'] += f
            if pdf_out is not None:
                dl = []
                m = []
                for index, sample in enumerate(samples):
                    dl.append(df.loc[df['sample_list'] == sample, 'quant'])
                    m.append(np.median(dl[index]))
                fig3 = plt.figure(figsize=(pdf_width, pdf_height))
                y_pos = np.arange(1, len(samples) + 1)
                plt.boxplot(dl,
                            flierprops=dict(marker='o', markerfacecolor='green', markersize=1, markeredgecolor='none'))
                plt.xticks(y_pos, samples, rotation=90)
                figs.append(fig3)

        with PdfPages(pdf_out) as pdf:
            for fig in figs:
                plt.figure(fig)
                pdf.savefig()
        return df
    else:
        raise TypeError("quant_table isn't pd.Dataframe")


def create_protein_list(preprocessed_data):
    if isinstance(preprocessed_data, pd.DataFrame):
        if any(pd.isna(preprocessed_data['protein_list'])):
            raise Exception("NA value in protein_list")
            
        if any(pd.isna(preprocessed_data['sample_list'])):
            raise Exception("NA value in sample_list")
            
        if any(pd.isna(preprocessed_data['id'])):
            raise Exception("NA value in id")
            
        if any(pd.isna(preprocessed_data['quant'])):
            raise Exception("NA value in quant") 
            
        proteins = preprocessed_data['protein_list'].unique()
        samples = preprocessed_data['sample_list'].unique()
        print("Create quantification list..")
        print("# rows = {0}, # columns = {1}".format(proteins.shape[0], samples.shape[0]))
        
        p_list = {}
        # progress display
        threes_display = 0
        filled = 0
        step = proteins.shape[0] / 20
        for i in range(proteins.shape[0]):
            if i >= threes_display-1:
                print('\r[{:}] {:.0%}'.format('#'*filled + ' '*(20-filled), i/proteins.shape[0]), end = '')
                threes_display+=step
                filled+=1
        
            tmp = preprocessed_data[preprocessed_data["protein_list"] == proteins[i]]
            if tmp.shape[0] > 0:
                dupl = tmp[['id', 'sample_list']].duplicated()
                if dupl.any():
                    for _sample, _id in tmp.loc[dupl, ['id', 'sample_list']]:
                        print("sample {0}; id {1} not unique.".format(_sample, _id))
                    raise Exception("Duplicate entry")
                else:
                    m = pd.DataFrame(columns=samples, index=tmp['id'].unique())
                    for j in tmp.index:
                        m.loc[tmp.at[j, 'id'], tmp.at[j, 'sample_list']] = tmp.at[j, 'quant']
                    p_list[proteins[i]] = m
        print('\r[{:}] 100% Complete.'.format('#'*20))
        
        return p_list
    else:
        raise TypeError("preprocessed_data isn't pd.Dataframe")


def maxLFQ(X):
    X = np.array(X, dtype=np.float64)
    
    # check value None in data
    if np.all(np.isnan(X)):
        return dict({"estimate": None, "annotation": "NA"})
    # check number row
    if X.shape[0] == 1:
        return dict({"estimate": X[0, ], "annotation": ""})
    
    N = X.shape[1]  # [row, col]
    cc = 0
    g = np.full(N, 0, dtype = np.uint)
    
    def spread(i):
        g[i] = cc
        for r in range(X.shape[0]):
            if not np.isnan(X[r, i]):
                for k in range(X.shape[1]):
                    if (not np.isnan(X[r, k])) and (g[k] == 0):
                        spread(k)
    
    def nanmedian2(s):
        if np.all(np.isnan(s)):
            return np.nan
        return np.nanmedian(s)

    # MaxLFQ
    def maxLFQdo(X):

        X = np.array(X)
        Ncol = X.shape[1]
        AtA = np.zeros((Ncol, Ncol))
        Atb = np.zeros(Ncol)

        for i in range(Ncol - 1):
            for j in range(i + 1, Ncol):
                r_i_j = nanmedian2(np.array(-X[:, i] + X[:, j]))
                if not np.isnan(r_i_j):
                    AtA[i, j] = AtA[j, i] = -1
                    AtA[i, i] = AtA[i, i] + 1
                    AtA[j, j] = AtA[j, j] + 1

                    Atb[i] = Atb[i] - r_i_j
                    Atb[j] = Atb[j] + r_i_j

        l = np.append(np.append(2*AtA, np.ones((Ncol, 1)), axis=1),
                      np.append(np.ones(Ncol), 0).reshape(1, Ncol+1),
                      axis=0)        

        r = np.append(2*Atb,
                      [np.nanmean(X)*Ncol],
                      axis=0).reshape((Ncol+1, 1))
        
        x = np.linalg.solve(l, r)
        
        return x.flatten()[:Ncol]

    for i in range(N):
        if g[i] == 0:
            cc += 1
            spread(i)

    w = np.full(N, np.nan)
    for i in range(cc):
        ind = np.array(g == (i + 1))
        if sum(ind) == 1:            
            w[ind] = nanmedian2(np.array(X[:,ind]))            
        else:
            w[ind] = maxLFQdo(X[:,ind])        

    if np.all(np.isnan(w)):
        return dict({"estimate": w, "annotation": "NA"})
    else:
        g_quantified_samples = g[~np.isnan(w)]        
        if np.all([g_quantified_samples[0] == x for x in g_quantified_samples]):                        
            return dict({"estimate": w, "annotation": ""})
        else:        
            g[np.isnan(w)] = 0
            return dict({"estimate": w, "annotation": ";".join(['NA' if (x == 0) else str(x) for x in g])})

          
def create_protein_table(protein_list, method = "maxlfq"):
    if not isinstance(protein_list, dict):
        raise TypeError("Only dict are allowed")

    if len(protein_list) == 0:
        return None

    tab = pd.DataFrame(None, columns=list(protein_list.values())[0].columns, index=list(protein_list))
    annotation = pd.DataFrame('', columns= ['MaxLFQ_annotation'], index=list(protein_list))
        
    # progress display
    nrow = tab.shape[0]
    threes_display = 0
    filled = 0
    step = nrow / 20
    for i in range(nrow):
        if i >= threes_display-1:
            print('\r[{:}] {:.0%}'.format('#'*filled + ' '*(20-filled), i/nrow), end = '')
            threes_display += step
            filled += 1
    
        if method == "maxlfq":
            out = maxLFQ(list(protein_list.values())[i])
    
        else:
            raise Exception("Unknown method: ", method)

        tab.iloc[i, :] = out['estimate']
        annotation.iloc[i, 0] = out['annotation']

    print('\r[{:}] 100% Complete.'.format('#'*20))
    
    return dict({"estimate": tab, "annotation": annotation})


def plot_protein(X, 
                 main = "", 
                 col = [], 
                 split = 0.6):
    #Input X is a matrix
    if col == []:
        for i in range(1, X.shape[0]):
            col.append(i)

    plt.close('all')
    X = X.to_numpy()
    a = []
    j = 0
    x = np.arange(len(X))
    ys = [i+x+(i*x)**2 for i in range(len(X))]
    colors = cm.rainbow(np.linspace(0, 1, len(ys)))
    for y, c in zip(ys, colors):
        a.append(c)
    fig4 = plt.figure()
    plt.ylabel('Intensity')
    plt.xlabel('Sample')
    plt.title(main)
    plt.xticks(np.linspace(1,24,24))
    plt.yticks(np.linspace(0,int(np.nanmax(X)),9))
    for i in range(len(X)):
        plt.plot(np.linspace(1,24,24),X[i],color = a[j],marker = 'o')
        j += 1
    plt.show() 


def extract_annotation(protein_ids, 
                       quant_table,
                       primary_id_col = "PG.ProteinGroups",
                       annotation_columns = None):
    
    #concatenate input columns
    all_columns = np.concatenate(([primary_id_col], annotation_columns), axis = 0) 

    #get column name of input data 
    colnames, index = [], []
    for a, b in quant_table.items():
        colnames.append(a)

    #check all_columns in data ?
    for i in range(len(all_columns)):
        if(all_columns[i] in colnames):
            index.append(colnames.index(all_columns[i]))
        else:
            index.append(None)
    
    #if exist value NaN then break
    if -1 in index:
        raise Exception("The input table has no column")

    for i in range(len(protein_ids)):
        if(protein_ids[i] in quant_table[primary_id_col]):
            index.append(quant_table[primary_id_col].index(protein_ids[i]))
        else:
            index.append(None)

    if -1 in index:
        raise Exception("Cannot find")

    #get rows index of column in all_columns    
    ind = [quant_table[primary_id_col].index(x) for x in protein_ids]
    tab = quant_table[np.concatenate(([primary_id_col], annotation_columns), axis = 0)].iloc[ind]
    
    return(tab)


def create_site_key(protein_ids, ptm_locations):
    n = len(protein_ids)
    
    all_sites = [''] * n   
    first_sites = [''] * n
    
    for i in range(n):
        
        pp = protein_ids[i].split(';')
        ss = ptm_locations[i].split(';')
                
        if len(pp) > 0 and len(pp) == len(ss):
            a = [''] * len(pp)
            for j in range(len(pp)):
                p = pp[j]
                s = ss[j]
                
                s2 = re.findall(r'\(.*?\)', s) 
                
                b = [''] * len(s2) 
                
                for k in range(len(s2)):
                    r = re.findall( '[STY][0-9]+', s2[k]) 
                    if len(r) > 0:
                        tmp = [p + '_' + rr + '_M' + str(min(len(r), 3)) for rr in r]
                        b[k] = ';'.join(tmp)                                    
                a[j] = ';'.join(b)
                                    
            first = True
            for j in range(len(pp)):
                if pp[j] != '' and a[j] != '':
                    if first:
                        all_sites[i] = a[j]
                        first_sites[i] = a[j]
                        first = False
                    else:
                        all_sites[i] = all_sites[i] + ';' + a[j]
                        first_sites[i] = first_sites[i] + ';' + a[j]
            
    unique_sites = {y for x in all_sites for y in x.split(';')}
    
    unique_first_sites = {y for x in first_sites for y in x.split(';')}

    return(all_sites, unique_sites, first_sites, unique_first_sites)


def create_site_report_longformat(tab, keys):
    
    ids = [y for x in keys for y in x.split(';')]
    
    tmp2 = [len(x.split(';')) for x in keys]
    r = [i for i in range(len(tmp2)) for j in range(tmp2[i])]
        
    tab2 = tab.reset_index(drop = True)
    return(tab2.iloc[r].assign(site = ids))


def fast_MaxLFQ(norm_data):
    p_names = norm_data['protein_list'].unique()
    s_names = norm_data['sample_list'].unique()

    protein_to_index = {protein: index+1 for index, protein in enumerate(p_names)}
    sample_to_index = {protein: index+1 for index, protein in enumerate(s_names)}
    id_to_index = {protein: index+1 for index, protein in enumerate(norm_data['id'].unique())}

    import _msproteomics

    res = _msproteomics.iq_maxLFQ(np.asarray(norm_data['protein_list'].map(protein_to_index).values.tolist()), 
                                  np.asarray(norm_data['id'].map(id_to_index).values.tolist()), 
                                  np.asarray(norm_data['sample_list'].map(sample_to_index).values.tolist()), 
                                  np.asarray(norm_data['quant'].values.tolist()))

    return {'estimate' : pd.DataFrame(res['estimate'].reshape(len(res['_r_names']), len(res['_c_names'])), 
                                      columns = s_names[res['_c_names'] - 1],
                                      index = p_names[res['_r_names'] - 1]),
            'annotation' : pd.DataFrame(res['annotation'], columns = ['MaxLFQ_annotation'], index = p_names[res['_r_names'] - 1])}


def create_report_wideformat(report_lf,
                             sample_id_col = "R.FileName",
                             intensity_col = "log2_intensity",
                             primary_id = 'site',
                             secondary_id_cols = ["EG.PrecursorId", "EG.Library", "FG.Charge", "F.FrgIon", "F.Charge", "F.FrgLossType"],
                             annotation_cols = ["PG.Organisms"],
                             method = "maxlfq",
                             check = False):
    
    print("Concatenating secondary ids...")
    second_id = report_lf[secondary_id_cols[0]].astype(str)
    for col in range(1, len(secondary_id_cols)):
        second_id += ':' + report_lf[secondary_id_cols[col]].astype(str)

    aa = pd.DataFrame({'protein_list': report_lf[primary_id],
                       'sample_list': report_lf[sample_id_col],
                       'quant': report_lf[intensity_col],
                       'id': second_id})
    aa.dropna(axis=0)    
    
    if method == 'maxlfq':        

        result = fast_MaxLFQ(aa)    

        if check:
            
            print("checking pure Python implementation...")

            p_list = create_protein_list(aa)        
            res = create_protein_table(p_list, method = 'maxlfq')        
            
            # check difference here
            diff = np.double(result['estimate']) - np.double(res['estimate'])
            print('Max difference  =', np.max(diff, initial = 0, where = ~np.isnan(diff)))
            print('Min difference  =', np.min(diff, initial = 0, where = ~np.isnan(diff)))
            print('Same NAs        =', np.all(np.isnan(np.double(result['estimate'])) == np.isnan(np.double(res['estimate']))))
            print('Same row names  =', np.all(result['estimate'].index == res['estimate'].index))            
            print('Same col names  =', np.all(result['estimate'].columns == res['estimate'].columns))
            print('Same annotation =', np.all(result['annotation'] == res['annotation']))

    elif method == 'sum':        
        p_list = create_protein_list(aa)
        result = {'estimate': pd.DataFrame(np.nan, columns=list(p_list.values())[0].columns, index = list(p_list))}
        for i in p_list:
            tmp = np.exp2(p_list[i].to_numpy().astype(float))
            tmp2 = np.nansum(tmp, axis = 0)
            result['estimate'].loc[i] = np.log2(tmp2, out = np.full_like(tmp2, np.nan), where = tmp2 > 0)

    elif method == 'max':
        p_list = create_protein_list(aa)
        result = {'estimate': pd.DataFrame(np.nan, columns=list(p_list.values())[0].columns, index = list(p_list))}

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            for i in p_list:
                tmp = p_list[i].to_numpy().astype(float)            
                result['estimate'].loc[i] = np.nanmax(tmp, axis = 0)

    else:
        raise Exception('Unknown method.')
         
    a = list(report_lf[primary_id])
    ind = [a.index(x) for x in list(result['estimate'].index)]

    tmp = report_lf[[primary_id] + annotation_cols].iloc[ind]        
    tmp = tmp.set_axis(result['estimate'].index, axis='index')
    if method == 'maxlfq':                
        ret = tmp.join(result['annotation']).join(result['estimate'])
    else:
        ret = tmp.join(result['estimate'])

    return ret


def normalize(report_lf, sample_id_col, intensity_col):
    
    s = set(report_lf[sample_id_col])
    print(len(s) , " samples:")
    for si in s:
        print(si)

    #int_log2 = np.log2(report_lf[intensity_col],  out = np.zeros_like(report_lf[intensity_col]), where = (report_lf[intensity_col] != 0))
    int_log2 = np.log2(report_lf[intensity_col].values, out = np.full_like(report_lf[intensity_col], np.nan), where = report_lf[intensity_col] > 0)

    print(sum(np.isnan(int_log2)), " NA(s).")
    
    m = {si: np.median(int_log2[report_lf[sample_id_col] == si]) for si in s}
        
    sm = 0
    for i in m:
        sm += m[i]
    sm /= len(m)  

    for i in m:
        m[i] = sm - m[i]
        
    for i in m:
        idx = report_lf[sample_id_col] == i
        int_log2[idx] = int_log2[idx] + m[i]

    return int_log2


def create_peptide_key(modified_sequence, 
                       regex_str = r'\[[^\\[]+\]',
                       target_modification = '[Phospho (STY)]',
                       modifications = None):
    import re

    if modifications is None:
        m = {i for s in modified_sequence for i in re.findall(regex_str, s)}

        if target_modification in m:
            m.remove(target_modification)
            mods = [target_modification] + sorted(list(m))
        else :
            raise Exception('target modification not found.')
    else:
        mods = modifications

    ptm_key = [''] * len(modified_sequence)
    STY_count = [0] * len(modified_sequence)

    from collections import Counter

    for i, s in enumerate(modified_sequence):
        a = Counter(re.findall(regex_str, s))
        STY_count[i] = a[mods[0]]
        ptm_key[i] = re.sub(regex_str, '', s) + '_' + str(STY_count[i])
        for j in range(1, len(mods)):
            ptm_key[i] += '.' + str(a[mods[j]])

    return(ptm_key, STY_count)
