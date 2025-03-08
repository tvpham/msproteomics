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

import sys
import argparse
import numpy as np
import pandas as pd

import msproteomics

def main():
    
    parser = argparse.ArgumentParser(description='Create phosphosite and phosphopeptide report.', usage = 'sitereport inputfile [...]')
    
    parser.add_argument('file')
    parser.add_argument('-tool', '--tool', default = 'generic', choices=['generic', 'diann', 'sn'], type = str, help = 'Processing tool, diann for DIA-NN and sn for Spectronaut')
    parser.add_argument('-output_site', '--output_site', default = 'phosphosite-report.txt', type = str, help = 'Filename for phosphosite report')
    parser.add_argument('-output_peptide', '--output_peptide', default = 'phosphopeptide-report.txt', type = str, help = 'Filename for phosphopeptide report')
    parser.add_argument('-sample_id_col', '--sample_id_col', default = None, type = str, help = 'The column containing sample id')
    parser.add_argument('-intensity_col', '--intensity_col', default = None, type = str, help = 'The column containing intensity')
    parser.add_argument('-secondary_id_cols', '--secondary_id_cols', default = None, nargs='*', type = str, help = 'Columns forming secondary ids')
    parser.add_argument('-annotation_cols', '--annotation_cols', default = None, nargs='*', type = str, help = 'Annotation columns')

    parser.add_argument('-protein_id_col', '--protein_id_col', default = None, type = str, help = 'The column containing protein id')
    parser.add_argument('-site_id_col', '--site_id_col', default = None, type = str, help = 'The column containing site id')

    parser.add_argument('-site_filter_double_less', '--site_filter_double_less', default = None, nargs=2, action='append', type = str, help = 'Site filtering double less than or equal to')
    parser.add_argument('-site_filter_double_greater', '--site_filter_double_greater', default = None, nargs=2, action='append', type = str, help = 'Site filtering double greater than or equal to')
    parser.add_argument('-site_filter_string_equal', '--site_filter_string_equal', default = None, nargs=2, action='append', type = str, help = 'Site filtering string equal')
    parser.add_argument('-site_filter_string_not_equal', '--site_filter_string_not_equal', default = None, nargs=2, action='append', type = str, help = 'Site filtering string not equal')

    parser.add_argument('-modified_sequence_col', '--modified_sequence_col', default = None, type = str, help = 'The column containing modified sequences')
    parser.add_argument('-regex_str', '--regex_str', default = None, type = str, help = 'Regular expression to extract sites')
    parser.add_argument('-target_modification', '--target_modification', default = None, type = str, help = 'Target modification')
    parser.add_argument('-modifications', '--modifications', default = None, nargs='*', type = str, help = 'List of modifications')

    parser.add_argument('-peptide_filter_double_less', '--peptide_filter_double_less', default = None, nargs=2, action='append', type = str, help = 'Peptide filtering double less than or equal to')
    parser.add_argument('-peptide_filter_double_greater', '--peptide_filter_double_greater', default = None, nargs=2, action='append', type = str, help = 'Peptide filtering double greater than or equal to')
    parser.add_argument('-peptide_filter_string_equal', '--peptide_filter_string_equal', default = None, nargs=2, action='append', type = str, help = 'Peptide filtering string equal')
    parser.add_argument('-peptide_filter_string_not_equal', '--peptide_filter_string_not_equal', default = None, nargs=2, action='append', type = str, help = 'Peptide filtering string not equal')

    parser.add_argument('-normalize', '--normalize', default = 'none', type = str, choices=['none', 'median'], help = 'Normalization method')
    parser.add_argument('-quant_method', '--quant_method', default = 'maxlfq', type = str, choices=['sum', 'max', 'maxlfq'], help = 'Quantification method')

    parser.add_argument('-check', '--check', default = False, action = 'store_true')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        print('Missing input filename. Try phosphoreport -h.')
        sys.exit(0)
    
    inputfile = sys.argv[1]

    if args.tool == "generic" or args.tool == "diann":
        if args.sample_id_col is None:
            args.sample_id_col = 'Run'
        if args.intensity_col is None:
            args.intensity_col = 'Fragment.Intensity'
        if args.secondary_id_cols is None:
            args.secondary_id_cols = ['Precursor.Id', 'Fragment.Rel.Id']
        if args.annotation_cols is None:
            args.annotation_cols = ['Fasta.Files']

        # site paramters
        if args.protein_id_col is None:
            args.protein_id_col = 'Protein.Ids'
        if args.site_id_col is None:
            args.site_id_col = 'Phospho.Site.Specs'

        if args.site_filter_double_less is None:
            args.site_filter_double_less = [['Global.Q.Value', '0.01']]
        elif args.site_filter_double_less == [['none', 'none']]:
            args.site_filter_double_less = None

        if args.site_filter_double_greater is None:
            args.site_filter_double_greater = [['PTM.Site.Confidence', '0.01']]
        elif args.site_filter_double_greater == [['none', 'none']]:
            args.site_filter_double_greater = None

        # peptide paramters
        if args.modified_sequence_col is None:
            args.modified_sequence_col = 'Modified.Sequence'
        if args.regex_str is None:
            args.regex_str = '\\(UniMod:[0-9]*\\)'
        if args.target_modification is None:
            args.target_modification = '(UniMod:21)'                

        if args.peptide_filter_double_less is None:
            args.peptide_filter_double_less = [['Global.Q.Value', '0.01']]
        elif args.peptide_filter_double_less == [['none', 'none']]:
            args.peptide_filter_double_less = None
    
    elif args.tool == 'sn':

        if args.sample_id_col is None:
            args.sample_id_col = 'R.FileName'
        if args.intensity_col is None:
            args.intensity_col = 'F.PeakArea'
        if args.secondary_id_cols is None:
            args.secondary_id_cols = ['EG.PrecursorId', 'EG.Library', 'FG.Charge', 'F.FrgIon', 'F.Charge', 'F.FrgLossType']
        if args.annotation_cols is None:
            args.annotation_cols = ['PG.Organisms']

        # site paramters
        if args.protein_id_col is None:
            args.protein_id_col = 'PG.ProteinAccessions'
        if args.site_id_col is None:
            args.site_id_col = 'EG.ProteinPTMLocations'
        
        if args.site_filter_double_greater is None:
            args.site_filter_double_greater = [['EG.PTMAssayProbability', '0.75']]
        elif args.site_filter_double_greater == [['none', 'none']]:
            args.site_filter_double_greater = None

        # peptide paramters
        if args.modified_sequence_col is None:
            args.modified_sequence_col = 'EG.ModifiedSequence'
            
        if args.regex_str is None:
            args.regex_str = '\[[^\\[]+\]'
        if args.target_modification is None:
            args.target_modification = '[Phospho (STY)]'
    
    else:
        raise Exception('Unsupported tool.')

    print('\nGeneral setting:')
    print('  input data        = ', inputfile)
    print('  processing tool   = ', args.tool)
    print('  sample_id_col     = ', args.sample_id_col)
    print('  intensity_col     = ', args.intensity_col)
    print('  secondary_id_cols = ', args.secondary_id_cols)        
    print('  annotation_cols   = ', args.annotation_cols)
    print('  normalize         = ', args.normalize)    
    print('  quant_method      = ', args.quant_method)                
    
    print('\nPhosphosite setting:')
    print('  protein_id_col               = ', args.protein_id_col) 
    print('  site_id_col                  = ', args.site_id_col) 
    print('  site_filter_double_less      = ', args.site_filter_double_less) 
    print('  site_filter_double_greater   = ', args.site_filter_double_greater) 
    print('  site_filter_string_equal     = ', args.site_filter_string_equal) 
    print('  site_filter_string_not_equal = ', args.site_filter_string_not_equal) 
    print('  output_site                  = ', args.output_site) 

    print('\nPhosphopeptide setting:')
    print('  modified_sequence_col           = ', args.modified_sequence_col) 
    print('  regex_str                       = ', args.regex_str) 
    if args.modifications is None:
        print('  target_modification             = ', args.target_modification) 
    else:
        print('  modifications                   = ', args.modifications) 

    print('  peptide_filter_double_less      = ', args.peptide_filter_double_less) 
    print('  peptide_filter_string_not_equal = ', args.peptide_filter_string_not_equal) 
    print('  output_peptide                  = ', args.output_peptide) 
        
    print('\nLoading data file:', inputfile)

    required_cols = args.secondary_id_cols + [args.sample_id_col] + [args.intensity_col] + args.annotation_cols + [args.protein_id_col] + [args.site_id_col] + [args.modified_sequence_col]
    
    if args.site_filter_double_less is not None:
        for c in args.site_filter_double_less:
            if c[0] not in required_cols:
                required_cols.append(c[0])
    if args.site_filter_double_greater is not None:                
        for c in args.site_filter_double_greater:
            if c[0] not in required_cols:
                required_cols.append(c[0])
    if args.site_filter_string_equal is not None:
        for c in args.site_filter_string_equal:
            if c[0] not in required_cols:
                required_cols.append(c[0])
    if args.site_filter_string_not_equal is not None:
        for c in args.site_filter_string_not_equal:
            if c[0] not in required_cols:
                required_cols.append(c[0])

    if args.peptide_filter_double_less is not None:
        for c in args.peptide_filter_double_less:
            if c[0] not in required_cols:
                required_cols.append(c[0])
    if args.peptide_filter_double_greater is not None:                
        for c in args.peptide_filter_double_greater:
            if c[0] not in required_cols:
                required_cols.append(c[0])
    if args.peptide_filter_string_equal is not None:
        for c in args.peptide_filter_string_equal:
            if c[0] not in required_cols:
                required_cols.append(c[0])
    if args.peptide_filter_string_not_equal is not None:
        for c in args.peptide_filter_string_not_equal:
            if c[0] not in required_cols:
                required_cols.append(c[0])


    d = pd.read_csv(inputfile, sep = '\t', usecols = required_cols, dtype = {args.intensity_col:'float64'})

    print(d.shape[0], 'rows x', d.shape[1], 'columns read')        

    d = d.loc[d[args.intensity_col] > 0]
    d.reset_index()

    print(d.shape[0], 'rows after filtering out <= 0 intensities.')    

    if args.normalize == 'median':        
        log2_intensity = msproteomics.normalize(d, args.sample_id_col, args.intensity_col)        
    else:            
        log2_intensity = np.log2(d[args.intensity_col].values, out = np.full_like(d[args.intensity_col], np.nan), where = d[args.intensity_col] > 0)        
    
    if sum(np.isnan(log2_intensity)) > 0:        
        print("Warning: missing values in the intensity column.")

    d = d.assign(log2_intensity = log2_intensity)    

    #-- site filtering
    index = np.full(d.shape[0], True)

    if args.site_filter_double_less is not None:
        for f in args.site_filter_double_less:
            col = f[0]            
            val = float(f[1])            
            index = index & (~d[col].isna()) & (d[col].astype(np.float64) <= val)

    if args.site_filter_double_greater is not None:
        for f in args.site_filter_double_greater:
            col = f[0]            
            val = float(f[1])
            index = index & (~d[col].isna()) & (d[col].astype(np.float64) >= val)

    if args.site_filter_string_not_equal is not None:
        for f in args.site_filter_string_not_equal:
            col = f[0]
            val = f[1]
            index = index & (~d[col].isna()) & (d[col].astype(str) != val)

    if args.site_filter_string_equal is not None:
        for f in args.site_filter_string_equal:
            col = f[0]
            val = f[1]
            index = index & (~d[col].isna()) & (d[col].astype(str) == val)

    d_filtered = d.loc[index]
    
    print('{0} rows after other filtering(s)\n'.format(d_filtered.shape[0]))

    protein_ids = d_filtered[args.protein_id_col].fillna('').tolist()
    
    ptm_locations = d_filtered[args.site_id_col].fillna('').tolist()

    print('Creating site identifiers')
    (all_sites, unique_sites, first_sites, unique_first_sites) = msproteomics.create_site_key(protein_ids, ptm_locations)
    d_filtered = msproteomics.create_site_report_longformat(d_filtered.assign(all_sites = all_sites), first_sites)
    d_filtered = d_filtered.loc[d_filtered['site'] != '']

    print('Creating site report')
    result = msproteomics.create_report_wideformat(d_filtered, method = args.quant_method,
                                            sample_id_col = args.sample_id_col,
                                            intensity_col = "log2_intensity", 
                                            secondary_id_cols = args.secondary_id_cols, annotation_cols = args.annotation_cols + ['all_sites'],
                                            check = args.check)

    result.to_csv(args.output_site, sep = '\t', index = False)


    #-- peptide filtering
    index = np.full(d.shape[0], True)

    if args.peptide_filter_double_less is not None:
        for f in args.peptide_filter_double_less:
            col = f[0]
            val = float(f[1])
            index = index & (~d[col].isna()) & (d[col].astype(np.float64) <= val)

    if args.peptide_filter_double_greater is not None:
        for f in args.peptide_filter_double_greater:
            col = f[0]
            val = float(f[1])
            index = index & (~d[col].isna()) & (d[col].astype(np.float64) >= val)

    if args.peptide_filter_string_not_equal is not None:
        for f in args.peptide_filter_string_not_equal:
            col = f[0]
            val = f[1]
            index = index & (~d[col].isna()) & (d[col].astype(str) != val)

    if args.peptide_filter_string_equal is not None:
        for f in args.peptide_filter_string_equal:
            col = f[0]
            val = f[1]
            index = index & (~d[col].isna()) & (d[col].astype(str) == val)

    d_filtered = d.loc[index]        
    print('\n{0} rows after peptide filtering'.format(d_filtered.shape[0]))

    print('Creating peptide identifiers')
    if args.modifications is None:
        (ptm_key, STY_count) = msproteomics.create_peptide_key(d_filtered[args.modified_sequence_col].tolist(),
                                                               regex_str = args.regex_str,
                                                               target_modification = args.target_modification)
    else:                                                               
        (ptm_key, STY_count) = msproteomics.create_peptide_key(d_filtered[args.modified_sequence_col].tolist(),
                                                               regex_str = args.regex_str,
                                                               modifications = args.modifications)
    print('Creating peptide report')
    result = msproteomics.create_report_wideformat(d_filtered.assign(ptm_key = ptm_key, STY_count = STY_count),
                                                    sample_id_col = args.sample_id_col,
                                                    intensity_col = "log2_intensity",
                                                    primary_id = 'ptm_key', 
                                                    secondary_id_cols = args.secondary_id_cols,
                                                    annotation_cols = ['STY_count'] + args.annotation_cols,
                                                    method = args.quant_method,
                                                    check = args.check)
    result.to_csv(args.output_peptide, sep = '\t', index = False)
