import pandas
import numpy as np
from datetime import datetime
import re
import os

'''
This script performs the main quantile filtering function of DIA-SIFT. 

INPUT:
- [cwd] folder containing processed .csv files (from merge_and_preprocess_PLGS_outputs.py)
    "[SampleName]_peptideMergedAndPreprocessed.csv"
    "[SampleName]_fragmentMergedAndPreprocessed.csv"
    
OUTPUT:
- filtered & unfiltered .csv files in the same directory
    "[SampleName]_[iqrfactor]_[minnumobs]_unfiltered.csv"
    "[SampleName]_[iqrfactor]_[minnumobs]_DIA-SIFTed.csv"
    "[SampleName]_[iqrfactor]_[minnumobs]_unfiltered_protein.csv"
    "[SampleName]_[iqrfactor]_[minnumobs]_DIA-SIFTed_protein.csv"
'''

# SET INPUTS HERE:
iqrfactor = 0.2  # 0.2 recommended
minnumobs = 3  # 3 recommended
cwd = r'C:\Admin\Desktop\Sarah\DIA-SIFT_revisions\PLGS_outputs'

processed_suffix = 'MergedAndPreprocessed.csv'
f_suffix = '_fragmentMergedAndPreprocessed.csv'
p_suffix = '_peptideMergedAndPreprocessed.csv'

# regex for sample names, don't modify these
match = '(.+)_.+' + processed_suffix
filematch = '(.+)' + processed_suffix
f_filematch = '(.+)' + f_suffix
p_filematch = '(.+)' + p_suffix

start1 = datetime.now()

def overall(finalPeptide, finalFragment, iqrfactor, minnumobs):
    print('    fragment file: ', finalFragment)
    print('    peptide file: ', finalPeptide)
    samplename = re.sub('_fragmentMergedAndPreprocessed.csv', '', finalFragment)
    outputName = samplename

    # add 'ID #' column to merged peptide & fragment files
    fullpeptdf = pandas.read_csv(finalPeptide)
    fullfragdf = pandas.read_csv(finalFragment)
    fullpeptdf['ID #'] = fullpeptdf.index
    fullpeptdf['ID #'] = 'P_' + fullpeptdf['ID #'].astype(str)
    fullfragdf['ID #'] = fullfragdf.index
    fullfragdf['ID #'] = 'F_' + fullfragdf['ID #'].astype(str)

    ms1df = fullpeptdf
    ms2df = fullfragdf

    ms2df = ms2df.rename(columns={'fragment.str': 'peptide.mhp', 'L/H Ratios': 'precursor.pairLtoHRatio'})
    ms1df['Type'] = 'P'
    ms2df['Type'] = 'F'

    ms1ms2_combined = pandas.concat([ms1df, ms2df], axis=0)
    ms1ms2_combined.index.name = 'protein.Entry'

    def q1(pivot):
        Q1 = np.percentile(pivot, 25)
        return Q1

    def q3(pivot):
        Q3 = np.percentile(pivot, 75)
        return Q3

    def iqr(pivot):
        get_IQR = np.percentile(pivot, 75) - np.percentile(pivot, 25)
        return get_IQR

    def lower(pivot):
        get_lower = np.percentile(pivot, 25) - (iqrfactor*(np.percentile(pivot, 75) - np.percentile(pivot, 25)))
        return get_lower

    def upper(pivot):
        get_upper = np.percentile(pivot, 75) + (iqrfactor*(np.percentile(pivot, 75) - np.percentile(pivot, 25)))
        return get_upper

    def cv(pivot):
        get_cv = np.std(pivot) / np.mean(pivot)
        return get_cv

    def num_obs(pivot):
        get_num_obs = len(pivot)
        return get_num_obs

    # use pivot tables to get [mean, median, Q1, Q3] for all observed ratios of each protein
    ms1ms2_combined = ms1ms2_combined.reset_index(drop=True)
    ms1ms2_pivot = pandas.pivot_table(ms1ms2_combined, values='precursor.pairLtoHRatio', index=['protein.Entry'], aggfunc=[np.median, np.mean, np.std, cv, num_obs, q1, q3, iqr, lower, upper])
    ms1ms2_pivot.columns = ms1ms2_pivot.columns.droplevel(1)  # remove multiindex
    ms1ms2_pivot = ms1ms2_pivot.reset_index()  # remove multiindex (again)
    ms1ms2_merge = pandas.merge(ms1ms2_combined, ms1ms2_pivot, how='left', on='protein.Entry')

    # perform quartile tests within each row, ** maintain proteins w/ two observations & CV <= 20% **
    def oqr_test(ms12_value, upper_value, lower_value, num_obs_value, cv_value):
        if 0 <= ms12_value <= upper_value and ms12_value >= lower_value and num_obs_value >= minnumobs:
            ms12_result = 1
        elif num_obs_value == 2 and cv_value <= 0.2:  # keep proteins with only two measurements if they seem to agree
            ms12_result = 1
        else:
            ms12_result = 0
        return ms12_result

    ms1ms2_merge['overall quartile test?'] = ""
    column = ms1ms2_merge.columns.get_loc('overall quartile test?')
    row = 0
    for value in ms1ms2_merge['precursor.pairLtoHRatio']:
        ms12_upper = ms1ms2_merge.iloc[row]['upper']
        ms12_lower = ms1ms2_merge.iloc[row]['lower']
        ms12_num_obs = ms1ms2_merge.iloc[row]['num_obs']
        ms12_cv = ms1ms2_merge.iloc[row]['cv']
        result = oqr_test(value, ms12_upper, ms12_lower, ms12_num_obs, ms12_cv)
        ms1ms2_merge.iloc[row, column] = result
        row += 1

    # make new dataframe with events that pass the overall quartile filter
    filteredDF = ms1ms2_merge.loc[ms1ms2_merge['overall quartile test?'] == 1]

    unfiltered_outputName = samplename + "_unfiltered.csv"
    filtered_outputName = samplename + "_DIA-SIFTed.csv"
    unfiltered_protein_outputName = samplename + "_unfiltered_protein.csv"
    filtered_protein_outputName = samplename + "_DIA-SIFTed_protein.csv"

    unfiltered_protein = ms1ms2_merge[['protein.Accession', 'protein.Description', 'protein.Entry',
                                      'median', 'mean', 'std', 'cv', 'num_obs']]
    unfiltered_protein = unfiltered_protein.drop_duplicates(['protein.Entry'])
    unfiltered_protein = unfiltered_protein[unfiltered_protein.num_obs >= 2]
    protein_info = unfiltered_protein[["protein.Accession", "protein.Description", "protein.Entry"]]

    clean_filteredDF = filteredDF.drop(columns=['median', 'mean', 'std', 'cv', 'num_obs', 'q1', 'q3', 'iqr', 'lower', 'upper'])
    filtered_pivot = pandas.pivot_table(clean_filteredDF, values='precursor.pairLtoHRatio', index=['protein.Entry'], aggfunc=[np.median, np.mean, np.std, cv, num_obs])
    filtered_pivot.columns = filtered_pivot.columns.droplevel(1)  # remove multiindex
    filtered_pivot = filtered_pivot.reset_index()  # remove multiindex (again)
    filtered_pivot = pandas.merge(protein_info, filtered_pivot, on="protein.Entry", how="inner")

    ms1ms2_merge.to_csv(unfiltered_outputName, index=False)
    filteredDF.to_csv(filtered_outputName, index=False)
    unfiltered_protein.to_csv(unfiltered_protein_outputName, index=False)
    filtered_pivot.to_csv(filtered_protein_outputName, index=False)

    print('Done with ' + outputName + ' !  ' + str(datetime.now() - start1) + "\n")


def perform_filtering(workingdir):
    print(workingdir)
    preprocessedfiles = [f for f in os.listdir(workingdir) if re.search(filematch, f)]
    n = 1
    samplelist = []
    for f in preprocessedfiles:
        samplename = re.search(match, f).group(1)
        if samplename:
            samplelist.append(samplename)
    samplelist = set(samplelist)
    for sample in samplelist:
        print('Filtering sample {} of {} ... '.format(n, len(samplelist)))
        files = [f for f in preprocessedfiles if re.search(sample, f)]
        fragfile = os.path.join(workingdir, files[0])
        pepfile = os.path.join(workingdir, files[1])
        overall(pepfile, fragfile, iqrfactor, minnumobs)
        n += 1

perform_filtering(cwd)

