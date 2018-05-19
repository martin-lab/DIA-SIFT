import pandas
import os
import re

'''
This script pre-processes ProteinLynx Global Server (PLGS) output .csv files for SILAC quantile filtering (DIA-SIFT).

INPUTS (each of these must be defined below):
- [cwd] folder containing peptide and fragment files from 3 replicate injections of each sample
- [prefix] include only the first part of the filenames right up until a unique sample identifier
- [p_suffix] include only the part of the peptide filenames immediately following the unique sample identifier
- [f_suffix] include only the part of the fragment filenames immediately following the unique sample identifier

EXAMPLE INPUTS:
SILAC_HeLa_SampleA_replicate001_IA_final_peptide.csv
SILAC_HeLa_SampleA_replicate001_IA_final_fragment.csv

[SILAC_HeLa_]SampleA_replicate001[_IA_final_peptide.csv]
      ^                                     ^
    prefix                               p_suffix
    
[SILAC_HeLa_]SampleA_replicate001[_IA_final_fragment.csv]
      ^                                     ^
    prefix                               f_suffix
    
prefix = 'SILAC_HeLa_'
p_suffix = '_IA_final_peptide.csv'
f_suffix = '_IA_final_fragment.csv'
(the unique sample identifier would be 'SampleA_replicate001')


OUTPUTS:
- merged and pre-processed .csv files for peptides & fragments:
    "[SampleName]_peptideMergedAndPreprocessed.csv"
    "[SampleName]_fragmentMergedAndPreprocessed.csv"
'''

# SET INPUTS HERE (see above for details):
cwd = r'C:\Admin\Desktop\Sarah\DIA-SIFT_revisions\PLGS_outputs'  # leave a lowercase 'r' before the working directory string
prefix = '293T_SILAC_' # this should be common to all files in the project
p_suffix = '.+_IA_final_peptide.csv'  # make sure to use re to match different injection numbers!
f_suffix = '.+_IA_final_fragment.csv'  # make sure to use re to match different injection numbers!

f_filematch = prefix + '(.+)' + f_suffix
p_filematch = prefix + '(.+)' + p_suffix

def computeHeavyorLight(string):
    if "SILAC" in string:
        return "Heavy"
    else:
        return "Light"


def checkMatchType(string1, string2):
    if string1 == string2:
        return "Yes"
    elif string1 == "PepFrag1" and string2 == "PepFrag2":
        return "Yes"
    elif string1 == "PepFrag2" and string2 == "PepFrag1":
        return "Yes"
    else:
        return "No"


def computeSILACRatios(light, heavy):
    return light / heavy


def computeSILAC(inputDF):
    initialDF = inputDF[['protein.Entry', 'protein.dataBaseType', 'peptide.matchType',
                         'peptide.modification', 'peptide.seq', 'fragment.fragmentType', 'fragment.fragInd',
                         'Neutral.LossType', 'fragment.str', 'fragment.seq', 'fragment.fragSite', 'product.inten',
                         'fragmentProduct.deltaMhpPPM', 'peptidePrecursor.deltaMhpPPM']]  # deletes unnecessary columns

    initialDF = initialDF[initialDF['fragment.fragmentType'] != "b"]  # removes b ions
    initialDF = initialDF[initialDF['protein.dataBaseType'] != "Reverse"]  # removes Reverse fragments
    initialDF = initialDF[initialDF['Neutral.LossType'] == "None"]  # keep fragments with no loss type
    initialDF = initialDF[initialDF['fragmentProduct.deltaMhpPPM'] < 20]
    initialDF = initialDF[initialDF['fragmentProduct.deltaMhpPPM'] > -20]
    initialDF = initialDF[initialDF['peptidePrecursor.deltaMhpPPM'] < 10]  # remove fragments outside of ppm range
    initialDF = initialDF[initialDF['peptidePrecursor.deltaMhpPPM'] > -10]  # remove fragments outside of ppm range

    initialDF['SILAC'] = initialDF.apply(lambda row: computeHeavyorLight(row['peptide.modification']), axis=1)

    HeavyDF = initialDF[initialDF['SILAC'] == "Heavy"]
    LightDF = initialDF[initialDF['SILAC'] == "Light"]
    SILACDF = pandas.merge(HeavyDF, LightDF, on=['protein.Entry', 'peptide.seq', 'fragment.str', 'fragment.seq'],
                           how='outer')  # merge Light and Heavy sheets on 4 parameters

    SILACDF['L/H Ratios'] = SILACDF.apply(
        lambda row: computeSILACRatios(row['product.inten_y'], row['product.inten_x']), axis=1)
    return SILACDF  # calculate ratios based on intensities


def processFinalFragment(scout1, scout2, scout3):
    initialDF3 = pandas.read_csv(scout3, encoding="latin-1")
    initialDF1 = pandas.read_csv(scout1, encoding="latin-1")
    initialDF2 = pandas.read_csv(scout2, encoding="latin-1")

    SILACDF1 = computeSILAC(initialDF1)
    SILACDF2 = computeSILAC(initialDF2)
    SILACDF3 = computeSILAC(initialDF3)

    # combine 3 separate DF's
    SILACDF4 = SILACDF1.append(SILACDF2)
    SILACDF_merged = SILACDF4.append(SILACDF3)  # three Scouts
    SILACDF_merged['matchType Match?'] = SILACDF_merged.apply(
        lambda row: checkMatchType(row['peptide.matchType_x'], row['peptide.matchType_y']), axis=1)
    SILACDF_final = SILACDF_merged[SILACDF_merged['matchType Match?'] == "Yes"]

    filename = re.findall(f_filematch, scout1)[0]
    outputFileName = os.path.join(cwd, filename + "_fragmentMergedAndPreprocessed.csv")
    SILACDF_final.to_csv(outputFileName)


def processFinalPeptideDF(finalpeptideDF, finalpeptideFileName):  # outputs excel file with 8 sheets (ratios, potential uniques, pivot table)
    initialDF = finalpeptideDF
    initialDF = initialDF[initialDF['protein.dataBaseType'] != "Reverse"]

    # only use peptides with mass error < 10 ppm
    initialDF = initialDF[initialDF['peptidePrecursor.deltaMhpPPM'] < 10]
    initialDF = initialDF[initialDF['peptidePrecursor.deltaMhpPPM'] > -10]
    ratiosDF = initialDF[initialDF['precursor.pairLtoHRatio'].notnull()]  # ratios

    filename = os.path.join(cwd, finalpeptideFileName + "_peptideMergedAndPreprocessed.csv")
    ratiosDF.to_csv(filename)


def process_fragment_reps(workingdir):
    fragfiles = [f for f in os.listdir(workingdir) if re.search(f_filematch, f)]
    n = 1
    samplelist =[]
    for f in fragfiles:
        samplename = re.search(f_filematch, f).group(1)
        if samplename:
            samplelist.append(samplename)
    samplelist = set(samplelist)

    for sample in samplelist:
        print('Processing fragment files for sample {} of {} ... '.format(n, len(samplelist))),
        reps = [f for f in fragfiles if re.search(sample, f)]
        rep1 = os.path.join(workingdir, reps[0])
        rep2 = os.path.join(workingdir, reps[1])
        rep3 = os.path.join(workingdir, reps[2])
        processFinalFragment(rep1, rep2, rep3)
        print('Done processing fragment files for sample: ', sample)
        n += 1


def process_peptide_reps(workingdir):
    pepfiles = [f for f in os.listdir(workingdir) if re.search(p_filematch, f)]
    n = 1
    samplelist = []
    for f in pepfiles:
        samplename = re.search(p_filematch, f).group(1)
        if samplename:
            samplelist.append(samplename)
    samplelist = set(samplelist)

    for name in samplelist:
        print('Processing peptide files for sample {} of {} ... '.format(n, len(samplelist))),
        reps = [f for f in pepfiles if re.search(name, f)]
        rep1 = os.path.join(workingdir, reps[0])
        rep2 = os.path.join(workingdir, reps[1])
        rep3 = os.path.join(workingdir, reps[2])

        rep1df = pandas.read_csv(rep1, header=0, encoding="latin-1")
        rep2df = pandas.read_csv(rep2, header=0, encoding="latin-1")
        rep3df = pandas.read_csv(rep3, header=0, encoding="latin-1")

        mergedPepDF = rep1df.append(rep2df, ignore_index=True)
        mergedPepDF = mergedPepDF.append(rep3df, ignore_index=True)

        outfilename = os.path.join(workingdir, name)
        processFinalPeptideDF(mergedPepDF, outfilename)
        print('Done with sample: ', name),
        n += 1


process_peptide_reps(cwd)
process_fragment_reps(cwd)










