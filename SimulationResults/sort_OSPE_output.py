import glob
import numpy as np

"""
the purpose of this code is to sort all of the OSPE output files based on their end state.
The OSPE output filenames are sorted into 3 lists:
    zeros : systems that didn't survive (planet spiraled into star)
    ones_and_twos : systems that survived (although some may not be finished running)
    negatives : negatives : files that produce an error
    
"""

# ===================================================

def get_endstate(fname):
    """
    gets the final state of an OSPE run
    :param fname : filename of output file, e.g. 'output/files/output_1.txt' , str
    :return: endstate (-1,0,1,2) , int
    """
    lines = open(fname).readlines() 
    endstate = int(lines[-1].split('\t')[0])
    return endstate

# ===================================================

def get_values(fname):
    """
    extracts system parameters for initial and final conditions
    :param fname : filename of output file, e.g. 'output/files/output_1.txt' , str
    :return: initial and final parameters (a_i, a_f, e_i, e_f, i_i, i_f, m3)  , tuple
    """
    # open file
    lines = open(fname).readlines()
    # get initial and final conditions
    initial_values = np.array([float(i) for i in lines[0].split('\t')[:-1]])
    final_values = np.array([float(i) for i in lines[-1].split('\t')[:-1]])
    # store values 
    a_i = initial_values[7]; a_f = final_values[7]
    e_i = initial_values[3]; e_f = final_values[3]
    i_i = initial_values[10]; i_f = final_values[10]
    m3 = initial_values[26]
    return (a_i, a_f, e_i, e_f, i_i, i_f, m3)  
    
# ===================================================

# sort all files by end state and print results
zeros = [] # didn't survive (spiraled into star)
ones_and_twos = [] # survived
negatives = [] # error

output_file_dir = input('Output file directory: ')
fnames = glob.glob(output_file_dir)
for name in fnames:
    s = get_endstate(name)
    if s == 0:
        zeros.append(name)
    if s == 1 or s == 2:
        ones_and_twos.append(name)
    if s == -1:
        negatives.append(name)

print(len(zeros), '/', len(fnames), 'didn\'t survive')
print(len(ones_and_twos), '/', len(fnames), 'survived')
print(len(negatives), '/', len(fnames), 'produced error')

# next get a list of all the file names that need to be re-run
for fname in ones_and_twos:
    lines = open(fname).readlines()
    # get final time
    final_values = np.array([float(i) for i in lines[-1].split('\t')[:-1]])
    final_time = final_values[2]
    if final_time < 7.5e9:
        print(fname, final_time)

# ===================================================

# store all files that satisfy end conditions of KOI367
def get_relevant_systems():
    """
    compiles list of all files that have evolved to a similar orbital configuration as KOI-367
     :return: list of filenames for output file name strings (list)
    """
    a_subsat = 0.197; a_subsat_err = 0.021
    e_subsat = 0.836; e_subsat_err = 0.013
    relevant_systems = []

    for name in ones_and_twos:
        values = get_values(name)
        af = values[1]; ef = values[4]
        if (a_subsat/3. < af) & (af < a_subsat*3):
            if (e_subsat-e_subsat_err < ef) & (ef < e_subsat+e_subsat_err):
                relevant_systems.append(name)

    return relevant_systems
