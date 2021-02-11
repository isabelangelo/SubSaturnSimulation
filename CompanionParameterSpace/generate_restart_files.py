from generate_inicon_files import *
import numpy as np
import glob

output_dir = '../SimulationResults/Run5/'
inicon_dir = '../SimulationResults/Run5/IniCon'

AU = 1.496E13       # astronomical unit, cm 
Rsun = 6.96e10    # radius of sun, cm 
pi = np.pi          # pi, --

# need unfinished numbers and unfinished filenames
loop_number = int(input('previous loop number: '))
loop_str = '.loop'+str(loop_number)
unfinished_filenames = sorted(glob.glob(output_dir+'Output/output_*.txt'+loop_str))
unfinished_numbers = [fname[40:-10] for fname in unfinished_filenames]


print(unfinished_numbers)
print(len(unfinished_numbers))

## temporary
#unfinished_filenames = ['./output_15.txt']
#unfinished_numbers = [15]


for i in range(len(unfinished_filenames)):
    
    filename = unfinished_filenames[i]
    filenumber = unfinished_numbers[i]

    # read in parameters from output file
    last_line = open(filename).readlines()[-1].split('\t')[:-1]
    fp = [float(i) for i in last_line]
    age =  6.31*1.2*1000 # time for simulation to run to
    end_time = fp[2] # start time for simulation
    e1,e2,a1,i,beta = fp[3], fp[4], fp[7], fp[10],fp[13]
    g1, g2 = np.rad2deg(fp[5]), np.rad2deg(fp[6])
    m1,m2,a2,m3 = fp[21],fp[23],fp[25],fp[26]
    R1, R2 = fp[22]*AU/Rsun, fp[24]*AU/Rsun # convert from AU to Rsun
    beta2,gamma,gamma2 = fp[32], fp[33],fp[34]
    
    # compute spins from output
    spin1h, spin1e, spin1q = fp[11], fp[15], fp[16]
    spin2h, spin2e, spin2q = fp[19], fp[17], fp[18]
    spin1P_radyr = np.sqrt(spin1h**2.+spin1e**2.+spin1q**2.)
    spin2P_radyr = np.sqrt(spin2h**2.+spin2e**2.+spin2q**2.)
    spin1P = 2.*pi*365.25/spin1P_radyr
    spin2P = 2.*pi*365.25/spin2P_radyr
    
    final_conditions = np.array([m1,m2,m3,R1,R2,spin1P,spin2P,beta,beta2,gamma,gamma2,a1,a2,e1,e2,g1,g2,i,end_time,age])
    
    # write new triple.in file
    file_end = str(filenumber)
#    print(file_end)
    create_file(final_conditions, inicon_dir, file_end, start_time = end_time, restart=True)
    print(file_end)
