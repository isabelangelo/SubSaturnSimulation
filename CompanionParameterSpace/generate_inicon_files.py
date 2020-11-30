from generate_inicon import *

# ==========================================================================
"""
Functions for generating initial conditions files in OSPE

The create_file function reads in initial conditions generated in generate_inicon 
and writes to triple.in file
"""
# ==========================================================================

# flags for OSPE
quadrupole = 0
octupole = 0
GR = 0
TF = 0
ML1 = 0
MB1 = 0
ML2 = 1
MB2 = 1
SSE1 = 0
SSE2 = 1
SSE3 = 1
tMSMyr = 0

# ==========================================================================

def create_file(initial_conditions, file_dir, file_end, start_time = 1e-16, restart=False):
    """
    Function to create a triple.in file with input initial conditions
    :param initial_conditions: array of initial conditions from generate_inicon
    :param file_dir: directory to store triple.in files in
    :param file_end: file suffix (i.e. '1' to write triple.in1)
    :param start_time: system age to start simulations on
    :param restart: flag set to True if code is being re-started on Hoffman
    :return: writes a triple.inX file and saves to file_dir
    """
    # open new file
    filename = file_dir + '/' + 'triple.in' + file_end
    out = open(filename,'w')
    
    # write header
    out.write("#########################################################################\n")
    out.write("## Read in the initial orbital parameters of the triple system         ##\n")
    out.write("## Mass in M_sol, semimajor axis in AU                                 ##\n")
    out.write("#########################################################################\n")
    out.write("\n")
    
    if restart==False:
        out.write("__m1_______m2___________m3________R1(Rsun)____R2(Rsun)____Spin1P(day)__Spin2P(day)_ _beta(s_O_deg)__beta2(s_O_deg)__gamma(s_O_deg)__gamma2(s_O_deg)__a1______a2____e1____e2____g1(deg)__g2(deg)__i(deg)__age(Myr)___ :::\n")
    
    elif restart==True:
        out.write("__m1_______m2___________m3________R1(Rsun)____R2(Rsun)____Spin1P(day)__Spin2P(day)_    _beta(s_O_deg)__beta2(s_O_deg)__gamma(s_O_deg)__gamma2(s_O_deg)__a1______a2____e1____e2____g1(deg)__g2    (deg)__i(deg)__t_i(yr)__age(Myr)___ :::\n")
        
    out.write("\n")
    
    # add initial conditions
    out.write('  '.join([str(i) for i in initial_conditions])+'\n')
    
    # add header for flags
    out.write("\n")
    out.write("#########################################################################\n")
    out.write("## Flags                                                               ##\n")
    
    # quadrupole
    out.write("#########################################################################\n")
    out.write("\n")
    out.write("__QUADRUPOLE(0_yes_1_no)___ :::\n")
    out.write(str(quadrupole)+'\n') # for quadrupole?
    out.write("\n")
    
    # octupole
    out.write("#########################################################################\n")
    out.write("__OCTUPOLE(0_yes_1_no)___ :::\n")
    out.write("\n")
    out.write(str(octupole)+'\n') # for octupole?
    out.write("\n")
    
    # GR
    out.write("#########################################################################\n")
    out.write("__GR(0_yes_1_no)___ :::\n")
    out.write("\n")
    out.write(str(GR)+'\n') # for GR?
    out.write("\n")
    
    # tidal forces
    out.write("#########################################################################\n")
    out.write("__TF(0_yes_1_no)___ :::\n")
    out.write("\n")
    out.write(str(TF)+'\n') #for TF?
    out.write("\n")
    
    out.write("########################################################################\n")
    out.write("__ML1(0_yes_1_no)___ :::\n")
    out.write("\n")
    out.write(str(ML1)+'\n') # for ML1?
    out.write("\n")
    
    out.write("########################################################################\n")
    out.write("__MB1(0_yes_1_no)___ :::\n")
    out.write("\n")
    out.write(str(MB1)+'\n') # for MB1?
    out.write("\n")
    
    out.write("########################################################################\n")
    out.write("__ML2(0_yes_1_no)___ :::\n")
    out.write("\n")
    out.write(str(ML2)+'\n') # for ML2?
    out.write("\n")
    
    out.write("########################################################################\n")
    out.write("__MB2(0_yes_1_no)___ :::\n")
    out.write("\n")
    out.write(str(MB2)+'\n') # for MB2?
    out.write("\n")
    
    out.write("########################################################################\n")
    out.write("__SSE1(0_yes_1_no)___ :::\n")
    out.write("\n")
    out.write(str(SSE1)+'\n') # for SSE1?
    out.write("\n")
    
    out.write("########################################################################\n")
    out.write("__SSE2(0_yes_1_no)___ :::\n")
    out.write("\n")
    out.write(str(SSE2)+'\n') # for SSE2?
    out.write("\n")
    
    out.write("########################################################################\n")
    out.write("__SSE3(0_yes_1_no)___ :::\n")
    out.write("\n")
    out.write(str(SSE3)+'\n') # for SSE3?
    out.write("\n")
    
    out.write("########################################################################\n")
    out.write("__tMSMyr(0_yes_1_no)___ :::\n")
    out.write("\n")
    out.write(str(tMSMyr)+'\n') # for tMSMyr?
    out.write("\n")
    
    # add control parameters section
    out.write("#########################################################################\n")
    out.write("## Control Parameters                                                  ##\n")
    out.write("#########################################################################\n")
    out.write("\n")
    out.write("__eps___ :::\n")
    out.write(str(start_time)+"\n")
 #   out.write("1e-16\n")
    
    out.close()

# ==========================================================================

def generate_set_of_files(n, file_dir):
    """
    Creates a set of triple.in files by iteratively running create_file()
    :param n: number of triple.in files to generate (will create triple.in1-n)
    :param file_dir: directory to store triple.in files in
    :return: writes a n triple.inX file and saves to file_dir
    """
    for j in range(1,n+1):
        inicons = initial_conditions()
        file_end = str(j)
        print(initial_conditions)
        create_file(inicons, file_dir,file_end)

# ==========================================================================

def generate_restart_file(filename):

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
    file_end = str(filename[7:-4]) # store file number
    file_end = str([int(n) for n in filename.split() if n.isdigit()][0])
    create_file(final_conditions, 'unfinished_inicons', file_end, start_time = end_time, restart=True)








