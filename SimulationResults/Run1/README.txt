This code processes OSPE output in the following ways:

    1. figures out which OSPE runs are unfinished and prepares re-runs
    2. performs large-scale data analysis on finished runs
    
The process to restart OSPE on unfinished runs is as follows:
 1 run unfinished_runs.py : generates new triple.in files, moves old outputs
 2 *run prepare_rerun.sh in hoffman : moves firstrun inicons and output files to make room for new ones
 3 run scp_unfinished_inicon.sh : moves new inicons to hoffman (you will need to untar after you scp and before you run restart code!!)
 4 *run RunNewOSPE_unfinished_runs.prl : runs OSPE on new computer (be sur to run the one in the "Restart" directory!
 5 scp new output files here to prepare for data analysis
 6 run plot_OSPE_results : plots histogram results of run
 * = takes in the file numbers as input