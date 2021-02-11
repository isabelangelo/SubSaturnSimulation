import pandas as pd
import os
import matplotlib.pyplot as plt

# set of models to use
run_number = 6

# store values for Saturn
a_subsat = 0.1974
e_subsat = 0.816

class Model(object):
    """
    Model object for OSPE output
    
    :param n_model: model number, i.e. 1 for triple.in1, output_1.txt (int)
    
    :attr: # list attributes here
    """
    def __init__(self, n_model, sort=True):

        # get data from separate output files
        self.n_model=n_model
        dataframes = []
        for i in [1,2]:
            output_filepath = 'Run'+str(run_number)+'/Output/output_'+str(n_model)+'.txt.loop'+str(i)
            if os.path.isfile(output_filepath):
                data = pd.read_csv(output_filepath,sep='\t')
                dataframes.append(data)

        final_output_filepath = 'Run'+str(run_number)+'/Output/output_'+str(n_model)+'.txt'
        final_data = pd.read_csv(final_output_filepath,sep='\t')
        dataframes.append(final_data)
        
        # combine into one dataframe, sort by t
        df = pd.concat(dataframes, ignore_index=True).sort_values(' t ')
        
        # remove rows with 'evaluating at' print statement
        indexes_to_drop = []
        for index, row in df.iterrows():
            current_row = row['sur ']
            if isinstance(current_row, str) and 'evaluating' in current_row:
                indexes_to_drop.append(index)
        self.data=df.drop(indexes_to_drop)
            
        # store important values
        self.time = self.data[' t '].to_numpy()
        self.e1 = self.data[' e1 '].to_numpy(dtype='float32')
        self.e2 = self.data[' e2 '].to_numpy(dtype='float32')
        self.i1 = self.data[' i1 '].to_numpy()
        self.i2 = self.data[' i2 '].to_numpy()
        self.a1 = self.data[' a1 '].to_numpy()
        self.a2 = self.data[' a2 '].to_numpy()
        self.beta1 = self.data[' beta '].to_numpy()
        self.beta2 = self.data[' beta2 '].to_numpy()
        
        # store endstates
        self.sur = self.data['sur '].to_numpy()[-1]
        self.sur2 = self.data[' sur2 '].to_numpy()[-1]
    
    # create plot of model timeseries
    def plot_timeseries(self,m='.'):
        # make plot
        fig,ax = plt.subplots(4,2,figsize=(10,10))
        ax[0,0].plot(self.time/1e9,self.e1,'b', marker=m,linewidth=0.7)
        ax[0,0].set_ylabel('eccentricity',fontsize=15,family='serif')
        ax[0,0].set_title('Kepler-1565b',fontsize=15,family='serif')
        ax[0,1].plot(self.time/1e9,self.e2,'b', marker=m,linewidth=0.7)
        ax[0,1].set_title('Outer Companion',fontsize=15,family='serif')
        ax[1,0].plot(self.time/1e9,self.a1,'r', marker=m,linewidth=0.7)
        ax[1,0].set_ylabel('semi-major \n axis (au)',fontsize=15,family='serif')
        ax[1,1].plot(self.time/1e9,self.a2,'r', marker=m,linewidth=0.7)
        ax[2,0].plot(self.time/1e9,self.i1,'c', marker=m,linewidth=0.7)
        ax[2,0].set_ylabel('inclination (deg)',fontsize=15,family='serif')
        ax[2,1].plot(self.time/1e9,self.i2,'c', marker=m,linewidth=0.7)
        ax[3,0].plot(self.time/1e9,self.beta1,'k', marker=m,linewidth=0.7)
        ax[3,0].set_ylabel('$\lambda$ (deg)',fontsize=15,family='serif')
        ax[3,0].set_xlabel('time [Gyr]',fontsize=15,family='serif')
        ax[3,1].plot(self.time/1e9,self.beta2,'k', marker=m,linewidth=0.7)
        ax[3,1].set_xlabel('time [Gyr]',fontsize=15,family='serif')
        
        # remove space between panels
        plt.subplots_adjust(hspace=0)
        plt.subplots_adjust(wspace=0)
        
        # remove offset labels for sma
        ax[1,0].ticklabel_format(useOffset=False)
        ax[1,1].ticklabel_format(useOffset=False)
        
        # move outer planet ticks to the right
        for ax_i in ax[:,1]:
            ax_i.yaxis.tick_right()
            
#        for axes in ax:
#            for ax_i in axes:
#                ax_i.set_xlim(0,2.5)
            
               
    
        # set figure title
        fig.suptitle('Model '+str(self.n_model), fontsize=15,family='serif')
        
        