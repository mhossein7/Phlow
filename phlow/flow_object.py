import contextlib
import io
import numpy as np 
from matplotlib import pyplot as plt
import FlowCal
from matplotlib.legend_handler import HandlerTuple 
import matplotlib.patches as mpatches


with contextlib.redirect_stderr(io.StringIO()):
    import seaborn as sns




def sample_seaborn_palette(
    palette,
    x,
    min_frac=0.15,
    max_frac=0.95,
):
    """
    Sample x colors evenly from a seaborn palette, avoiding very light/dark extremes.

    Parameters
    ----------
    palette : str or sequence
        Seaborn palette name or palette object.
    x : int
        Number of colors to sample.
    min_frac : float, optional
        Lower bound of palette range (0–1). Default avoids near-white.
    max_frac : float, optional
        Upper bound of palette range (0–1). Default avoids near-black.
    """

    if x <= 0:
        raise ValueError("x must be a positive integer")
    if not 0 <= min_frac < max_frac <= 1:
        raise ValueError("min_frac and max_frac must satisfy 0 ≤ min < max ≤ 1")

    cmap = sns.color_palette(palette, as_cmap=True)

    values = np.linspace(min_frac, max_frac, x)
    colors = [cmap(v)[:3] for v in values]

    return colors


def rowsCols(a):
    if len(a.shape) > 1:
        rows = a.shape[0]
        cols = a.shape[1]
    else:
        rows = 1
        cols = a.shape[0]
    return (rows, cols)

# Blahut-Arimoto alghorithm
def compute_capacity(p):
    # Check for negative entries
    if np.any(p < 0):
        print('Error: some entry in the input matrix is negative')
        return 0
    
    # Check for zero columns
    if np.any(np.sum(p, axis=0) == 0):
        print('Error: there is a zero column in the input matrix')
        return 0
    
    # Check for zero rows
    row_sum = np.sum(p, axis=1)
    if np.any(row_sum == 0):
        print('Error: there is a zero row in the input matrix')
        return 0
    else:
        p = np.diag(1 / row_sum) @ p  # Normalize rows to sum to 1
    
    m, n = p.shape
    
    # Initial distribution
    r = np.ones(m) / m
    q = np.zeros((m, n))
    error_tolerance = 1e-7 / m
    
    # Normalize rows of p
    for i in range(m):
        p[i, :] /= np.sum(p[i, :])
    
    for _ in range(100000):
        for j in range(n):
            q[:, j] = r * p[:, j]
            q[:, j] /= np.sum(q[:, j])
        
        r1 = np.prod(q ** p, axis=1)
        r1 /= np.sum(r1)
        
        if np.linalg.norm(r1 - r) < error_tolerance:
            break
        r = r1
    
    # Compute capacity
    C = 0
    for i in range(m):
        for j in range(n):
            if r[i] > 0 and q[i, j] > 0:
                C += r[i] * p[i, j] * np.log(q[i, j] / r[i])
    
    return C/np.log(2)


def log_discretizer(X,n,t_low,t_high):
    '''
    X: data in a m*d dimension (m: # of conditions, d: # of data points)
    n: number of bins for discretizing the data
    t_low: lower threshold for discretization
    t_high: higher threshold for discretization
    '''
    m,d = rowsCols(X)
    p = np.zeros((m,n-1))
    
    edges = np.logspace(t_low,t_high,num=n)
    binned_data = np.digitize(X,edges)
    
    for i in range(m):
        for j in range(n-1):
            p[i,j] = np.sum(binned_data[i]==j+1)/d
    
    non_zero_cols = ~np.all(p == 0, axis=0)
    p_nz = p[:,non_zero_cols]

    return p_nz

def MI_calculator(Data,num_cell,num_inputs,bins=[5,50],t_low=0.1,t_high=4):
    
    n = np.arange(bins[0],bins[1])
    C = np.zeros(len(n))
    
    for i, ni in enumerate(n):
        if isinstance(num_cell,(list,np.ndarray)):
            C_k = np.zeros(len(num_cell))
            for k,num_k in enumerate(num_cell):
                X = np.zeros((num_inputs,num_k))
                for j in range(num_inputs):
                    X[j,:] = Data[f'Cond{j}'][:num_k]
                C_k[k] = compute_capacity(log_discretizer(X,ni,t_low,t_high))
            x = 1/num_cell.astype(float)
            coeffs = np.polyfit(x,C_k,1)
            C[i] = coeffs[-1]
            
        else:        
            X = np.zeros((num_inputs,num_cell))
            for j in range(num_inputs):
                    X[j,:] = Data[f'Cond{j+1}'][:num_cell]
            C[i] = compute_capacity(log_discretizer(X,ni,t_low,t_high))

    C_unbiased = np.mean(C[20:])
    return C_unbiased



class flow_unit():
    def __init__(self,strain_name,details):
        self.name = strain_name
        self.detailed_name = details
        self.paths = {}
        self.multi_gain = False
        self.gains= {}
        self.num_conds = 4
        self.cond_labels = ['0','21','52','208']
        self.gains_for_conds = {f'Cond{i}':None for i in range(1,self.num_conds+1)}
        self.raw_gfp_data = {}
        self.raw_rfp_data = {}
        self.Data_gfp = {f'Cond{i}':[] for i in range(1,self.num_conds+1)}
        self.Data_rfp = {f'Cond{i}':[] for i in range(1,self.num_conds+1)}
        self.Data_corrected = {f'Cond{i}':[] for i in range(1,self.num_conds+1)}
        self.istrips = False
        self.log = True
        self.green_palette =  sample_seaborn_palette('Greens',self.num_conds)
        self.red_palette = sample_seaborn_palette('Reds',self.num_conds)
        self.blue_palette = sample_seaborn_palette('Blues',self.num_conds)
    
    
    
    def set_paths(self,paths):
        if isinstance(paths,list) and len(paths)>1:
            self.multi_gain = True
            for i in range(len(paths)):
                self.paths[f'{i+1}'] = paths[i]
        else:
            self.paths['1'] = paths[0] if isinstance(paths,list) else paths
    
    def set_gain(self,gains):
        if not self.paths:
            raise Exception('Please set paths to files first.')
        else:
            if isinstance(gains,list) and len(gains)>1:
                if len(gains)!=len(self.paths):
                    raise Exception('The number of gains provided does not match the number of file paths')
                else:
                    for i in range(len(gains)):
                        self.gains[f'{i+1}'] = gains[i]
            else:
                self.gains['1'] = gains[0] if isinstance(gains,list) else gains
                self.gains_for_conds = {f'Cond{i}':gains[0] if isinstance(gains,list) else gains for i in range(1,self.num_conds+1)}
    
    def set_number_of_conditions(self,cond):
        self.num_conds = cond
        
        self.gains_for_conds = {f'Cond{i}':None for i in range(1,self.num_conds+1)}
        self.Data_gfp = {f'Cond{i}':[] for i in range(1,self.num_conds+1)}
        self.Data_rfp = {f'Cond{i}':[] for i in range(1,self.num_conds+1)}
        self.Data_corrected = {f'Cond{i}':[] for i in range(1,self.num_conds+1)}
        
        self.green_palette =  sample_seaborn_palette('Greens',self.num_conds)
        self.red_palette = sample_seaborn_palette('Reds',self.num_conds)
        self.blue_palette = sample_seaborn_palette('Blues',self.num_conds)
    
    def set_gains_for_conds(self,gains_list):
        if not self.multi_gain:
            raise Exception('Please set paths to files first')
        if not self.gains:
            raise Exception('Plese set the gains first.')
        if len(gains_list)!= self.num_conds:
            raise Exception('The number of provided gains for each condition should match the number of total conditions')
        else:
            for i in range(1,self.num_conds+1):
                if gains_list[i-1] in self.gains.values():
                    self.gains_for_conds[f'Cond{i}'] = gains_list[i-1]
                else: 
                    raise Exception('At least one of the provided gains do not belong to set gains.')
    
    def set_triplicates(self,istrips):
        self.istrips = istrips
    
    def set_log(self,islog):
        self.log = islog
    
    def read_data(self):
        if not self.gains or not self.paths:
            raise Exception('Please set paths to files and/or set gains')
        
        if self.istrips:
            file_inds = np.arange(1,self.num_conds*3+1)
        else: file_inds = np.arange(1,self.num_conds+1)
        
        if len(self.gains.values())>1:
            for i in range(len(self.gains.values())):
                self.raw_gfp_data[f'{i+1}'] = [[] for _ in range(self.num_conds)]
                self.raw_rfp_data[f'{i+1}'] = [[] for _ in range(self.num_conds)]
        else: 
            self.raw_gfp_data['1'] = [[] for _ in range(self.num_conds)]
            self.raw_rfp_data['1'] = [[] for _ in range(self.num_conds)]
        
        for key,path in self.paths.items():
            gfp = []
            rfp = []
            for i in file_inds:
                pre_sample = FlowCal.io.FCSData(path+f'{i}'+'.fcs')
                sample = FlowCal.gate.high_low(pre_sample, channels=['FSC-HLin', 'SSC-HLin'], low=25)
                sample = FlowCal.gate.high_low(sample, channels=['GRN-B-HLin','RED-B-HLin'], low=1)
                gfps = sample[:,'GRN-B-HLin']
                rfps = sample[:,'RED-B-HLin']
                gfp.append(gfps)
                rfp.append(rfps)
                
            if not self.istrips:
                for i in range(len(gfp)):
                    self.raw_gfp_data[key][i] = gfp[i]
                    self.raw_rfp_data[key][i] = rfp[i]
            else:
                for i in range(self.num_conds):
                    self.raw_gfp_data[key][i] = np.hstack(([gfp[j] for j in range(i*3,(i+1)*3)]))
                    self.raw_rfp_data[key][i] = np.hstack(([rfp[j] for j in range(i*3,(i+1)*3)]))
        
    def compile_data(self,correction = True):
        '''
        Here, the data stored in self.raw_data will be compiled into self.Data based on the gains the user has chosen for each condition
        '''
        for i , (key , val) in enumerate(self.gains_for_conds.items()):
            grp = [k for k,v in self.gains.items() if v==val][0]
            if self.log:
                self.Data_gfp[key] = np.log10(np.array(self.raw_gfp_data[grp][i]) * 8/val)
                self.Data_rfp[key] = np.log10(np.array(self.raw_rfp_data[grp][i]))
            else:    
                self.Data_gfp[key] = np.array(self.raw_gfp_data[grp][i]) * 8/val
                self.Data_rfp[key] = np.array(self.raw_rfp_data[grp][i]) 
        
        if correction:
            tlog = self.log
            self.set_log(False)
            self.compile_data(False)
            for i , (key , val) in enumerate(self.gains_for_conds.items()):
                self.Data_corrected[key] = np.log10(self.Data_gfp[key]/self.Data_rfp[key]) if tlog else self.Data_gfp[key]/self.Data_rfp[key]
            if tlog:
                self.set_log(True)
                self.compile_data(False)
    
    def compute_pop_metrics(self):
        if all(len(v) == 0 for v in self.Data_gfp.values()): raise Exception('Please load and compile data first')
        tlog = True
        if not self.log:
            tlog = False
            self.log = True
            self.compile_data()
        
        GFP_metrics = {'Mean':[],'Var':[],'CV':[],'r':[]}
        RFP_metrics = {'Mean':[],'Var':[],'CV':[],'r':[]}
        Corrected_metrics = {'Mean':[],'Var':[],'CV':[]}
        for i in range(self.num_conds):
            gdata = self.Data_gfp[f'Cond{i+1}']
            rdata = self.Data_rfp[f'Cond{i+1}']
            cdata = self.Data_corrected[f'Cond{i+1}']
            GFP_metrics['Mean'].append(np.mean(gdata))
            GFP_metrics['Var'].append(np.var(gdata))
            GFP_metrics['CV'].append(np.std(gdata)/np.mean(gdata))
            GFP_metrics['r'].append(np.corrcoef(rdata,gdata))
            RFP_metrics['Mean'].append(np.mean(rdata))
            RFP_metrics['Var'].append(np.var(rdata))
            RFP_metrics['CV'].append(np.std(rdata)/np.mean(rdata))
            RFP_metrics['r'].append(np.corrcoef(rdata,gdata))
            Corrected_metrics['Mean'].append(np.mean(cdata))
            Corrected_metrics['Var'].append(np.var(cdata))
            Corrected_metrics['CV'].append(np.std(cdata)/np.mean(cdata))
        
        if not tlog:
            self.log = False
            self.compile_data()
        
        return (RFP_metrics,GFP_metrics,Corrected_metrics)
                
                
    def plot_histogram(self,axes = None,is_legend=False,is_norm = True):
        if all(len(v) == 0 for v in self.Data_gfp.values()): raise Exception('Please load and compile data first')
        if axes is not None and is_norm and len(axes)!=3: raise Exception('The input should have 3 axes')
        if axes is not None and not is_norm and len(axes)!=2: raise Exception('The input should have 2 axes')
        if axes is None: _ , axes = plt.subplots(3 if is_norm else 2)
        for i,data in enumerate(self.Data_rfp.values()):
            sns.kdeplot(data,ax=axes[0],color = self.red_palette[i],alpha = 0.3,fill=True,linewidth=1.5,common_norm=False,legend = False)
        for i,data in enumerate(self.Data_gfp.values()):
            sns.kdeplot(data,ax=axes[1],color = self.green_palette[i],alpha = 0.3,fill=True,linewidth=1.5,common_norm=False,legend = False)
        if is_norm:
                for i,data in enumerate(self.Data_corrected.values()):
                    sns.kdeplot(data,ax=axes[2],color = self.blue_palette[i],alpha = 0.3,fill=True,linewidth=1.5,common_norm=False,legend = False)
        axes[0].set_ylim([0,3])
        axes[1].set_ylim([0,3])
        if is_norm: axes[2].set_ylim([0,3])
        sns.despine()
        
        if is_legend:
            conditions = self.cond_labels
            # Create Patch objects manually
            # This ensures we have clean handles regardless of matplotlib/seaborn versions
            handles_rfp = [mpatches.Patch(facecolor=c, edgecolor=c, alpha=0.3, label=l) 
                        for c, l in zip(self.red_palette, conditions)]
            handles_gfp = [mpatches.Patch(facecolor=c, edgecolor=c, alpha=0.3, label=l) 
                        for c, l in zip(self.green_palette, conditions)]
            handles_bfp = [mpatches.Patch(facecolor=c, edgecolor=c, alpha=0.3, label=l) 
                        for c, l in zip(self.blue_palette, conditions)]

            # Zip them together for the side-by-side effect
            combined_handles = list(zip(handles_rfp, handles_gfp,handles_bfp))

            # Create the Legend
            legend = axes[0].legend(
                combined_handles, 
                conditions, # Labels
                loc="upper left", 
                bbox_to_anchor=(1, 1), 
                title=r'Green light intensity'+'\n'+ r'($\mu$W/cm$^{2}$)',
                alignment = 'center',
                handler_map={tuple: HandlerTuple(ndivide=None)}, 
                frameon=True,
                facecolor='white',
                framealpha=1.0
            )
            title = legend.get_title()
            title.set_horizontalalignment('center')  # Center-align the title text
            title.set_verticalalignment('bottom')    # Align the title at the top of the legend
            
        
        return axes
    
    def plot_scatter(self,axes=None,isstd = False,is_norm = True,marker = 'o'):
        if all(len(v) == 0 for v in self.Data_gfp.values()): raise Exception('Please load and compile data first')
        if axes is not None and is_norm and len(axes)!=3: raise Exception('The input should be a vector of three axes')
        if axes is not None and not is_norm and len(axes)!=2: raise Exception('The input should have 2 axes')
        if axes is None: _ , axes = plt.subplots(3 if is_norm else 2)
        
        RFP_mets , GFP_mets, Norm_mets = self.compute_pop_metrics()
        rfp_means = RFP_mets['Mean']
        gfp_means = GFP_mets['Mean']
        norm_means = Norm_mets['Mean']
        axes[0].scatter(np.arange(self.num_conds),rfp_means,color='red',marker=marker)
        axes[1].scatter(np.arange(self.num_conds),gfp_means,color = 'green',marker=marker)
        axes[0].set_xticks(np.arange(self.num_conds))
        axes[1].set_xticks(np.arange(self.num_conds))
        axes[0].set_xticklabels(self.cond_labels)
        axes[1].set_xticklabels(self.cond_labels)
        if is_norm: axes[2].scatter(rfp_means,gfp_means,color = 'blue',marker=marker)
        return axes
            
    
    def compute_CC(self,num_cells=None,bins=[5,50],t_low=[0.1,0.1,-1],t_high=[3,6,5]):
        if all(len(v) == 0 for v in self.Data_gfp.values()): raise Exception('Please load and compile data first')
        if self.log : 
        
            self.set_log(False)
            self.compile_data()
            if num_cells is None:
                num_cells = min([len(v) for v in self.Data_rfp.values()])
                
            CC_rfp = MI_calculator(self.Data_rfp,num_cells,self.num_conds,bins,t_low[0],t_high[0])
            CC_gfp = MI_calculator(self.Data_gfp,num_cells,self.num_conds,bins,t_low[1],t_high[1])
            CC_corrected = MI_calculator(self.Data_corrected,num_cells,self.num_conds,bins,t_low[2],t_high[2])
            self.set_log(True)
            self.compile_data()
        else:
    
            CC_rfp = MI_calculator(self.Data_rfp,num_cells,self.num_conds,bins,t_low[0],t_high[0])
            CC_gfp = MI_calculator(self.Data_gfp,num_cells,self.num_conds,bins,t_low[1],t_high[1])
            CC_corrected = MI_calculator(self.Data_corrected,num_cells,self.num_conds,bins,t_low[2],t_high[2])
        return [CC_rfp,CC_gfp,CC_corrected]
        
        
