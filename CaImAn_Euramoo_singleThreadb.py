import sys
import numpy as np
from scipy.io import savemat
import caiman as cm
from caiman.components_evaluation import evaluate_components
import caiman.source_extraction.cnmf as cnmf
from caiman.source_extraction.cnmf.utilities import extract_DF_F
fnames=[sys.argv[1]]	#name of the movie
final_frate=int(sys.argv[2]) # frame rate in Hz
K=int(sys.argv[3]) # number of neurons expected per patch, that seems to work well
n_processes = 5 # if using the intel nodes
single_thread=True   
dview=None
#is_patches=True
gSig=[3,3] # expected half size of neurons, works for nuclear GCaMP
merge_thresh=0.9 # merging threshold, max correlation allowed
p=2 #order of the autoregressive system
downsample_factor=1 # use .2 or .1 if file is large and you want a quick answer
final_frate=final_frate*downsample_factor
init_method = 'greedy_roi'
alpha_snmf=None #10e2  # this controls sparsity
idx_xy=None
base_name=fnames[0]
fname_new=cm.save_memmap_each(fnames, dview=dview,base_name=base_name, resize_fact=(1, 1, downsample_factor), remove_init=0,idx_xy=idx_xy )
#fname_new=cm.save_memmap_join(fname_new,base_name='Yr', n_chunks=n_processes, dview=dview)
fname_new = fname_new[0]
Yr,dims,T=cm.load_memmap(fname_new)
d1, d2 = dims
images = np.reshape(Yr.T, [T] + list(dims), order='F')
Y = np.reshape(Yr, dims + (T,), order='F')
nb_back=2
rf = 100  # half-size of the patches in pixels. rf=25, patches are 50x50
stride = 10  # amounpl.it of overlap between the patches in pixels
K = K/10  # number of neurons expected per patch
p = 1  # order of the autoregressive system
save_results = False
#%% RUN ALGORITHM ON PATCHES
cnm = cnmf.CNMF(n_processes, k=K, gSig=gSig, merge_thresh=merge_thresh, p=0, dview=dview, Ain=None, rf=rf, stride=stride, memory_fact=1, method_init=init_method, alpha_snmf=alpha_snmf, only_init_patch=True, gnb=1,method_deconvolution='oasis')
cnm = cnm.fit(images)
A_tot = cnm.A
C_tot = cnm.C
YrA = cnm.YrA
bl = cnm.b
f_tot = cnm.f
sn = cnm.sn
S= cnm.S

traces=C_tot+YrA
	
	
tB = np.minimum(-2,np.floor(-5./30*final_frate))
tA = np.maximum(5,np.ceil(25./30*final_frate))
fitness_raw, fitness_delta, erfc_raw,erfc_delta,r_values,num_significant_samples = evaluate_components(Y,traces,A_tot,C_tot,bl,f_tot, final_frate, remove_baseline = True, N = 5, robust_std = False, Athresh = 0.1, Npeaks = 5, thresh_C = 0.3)
idx_components_r=np.where(r_values>=.5)[0]
idx_components_raw=np.where(fitness_raw<-40)[0]        
idx_components_delta=np.where(fitness_delta<-20)[0]
idx_components=np.union1d(idx_components_r,idx_components_raw)
idx_components=np.union1d(idx_components,idx_components_delta)  	
idx_components_bad=np.setdiff1d(range(len(traces)),idx_components)	
	
savemat(fnames[0][:-4]+'_output_analysis_patches.mat',mdict={'ROIs':A_tot,'DenoisedTraces':C_tot,'Baseline':bl, 'Noise':YrA, 'Spikes': S, 'idx_components':idx_components})
A_tot = A_tot.tocsc()[:, idx_components]
C_tot = C_tot[idx_components]

cnm = cnmf.CNMF(n_processes, k=A_tot.shape, gSig=gSig, merge_thresh=merge_thresh, p=p, dview=dview, Ain=A_tot, Cin=C_tot, f_in=f_tot, rf=None, stride=None, method_deconvolution='oasis')
cnm = cnm.fit(images)

A, C, b, f, YrA, sn, S = cnm.A, cnm.C, cnm.b, cnm.f, cnm.YrA, cnm.sn cnm.S
traces = C + YrA
fitness_raw, fitness_delta, erfc_raw, erfc_delta, r_values, significant_samples = evaluate_components(Y, traces, A, C, b, f, final_frate, remove_baseline=True, N=5, robust_std=False, Athresh=0.1, Npeaks=5,  thresh_C=0.3)

idx_components_r = np.where(r_values >= .95)[0]
idx_components_raw = np.where(fitness_raw < -100)[0]
idx_components_delta = np.where(fitness_delta < -100)[0]
idx_components = np.union1d(idx_components_r, idx_components_raw)
idx_components = np.union1d(idx_components, idx_components_delta)
#idx_blobs = np.intersect1d(idx_components, idx_blobs)
idx_components_bad = np.setdiff1d(list(range(len(traces))), idx_components)
C_dff = extract_DF_F(Yr, A.tocsc()[:, idx_components], C[idx_components, :], cnm.bl[idx_components], quantileMin = 8, frames_window = 200, dview = dview)

savemat(fnames[0][:-4]+'_output_analysis_matlab.mat',mdict={'ROIs':A,'DenoisedTraces':C,'Baseline':b, 'Noise':YrA, 'Spikes': S, 'idx_components':idx_components, 'DFF':C_dff})