import os,sys
import numpy as np
from msmbuilder import tpt
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
from sklearn.externals import joblib

ids = range(6383,6391)

for p_id in ids:

    msm = joblib.load('MSMs-{pid}-macro40/MSMs-{pid}-macro40.pkl'.format(pid=p_id))
    committors = tpt.committors(sources=[4],sinks=[13],msm=msm)
    net_fluxes = tpt.net_fluxes(sources=[4],sinks=[13],msm = msm,for_committors=committors)
    top_10_paths_sub,top_10_fluxes_sub = tpt.paths(sources=[4],sinks=[13],net_flux=net_fluxes,num_paths=10,remove_path='subtract')
    top_10_paths_bot,top_10_fluxes_bot = tpt.paths(sources=[4],sinks=[13],net_flux=net_fluxes,num_paths=10,remove_path='bottleneck')
    output_dir = "Data-{pid}-macro40".format(pid=p_id)
    print top_10_paths_sub,top_10_fluxes_sub
    print top_10_paths_bot,top_10_fluxes_bot

    fn_committors = os.path.join(output_dir,"committors.npy")
    fn_n_fluxes = os.path.join(output_dir,"net_Fluxes.mtx")
    fn_10_paths_sub = os.path.join(output_dir,"top_paths_substract.npy")
    fn_10_fluxes_sub = os.path.join(output_dir,"top_fluxes_substract.npy")
    fn_10_paths_bot = os.path.join(output_dir,"top_paths_bottleneck.npy")
    fn_10_fluxes_bot = os.path.join(output_dir,"top_fluxes_bottleneck.npy")

    np.save(fn_committors,committors)
    mmwrite(fn_n_fluxes,coo_matrix(net_fluxes))
    np.save(fn_10_paths_sub,top_10_paths_sub)
    np.save(fn_10_fluxes_sub,top_10_fluxes_sub)
    np.save(fn_10_paths_bot,top_10_paths_bot)
    np.save(fn_10_fluxes_bot,top_10_fluxes_bot)

