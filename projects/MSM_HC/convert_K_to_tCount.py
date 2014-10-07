from scipy.io import mmread,mmwrite
from scipy import sparse

c = mmread('tCounts.MSM.mtx')
k = mmread('K_lifsonroig.gamma_0.12.mtx')

for row in range(c.shape[0]):
    k[row] = k[row]*c.data[c.row==row].sum()

k_sparse = sparse.coo_matrix(k)

mmwrite('tCounts.MaxCal.mtx',k_sparse)

