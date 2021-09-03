import numpy as np


def vlr(x, y):
    return np.var(x - y, ddof=1) # - 2 * np.cov(x, y, ddof=1)[1, 0] + np.var(x, ddof=1) + np.var(y, ddof=1) 


def phi(x, y):
    return vlr(x, y) / np.var(y, ddof=1)


def rho(x,y):
    return 1 - (vlr(x, y) / (np.var(x, ddof=1) + np.var(y, ddof=1)))


def phs(x, y):
    r = rho(x, y)
    return (1 - r) / (1 + r)


def get_prop_matrix(mat, func):
    r, c = mat.shape
    corr_mat = np.zeros((c, c))
    for i in range(c):
        for j in range(c):
            corr_mat[i, j] = func(mat[:, i], mat[:, j])
    return corr_mat


def pairwise_vlr(mat):
    cov = np.cov(mat.T, ddof=1)
    diagonal = np.diagonal(cov)
    return -2 * cov + diagonal[:, np.newaxis] + diagonal


def pairwise_phi(mat):
    return pairwise_vlr(mat) / np.var(mat, axis=0, ddof=1)


def pairwise_rho(mat):
    variances = np.var(mat, axis=0, ddof=1)
    return 1 - (pairwise_vlr(mat) / np.add.outer(variances, variances))


def pairwise_phs(mat):
    r = pairwise_rho(mat)
    return (1 - r) / (1 + r)


# from skbio.stats.composition import clr
# counts = np.array([[12,2,3,4],[5,6,7,8],[9,11,12,13]])# print(get_prop_matrix(counts, vlr, np.log))
# counts = clr(counts)

# print(get_prop_matrix(counts, vlr))
# print(get_prop_matrix(counts, phi))
# print(get_prop_matrix(counts, phs))
# print(get_prop_matrix(counts, rho))

# print(pairwise_vlr(counts))
# print(pairwise_phi(counts))
# print(pairwise_phs(counts))
# print(pairwise_rho(counts))

# # compare to each other 
# print(np.isclose(get_prop_matrix(counts, vlr), pairwise_vlr(counts)).all())
# print(np.isclose(get_prop_matrix(counts, phi), pairwise_phi(counts)).all())
# print(np.isclose(get_prop_matrix(counts, phs), pairwise_phs(counts)).all())
# print(np.isclose(get_prop_matrix(counts, rho), pairwise_rho(counts)).all())

# # compare to propr
# from rpy2.robjects.packages import importr
# from rpy2.robjects import numpy2ri
# propr = importr("propr")
# numpy2ri.activate()

# print(np.isclose(propr.lr2vlr(counts), pairwise_vlr(counts)).all())
# print(np.isclose(propr.lr2phi(counts), pairwise_phi(counts)).all())
# print(np.isclose(propr.lr2phs(counts), pairwise_phs(counts)).all())
# print(np.isclose(propr.lr2rho(counts), pairwise_rho(counts)).all())
