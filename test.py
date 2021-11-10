
import numpy as np
import scipy.sparse as sp
import time

from mbest import *

if __name__ == '__main__':
	
	in_vals = np.random.randn(100,3) + 10
	out_vals = np.random.randn(100,3) + 12

	N = in_vals.shape[0]
	M = out_vals.shape[0]
	d = N*(M+1)

	rad = 1

	st = time.time()

	A, Aeq, b, beq = build_initial_constr(N, M)

	f = build_C(in_vals, out_vals, rad)

	K = 50
	x, xv = MBest(f, A, b, Aeq, beq, K)

	rt = time.time() - st

	print('Finished in ', np.round(rt, 2), 'seconds.')
	print(xv)


