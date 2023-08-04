import numpy as np

def generate_symmetric_matrix(size):
    if size < 2:
        raise ValueError("Matrix size must be at least 2x2")
    
    matrix = np.zeros((size, size))
    
    for i in range(size):
        for j in range(i, size):
            if i == j:
                matrix[i, j] = np.random.uniform(0, 100)  # Positive diagonal elements
            else:
                temp_mean=matrix[i, i]/size
                temp_std_dev=temp_mean*np.random.uniform(0,2)
                value = np.random.normal(temp_mean, temp_std_dev)  # Negative non-diagonal elements
                matrix[i, j] = -abs(value)
                matrix[j, i] = -abs(value)
    
    # Ensure row and column sums are zero
    row_sums = np.sum(matrix, axis=1)
    
    for i in range(size):
        matrix[i, i] -= row_sums[i]
    
    
    return matrix

def floating_net_reduction(matrix, n_signal):
	n = matrix.shape[0]
	if n<n_signal:
		raise ValueError("n_signal should be smaller than the size of the matrix.")
	A=matrix[:n_signal,:n_signal]
	X=matrix[:n_signal,n_signal:]
	Y=matrix[n_signal:,:n_signal]
	Z=matrix[n_signal:,n_signal:]
	
	Z_inv=np.linalg.inv(Z)
	C = A - np.dot(np.dot(X, Z_inv), Y)
	return C
	
def using_neumann_series(matrix, n_signal, max_iterations=1000, tolerance=1e-5):
	n = matrix.shape[0]
	if n<n_signal:
		raise ValueError("n_signal should be smaller than the size of the matrix.")
	A=matrix[:n_signal,:n_signal]
	X=matrix[:n_signal,n_signal:]
	Y=matrix[n_signal:,:n_signal]
	Z=matrix[n_signal:,n_signal:]
	
	T=np.diag(1.0/np.diag(Z))
	I=np.eye(n-n_signal)
	K=I-np.dot(Z, T)
	result=I
	term = K
	iteration = 1
	while iteration <= max_iterations and np.linalg.norm(term, 'fro') >= tolerance:
		result += term
		term = np.dot(term, K)
		iteration += 1
	Z_inv=np.dot(T, result)
	C = A - np.dot(np.dot(X, Z_inv), Y)
	return C, iteration
	
	
	
