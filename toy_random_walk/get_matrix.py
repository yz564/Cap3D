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

def generate_special_symmetric_matrix(size):
    # Generate a matrix with random off-diagonal elements
    matrix = np.random.uniform(low=-0.01, high=0, size=(size, size))
    matrix = np.tril(matrix, k=-1)  # Ensure lower triangular part is non-positive

    # Make the matrix symmetric by copying lower triangle to upper triangle
    matrix = matrix + matrix.T

    # Add a random number of large absolute values in each row
    for i in range(size):
        num_large_elements = np.random.randint(1, 11)  # Random number from 1 to 10
        for _ in range(num_large_elements):
            col = np.random.randint(0, size)
            while col == i:
                col = np.random.randint(0, size)
            large_value = np.random.uniform(low=-0.1, high=-10)
            matrix[i, col] = matrix[col, i] = large_value

        
    # Set diagonal elements as negative sum of remaining elements
    for i in range(size):
        diag_element = -np.sum(matrix[i]) + matrix[i, i]
        matrix[i, i] = diag_element


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
	
	
	
