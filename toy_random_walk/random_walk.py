import numpy as np
import concurrent.futures

def get_transiant_probability_matrix(matrix):
	rows, cols = matrix.shape
	p_matrix=np.copy(matrix)
	for col in range(cols):
		t=-matrix[col,col]
		if t !=0:
			p_matrix[:,col]/=t
			p_matrix[col,col]=0.0
	return np.transpose(p_matrix)

def single_random_walk(P, T, initial_id, n_signal,N=100):
	E=0
	current_id=initial_id
	all_ids=np.arange(len(T))
	next_id=np.random.choice(all_ids,p=P[current_id])
	if next_id<n_signal:
		return E,current_id
	first_floating_id=next_id
	last_floating_id=next_id
	current_id=next_id
	#steps=0
	while current_id>=n_signal:
		last_floating_id=current_id
		current_id=np.random.choice(all_ids,p=P[current_id])
	#if steps==N:
	#	return single_random_walk(P, T, initial_id, n_signal)
	E=(T[initial_id]/T[first_floating_id])*T[last_floating_id]
	return E, current_id
	

def random_walk(matrix, n_signal, N=10000):
	B = np.zeros((n_signal, n_signal))
	A=matrix[:n_signal,:n_signal]
	P=get_transiant_probability_matrix(matrix)
	T=np.diag(matrix)
	for i in range(n_signal):
		iteration = 1
		E=np.zeros(n_signal)
		while iteration <= N:
			value,j=single_random_walk(P, T, i, n_signal)
			E[j]+=value
			iteration+=1
		B[i,:]=E/iteration
	return A-B
