import numpy as np
import math 
from random import random
import time

n = 15

# in all vectors and loops the element in position 0 is ignored
# we use indices from 1 to n, not from 0 to n-1

states = np.zeros(n+1)

for i in range(1,6):
	states[i] = 1.

print("states",states)

# preprocessing 
w = np.zeros(n+2)
r = np.zeros(n+2)
c = np.zeros(n+2)
const_c = 0
for i in reversed(range(1,n+1)):
	w[i] = w[i+1] + states[i]
	c[i] = (n-i+1)*states[i] - w[i]
	r[i] = r[i+1] + c[i]
	const_c += c[i]

print("w",w)
print("c",c)
print("r",r)
print("const_c",const_c)

# start sampling 
m = 10000000

pairs_linear = dict()
pairs_binary = dict()
source_linear = dict()
source_binary = dict()
target_linear = dict()
target_binary = dict()


print("linear sampling...")

# linear sampling 
start = time.time()
for i in range(m):
	s = 0
	t = 0
	rand_val = random()
	for j in range(1,n+1):
		k = (const_c-r[j+1])/const_c
		if rand_val <= k:
			s = j
			break
	rand_val = random()
	for j in range(s,n+1):
		k = ((j-s+1)*states[s] - w[s] + w[j+1])/((n-s+1)*states[s] - w[s])
		if rand_val <= k:
			t = j
			break
	if (s,t) in pairs_linear:
		pairs_linear[(s,t)] = pairs_linear[(s,t)] + 1./m
	else:
		pairs_linear[(s,t)] = 1./m
	if s in source_linear:
		source_linear[s] = source_linear[s] + 1./m
	else:
		source_linear[s] = 1./m
	if t in target_linear:
		target_linear[t] = target_linear[t] + 1./m
	else:
		target_linear[t] = 1./m
	#print(rand_val)


end = time.time()

print()
print("source_linear\n",source_linear)
print()
print("target_linear\n",target_linear)
print()
print("pairs_linear\n",pairs_linear)
print()
print("done in",end - start)


print("binary sampling...")

# binary sampling 
start = time.time()
for i in range(m):
	s = 0
	t = 0

	# select source 
	a = 1 
	b = n 
	d = math.floor((a+b)/2.)
	rand_val = random()

	while a <= b:
		k = (const_c-r[d+1])/const_c
		if rand_val <= k:
			b = d-1
		else:
			a = d+1
		d = math.floor((a+b)/2.)
	s = b+1

	# select target  
	a = s
	b = n 
	d = math.floor((a+b)/2.)
	rand_val = random()
	while a <= b:
		k = ((d-s+1)*states[s] - w[s] + w[d+1])/((n-s+1)*states[s] - w[s])
		if rand_val <= k:
			b = d-1
		else:
			a = d+1
		d = math.floor((a+b)/2.)
	t = b+1

	if (s,t) in pairs_binary:
		pairs_binary[(s,t)] = pairs_binary[(s,t)] + 1./m
	else:
		pairs_binary[(s,t)] = 1./m
	if s in source_binary:
		source_binary[s] = source_binary[s] + 1./m
	else:
		source_binary[s] = 1./m
	if t in target_binary:
		target_binary[t] = target_binary[t] + 1./m
	else:
		target_binary[t] = 1./m
	#print(rand_val)

end = time.time()

print()
print("source_binary\n",source_binary)
print()
print("target_binary\n",target_binary)
print()
print("pairs_binary\n",pairs_binary)
print()
print("done in",end - start)


for (s,t) in pairs_linear:
	max_diff = 0.
	if (s,t) in pairs_binary:
		diff = abs( pairs_linear[(s,t)] - pairs_binary[(s,t)] )
		max_diff = max(max_diff , diff)
print("max diff of pairs sampling is",max_diff)
