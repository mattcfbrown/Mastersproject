#Testing something quickly

num_genes = 25
for i in range(num_genes):
    for j in range(i,num_genes):
        p = 1

y = ['a','b','c','d','e']
z = ['a','f','b','d','q']

agreed_no = set(y).intersection(set(z))
agreed_no = list(agreed_no)
sim_found_no = [x for x in y if x not in agreed_no]
orig_found_no = [x for x in z if x not in agreed_no]

print(sim_found_no)
print(orig_found_no)


