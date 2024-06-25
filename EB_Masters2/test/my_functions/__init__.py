
def max_prior(genes,prior_info):
    gene1 = genes[0]
    gene2 = genes[1]
    return max(prior_info[int(gene1[1:]),int(gene2[1:])],
               prior_info[int(gene2[1:]),int(gene1[1:])]) *2.2

def search(gene_list,list1,list2):
    found = []
    for genes in gene_list:
        gene1 = str(genes[0])
        gene2 = str(genes[1])
        search = gene1 + '----' + gene2
        if search in list1:
            found.append(1)
        elif search in list2:
            found.append(2)
        else:
            found.append(3)
    
    return found

def index(connection_type,x):
    agreed_index = [i for i, value in enumerate(connection_type) if value == x]
    return agreed_index

def posterior(priors,data,index):
    return list(1 - (priors)*data[0][index]/data[1][index])