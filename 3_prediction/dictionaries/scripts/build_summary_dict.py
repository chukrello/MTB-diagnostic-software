dicts = ['../Walker_dictionary.txt', '../Coll_dictionary.txt', '../Farhat_dictionary.txt', '../Tbdb_dictionary.txt']

sum_dict = []

for dict in dicts:
	sum_dict +=  [line[:-1] for line in open(dict).readlines()]

sum_dict = sorted(list(set(sum_dict)))

file = open('Summary_dictionary.txt', 'w')

for el in sum_dict:
	file.write(el + '\n')

file.close()