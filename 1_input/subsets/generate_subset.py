#farhat

Author_set = [line[:-1] for line in open('Walker_subset.txt').readlines()]
Subset_for_author_test = [line[:-1] for line in open('Full_subset.txt').readlines() if line[:-1] not in Author_set]

file = open('Subset_for_Walker_test.txt', 'w')

for line in Subset_for_author_test:
    file.write(line + '\n')
file.close()