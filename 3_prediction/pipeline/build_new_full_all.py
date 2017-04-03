full = [line[:-1].split('\t') for line in open('full_all.txt').readlines()]
dict = [line[:-1].split('\t')[0] for line in open('new_dictionary.txt').readlines()]

new_full = []

for element in full:
    new_full.append([element[0]])
    print(element)
    for mut in element[1:]:
        if mut != '' and mut in dict:
            new_full[-1].append(mut)

file = open('new_full_4.txt', 'w')

for line in new_full:
    for element in line[:-1]:
        file.write(element + '\t')
    file.write(line[-1] + '\n')

file.close()