# Script shows differences of dictionaries

dictionary1_name = '../Coll.txt'
dictionary2_name = '../Walker.txt'
drug_name = 'INH'

dict1 = [line for line in open(dictionary1_name).readlines() if drug_name in line]
dict2 = [line for line in open(dictionary2_name).readlines() if drug_name in line]

differences_1 = []
differences_2 = []

for el in dict1:
	if el not in dict2:
		differences_1.append(el)

for el in dict2:
	if el not in dict1:
		differences_2.append(el)

#print('Difference first dict')
#print(differences_1)

#print('Difference second dict')
#print(differences_2)

for el in differences_2:
		print(el[:-1])