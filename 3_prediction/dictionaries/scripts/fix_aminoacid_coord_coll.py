# Coll has coord shift in gyrB for 39 aminoacids

lines = [line for line in open('../Tbdb.txt').readlines()]

new_lines = []

for i in range(len(lines)):
	if 'gyrB' not in lines[i]:
		new_lines.append(lines[i])
	else:
		coord = int(lines[i].split(' ')[1][1:-1])
		coord += 39
		split = lines[i].split(' ')
		new_lines.append(split[0]+ ' ' + split[1][0]+str(coord)+split[1][-1] + ' ' + split[2])

file = open('Tbdb.txt', 'w')

for el in new_lines:
	file.write(el)

file.close()