lines = [line for line in open('../Coll.txt').readlines()]

new_lines = []

for line in lines:
	number = int(line.split(' ')[1][1:-1])

	if number < 4000:
		new_lines.append(line)
	else:
		print(line)

file = open('Coll_new', 'w')

for el in new_lines:
	file.write(el)

file.close()