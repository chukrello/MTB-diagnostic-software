# Fixes Pankhurst dictionary
# Delete \r

lines = [line for line in open('Pankhurst_old.txt').readlines()]
new_lines = []

for i in range(len(lines)):
	tmp = ''
	for letter in lines[i]:
		if letter != '\r':
			tmp += letter

	new_lines.append(tmp)

file = open('Pankhurst.txt', 'w')

for line in new_lines:
	file.write(line)

file.close()
