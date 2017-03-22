# Fixes Tbdb.txt file
# There was capital letter for some reason

lines = [line[0].lower() + line[1:-3] + '\n' for line in open('Tbdb_old.txt').readlines()]

file = open('Tbdb.txt', 'w')

for line in lines:
	file.write(line)

file.close()
