ins = [line[:-1].split('\t') for line in open('dels.txt').readlines()]

for el in ins:
    ref = el[2]
    alt = el[3]

    delta = len(ref) - len(alt)

    if alt[0] + ref[1:delta+1] + alt[1:] != ref:
        print(el)