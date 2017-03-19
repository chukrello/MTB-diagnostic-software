ins = [line[:-1].split('\t') for line in open('ins.txt').readlines()]

for el in ins:
    ref = el[2]
    alt = el[3]

    delta = len(alt) - len(ref)

    if ref[0] + alt[1:delta+1] + ref[1:] != alt:
        print(el)