string = 'MGKNEARRSALAPDHGTVVCDPLRRLNRMHATPEESIRIVAAQKKKAQDEYGAASITILEGLEAVRKRPGMYIGSTGERGLHHLIWEVVDNAVDEAMAGYATTVNVVLLEDGGVEVADDGRGIPVATHASGIPTVDVVMTQLHAGGKFDSDAYAISGGLHGVGVSVVNALSTRLEVEIKRDGYEWSQVYEKSEPLGLKQGAPTKKTGSTVRFWADPAVFETTEYDFETVARRLQEMAFLNKGLTINLTDERVTQDEVVDEVVSDVAEAPKSASERAAESTAPHKVKSRTFHYPGGLVDFVKHINRTKNAIHSSIVDFSGKGTGHEVEIAMQWNAGYSESVHTFANTINTHEGGTHEEGFRSALTSVVNKYAKDRKLLKDKDPNLTGDDIREGLAAVISVKVSEPQFEGQTKTKLGNTEVKSFVQKVCNEQLTHWFEANPTDAKVVVNKAVSSAQARIAARKARELVRRKSATDIGGLPGKLADCRSTDPRKSELYVVEGDSAGGSAKSGRDSMFQAILPLRGKIINVEKARIDRVLKNTEVQAIITALGTGIHDEFDIGKLRYHKIVLMADADVDGQHISTLLLTLLFRFMRPLIENGHVFLAQPPLYKLKWQRSDPEFAYSDRERDGLLEAGLKAGKKINKEDGIQRYKGLGEMDAKELWETTMDPSVRVLRQVTLDDAAAADELFSILMGEDVDARRSFITRNAKDVRFLDV*'

for i in range(len(string)):
	print(i)
	if string[i] == 'G':
		if string[i+24] == 'D':
			if string[i+29] == 'N':
				if string[i+31] == 'E':
					if string[i+34] == 'D':
						if string[i+68] == 'Q':
							print(i)
