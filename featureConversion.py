#package used for data analysis tools
import pandas as pd 


#Calculate Melting Temperature of an input DNA strand according to following formula:
#input dna is of type String
#Formula Source: http://www.biophp.org/minitools/melting_temperature/demo.php?formula=basic
def calcMeltingTemp (dna):
	numA = 0
	numG = 0
	numT = 0
	numC = 0
	for nucleobase in dna:
		if nucleobase == 'A':
			numA += 1
		if nucleobase == 'G':
			numG += 1 
		if nucleobase == 'C':
		 	numC += 1 
		if nucleobase == 'T':
			numT += 1
	
	meltingTemp = 41 * ((numG + numC - 16.4) / (numA + numT + numG + numC)) + 64.9					
	return meltingTemp


#Calculate GC content of an input DNA strand
def gccontent (dna):
	numA = 0
	numG = 0
	numT = 0
	numC = 0
	for nucleobase in dna:
		if nucleobase == 'A':
			numA += 1
		if nucleobase == 'G':
			numG += 1 
		if nucleobase == 'C':
		 	numC += 1 
		if nucleobase == 'T':
			numT += 1
	gcpercent = float(numG + numC) / float(numC + numG + numT + numA)
	return gcpercent		


#Calculate molecular weight of an input DNA strand
#Formula Source: https://www.thermofisher.com/us/en/home/references/ambion-tech-support/rna-tools-and-calculators/dna-and-rna-molecular-weights-and-conversions.html
def molweight (dna):
	numA = 0
	numG = 0
	numT = 0
	numC = 0
	for nucleobase in dna:
		if nucleobase == 'A':
			numA += 1
		if nucleobase == 'G':
			numG += 1 
		if nucleobase == 'C':
		 	numC += 1 
		if nucleobase == 'T':
			numT += 1
	molecularweight = (numA * 313.2) + (numT * 304.2) + (numC * 289.2) + (numG * 329.2) + 79.0
	return molecularweight		


#Calculate extinction coefficient of an input DNA strand
#Formula Source: http://www.atdbio.com/content/1/Ultraviolet-absorbance-of-oligonucleotides
def extinctcoeff (dna):
	numA = 0
	numG = 0
	numT = 0
	numC = 0
	for nucleobase in dna:
		if nucleobase == 'A':
			numA += 1
		if nucleobase == 'G':
			numG += 1 
		if nucleobase == 'C':
		 	numC += 1 
		if nucleobase == 'T':
			numT += 1
	extinctioncoefficient = ((numA * 15.4) + (numC * 7.4) + (numG * 11.5) + (numT * 8.7)) * 0.9 * 1000
	return extinctioncoefficient		


#creates an empty list for each DNA strand characteristic
dnaStrands = list()
meltingTemps = list()
gccontents = list()
molweights = list()
extinctcoeffs = list()		


#Reads .txt file containing DNA strands
dnaFile = open('Bm_space_320_guides.txt')
dnaRows = dnaFile.readlines()

#Goes through each row to analyze each DNA strand 
for line in dnaRows:
	#line = 1 DNA strand

	#add dna strands to dnaStrands list
	dnaStrands.append(line)
	#adds melting temperature of one dna strand to melting temperatures list
	meltingTemps.append(calcMeltingTemp(line))
	#adds gc content of one dna strand to gc contents list
	gccontents.append(gccontent(line))
	#adds molecular weight of one dna strand to molecular weights list
	molweights.append(molweight(line))
	#adds extinction coefficient of one dna strand to extinction coefficients list
	extinctcoeffs.append(extinctcoeff(line))


#Creates dataframe with dictionaries for each DNA characteristic
df = pd.DataFrame({'DNA Strand': dnaStrands,
				  'Melting Temperature (Celsius)': meltingTemps,
				  'GC Content (%)': gccontents,
				  'Molecular Weight (g/mole)': molweights,
				  'Extinction Coefficient (M^-1 cm^-1)': extinctcoeffs})
#Converts dataframe to csv file
df.to_csv('features.csv')




	





