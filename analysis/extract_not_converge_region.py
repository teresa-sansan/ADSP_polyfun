import pandas as pd
df = pd.read_csv('not_converge.txt', sep = '.', names=["chr","start","end"])

def region(type):
	chr = df[df.chr == type]    
	pos =  list(chr.start.append(chr.end).astype(int))
	pos.sort()
	region = []
	for i in range(len(pos)):
		if(pos.count(pos[i])>1):
			region.append('-')
		else:
			region.append(pos[i])

	return(region)

print("chr1, ", region("max_snp_3_chr1"),'\n')
print("chr2, ", region("max_snp_3_chr2"),'\n')   
print("chr3, ", region("max_snp_3_chr3"),'\n')   
