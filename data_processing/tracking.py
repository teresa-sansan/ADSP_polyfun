## This script is a tracking record of the process of fixing convergence issue.
## 1 MB window size with 0.5 MB overlap


table = pd.read_table('/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/max_snp_10/IBSS_not_converge_list.txt',sep = '.', names=["CHR","POS", else]
table = table[1:]
table.index -=1 
new = pd.DataFrame(table.POS.str.split('_').tolist(), columns=["START","END"]) ##extract last 1MB 
new["CHR"] = table.CHR.str.replace('chr','') 
new = new.astype(int)

new['block_1'] = new.START -500000
new['block_2'] = new.START + 500000   ## a.k.a new.end - 500000
new['block_3'] = new.END + 500000  

new1 = pd.concat([new.block_1, new.block_2, new.CHR], axis=1,keys=["START","END","CHR"]).sort_values(['CHR',"START"])
new2 = pd.concat([new.START, new.END, new.CHR], axis=1,keys=["START","END","CHR"]).sort_values(['CHR',"START"])
new3 = pd.concat([new.block_2, new.block_3, new.CHR], axis=1,keys=["START","END","CHR"]).sort_values(['CHR',"START"])
 



df = pd.concat([new1,new2,new3]) 
df = df.astype(int)
df = df[df.START > 0]

df = df.sort_values(['CHR',"START"]).reset_index().drop(columns=["index"])
