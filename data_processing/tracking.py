table = pd.read_table('/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/max_snp_10/IBSS_not_converge_list.txt',sep = '.', names=["CHR","POS", else]
table = table[1:]
table.index -=1 
new = pd.DataFrame(table.POS.str.split('_').tolist(), columns=["START","MID"])
new["CHR"] = table.CHR.str.replace('chr','') 

new1 = pd.concat([new.START, new.MID, new.CHR], axis=1,keys=["START","END","CHR"]).sort_values(['CHR',"START"])
new2 = pd.concat([new.MID, new.END, new.CHR], axis=1,keys=["START","END","CHR"]).sort_values(['CHR',"START"])

df = pd.concat([new1,new2]) 

df.astype(int)

df = df.sort_values(['CHR',"START"]).reset_index().drop(columns=["index"])
