import pandas as pd

#path = '/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/polypred/'
path = '/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/polypred/PLINKupdate_max'

prs1 = pd.read_csv(path+'1.prs', sep = '\t' )
prs3 = pd.read_csv(path+'3.prs', sep = '\t') 
prs5 = pd.read_csv(path+'5.prs', sep = '\t')
prs7 = pd.read_csv(path+'7.prs', sep = '\t')
prs10 =pd.read_csv(path+'10.prs', sep = '\t')

prs = pd.DataFrame({'PRS1':prs1.PRS,'PRS3':prs3.PRS, 'PRS5':prs5.PRS,'PRS7':prs7.PRS, 'PRS10':prs10.PRS})
prs['IID'] = prs1.IID
merged = prs1.drop(columns = 'PRS').merge(prs)

merged.to_csv('/gpfs/commons/home/tlin/output/prs/updated_plink_prs.tsv', sep = '\t', index = False)