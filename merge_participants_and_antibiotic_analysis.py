import pandas as pd

a =pd.read_csv('insight_export8_antibiotics.csv')

data = np.loadtxt('bin_shift_data')

pats = final_spreadsheet['pat id'].unique()
pats= np.array(pats[~np.isnan(pats)])

pats = pats[pats!=4635].shape

antibio = pd.read_csv('antibiotics.csv')
antibio_pats = np.array(antibio['Participant ID'])
temp = pd.merge(antibio_pats, pats, how = 'inner')

pats = np.sort(pats)
antibio_pats = np.sort(antibio_pats)

pat2 = np.intersect1d(pats,antibio_pats) #there are 75 participants on antibiotics
antibio_pats = antibio[antibio['Participant ID'].isin(pat2)]
antibio_pats.to_csv('antibio_pats.csv')

#get rid of multiple entries manually
pats =pd.read_csv('antibio_pats.csv')
pats.columns = ['number', 'pat id', 'Centre', 'Antibiotic', 'GA given weeks', 
		'GA given days', 'no. of days', 'indication', 'indication other', 'antibiotic_list']

f = final_spreadsheet[['pat id', 'clinical group']]
f = f.drop_duplicates()
f = f[f['pat id']!= 4635]

temp = pd.merge(f,pats, how = 'left', on ='pat id')

temp = temp.dropna(axis =0, thresh =2 )

a = temp['clinical group']
b = [4 if i >4 else i for i in a]
b = np.array(b).astype(int)
temp['clinical group'] = b
temp.groupby('clinical group')['Antibiotic'].value_counts()

t = temp[temp['indication']=='Urinary tract infection']
