import pandas as pd
import numpy as np
import sklearn.decomposition as decomp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re
#working directory
#/home/louic/Desktop/vaginal_samples/data_summer2016_vikash_data/vag_project/
#load and split the nmr title data file
pat_metadata = []
all_dat_dict = []

tit = pd.read_csv('titles.csv', sep ='\t')
clinical_groups = pd.read_csv('unblinded.csv')
dataset = pd.read_csv('Insight_metabolome_nmr -PSRT 11.11.16.csv')

cg = dataset[['Pt. ID', 'Clinical group', 'HVS cat', 'Pulled barcode']]

folders =['160627_TK_vagswab','160628_TK_vagswab', '160728_TK_vagswab_A','160629_TK_vagswab',
	  '160729_TK_vagswab','160730_TK_vagswab', '160731_TK_vagswab', '160801_TK_vagswab',
	  '160728_TK_vagswab_B', '160802_TK_vagswab_rerun']

#dictionary comprehension to split the nmr title list into different sets that can be 
#worked with. Make it a dataframe

nmr_data_dict = { i: tit[tit['folder']==i] for i in folders}  
for i in folders:
  nmr_data_dict[i] = pd.DataFrame(nmr_data_dict[i])

  #load and split the patient metadatafile
pat_metadata = []

x   = pd.read_csv('Vaginal samples blinded for NMR-Day 1-27-06-16.csv')
x2  = pd.read_csv('Vaginal samples blinded for NMR-Day 2-28-06-16.csv')
x3  = pd.read_csv('Vaginal samples blinded for NMR-Day 3-28-07-16.csv')
x4  = pd.read_csv('Vaginal samples blinded for NMR-Day 3-29-06-16.csv')
x5  = pd.read_csv('Vaginal samples blinded for NMR-Day 4-29-07-16.csv')
x6  = pd.read_csv('Vaginal samples blinded for NMR-Day 5-30-07-16.csv')
x7  = pd.read_csv('Vaginal samples blinded for NMR-Day 6-31-07-16.csv')
x8  = pd.read_csv('Vaginal samples blinded for NMR-Day 7-01-08-16.csv')
x9  = pd.read_csv('Vaginal samples blinded for NMR-Day 8-28-07-16.csv')
x10 = pd.read_csv('Vaginal samples blinded for NMR-Rerun-02-08-16.csv')
file_list= [x,x2,x3,x4,x5,x6,x7,x8,x9,x10]
##merge the row and well number text in the loaded files as they've been split
##this will facilitate merging the nmr metadata with the patient metadata
##first have to convert numbers to string can then concatenate
##making dictionary with all the patient metadata divided by date
pat_metadata= { folders[i]: file_list[i] for i in range(len(folders))}

for i in folders:
  #i=0
  #numpy doesnt allow NaN in interger arrays. so to work around, I convert them to float
  #then string, then chop off the ends of the string to make it look like nmr_data_dict['well']
  pat_metadata[i]['Destination 96-well plate Column'] = pd.to_numeric(pat_metadata[i]['Destination 96-well plate Column'] , errors = 'coerce')
  pat_metadata[i]['Destination 96-well plate Column'] = pat_metadata[i]['Destination 96-well plate Column'].astype(float)
  pat_metadata[i]['Destination 96-well plate Column'] = pat_metadata[i]['Destination 96-well plate Column'].astype(str)  
  pat_metadata[i]['Destination 96-well plate Column'] = pat_metadata[i]['Destination 96-well plate Column'].map(lambda x: x[:-2])
  pat_metadata[i]['Destination 96-well plate Column'] = pat_metadata[i]['Destination 96-well plate Column'].apply(str)
  pat_metadata[i]['well'] = pd.Series(
			      pat_metadata[i]['Destination 96-well plate Row'].str.cat(
			      pat_metadata[i]['Destination 96-well plate Column']))
  #make it a dataframe
  pat_metadata[i] = pd.DataFrame(pat_metadata[i])


#now merge individual metadata sets, do outer merge to not miss out on anything,
#can later elminate rows that lack nmr metadata
final_df = {}
labels = pat_metadata.keys()
for i in labels:
  final_df[i] = pd.merge(nmr_data_dict[i], pat_metadata[i], how= 'outer', on = 'well')
final_df  = pd.concat(final_df)  

nmr_final_df = final_df[final_df['nmr_identifier'].isnull() == False]
print final_df.shape
print nmr_final_df.shape

nmr_list = ['HVS sample barcode', 'NMR File name', 'nmr_identifier', 'well']
df = nmr_final_df[nmr_list]
df.columns = ['barcode', 'nmr file name', 'nmr identifier', 'well']
cg.columns = ['pat id', 'clinical group', 'hvs cat', 'barcode']
#sort
df = df.sort_values(by  =  'barcode')
cg = cg.sort_values(by = 'barcode')
a = pd.merge(cg,df,how = 'outer', on = 'barcode')
a = a[a['nmr identifier'].isnull() == False]
a = a.sort_values(by = 'nmr identifier')
#all the ones with wells == 10 are cut off with data, have to figure out what happened
#turns out the regex function was cutting of the 0 from 10, so any well number 10 was deleted
#tried again with pd.applymap for elementwise operation removing the end of the string [:-2]
#applymap only works on dataframes. have to use map to series
#vstack all the dictionaries to make one final spreadsheet

#end up with dataframe of size 717, while we have size 715 for the nmr data
#inspection shows there are two duplicates at 593 and 594
final_spreadsheet = a.drop_duplicates()
final_spreadsheet = final_spreadsheet.reset_index(drop = True)
final_spreadsheet.to_csv('final_spreadsheet.csv')


#now to remove data that isn't going to be used in the analysis due to errors in nmr runs

def plot_spects(ppm,data,cx,cy):
    countx = 0
    county = 0
    for i in data:
        plt.plot(ppm+countx, i+county)
        countx+=cx
        county+=cy
    plt.gca().invert_xaxis()
    plt.show()
    

def scaletointegral(ydata):
        integral = ydata.sum(1) / 100
        newydata = np.transpose(ydata.T / integral)
        return np.array(newydata)
        
def scalepqn(ydata):
      # 1: scale to integral
      newydata = scaletointegral(ydata)
      # 2: reference spectrum is median of all samples
      yref = np.median(ydata, 0)
      # 3: quotient of test spectra with ref spectra
      yquot = newydata / yref
      # 4: median of those quotients
      ymed = np.median(yquot, 0)
      # 5: divide by this median
      newydata = newydata / ymed
      return np.array(newydata)
      
def do_pca(data):
  pca = decomp.PCA(n_components  = 2)
  pca.fit(data)
  pca_score = pca.explained_variance_ratio_
  V = pca.components_
  X = pca.transform(data)
  
  return X, V, pca_score
  
#do a PCA to identify outliers (bad phases, shim etc..)  

data = np.loadtxt('bin_shift_data')
data_to_remove = list(np.loadtxt('to_remove_before_analysis.csv'))
metadata = final_spreadsheet.drop(final_spreadsheet.index[data_to_remove])
metadata = metadata.reset_index( drop = True)
ppm =data[0]
data = data[1:]
data = np.delete(data, data_to_remove, 0)

#remove regions of EDTA
regions_to_remove =np.array( [ [1.17, 1.20], [1.14, 1.16], [2.70, 2.72], [2.55, 2.59], [3.05, 3.26],
		      [3.65, 3.90], [3.17, 3.39], [7.65, 7.70], [8.08, 8.12]])
n,m = regions_to_remove.shape
regions_to_remove = regions_to_remove.flatten()
idx_to_remove = [np.argmin(np.abs(ppm-i)) for i in regions_to_remove]
idx_to_remove = np.reshape(idx_to_remove,[n,m])
idx_to_remove.sort()
for i in idx_to_remove:
  data = np.delete(data,range(i[0], i[1]), axis  =1)
  ppm = np.delete(ppm, range(i[0], i[1]), axis =0)

       
#normalise and scale data
data = scalepqn(data)
data = data-data.mean(axis=0)
data = np.nan_to_num(data)


#PCA of run order within groups
test = np.unique(metadata['nmr file name'])
for i in test:
  X=[]
  nmr_data = data[np.where(metadata['nmr file name']==i)]
  X = do_pca(nmr_data)
  explained_variance = X[2]*100
  X  = X[0]
  n = len(X)
  m = range(n)
  X = pd.DataFrame((X))
  X.columns = ['PC1','PC2']
  X['color'] = m
  colors_to_use = sns.cubehelix_palette(n, light = 0.8)
  #sns.palplot(colors_to_use);plt.show()

  ax = sns.lmplot('PC1', 'PC2', data =X, fit_reg = False, hue = 'color', palette = colors_to_use ,scatter_kws = {'s':150, 'edgecolor' :'black'}, legend = False,
		  size = 7, aspect = 1.5
		  )
  ax.set(title = i,
	xlabel = "PC1 (%.2f%% explained variance)" %(explained_variance[0]), 
	ylabel = "PC2 (%.2f%% explained variance)" %(explained_variance[1])
	)
  plt.tight_layout()
  plt.show()
#plt.scatter(X[:,0], X[:,1])
#for i, txt in enumerate(labels):
  #plt.annotate(txt, (X[i,0], X[i,1]))
#plt.show()

#plot pca according to run day
X = do_pca(data)
explained_variance = X[2]*100
loadings = X[1]
X = X[0]
pcadat = pd.DataFrame((X))
pcadat.columns =['PC1', 'PC2']
pcadat['nmr file name'] = metadata['nmr file name']
#loadings plot
loadings = pd.DataFrame(loadings).T
loadings['ppm'] = ppm
loadings.columns =['PC1', 'PC2', 'ppm']

sns.set(font_scale = 1.5)
colors = sns.color_palette("Paired", 10)

ax = sns.lmplot('PC1', 'PC2', data = pcadat, hue = 'nmr file name', fit_reg = False, palette =colors,
	    scatter_kws ={"s":75} )
#for i, txt in enumerate(pcadat['names']):
 # plt.annotate(i, (pcadat['PC1'].ix[i], pcadat['PC2'].ix[i] ))   
ax.set(title = "PCA Scores Plot",
       xlabel = "PC1 (%.2f%% explained variance)" %(explained_variance[0]), 
       ylabel = "PC2 (%.2f%% explained variance)" %(explained_variance[1])
       )
plt.show()


ax = sns.lmplot('PC1', 'PC2', data = loadings, fit_reg = False, palette = colors, scatter_kws = {'s':50})
for i, txt in enumerate(loadings['ppm']):
  
  plt.annotate("%.3f"%txt, (loadings['PC1'].ix[i], loadings['PC2'].ix[i]))
ax.set(title = "PCA Loadings",
       xlabel = "PC1 (%.2f%% explained variance)" %(explained_variance[0]), 
       ylabel = "PC2 (%.2f%% explained variance)" %(explained_variance[1])
       )
plt.show()
  
#merge the pcadat data with the clinical data and make pca plots

pcadat = pd.concat([pcadat,metadata],axis = 1)
pcadat = test.drop('index', axis =1)
def pca_scores_and_density(data, x ='PC1', y = 'PC2', hue = 'group'):
  g = sns.PairGrid(data, vars = [x,y], hue = hue)
  g = g.map_diag(sns.kdeplot)
  g = g.map_offdiag(plt.scatter)
  g = g.add_legend()
  plt.show()
g = sns.PairGrid(pcadat, vars = ['PC1', 'PC2'], hue = 'clinical group')
g = g.map_diag(sns.kdeplot)
g = g.map_offdiag(plt.scatter)
g = g.add_legend()
plt.show()
g = sns.PairGrid(pcadat, vars = ['PC1', 'PC2'], hue = 'hvs cat')
g = g.map_diag(sns.kdeplot)
g = g.map_offdiag(plt.scatter)
g = g.add_legend()
plt.show()

def pcascores(data, explained_variance, group):
  sns.set(font_scale = 1.5)
  colors = sns.color_palette('dark')
  ax = sns.lmplot('PC1', 'PC2', data = data, hue = group, fit_reg = False, palette =colors,
	      scatter_kws ={"s":75, "alpha":0.4} )
  d2 = data.groupby(group)
  d2 = d2.mean()
  ax.axes[0][0].scatter(x = 'PC1', y = 'PC2', data = d2, s =250, color = 'grey')
  for i in d2.index:
    ax.axes[0][0].annotate("%i"%i, xy = (d2.ix[i]['PC1'], d2.ix[i]['PC2']))
  #for i, txt in enumerate(pcadat['names']):
  # plt.annotate(i, (pcadat['PC1'].ix[i], pcadat['PC2'].ix[i] ))   
  ax.set(title = "PCA Scores Plot",
	xlabel = "PC1 (%.2f%% explained variance)" %(explained_variance[0]), 
	ylabel = "PC2 (%.2f%% explained variance)" %(explained_variance[1])
	)
	
	
  plt.show()
		      #remove outliers
#outliers =pd.read_csv('outliers.csv', header = None)
#outliers = np.array(outliers)
#outliers = np.reshape(outliers, len(outliers))
#outlier_data = data[outliers]
#norm_data = np.delete(data,outliers,axis = 0)
#for i in outlier_data:
  #plt.plot(ppm, i, 'r')
#for i in norm_data:
  #plt.plot(ppm, i, 'b')
#plt.gca().invert_xaxis()
#plt.show()


#data2 =np.delete(data, to_remove,0)
#file_names2 = np.delete(file_names, to_remove, 0)
#g2 = np.delete(g,to_remove,0)
#X2 = do_pca(data2)

#loadings = X2[1]
#explained_variance = X2[2]*100
#X2 = X2[0]
#pat = mpl.patches.Patch (color=g2 , label = labels)
#pat = mpl.cm.hot(g2)
#pat2 = []
#for i,j in enumerate(pat):
    #pat2.append(mpl.patches.Patch(color = j, label=  file_names2[i]))

    
#pcadat = pd.DataFrame((X2))
#pcadat.columns =['PC1', 'PC2']
#pcadat['colors'] =g2
#pcadat['names'] = file_names2
##loadings plot
#loadings = pd.DataFrame(loadings).T
#loadings['ppm'] = ppm
#loadings.columns =['PC1', 'PC2', 'ppm']

#sns.set(font_scale = 1.5)
#colors = sns.color_palette("Paired", 10)

#ax = sns.lmplot('PC1', 'PC2', data = pcadat, hue = 'names', fit_reg = False, palette =colors,
	    #scatter_kws ={"s":50} )
#for i, txt in enumerate(pcadat['names']):
  #plt.annotate(i, (pcadat['PC1'].ix[i], pcadat['PC2'].ix[i] ))   
#ax.set(title = "PCA Scores Plot",
       #xlabel = "PC1 (%.2f%% explained variance)" %(explained_variance[0]), 
       #ylabel = "PC2 (%.2f%% explained variance)" %(explained_variance[1])
       #)
#plt.show()

#ax = sns.lmplot('PC1', 'PC2', data = loadings, fit_reg = False, palette = colors, scatter_kws = {'s':50})
#for i, txt in enumerate(loadings['ppm']):
  
  #plt.annotate("%.3f"%txt, (loadings['PC1'].ix[i], loadings['PC2'].ix[i]))
#ax.set(title = "PCA Loadings",
       #xlabel = "PC1 (%.2f%% explained variance)" %(explained_variance[0]), 
       #ylabel = "PC2 (%.2f%% explained variance)" %(explained_variance[1])
       #)
#plt.show()

count = 0
for i in y2:
  if i==0:
    count = i
  elif i>0:
    
    if i-count ==1:
      count =i#
      next

    else:
      print "%i - %i = %i" %(i,count, i-count)
      count = i
