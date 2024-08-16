#!/usr/bin/env python
# coding: utf-8

# In[1140]:


import os
import pandas as pd
import csv
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.stats


# In[1266]:


# declare CSV variable names
pd.options.display.max_rows = 9999
pd.options.display.max_columns = 9999
model = "T:\\Lab member files\\Connor\\Excel\\DepMap\\Model.csv"
osm = "T:\\Lab member files\\Connor\\Excel\\DepMap\\OmicsSomaticMutations.csv"
effect = "T:\\Lab member files\\Connor\\Excel\\DepMap\\CRISPRGeneEffect.csv"


# In[1671]:


# convert CSVs to pandas dataframes
model_df = pd.read_csv(model, header=0, index_col=0)
osm_df = pd.read_csv(osm, header=0)
effect_df = pd.read_csv(effect, header=0, index_col=0)


# In[1672]:


# filter model.csv to only keep indices which contain cancer type declared
cancer = "Non-Small Cell Lung Cancer"
filter_model = model_df[model_df["OncotreePrimaryDisease"] == cancer]


# In[1673]:


# filter the OmicsSomaticMutations.csv to only include model ids that have declared cancer type
osm_filter = pd.merge(filter_model, osm_df, on=['ModelID'], how='inner')


# In[1674]:


# filter that dataframe to only include the gene and mutation declared
gene = 'KRAS'
mutation = "p.G12C"
osm_gene_filter = osm_filter[(osm_filter['HugoSymbol'] == gene)]
osm_mutant_filter = osm_gene_filter[(osm_gene_filter['ProteinChange'] == mutation)]


# In[1675]:


# filter such that another dataframe is created for all modelIDs that have a mutation in another gene
comutation_gene = "STK11"
osm_filter = pd.merge(filter_model, osm_df, on=['ModelID'], how='inner')
osm_comutant = osm_filter[(osm_filter['HugoSymbol'] == comutation_gene)]


# In[1676]:


# compare the two dataframes. The models that have a mutation for both the gene and comutant gene will be in a comutant dataframe.
# the models that only have the gene mutation will be under no_comutation

comutant = osm_mutant_filter.merge(osm_comutant, on=['ModelID'], how='inner')
no_comutant = osm_mutant_filter.merge(osm_comutant, on=['ModelID'], how='left')
no_comutant_filter = no_comutant.merge(comutant, on = ['ModelID'], how = 'left')


# In[1677]:


# merge the left and inner parts to make a new merged dataframe with duplicates of the comutants that can be dropped
no_comutant = pd.concat([no_comutant_dropped, comutant_filter])
no_comutant = no_comutant.drop_duplicates(keep =  False)
# filter the dataframes just created to only include the model ids of the ones that match the parameters
comutant_filter = comutant.filter(['ModelID'])
no_comutant_dropped = no_comutant.filter(['ModelID'])


# In[1678]:


# using CRISPRGeneEffect.csv, filter into two dataframes for mutant and comutant and average the values over the row
comutant_effect = comutant_filter.merge(effect_df, on=['ModelID'], how='inner')
comutant_effect = comutant_effect.set_index('ModelID')

comutant_effect = comutant_effect.drop_duplicates(keep =  False)
no_comutant_effect = no_comutant.merge(effect_df, on=['ModelID'], how='inner')
no_comutant_effect = no_comutant_effect.set_index('ModelID')
no_comutant_effect = no_comutant_effect.drop_duplicates(keep =  False)


# In[1679]:


# make both dataframes the same shape so the t-test can be run
shape_no_comutant = no_comutant_effect.shape 
newlen = (no_comutant_effect.shape[0] - comutant_effect.shape[0])
drophere = newlen
no_comutant_effect = no_comutant_effect.drop(no_comutant_effect.index[drophere:])


# In[1680]:


# t-test to get p-values
p_val = scipy.stats.ttest_ind(comutant_effect, no_comutant_effect, axis = 1, nan_policy = 'omit')
pvalue = p_val.pvalue


# In[1681]:


# convert pvalues to be used for yaxis
yaxis = -np.log(pvalue)
yaxis = pd.DataFrame(yaxis)
yaxis.rename(columns = {0 : '-logpvalue'}, inplace =True)


# In[1682]:


# convert p values into dataframe
pvalue = pd.DataFrame(pvalue)
pvalue.rename(columns = {0 : 'pvalue'}, inplace= True)
print(pvalue)


# In[1683]:


# find the mean of the dataframes
comutant_mean = comutant_effect.mean(axis = 1, numeric_only = True)
comutant_mean = comutant_mean.reset_index(drop = True)
no_comutant_mean = no_comutant_effect.mean(axis = 1, numeric_only = True)
no_comutant_mean = no_comutant_mean.reset_index(drop = True)


# In[1684]:


# calculate fold change
fold_change = no_comutant_mean.divide(comutant_mean)
print(fold_change)


# In[1685]:


xaxis = np.emath.logn(2, fold_change)
xaxis = pd.DataFrame(xaxis)
xaxis.rename(columns = {0 : 'log2FC'}, inplace= True)


# In[1686]:


# convert fold change series into dataframe
fold_change = fold_change.to_frame(name='fold_change')
print(fold_change)


# In[1687]:


volcano = pd.concat([fold_change, pvalue, xaxis, yaxis], axis =1)
print(volcano)


# In[1688]:


volcano.plot(x = 'log2FC', y = '-logpvalue', kind = 'scatter')


# In[ ]:




