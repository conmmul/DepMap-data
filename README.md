**Instructions**


**Create .csv file variables**

* When declaring the file name variables, be sure to edit the first column in the CRISPRGeneEffect.csv to have ‘ModelID’ as the name for the first column and save before passing them into a dataframe.

**Declaring variables to sort the .csv files**
* To change the type of cancer for the first mutation, the variable ‘cancer’ must be declared to match exactly what is written in the model.csv file. This will filter to only the model ids that have that type of cancer.
* To change the gene and the exact protein change in the first mutant, the variable ‘gene’ must match exactly what is written in the OmicsSomaticMutation.csv file under the “HugoSymbol” column. The gene change must match exactly what is written in the OmicsSomaticMutation.csv file under the “ProteinChange” column. 
* For the second gene that is being compared, the same steps should be taken to change the type of cancer as well as the gene that is being compared.

**Co-mutation and gene change boolean expressions**

* In order to be able to compare a mutation in one gene to another co-mutated gene with no specified protein change, you would want to make ‘is_gene_change’ = False and ‘comutation’ = True. For example if I wanted to compare NSCLC KRAS G12C mutants to NSCLC KRAS G12C / SMARCA1 co-mutants, the code would look like this:

``` Ruby
# filter model.csv to only keep indices which contain cancer type declared
cancer = "Non-Small Cell Lung Cancer"
filter_model = model_df[model_df["OncotreePrimaryDisease"] == cancer]
# filter the OmicsSomaticMutations.csv to only include model ids that have declared cancer type
osm_filter = pd.merge(filter_model, osm_df, on=['ModelID'], how='inner')
# filter that dataframe to only include the gene and mutation declared
gene = 'KRAS'
mutation = "p.G12C"
osm_gene_filter = osm_filter[(osm_filter['HugoSymbol'] == gene)]
osm_mutant_filter = osm_gene_filter[(osm_gene_filter['ProteinChange'] == mutation)]

--------------------------------------------------------------------------------------------------------

# filter such that another dataframe is created for all modelIDs that have a mutation in another gene
# if there is a gene change in the second protein 'is_gene_change' should = True
cancer_2 = "Non-Small Cell Lung Cancer"
filter_model_2 = model_df[model_df["OncotreePrimaryDisease"] == cancer_2]
comutation_gene = "SMARCA4"
osm_filter_2 = pd.merge(filter_model_2, osm_df, on=['ModelID'], how='inner')
osm_comutant = osm_filter_2[(osm_filter_2['HugoSymbol'] == comutation_gene)]

is_gene_change = False

if is_gene_change == False:
    pass
else:
    comutation_change = "p.Q61H"
    osm_comutant_filter = osm_comutant[(osm_comutant['ProteinChange'] == comutation_change)]

--------------------------------------------------------------------------------------------------------

# compare the two dataframes. The models that have a mutation for both the gene and comutant gene will be in a comutant dataframe.
# the models that only have the gene mutation will be under no_comutation
# if there is a comutation, comutation = True, if not comutation = False
comutation = True

if comutation == True:
    comutant = osm_mutant_filter.merge(osm_comutant, on=['ModelID'], how='inner')
    no_comutant = osm_mutant_filter.merge(osm_comutant, on=['ModelID'], how='left')

else:
    comutant = osm_comutant_filter
    no_comutant = osm_mutant_filter

```

* In order to be able to compare a mutation in one gene to another non-co-mutated gene with a specific protein change, you would want to make ‘is_gene_change’ = True and ‘comutation’ = False. This can be helpful if you want to compare, for example, NSCLC KRAS G12C mutants to NSCLC KRAS G12V mutants. 


**Making the plot**

* The rest of the code after filtering the dataframes is done automatically, however, the plot will likely need to be edited in order to look clean. Since the -log(p-value) varies drastically over each comparison, the values here may require editing:
``` Ruby
 if volcano.iloc[i].nlogpvalue > 1.25 and (volcano.iloc[i].log2FC > 5 or volcano.iloc[i].log2FC < 5):

```
* In this instance the 1.25 and +/- 5 can be raised or lowered to change the amount of values that have the gene names labeled.

* The plot title can also be changed by altering the string in plt.title (there is a comment indicating where this change should be done) 
	
	
