import pandas as pd
nom=pd.read_csv('csv.csv')
up=pd.read_csv('up_csv.csv')
dn=pd.read_csv('down_csv.csv')
output_df=up
output_df[' formula']='2*('+up[' formula']+'-'+nom[' formula']+')+'+nom[' formula']
output_df[' jetFlavor']='1'
output_df.to_csv('outup_csv.csv',index=False)
# print('2*'+(nom[' formula']+'-'+dn[' formula'])+'+'+nom[' formula'])
# print(type(nom[' formula']))
# print(nom[' formula'])
# print(up[' formula'])
