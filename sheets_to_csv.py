import pandas as pd
sheet = 0
max_sheet = 11
while sheet<=max_sheet-1:
	nik = pd.read_excel('RMSE_results.xlsx', sheetname = sheet)
	out_file = 'sheet'+str(sheet)+'.csv'
	nik.to_csv(out_file,header=None,index=True)
	sheet = sheet+1

