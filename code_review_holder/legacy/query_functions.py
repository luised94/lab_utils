# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 00:19:20 2020

@author: liusm
"""
# import os

# file = os.chdir(r"C:\Users\liusm\Desktop\ARS_and_Gene_Analysis\YeastMine")

# from query_functions import query_to_df


# from intermine.webservice import Service 
import pandas as pd
import numpy as np
from datetime import date
#Takes query from intermine library, views_list [] and [] list of columns and outputs and excel file by converting query into dataframe
#and outputting as excel file (.xlsx)
def query_to_df(query, views_list, column_names, To_excel=1, Filename=str(date.today())+'query.xlsx'):
    if len(query.rows()) > 0:
        number_of_views = len(views_list)
        
        data = np.zeros((len(query.rows()),  number_of_views))
        df = pd.DataFrame(data, columns = views_list)
        
        row_count = 0
        
        for row in query.rows():    
        
            column_count = 0
        
            for view in views_list:
                if column_count >= number_of_views+1:
                    break
                else:
                    string_holder = str(row[view]).replace(',', '_' )
                    df.iloc[row_count, column_count] = string_holder.replace('\n', '')
                    column_count += 1
        
            row_count += 1
        df.columns = column_names 
        if To_excel:
            filename = Filename 
            file = open(filename, 'w')
            df.to_excel(file, index=False)
            file.close()
        else:
            return df
    else:
        print("No genes in query")