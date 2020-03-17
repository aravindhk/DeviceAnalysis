"""
Created on Tue March 5 03:14:00 2020

@author: Aravindh Kumar
"""
""" 
Python3 code to rename auto-cascade files in a directory or folder 
to make it easier to parse for analysis
"""
# importing os module 
import os 
#import numpy as np
#import pandas as pd

# Function to rename multiple files 
def main():
    """ 
    Input - IDVG_Die_X5_Y-1_Dev3_Rep1.csv
    Output - 1-5-CH3-G.csv
    """
    current_dir = os.getcwd()
    for filename in os.listdir(current_dir):     
        index1 = filename.find('IDVG_Die_X')
        index2 = filename.find('_Y-')
        index3 = filename.find('_Dev')
        index4 = filename.find('.')
#        print("Index1 ="+str(index1)+" Index2 ="+str(index2)+" Index1 ="+str(index3))        
        if index1 != -1 and index2 != -1 and index3 != -1:
            col_index = int(filename[index1+10])
            row_index = int(filename[index2+3])
            ch_index = int(filename[index3+4])
#            print("Row index ="+str(row_index)+" Col index ="+str(col_index)+" Ch index ="+str(ch_index))
            new_filename = str(row_index)+"-"+str(col_index)+"-CH"+str(ch_index)+"-G"+filename[index4:]
#            print(new_filename)        
            src = current_dir +'\\'+ filename
            dst = current_dir +'\\'+ new_filename    
    		# rename() function will 
    		# rename all the files 
            os.rename(src, dst)
            
# Driver Code                        
if __name__ == '__main__': 
	
	# Calling main() function 
    main() 