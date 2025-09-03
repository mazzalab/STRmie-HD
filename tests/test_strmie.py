import pandas as pd 
import subprocess
import ast
import numpy as np


output_Complete_Pipeline="tests/output_file/Final_report.xlsx"
output_Index_Calculation="tests/output_file/indices_calculation.xlsx"

output_expected_Complete_Pipeline="/tests/expected_results/Final_report_expected.xlsx"
output_expected_Index_Calculation="/tests/expected_results/indices_calculation_expected.xlsx"


def test_complete_pipeline():	

 subprocess.call(f"strmie --mode Complete_Pipeline -f tests/input_file/ -o tests/output_file/",shell=True)

 df_expected=pd.read_excel(output_expected_Complete_Pipeline)
 df=pd.read_excel(output_Complete_Pipeline)

 # sort by sample 
 df = df.sort_values(by="Sample").reset_index(drop=True)
 df_expected = df_expected.sort_values(by="Sample").reset_index(drop=True)


 #v=(df==df_expected).all()
 #assert all(v)==True
 # confronto robusto con tolleranza sui float
 pd.testing.assert_frame_equal(df, df_expected, check_dtype=False, atol=1e-6, rtol=1e-6)

def test_Index_Calculation():	

 subprocess.call(f"--mode Index_Calculation -f tests/input_file/ -o tests/output_file/ -p tests/input_file/CAG_data_for_recalculating_indices.xlsx",shell=True)

 df_expected=pd.read_excel(output_expected_Index_Calculation)
 df=pd.read_excel(output_Index_Calculation)
 
  # sort by sample 
 df = df.sort_values(by="Sample").reset_index(drop=True)
 df_expected = df_expected.sort_values(by="Sample").reset_index(drop=True)

 #v=(df==df_expected).all()
 #assert all(v)==True
 # confronto robusto con tolleranza per numeri float
 pd.testing.assert_frame_equal(df, df_expected, check_dtype=False, atol=1e-6, rtol=1e-6)
