import pandas as pd 
import subprocess
import ast
import numpy as np


output_Complete_Pipeline="pytest_STRmie/output_file/Final_report.xlsx"
output_Index_Calculation="pytest_STRmie/output_file/indices_calculation.xlsx"

output_expected_Complete_Pipeline="pytest_STRmie/expected_results/Final_report_expected.xlsx"
output_expected_Index_Calculation="pytest_STRmie/expected_results/indices_calculation_expected.xlsx"


def test_complete_pipeline():	

 subprocess.call(f"strmie --mode Complete_Pipeline -f pytest_STRmie/input_file/ -o pytest_STRmie/output_file/",shell=True)

 df_expected=pd.read_excel(output_expected_Complete_Pipeline)
 df=pd.read_excel(output_Complete_Pipeline)

 #v=(df==df_expected).all()
 #assert all(v)==True
 # confronto robusto con tolleranza sui float
 pd.testing.assert_frame_equal(df, df_expected, check_dtype=False, atol=1e-6, rtol=1e-6)

def test_Index_Calculation():	

 subprocess.call(f"--mode Index_Calculation -f pytest_STRmie/input_file/ -o pytest_STRmie/output_file/ -p pytest_STRmie/input_file/CAG_data_for_recalculating_indices.xlsx",shell=True)

 df_expected=pd.read_excel(output_expected_Index_Calculation)
 df=pd.read_excel(output_Index_Calculation)

 #v=(df==df_expected).all()
 #assert all(v)==True
 # confronto robusto con tolleranza per numeri float
 pd.testing.assert_frame_equal(df, df_expected, check_dtype=False, atol=1e-6, rtol=1e-6)
