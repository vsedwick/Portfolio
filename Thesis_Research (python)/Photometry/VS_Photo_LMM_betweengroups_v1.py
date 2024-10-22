import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.formula.api import mixedlm
import os
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multicomp import pairwise_tukeyhsd, MultiComparison
from scipy.stats import ttest_ind
import numpy as np
# import rpy2.robjects as ro
# from rpy2.robjects import pandas2ri
# from rpy2.robjects.packages import importr

#per trace
# formula = 'Outcome ~ Fixed_Effect_1 + Fixed_Effect_2 + ... + (1 | Random_Effect)'
file = r"E:\Photometry-Fall2022\Final Analysis\New_Pup_SORT\peri-zscore\Pups\cleaning-IndV_ID_Pups_during_Individual Z-score_summary.csv_compiled.xlsx"
# mouse = {'3': ['3', '23', '33'], '7':['7', '27', '37'], '10': ['10', '210', '310'], '11': ['11', '211', '311']}
mouse = {'3': ['3', '23', '33','43'], '7':['7', '27', '37'], '10': ['10', '210', '310'], '11': ['11', '211', '311'], 
         '1': ['1', '21', '31'], '6': ['6', '26', '36'], '8': ['8', '28', '38'], '16': ['156', '16', '216', '316', '416'],
         '5': ['5', '25', '35', '45', '55'], '2': ['2'], '12': ['12', '212', '312', '412','512'], '4': ['4','24','34','44','54'], 
         '13': ['13','213','313','413'], '82': ['82', '282', '382','482'], '9': ['9', '29', '39','49']}
# mouse = {'1': ['1', '21', '31'], '6': ['6', '26', '36'], '8': ['8', '28', '38']}

def main(file, mouse):
    #TRACE IDs
    traces = {trace_id: mouse_id for mouse_id, trace_ids in mouse.items() for trace_id in trace_ids}
    file_folder = file.split('\\'); file_folder = '\\'.join(file_folder[0:-1])
    full_file = pd.ExcelFile(file)
    _ = make_folder('Statistics', file_folder)
    output = file.split('\\')
    output = output[-1].replace('.xlsx','')    
    save_folder = make_folder(output, _)
    for behavior in full_file.sheet_names:
        print(behavior.upper())
        df = pd.read_excel(file, sheet_name=behavior)
        if df.shape[0] <= 1:
            continue
        else:
            # Transform the data into a long format
            data_long = []
            data_long_noref = []
            data_noref_2pairs = []

            # Iterate over the columns by pairs
            for i in range(0, len(df.columns), 2):
                id_col = df.columns[i]  # ID column for the test
                test_col = df.columns[i + 1]  # Corresponding test measurement column
                
                # Extract test name by removing 'ID_' prefix
                test_name = id_col.replace('ID_', '')

                # Iterate through each row of these columns
                for _, row in df[[id_col, test_col]].dropna().iterrows():
                    trace_id = str(int(row[id_col]))  # Convert to string to match dictionary keys
                    mouse_id = traces.get(trace_id)  # Map trace ID to mouse ID
                    if mouse_id and test_name!='Reference':
                        data_noref_2pairs.append({
                            'Mouse_ID': mouse_id,
                            'Trace_ID': trace_id,
                            'Stimulus': test_name,
                            'mean_zF': row[test_col]
                        })     
                    if mouse_id and test_name != 'Reference':  # Ensure the trace ID was found in the mapping
                        data_long_noref.append({
                            'Mouse_ID': mouse_id,
                            'Trace_ID': trace_id,
                            'Stimulus': test_name,
                            'mean_zF': row[test_col]
                        })                   
                    if mouse_id and 'Reference' in df.columns:  # Ensure the trace ID was found in the mapping
                        data_long.append({
                            'Mouse_ID': mouse_id,
                            'Trace_ID': trace_id,
                            'Stimulus': test_name,
                            'mean_zF': row[test_col]
                        })
            # Transform the data into a long format
            # Convert the list to a DataFrame

            df_long = pd.DataFrame(data_long)
            df_long_noref = pd.DataFrame(data_long_noref)
            data_noref_2pairs = pd.DataFrame(data_noref_2pairs)
            # Corrected model formula: Just the fixed effect (intercept in this case)
            # # Assuming 'Stimulus' is treated as a categorical variable and 'mean_zF' is your dependent variable
            if 'Reference' in df.columns:
                model_formula = 'mean_zF ~ C(Stimulus, Treatment(reference = "Reference"))'
                model = smf.mixedlm(model_formula, df_long, groups=df_long['Mouse_ID'])
                result = model.fit()
                # print(result.pvalues)
                # Save to Excel
                save_name = os.path.join(save_folder, 'LMM_results.xlsx')
                save_two(save_name, behavior, result)

                #MainEffects_LM
                results = nested_anova(df_long_noref)
                save_name = os.path.join(save_folder, 'ANOVA_LM(LMM).xlsx')
                save_results(results , behavior, save_name)

            # Identify all unique stimulus groups
            stimulus_groups = df_long_noref['Stimulus'].unique()
            pairs = []
            for i in stimulus_groups:
                for j in stimulus_groups: 
                    if j != i and (j,i) not in pairs: 
                        pairs.append((i,j))
            
            if len(stimulus_groups)>1:
                if len(data_noref_2pairs['Stimulus'].unique())==2 and len(pairs) == 1:
                    #Nested Unpaired T-test
                    group1 = data_noref_2pairs[data_noref_2pairs['Stimulus'] == data_noref_2pairs['Stimulus'].unique()[0]]['mean_zF']
                    group2 = data_noref_2pairs[data_noref_2pairs['Stimulus'] == data_noref_2pairs['Stimulus'].unique()[1]]['mean_zF']

                    t_statistic, p_value = ttest_ind(group1, group2)
                    tdf = pd.DataFrame({'Group': [str(stimulus_groups),str(stimulus_groups)],
                        'Significance': [p_value, p_value],
                        "T Stat": [t_statistic, t_statistic]})
                    save_name = save_name = os.path.join(save_folder, 'Unpaired nested t-test.xlsx')
                    save_results(tdf, behavior, save_name)
                    
                if len(stimulus_groups)>2:
                    # Convert the Tukey test results to a DataFrame
                    # Perform the pairwise comparison (adjusting the method as necessary for your specific context)
                    multicomp = MultiComparison(df_long_noref['mean_zF'], df_long_noref['Stimulus'])
                    # Adjust this to 'Bonferroni' or another method if you prefer
                    tukey_result = multicomp.tukeyhsd(alpha=0.05)
                    tukey_df2 = tukey_result.summary()
                    tukey_df2 = pd.DataFrame(tukey_df2.data[1:], columns=tukey_df2.data[0])
                    tukey_df2 = format_floats(tukey_df2, 6)
                    save_name = os.path.join(save_folder, 'Nested Tukey_results.xlsx')
                    save_results(tukey_df2, behavior, save_name)
                    if len(stimulus_groups)>3 and 'Pup' not in stimulus_groups:
                        #nested 2way
                        results, tukey = nested_2way(df, traces)
                        # Save the results to an Excel file
                        save_name = os.path.join(save_folder, '2way Nested_Tukey_results.xlsx')
                        save_results(tukey, behavior, save_name)
                        save_name = os.path.join(save_folder, 'Nested 2way ANOVA_results.xlsx')
                        save_two(save_name, behavior, results)
        summarize_by_mice(df, mouse, save_folder, behavior, output) 

print('Statistical summaries will be found in the statistics folder')

def save_results(df, sheet_name, save_name):    
    pd.set_option('display.float_format', lambda x: '%.6f' % x)
    # Convert columns to float as needed and format to 6 decimal places
    for col in df.columns:
        try:
            df[col] = df[col].astype(float).map('{:,.6f}'.format)
        except ValueError:
            # This column cannot be converted to float, likely a non-numeric column, so ignore
            pass
    # Create a new Excel workbook or open an existing one
    if os.path.exists(save_name):
        with pd.ExcelWriter(save_name, mode='a', engine='openpyxl', if_sheet_exists='replace') as writer:
            df.to_excel(writer, sheet_name=f'{sheet_name}', index=True)
    else:
        with pd.ExcelWriter(save_name, mode='w', engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name=f'{sheet_name}', index=True)

def format_floats(df, decimal_places):
    for col in df.columns:
        if df[col].dtype.kind in 'fc':  # f for float, c for complex numbers
            df[col] = df[col].apply(lambda x: f'{x:.{decimal_places}f}')
    return df

def make_folder(x, project_home):

    mode = 0o666
    j=os.listdir(project_home)
    if x not in j:
        behavior_path=os.path.join(project_home, x)
        os.mkdir(behavior_path, mode)
    else:
        behavior_paths=os.path.join(project_home, x)
    behavior_paths=os.path.join(project_home, x)
    return behavior_paths

def nested_anova(df_long):
    pd.set_option('display.float_format', lambda x: '%.6f' % x)
    # Step 2: Fit your models
    # Full model with the recoded 'Stimulus' (excluding reference category)
    full_model = smf.ols('mean_zF ~ C(Stimulus)', data=df_long).fit()

    # Reduced model without 'Stimulus_recode'
    reduced_model = smf.ols('mean_zF ~ 1', data=df_long).fit()

    # Step 3: Compare models with anova_lm
    return anova_lm(reduced_model, full_model)

def nested_2way(df, traces):
    # Initialize an empty DataFrame for the long-format data
    long_df = []
    # Iterate over the columns by pairs
    for i in range(0, len(df.columns), 2):
        id_col = df.columns[i]  # ID column for the test
        test_col = df.columns[i + 1]  # Corresponding test measurement column
        
        # Extract test name by removing 'ID_' prefix
        test_name = id_col.replace('ID_', '')
        if test_name == 'V. Males Attacking'or test_name == 'Reference':
            continue
        else:
            # Iterate through each row of these columns
            for _, row in df[[id_col, test_col]].dropna().iterrows():
                trace_id = str(int(row[id_col]))  # Convert to string to match dictionary keys
                mouse_id = traces.get(trace_id)  # Map trace ID to mouse ID
                groupings = test_name.split('.')
                condition = groupings[0]
                if ' male' in groupings[1].lower():
                    sex = 'Male'
                elif 'female' in groupings[1].lower():
                    sex = 'Female'
                if mouse_id and test_name!='Reference':
                    long_df.append(({
                        'Mouse_ID': mouse_id,
                        'Trace_ID': trace_id,
                        'Condition': condition,
                        'Sex': sex,
                        'mean_zF': row[test_col]
                    }))
    long_df = pd.DataFrame(long_df)

    model = mixedlm("mean_zF ~ Sex * Condition", data=long_df, groups=long_df["Mouse_ID"], re_formula="~Condition")
    result = model.fit()

    # Prepare data for Tukey HSD test - you might need to adjust this part
    # Here, I'm assuming you want to test the interaction between 'Sex' and 'Condition'
    long_df['Interaction'] = long_df['Sex'] + "_" + long_df['Condition']


    # Prepare for Tukey HSD test
    multicomp = MultiComparison(long_df['mean_zF'], long_df['Condition'] + "_" + long_df['Sex'])
    tukey_result = multicomp.tukeyhsd(alpha=0.05)
    tukey_df2 = pd.DataFrame(tukey_result.summary().data[1:], columns=tukey_result.summary().data[0])
    
    # Assuming `format_floats` is a function you've defined to format the float columns
    tukey_df2 = format_floats(tukey_df2, 6)  # You'll need to define this function
    
    return result, tukey_df2
    
def save_two(save_name, behavior, result):
    
    # Convert summary tables to DataFrames
    pd.set_option('display.float_format', lambda x: '%.6f' % x)
    df1 = result.summary().tables[0]
    if not isinstance(df1, pd.DataFrame):
        df1 = pd.DataFrame(df1)  # Convert to DataFrame if not already

    pval_save = result.pvalues.map(lambda x: f"{x:.6g}").to_frame('p-values')  # Formatting p-values
    # print(behavior, pval_save)
    df2 = result.summary().tables[1]
    if not isinstance(df2, pd.DataFrame):
        df2 = pd.DataFrame(df2)  # Convert to DataFrame if not already

   

    if os.path.exists(save_name):
        with pd.ExcelWriter(save_name, mode='a', engine='openpyxl', if_sheet_exists='replace') as writer:
            # Get the maximum row to avoid overwriting data in case of existing sheet
            if behavior in writer.book.sheetnames:
                startrow = writer.book[behavior].max_row
            else:
                startrow = 0

            # Write each DataFrame to the same sheet at different start rows
            df1.to_excel(writer, sheet_name=behavior, startrow=startrow, index=False)
            # pval_save.to_excel(writer, sheet_name=behavior, startrow=startrow + len(df1) + 2, index=False)
            pval_save.to_excel(writer, sheet_name=behavior, startrow=startrow + len(df1) + len(pval_save) + 4)
    else:
        with pd.ExcelWriter(save_name, mode='w', engine='openpyxl') as writer:
            df1.to_excel(writer, sheet_name = behavior, startrow = 0, index = False)
            # pval_save.to_excel(writer, sheet_name = behavior, startrow = 10, index = False)
            pval_save.to_excel(writer, sheet_name = behavior, startrow = 20)

def summarize_by_mice(df, mouse, save_folder, sheet_name, output):
    new_df = {}
    for i in range(0, len(df.columns), 2):
        id_col = df.columns[i]
        if 'Reference' in id_col:
            continue
        else:
            id_id = [str(int(i)) for i in df[id_col].dropna()]; unique_id = [*set(id_id)]
            test_col = df.columns[i + 1]
            test_name = id_col.replace('ID_', '')
            avg_values = []
            keys = []
            for key in mouse:
                new_idcol = [i in mouse[key] for i in id_id]
                test_col = df[test_name].dropna()
                new_testcol = test_col[new_idcol]
                if len(new_testcol) < 2:
                    continue
                else:
                    avg_values.append(np.average(new_testcol))
                    keys.append(key)

            new_df.update({
                f'ID_{test_name}': keys,
                f'{test_name}': avg_values
            })
    
    new_df = pd.DataFrame(padding(new_df))

    save_name = os.path.join(save_folder, f'BYMICE.xlsx')
    if os.path.exists(save_name):
        with pd.ExcelWriter(save_name, mode='a', engine='openpyxl', if_sheet_exists='replace') as writer:
            new_df.to_excel(writer, sheet_name=sheet_name, startrow=0)
    else:
        with pd.ExcelWriter(save_name, mode='w', engine='openpyxl') as writer:
            new_df.to_excel(writer, sheet_name=sheet_name, startrow=0, index=False)
def padding(data):
    # Find the maximum length of the lists in the dictionary
    max_length = max(len(lst) for lst in data.values())

    # Fill each list in the dictionary with NaNs to ensure they all have the same length
    for key in data.keys():
        data[key] += [np.nan] * (max_length - len(data[key]))

    return data

if __name__ == "__main__":  # avoids running main if this file is imported
    main(file, mouse)