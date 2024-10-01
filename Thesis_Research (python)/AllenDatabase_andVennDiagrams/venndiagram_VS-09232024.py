"""
This script was written by Victoria Sedwick as a student in Anita Aury's lab at Albert Einstein College of Medicine. 
The script is meant to be used in combination with Allen Brain Atlas single cell sequencing database. If you have any questions, contact sedwick.victoria@gmail.com

https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas
"""


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os
import yaml
import warnings
import sys

variables_config = r"config.yaml"

def main(variables_config):
    """
    This script is based on gene data from a dataset output that can be obtained from "[Jupyter notebook file]".

    In main, the genes are filtered by their classifications (e.g. 'class', supertype', sex')
    Each output will be placed in separate folders based on the information that is being filtered.
    """
    global gene1, gene2, gene3, combo
    warnings.filterwarnings("ignore")

        # Example usage
    config = load_config(variables_config)

    dataset = config['dataset']
    image_save_folder = config['image_save_folder']
    file_path = config['file_path']
    genes = config['genes']
    gene1 = genes['gene1']
    gene2 = genes['gene2']
    gene3 = genes['gene3']
    

    combo = f'{gene1}-{gene2}-{gene3}'


    region = config['venn_settings']['region']

    colors = config['venn_settings']['colors']
    save_formats = config['venn_settings']['save_formats']

        #create save_path
    _ = make_folder(dataset, image_save_folder)
    save_path = make_folder(combo, _)

    #Read file path
    expression_df = pd.read_csv(file_path)

    #using subclass to filter for MEA region and lhx6 transcription factor
    identify_region = expression_df['subclass'].str.contains('MEA' or 'Lhx6')

    #new dataframe only containing MEA and lhx6
    filtered_df = expression_df[identify_region]

    #split sexes
    female_df = filtered_df[filtered_df['donor_sex'].str.contains('F')]
    male_df = filtered_df[filtered_df['donor_sex'].str.contains('M')]

    #VENN DIAGRAMs for both sexes
    subclass_list, class_list, supertype = make_lists(filtered_df); subclass_folder = make_folder('Subclasses', save_path)
    class_folder = make_folder('Classes', save_path); supertype_folder = make_folder('Supertypes', save_path)
    total_folder = make_folder('All Cells', save_path); coloc_folder = make_folder('Colocalizations', total_folder)
    
    #VENN DIAGRAMs for Females
    subclass_list, class_list, supertype = make_lists(female_df)
    f_subclass_folder = make_folder('Female', subclass_folder); f_class_folder = make_folder('Female', class_folder)
    f_supertype_folder = make_folder('Female', supertype_folder); 

    #VENN DIAGRAMs for Males
    subclass_list, class_list, supertype = make_lists(male_df)
    m_subclass_folder = make_folder('Male', subclass_folder); m_class_folder = make_folder('Male', class_folder)
    m_supertype_folder = make_folder('Male', supertype_folder); 

    # Repeat for females and males
    process_groups(filtered_df, 'ALL', subclass_folder, class_folder, supertype_folder, total_folder, coloc_folder, region, colors, save_formats, save_path)
    process_groups(female_df, 'F', f_subclass_folder, f_class_folder, f_supertype_folder, total_folder, coloc_folder, region, colors, save_formats, save_path)
    process_groups(male_df, 'M', m_subclass_folder, m_class_folder, m_supertype_folder, total_folder, coloc_folder, region, colors, save_formats, save_path)


def make_lists(df):
    """ Extracts the names of unique cell clusters in each category
    """
    return df['subclass'].unique(), df['class'].unique(), df['supertype'].unique()

def process_groups(df, sex, subclass_folder, class_folder, supertype_folder, total_folder, coloc_folder, region, colors, save_formats, save_path):

    """
    Processes venndiagrams for all categories
    """

    print(f"Processing {sex}_{region}")
    subclass_list, class_list, supertype = make_lists(df);
    venndiagrams(df, subclass_list, 'subclass', fr"{subclass_folder}\{sex}_Subclass", save_formats)
    venndiagrams(df, class_list, 'class', fr"{class_folder}\{sex}", save_formats)
    venndiagrams(df, supertype, 'supertype', fr"{supertype_folder}\{sex}_Supertype", save_formats)
    total_region(df, f'{region}_{sex}', fr"{total_folder}\{sex}", save_formats)
    colocalizations(df, f'{region}_{sex}', fr"{coloc_folder}\{sex}", colors, save_formats)


def venndiagrams(df, item, name, file_save, save_formats):

    """
    This function determines the proportions for the venn diagrams vased on the  values in the dataset
    The final output will be the venn diagram image as both a vector and png.
    """
    for i, class_name in enumerate(item):
        plt.figure()
        # Subset data for each class
        subset_data = df[df[name] == class_name]
        
        # Determine sets for genes based on non-zero expression
        set_gene1 = set(subset_data[subset_data[gene1] > 0]['cell_barcode'])
        set_gene2 = set(subset_data[subset_data[gene2] > 0]['cell_barcode'])
        set_gene3 = set(subset_data[subset_data[gene3] > 0]['cell_barcode'])
        # set_slc17a6 = set(subset_data[subset_data['Slc17a6'] > 0]['cell_barcode'])
        
        # Plot Venn diagram
        venn3([set_gene1, set_gene2, set_gene3], (gene1, gene2, gene3))
        plt.title(f'{class_name}')

        plt.tight_layout()
        # plt.show()
        # Save in multiple formats if specified
        for fmt in save_formats:
            plt.savefig(f'{file_save}.{fmt}')
        plt.close()

def total_region(df, region, file_save, save_formats):
    """This function calculates the total number of cells that express the three defines genes with
    proportionate circles. """

    # Determine sets for genes based on non-zero expression
    set_gene1 = set(df[df[gene1] > 0]['cell_barcode'])
    set_gene2 = set(df[df[gene2] > 0]['cell_barcode'])
    set_gene3 = set(df[df[gene3] > 0]['cell_barcode'])
    # set_slc17a6 = set(subset_data[subset_data['Slc17a6'] > 0]['cell_barcode'])
    
    # Plot Venn diagram
    venn3([set_gene1, set_gene2, set_gene3], (gene1, gene2, gene3))
    plt.title(f'{region}')

    plt.tight_layout()
    # plt.show()
        # Save in multiple formats if specified
    for fmt in save_formats:
        plt.savefig(f'{file_save}.{fmt}')
    plt.close()


def colocalizations(df, region, file_save, colors, save_formats):

    """
    This function looks at the colocalization of two genes compared to one rather than looking at the total expression of all the genes together.
    i.e. gene 1 will be 100% and gene 2/3 will be a subset of gene1.
    """
    # Define a threshold for positive expression (here using a simple nonzero assumption)
    gene_list = [gene1, gene2, gene3]

    for i in gene_list:
        # Filter cells that are positive for Crhr2
        positive_df = df[df[i] > 0]
        gene_listx = gene_list

        gene_listx = [m for m in gene_list if m != i]
        print(f"    Checking {gene_listx[0]} and {gene_listx[1]} that colocalize with {i}")
        # Determine overlap with the other two genes
        main_gene = set(positive_df[positive_df[i] > 0]['cell_barcode'])
        main_2gene = set(positive_df[positive_df[gene_listx[0]] > 0]['cell_barcode'])
        main_3gene = set(positive_df[positive_df[gene_listx[1]] > 0]['cell_barcode'])

        labels = (i, f'{i} & {gene_listx[0]}', f'{i} & {gene_listx[1]}')
        arrays = [main_gene, main_2gene, main_3gene]
        # Plot Venn diagram
        v = venn3(arrays, labels)

        if labels[0] == gene1:
            # Set custom colors
            v.get_patch_by_id('100').set_color(colors[0])    # Circle for the first gene
            
            p = (len(arrays[0])-len(arrays[1] & arrays[2]) - (len(arrays[1]) - len(arrays[1] & arrays[2])) - (len(arrays[2]) - len(arrays[1] & arrays[2])))*100/len(arrays[0])
            v.get_label_by_id('100').set_text(f"{p:.2f}%")

            v.get_patch_by_id('110').set_color(colors[1])  # Circle for the second gene
            p = (len(arrays[1]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('110').set_text(f"{p:.2f}%")

            v.get_patch_by_id('101').set_color(colors[2])   # Circle for the third gene
            p = (len(arrays[2]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('101').set_text(f"{p:.2f}%")
        
            p = (len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('111').set_text(f"{p:.2f}%")
            # v.get_patch_by_id('111').set_color('yellow')   # 2nd and 3rd overlap
        
        elif labels[0] == gene2: 
            # Set custom colors
            v.get_patch_by_id('100').set_color(colors[0])    # Circle for the first gene
            p = (len(arrays[0])-len(arrays[1] & arrays[2]) - (len(arrays[1]) - len(arrays[1] & arrays[2])) - (len(arrays[2]) - len(arrays[1] & arrays[2])))*100/len(arrays[0])
            v.get_label_by_id('100').set_text(f"{p:.2f}%")

            v.get_patch_by_id('110').set_color(colors[1])  # Circle for the second gene
            p = (len(arrays[1]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('110').set_text(f"{p:.2f}%")

            v.get_patch_by_id('101').set_color(colors[2])   # Circle for the third gene
            p = (len(arrays[2]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('101').set_text(f"{p:.2f}%")
        
            p = (len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('111').set_text(f"{p:.2f}%")
        
        elif labels[0] == gene3:
            # Set custom colors
            v.get_patch_by_id('100').set_color(colors[0])    # Circle for the first gene
            p = (len(arrays[0])-len(arrays[1] & arrays[2]) - (len(arrays[1]) - len(arrays[1] & arrays[2])) - (len(arrays[2]) - len(arrays[1] & arrays[2])))*100/len(arrays[0])
            v.get_label_by_id('100').set_text(f"{p:.2f}%")

            v.get_patch_by_id('110').set_color(colors[1])  # Circle for the second gene
            p = (len(arrays[1]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('110').set_text(f"{p:.2f}%")
        
            v.get_patch_by_id('101').set_color(colors[2])   # Circle for the third gene
            p = (len(arrays[2]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('101').set_text(f"{p:.2f}%")
        
            p = (len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('111').set_text(f"{p:.2f}%")

        plt.title(f'{region}_{i}')

        plt.tight_layout()
        # plt.show()
        # Save in multiple formats if specified
        for fmt in save_formats:
            plt.savefig(f'{file_save}.{fmt}')
        plt.close()

def make_folder(x, project_home):
    """
    makes new folders and or sets the directory
    """

    mode = 0o666
    j=os.listdir(project_home)
    if x not in j:
        behavior_path=os.path.join(project_home, x)
        os.mkdir(behavior_path, mode)
    else:
        behavior_paths=os.path.join(project_home, x)
    behavior_paths=os.path.join(project_home, x)
    return behavior_paths


# Load the YAML configuration file
def load_config(file_path):
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

# def unload_config(config):
    


if __name__ == "__main__": 
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main(variables_config)

