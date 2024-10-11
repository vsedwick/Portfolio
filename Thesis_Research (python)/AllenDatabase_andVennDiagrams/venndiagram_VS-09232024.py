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

    Args:
        variables_config (dict): Pre-determined parameters that will be utilized for the analysis.
    """
    warnings.filterwarnings("ignore")

    #Load and validate configuration parameters
    config = Config_Variables(load_config(variables_config))

        #create save_path
    _ = make_folder(config.dataset, config.image_save_folder)
    save_path = make_folder(config.combo, _)

    #Read file path
    expression_df = pd.read_csv(config.dataset_path)

    #using subclass to filter for MEA region and lhx6 transcription factor
    identify_region = expression_df['subclass'].str.contains('MEA' or 'Lhx6')

    #new dataframe only containing MEA and lhx6
    filtered_df = expression_df[identify_region]

    #split sexes
    female_df = filtered_df[filtered_df['donor_sex'].str.contains('F')]
    male_df = filtered_df[filtered_df['donor_sex'].str.contains('M')]

    #VENN DIAGRAMs for both sexes
    subclass_folder = make_folder('Subclasses', save_path)
    class_folder = make_folder('Classes', save_path); supertype_folder = make_folder('Supertypes', save_path)
    total_folder = make_folder('All Cells', save_path); coloc_folder = make_folder('Colocalizations', total_folder)
    
    #VENN DIAGRAMs for Females
    f_subclass_folder = make_folder('Female', subclass_folder); f_class_folder = make_folder('Female', class_folder)
    f_supertype_folder = make_folder('Female', supertype_folder); 

    #VENN DIAGRAMs for Males
    m_subclass_folder = make_folder('Male', subclass_folder); m_class_folder = make_folder('Male', class_folder)
    m_supertype_folder = make_folder('Male', supertype_folder); 

    generate_and_save(filtered_df, 'ALL', [subclass_folder, class_folder, supertype_folder, total_folder, coloc_folder], config)
    generate_and_save(female_df, 'F', [f_subclass_folder, f_class_folder, f_supertype_folder, total_folder, coloc_folder], config)
    generate_and_save(male_df, 'M', [m_subclass_folder, m_class_folder, m_supertype_folder, total_folder, coloc_folder], config)


class Config_Variables:
    """
    A class used to call and validate variables from the configuration file.
    
    Attributes:
        dataset : str
            The name of the dataset used for the analysis.
        image_save_folder : str
            The directory where the images will be saved. Created if it does not exist.
        dataset_path : str
            The file path to the dataset containing single-cell sequencing data.
        gene1 : str
            The first gene used for generating Venn diagrams.
        gene2 : str
            The second gene used for generating Venn diagrams.
        gene3 : str
            The third gene used for generating Venn diagrams.
        combo : str
            A string concatenating the three gene names for labeling purposes.
        region : str
            The brain region used for labeling Venn diagrams.
        colors : list or tuple
            A list or tuple of exactly three color codes (Hex codes, named colors, or RGB tuples).
        save_formats : list
            A list of valid image formats in which the Venn diagrams will be saved.

    Methods:
        __init__(self, config):
            Initializes the configuration settings by validating all required parameters.
    """

    def __init__(self, config):
        """
        Initializes the Config_Variables class with settings from the configuration dictionary.

        Args:
            config (dict): A dictionary containing configuration settings for the dataset, folder paths, gene parameters,
                        color settings, and image save formats.

        Raises:
            ValueError: If the number of genes or colors is not exactly 3.
            ValueError: If any of the listed image save formats are not acceptable.
        """        
    
        self.dataset = config['dataset']
        image_save_folder = config['image_save_folder']
        if not os.path.exists(image_save_folder):
            partition = image_save_folder.split("\\")
            baby_folder = partition[-1]; parents = "\\".join(partition[0:-1])
            self.image_save_folder = make_folder(baby_folder, parents)
        else:
            self.image_save_folder = image_save_folder
        self.dataset_path = config['dataset_path']

        #check for gene paramter
        if len(config['genes']) != 3:
            raise ValueError("Must list 3 genes. This program currently does not support less or more.")
        else:
            self.gene1 = config['genes']['gene1']
            self.gene2 = config['genes']['gene2']
            self.gene3 = config['genes']['gene3']
        
        self.combo = f'{self.gene1}-{self.gene2}-{self.gene3}'

        #Brain region, will be used for labeling purposes
        self.region = config['venn_settings']['region']

        colors = config['venn_settings']['colors']
        if len(colors) != 3 or not isinstance(colors, (list, tuple)):
            raise ValueError("Must list 3 color codes. Hex Codes, named colors, and RGB tuples are accepted.")
        else:
            self.colors = colors

        save_formats = config['venn_settings']['save_formats']
        acceptable_formats = ['png', 'pdf', 'svg', 'eps', 'jpg', 'jpeg', 'tiff', 'tif', 'raw', 'ps', 'pgf']
        if not all(fmt in acceptable_formats for fmt in save_formats):
            raise ValueError('The listed format(s) is not acceptable for images.')
        else:
            self.save_formats = save_formats

def extract_cellcluster_names(df):
    """
    Extracts the unique names of cell clusters from different hierarchical levels.

    Args:
        df (pd.DataFrame): A DataFrame containing single-cell sequencing data loaded from the config file.

    Returns:
        Unique names of cell clusters at different levels.
            - df['subclass'].unique(): Unique subclass names.
            - df['class'].unique(): Unique class names.
            - df['supertype'].unique(): Unique supertype names.
    """
    return df['subclass'].unique(), df['class'].unique(), df['supertype'].unique()

def generate_and_save(df, sex, directory_list, config):

    """
    Processes venndiagrams for all categories

    Args:
        df (pd.DataFrame): A DataFrame filtered by region and sex (filtered from original file defined in the config)
        sex (string): The sex of the sample being processed. Primarily used for naming and labeling purposes.
        directory_list (list of str): A list of file paths where venn diagram images will be saved
        config (dict): Variables from configuration file that will be used to call region, save_formats, and colors
    
    Returns:
        None: The function saves venn diagrams to their respective directories.

    """


    print(f"Processing {sex}_{config.region}")
    subclass_list, class_list, supertype = extract_cellcluster_names(df);
    venndiagrams(df, subclass_list, 'subclass', fr"{directory_list[0]}\{sex}_Subclass", config)
    venndiagrams(df, class_list, 'class', fr"{directory_list[1]}\{sex}", config)
    venndiagrams(df, supertype, 'supertype', fr"{directory_list[2]}\{sex}_Supertype", config)
    total_region(df, f'{config.region}_{sex}', fr"{directory_list[3]}\{sex}", config)
    colocalizations(df, f'{config.region}_{sex}', fr"{directory_list[4]}\{sex}", config)

    print(f"Images for {sex} {config.region} genes are saved.")


def venndiagrams(df, unique_names, level, file_save, config):

    """
    This function determines the proportions for the venn diagrams vased on the  values in the dataset
    The final output will be the venn diagram image as both a vector and png.

    Args:
        df (pd.DataFrame): A DataFrame filtered by region and sex (filtered from original file defined in the config)
        unique_names (list of str): A list of unique cell cluster names at a specific hierarchichal level
        level (str): hierarchichal level (class, subclass, supertype). Used to further filter dataframe.
        file_save (str): Directory where generated images will be saved.
        config (dict): Variables from configuration file. Used to call gene names and file formats that images will be saved as.
    
    Returns:
        None. Generates and saves venn diragram plots to the specified directories
    """

    for i, class_name in enumerate(unique_names):
        plt.figure()
        # Subset data for each class
        subset_data = df[df[level] == class_name]
        
        # Determine sets for genes based on non-zero expression
        set_gene1 = set(subset_data[subset_data[config.gene1] > 0]['cell_barcode'])
        set_gene2 = set(subset_data[subset_data[config.gene2] > 0]['cell_barcode'])
        set_gene3 = set(subset_data[subset_data[config.gene3] > 0]['cell_barcode'])
        
        # Plot Venn diagram
        venn3([set_gene1, set_gene2, set_gene3], (config.gene1, config.gene2, config.gene3))
        plt.title(f'{class_name}')

        plt.tight_layout()

        # Save in multiple formats if specified
        for fmt in config.save_formats:
            plt.savefig(f'{file_save}.{fmt}')
        plt.close()

def total_region(df, region, file_save, config):
    """This function calculates the total number of cells that express the three defined genes with
    proportionate circles. 
    
    Args:
        df (pd.DataFrame): A DataFrame filtered by region and sex (filtered from original file defined in the config)
        region (str): A label that includes the sex and brain region of the data indexed
        file_save (str): The file directory where the images will be saved.
        config (dict): Variables from configuration file. Used to call gene names and file formats that images will be saved as.
    
    Returns: 
        None. Saves the total cell counts without considering clusters.
        """

    # Determine sets for genes based on non-zero expression
    set_gene1 = set(df[df[config.gene1] > 0]['cell_barcode'])
    set_gene2 = set(df[df[config.gene2] > 0]['cell_barcode'])
    set_gene3 = set(df[df[config.gene3] > 0]['cell_barcode'])
    # set_slc17a6 = set(subset_data[subset_data['Slc17a6'] > 0]['cell_barcode'])
    
    # Plot Venn diagram
    venn3([set_gene1, set_gene2, set_gene3], (config.gene1, config.gene2, config.gene3))
    plt.title(f'{region}')

    plt.tight_layout()
    # plt.show()
        # Save in multiple formats if specified
    for fmt in config.save_formats:
        plt.savefig(f'{file_save}.{fmt}')
    plt.close()


def colocalizations(df, region, file_save, config):

    """
    This function looks at the colocalization of two genes compared to one rather than looking at the total expression of all the genes together.
    i.e. gene 1 will be 100% and gene 2/3 will be a subset of gene1.
    
    Args:
        df:
        region:
        file_save:
        config: 

    Returns:
        None. Indexes cell counts iterated over a single gene to quantify and visualize colocalizations of the remaining genes.
    """
    # Define a threshold for positive expression (here using a simple nonzero assumption)
    gene_list = [config.gene1, config.gene2, config.gene3]

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

        if labels[0] == config.gene1:
            # Set custom colors
            v.get_patch_by_id('100').set_color(config.colors[0])    # Circle for the first gene
            
            p = (len(arrays[0])-len(arrays[1] & arrays[2]) - (len(arrays[1]) - len(arrays[1] & arrays[2])) - (len(arrays[2]) - len(arrays[1] & arrays[2])))*100/len(arrays[0])
            v.get_label_by_id('100').set_text(f"{p:.2f}%")

            v.get_patch_by_id('110').set_color(config.colors[1])  # Circle for the second gene
            p = (len(arrays[1]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('110').set_text(f"{p:.2f}%")

            v.get_patch_by_id('101').set_color(config.colors[2])   # Circle for the third gene
            p = (len(arrays[2]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('101').set_text(f"{p:.2f}%")
        
            p = (len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('111').set_text(f"{p:.2f}%")
            # v.get_patch_by_id('111').set_color('yellow')   # 2nd and 3rd overlap
        
        elif labels[0] == config.gene2: 
            # Set custom colors
            v.get_patch_by_id('100').set_color(config.colors[0])    # Circle for the first gene
            p = (len(arrays[0])-len(arrays[1] & arrays[2]) - (len(arrays[1]) - len(arrays[1] & arrays[2])) - (len(arrays[2]) - len(arrays[1] & arrays[2])))*100/len(arrays[0])
            v.get_label_by_id('100').set_text(f"{p:.2f}%")

            v.get_patch_by_id('110').set_color(config.colors[1])  # Circle for the second gene
            p = (len(arrays[1]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('110').set_text(f"{p:.2f}%")

            v.get_patch_by_id('101').set_color(config.colors[2])   # Circle for the third gene
            p = (len(arrays[2]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('101').set_text(f"{p:.2f}%")
        
            p = (len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('111').set_text(f"{p:.2f}%")
        
        elif labels[0] == config.gene3:
            # Set custom colors
            v.get_patch_by_id('100').set_color(config.colors[0])    # Circle for the first gene
            p = (len(arrays[0])-len(arrays[1] & arrays[2]) - (len(arrays[1]) - len(arrays[1] & arrays[2])) - (len(arrays[2]) - len(arrays[1] & arrays[2])))*100/len(arrays[0])
            v.get_label_by_id('100').set_text(f"{p:.2f}%")

            v.get_patch_by_id('110').set_color(config.colors[1])  # Circle for the second gene
            p = (len(arrays[1]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('110').set_text(f"{p:.2f}%")
        
            v.get_patch_by_id('101').set_color(config.colors[2])   # Circle for the third gene
            p = (len(arrays[2]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('101').set_text(f"{p:.2f}%")
        
            p = (len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('111').set_text(f"{p:.2f}%")

        plt.title(f'{region}_{i}')

        plt.tight_layout()
        # plt.show()
        # Save in multiple formats if specified
        for fmt in config.save_formats:
            plt.savefig(f'{file_save}.{fmt}')
        plt.close()

def make_folder(parent_directory, new_folder):
    """
    Creates new folders and or sets the directory. 

    Args:
        parent_directory (str): The parent directory for the new folder.
        new_folder (str): The name of the new fodler to be created.

    Returns:
        full_path (str): The new directory where the folder was created
    
    Raises: 
        FileNotFoundError: If the specified parent directory does not exist.
        PermissionError: If the directory is not writable or the folder cannot be created.

    """

    mode = 0o666

    full_path = os.path.join(parent_directory, new_folder)

    #check if parent directory exists
    if not os.path.exists(parent_directory):
        raise FileNotFoundError(f"Parent directory '{parent_directory}' does not exist.")
    
    #checks if user has permission to write to that directory
    if not os.access(parent_directory, os.W_OK):
        raise PermissionError(f"Write permission denied for directory '{parent_directory}.")
    
    #Creates the folder if it doesnt exists
    if not os.path.exists(full_path):
        try:
            os.mkdir(full_path, mode)
        except OSError:
            raise PermissionError(f"Failed to create directory {full_path}. Check permissions: {OSError}")


    return full_path


# Load the YAML configuration file
def load_config(variables_config):
    """ 
    Loads the configuration file. can be specified as an argument or defined in the script.
    
    Args: 
        variables_config (str): The file path that directs to the config file. 
        
    Returns:
        config (dict): A dictionary of the configuration file that contains the analysis parameters.

    """
    with open(variables_config, 'r') as file:
        config = yaml.safe_load(file)
    return config

# def unload_config(config):
if __name__ == "__main__": 
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main(variables_config)

