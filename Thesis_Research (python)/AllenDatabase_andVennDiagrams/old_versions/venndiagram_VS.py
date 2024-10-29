import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os
# from matplotlib_venn import venn4

combo = 'crfr2-vgat-vglut2'
dataset = 'STR_dataset'
image_save_folder = r"C:\Users\sedwi\Desktop\SingleCell_AllenAtlas\venndiagram_drafts_VS\Images"

gene1 = 'Crhr2'
gene2 = 'Slc32a1'
gene3 = ''


def main():
    global gene1, gene2, gene3, combo
    file_path =  r"C:\Users\sedwi\Desktop\SingleCell_AllenAtlas\STR_Region_Expression.csv"

    #create save_path
    save_path = os.path.join(image_save_folder, f'{dataset}\\{combo}')
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
    total_folder = make_folder('Combined', save_path); coloc_folder = make_folder('Colocalizations', total_folder)
    venndiagrams(filtered_df, subclass_list, 'subclass', fr"{subclass_folder}\Sublass")
    venndiagrams(filtered_df, class_list, 'class', fr"{class_folder}\Class")
    venndiagrams(filtered_df, supertype, 'supertype', fr"{supertype_folder}\Supertype")
    total_region(filtered_df, 'MeP', fr"{total_folder}\All")
    colocalizations(filtered_df, 'MeP', fr"{coloc_folder}\All")

    #VENN DIAGRAMs for Females
    subclass_list, class_list, supertype = make_lists(female_df)
    f_subclass_folder = make_folder('Female', subclass_folder); f_class_folder = make_folder('Female', class_folder)
    f_supertype_folder = make_folder('Female', supertype_folder); 
    venndiagrams(female_df, subclass_list, 'subclass', fr"{f_subclass_folder}\F_Sublass")
    venndiagrams(female_df, class_list, 'class', fr"{f_class_folder}\F_Class")
    venndiagrams(female_df, supertype, 'supertype', fr"{f_supertype_folder}\F_Supertype")
    total_region(female_df, 'F_MeP', fr"{total_folder}\F")
    colocalizations(female_df, 'F_MeP', fr"{coloc_folder}\F")

    #VENN DIAGRAMs for Males
    subclass_list, class_list, supertype = make_lists(male_df)
    m_subclass_folder = make_folder('Male', subclass_folder); m_class_folder = make_folder('Male', class_folder)
    m_supertype_folder = make_folder('Male', supertype_folder); 
    venndiagrams(male_df, subclass_list, 'subclass', fr"{m_subclass_folder}\M_Sublass")
    venndiagrams(male_df, class_list, 'class', fr"{m_class_folder}\M_Class")
    venndiagrams(male_df, supertype, 'supertype', fr"{m_supertype_folder}\M_Supertype")
    total_region(male_df, 'M_MeP', fr"{total_folder}\M")
    colocalizations(male_df, 'M_MeP', fr"{coloc_folder}\M")

def make_lists(df):
    return df['subclass'].unique(), df['class'].unique(), df['supertype'].unique()


def venndiagrams(df, item, name, file_save):
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
        plt.savefig(fr"{file_save}_{class_name}.png")
        plt.savefig(fr"{file_save}_{class_name}.svg")
        plt.close()

def total_region(df, region, file_save):

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
    plt.savefig(fr"{file_save}_{region}.png")
    plt.savefig(fr"{file_save}_{region}.svg")
    plt.close()


def colocalizations(df, region, file_save):
    # Define a threshold for positive expression (here using a simple nonzero assumption)
    gene_list = [gene1, gene2, gene3]
    print(gene_list)

    for i in gene_list:
        # Filter cells that are positive for Crhr2
        positive_df = df[df[i] > 0]
        gene_listx = gene_list

        gene_listx = [m for m in gene_list if m != i]
        print(gene_listx)
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
            v.get_patch_by_id('100').set_color((1, 0, 0, 0.5))    # Circle for the first gene
            
            p = (len(arrays[0])-len(arrays[1] & arrays[2]) - (len(arrays[1]) - len(arrays[1] & arrays[2])) - (len(arrays[2]) - len(arrays[1] & arrays[2])))*100/len(arrays[0])
            v.get_label_by_id('100').set_text(f"{p:.2f}%")

            v.get_patch_by_id('110').set_color((0, 1, 0, 0.5))  # Circle for the second gene
            p = (len(arrays[1]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('110').set_text(f"{p:.2f}%")

            v.get_patch_by_id('101').set_color((0, 0, 1, 0.5))   # Circle for the third gene
            p = (len(arrays[2]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('101').set_text(f"{p:.2f}%")
        
            p = (len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('111').set_text(f"{p:.2f}%")
            # v.get_patch_by_id('111').set_color('yellow')   # 2nd and 3rd overlap
        
        elif labels[0] == gene2: 
            # Set custom colors
            v.get_patch_by_id('100').set_color((0, 1, 0, 0.5))    # Circle for the first gene
            p = (len(arrays[0])-len(arrays[1] & arrays[2]) - (len(arrays[1]) - len(arrays[1] & arrays[2])) - (len(arrays[2]) - len(arrays[1] & arrays[2])))*100/len(arrays[0])
            v.get_label_by_id('100').set_text(f"{p:.2f}%")

            v.get_patch_by_id('110').set_color((1, 0, 0, 0.5))  # Circle for the second gene
            p = (len(arrays[1]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('110').set_text(f"{p:.2f}%")

            v.get_patch_by_id('101').set_color((0, 0, 1, 0.5))   # Circle for the third gene
            p = (len(arrays[2]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('101').set_text(f"{p:.2f}%")
        
            p = (len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('111').set_text(f"{p:.2f}%")
        
        elif labels[0] == gene3:
            # Set custom colors
            v.get_patch_by_id('100').set_color((0, 0, 1, 0.5))    # Circle for the first gene
            p = (len(arrays[0])-len(arrays[1] & arrays[2]) - (len(arrays[1]) - len(arrays[1] & arrays[2])) - (len(arrays[2]) - len(arrays[1] & arrays[2])))*100/len(arrays[0])
            v.get_label_by_id('100').set_text(f"{p:.2f}%")

            v.get_patch_by_id('110').set_color((1, 0, 0, 0.5))  # Circle for the second gene
            p = (len(arrays[1]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('110').set_text(f"{p:.2f}%")
        
            v.get_patch_by_id('101').set_color((0, 1, 0, 0.5))   # Circle for the third gene
            p = (len(arrays[2]) - len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('101').set_text(f"{p:.2f}%")
        
            p = (len(arrays[1] & arrays[2]))*100/len(arrays[0])
            v.get_label_by_id('111').set_text(f"{p:.2f}%")

        plt.title(f'{region}_{i}')

        plt.tight_layout()
        # plt.show()
        plt.savefig(fr"{file_save}_{region}_{i}.png")
        plt.savefig(fr"{file_save}_{region}_{i}.svg")
        plt.close()

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

main()