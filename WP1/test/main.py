# This is a sample Python script.
import fragment.utils
import anndata as ad
import scanpy

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    #print(fragment.utils.read_fragment_file("/home/leon/Documents/s.bed"))
    #print(fragment.utils.read_fragment_directory("/home/leon/Documents/"))
    #print(fragment.utils.combine_fragment_dictionaries(fragment.utils.read_fragment_directory("/home/leon/Documents/")))
    a = fragment.utils.create_h5ad_object(fragment.utils.combine_fragment_dictionaries(fragment.utils.read_fragment_directory("/home/leon/Documents/")))
    print(a.obs_names)
    print(a.obs)
    print(a.to_df())
    print(a.obs_keys())
    a.to_df().at["AAACGCAAGCAAACCCGAGATA", "Mean"] = 1
    print(a.to_df())

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
