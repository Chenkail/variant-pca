# Import libraries
# Built-in
import math

# 3rd party
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


class VariantCollection:
    def __init__(self):
        """Initialize data in VariantCollection object"""
        
        # Create empty dataframe
        self.dataframe = None
        # Initialize header row list
        self.header = ["Cell"]
        # Initialize list of rows
        self.rows = []
        # Other variables
        self.variant_appearances = ["Total appearances:"]
        self.variants_in_row = []
    
    # -------- Class interface -------- #
    
    def import_data(self, in_file, mode=""):
        """Import data from a single vcf file"""
        
        chrm_index = 0
        pos_index = 1
        ref_index = 2
        alt_index = 3

        with open(in_file) as file:
            table_header = file.readline()
            data = file.readlines()
        
        # Initialize data list and set first element to be the cell barcode
        data_list = [0] * len(self.header)
        data_list[0] = in_file.split("/")[-1].split(".")[0]
        variant_count = 0
        

        for line in data:
            variant_count += 1

            line = line.strip().split("\t")
            chrm = line[chrm_index]
            pos = line[pos_index]
            ref = line[ref_index]
            alt = line[alt_index]

            location = chrm + "-" + pos
            delta = ref + "->" + alt
            variant = location + ":" + delta
            
            # If variant already in list
            if variant in self.header:
                variant_index = self.header.index(variant)
                data_list[variant_index] = 1
                self.variant_appearances[variant_index] += 1
            else:
                self.header.append(variant)
                self.variant_appearances.append(1)
                data_list.append(1)
            
        self.rows.append(data_list)
        self.variants_in_row.append(variant_count)
        self.fill_blanks()

        print("Imported data from " + in_file)
    
    
    def mass_import(self, input_file):
        """Imports from all files listed in input file"""
        
        cell_list = self.cell_list(input_file)
        total = len(cell_list)
        print(str(total) + " filenames found.")
        
        imported = 0
        failed = 0
        
        for file in cell_list:
            try:
                self.import_data(file)
                imported += 1
            except FileNotFoundError:
                print("File not found: " + file)
                failed += 1

        print("Mass import from " + input_file + " complete. " + 
              str(imported) + "/" + str(total) + " files imported, " + 
              str(failed) + " files not found.")
        
    
    def export_data(self, file="variants.csv"):
        """Exports data into a text file"""
        
        self.dataframe.to_csv(file, index=False)
    
    
    def filter_variants(self, minscore=5, mode=None, data_filter="variant"):
        """
        Remove all variants which don't show up above a certain number of times
        
        NOTE: WILL CAUSE UNEXPECTED BEHAVIOR IF LABELS ARE NOT UNIQUE
        """
        
        if data_filter == "cell":
            to_delete = []
            subcount = [0] * len(self.header)
            # Throw errors if attempting to access first element (safety)
            subcount[0] = None
            
            if mode == "trim-bottom":
                # Trim bottom n% of cells by variant count
                sorted_counts = sorted(self.variants_in_row)
                min_index = int(math.ceil((minscore/100) * 
                                          len(self.variants_in_row)))
                min_count = sorted_counts[min_index]

            else:
                min_count = minscore
            
            # Find rows to remove and add row index to list
            for row in range(len(self.rows)):
                if self.variants_in_row[row] < min_count:
                    to_delete.append(row)


            # Adjust appearance count
            for row in to_delete:
                for i in range(1, len(self.header)):
                    if self.rows[row][i] == 1:
                        self.variant_appearances[i] -= 1
                    
            # Remove rows
            for row in range(len(self.rows) - 1, -1, -1):
                if row in to_delete:
                    del self.rows[row]
        
        else:
            # Variant mode
            if mode == "percent":
                # Variant must show up in certain % of cells
                min_count = len(self.rows)
                min_count = int(math.ceil((minscore/100) * cell_count))
            elif mode == "trim-bottom":
                # Variant must not be in bottom n% of all variants
                sorted_counts = sorted(self.variant_appearances)
                min_index = int(math.ceil((minscore/100) * 
                                          len(self.variant_appearances)))
                min_count = sorted_counts[min_index]
            else:
                min_count = minscore


            for i in range(len(self.variant_appearances) - 1, 0, -1):
                if self.variant_appearances[i] < min_count:
                    # Wipe variant from database
                    del self.header[i]
                    for row in range(len(self.rows)):
                        del self.rows[row][i]
                    del self.variant_appearances[i]
    
    # -------- Visualizations -------- #

    def plot_kmeans(self):
        snips = self.dataframe.iloc[:, 1:].values
        snips = StandardScaler().fit_transform(snips)

        # PCA
        pca = PCA(n_components=2)
        prcomp = pca.fit_transform(snips)
        prinDf = pd.DataFrame(data=prcomp, columns=['PC1', 'PC2'])
        finalDf = pd.concat([prinDf, self.dataframe[['Cell']]], axis = 1)

        # K-means
        kmeans = KMeans(n_clusters=2).fit(snips)

        # Plot using k-means colors
        plt.scatter(prinDf['PC1'], prinDf['PC2'], c= kmeans.labels_.astype(float), s=50, alpha=0.75)

    # -------- Data syncing -------- #

    def generate_dataframe(self):
        """Create pandas dataframe from lists"""

        self.dataframe = pd.DataFrame(self.rows, columns=self.header)


    def update_data(self):
        """Update self.* variables to match self.dataframe"""

        self.header = self.dataframe.columns.tolist()
        self.rows = self.dataframe.values.tolist()
        self.variant_appearances = ["Total appearances:"]
        self.variants_in_row = []
    
    # -------- Backend -------- #

    def cell_list(self, input_file):
        """Returns barcodes given list"""

        with open(input_file) as file:
            data = file.readlines()

        cells = []
        for line in data:
            line2 = line.strip()
            if line2 != "":
                cells.append(line2)

        return cells
    
    
    def fill_blanks(self):
        """Add zeros to empty spots in table"""
        
        width = 0
        for row in self.rows:
            width = max(width, len(row))
        
        for row in self.rows:
            if len(row) < width:
                zeros = [0] * (width-len(row))
                row += zeros

