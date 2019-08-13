def list_to_tab_string(data):
    """Convert a list of strings to a string with elements separated by tabs"""

    string = ""
    for text in data:
        string = string + text + "\t"
    # Remove tab at end
    string = string[:-1]
    return string

class VariantCollection:
    def __init__(self):
        self.clear_data()
    
    
    def clear_data(self):
        """Reset and initialize data in VariantCollection object"""
        
        # Initialize header row list
        self.header = ["Cell"]
        # Initialize list of rows
        self.rows = []
    
    
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

    
    def import_data(self, in_file):
        """Import data from a single vcf file"""
        
        chrm_index = 0
        pos_index = 1
        ref_index = 2
        alt_index = 3

        with open(in_file) as file:
            table_header = file.readline()
            data = file.readlines()
        
        # Initialize data list and set first element to be the cell barcode
        data_list = ["0"] * len(self.header)
        data_list[0] = in_file.split("/")[-1].split(".")[0]

        for line in data:
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
                data_list[self.header.index(variant)] = "1"
            else:
                self.header.append(variant)
                data_list.append("1")
            
        self.rows.append(data_list)
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

        print("Mass import from " + input_file + " complete.")
        print(str(imported) + "/" + str(total) + " files imported, " + 
              str(failed) + " files not found.")
        
    
    def fill_blanks(self):
        """Add zeros to empty spots in table"""
        
        width = 0
        for row in self.rows:
            width = max(width, len(row))
        
        for row in self.rows:
            if len(row) < width:
                zeros = ["0"] * (width-len(row))
                row += zeros
    
    
    def export_data(self, file="variants.txt"):
        """Exports data into a tab-separated text file"""
        
        with open(file, 'w') as output:
            # Add header
            header = list_to_tab_string(self.header)
            output.write(header + "\n")
            
            # Add data 
            for line in self.rows:
                data = list_to_tab_string(line)
                output.write(data + "\n")
