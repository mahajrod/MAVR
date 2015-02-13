__author__ = 'mahajrod'

import numpy as np


def parse_R_table_file(table_file, separator=","):
    with open(table_file, "r") as fd_table:
        column_names = [name[1:-1]for name in fd_table.readline().strip().split(separator)[1:]]
        row_names = []
        table = []
        for line in fd_table:
            if line == "\n":
                break
            temp_line = line.strip().split(separator)
            row_names.append(temp_line[0][1:-1])
            table.append([float(x) for x in temp_line[1:]])
    return column_names, row_names, np.array(table)


def get_indices_from_names(names_list, queue_names_list):
    result = []
    for entry in queue_names_list:
        result.append(names_list.index(entry))
    return sorted(result)


#def get_submatrix_by_column(table, column_list, row_list):
#    return table[row_list][column_list]


if __name__ == "__main__":
    column_names, row_names, table = parse_R_table_file("/home/mahajrod/genetics/MH_selection/trees/all_species/tree_dist.csv", separator=",")
    print(column_names)
    print(row_names)
    print(np.shape(table))
    indices = get_indices_from_names(row_names, ['AY442348.1', 'NC_001778.1', 'NC_020655.1', 'HM143934.1'])
    #print(indices)
    extracted_table = table[np.array(indices), :][:, np.array(indices)]
    print(extracted_table)
    print (extracted_table[np.triu_indices(np.shape(extracted_table)[0], 1)])
    #print(get_submatrix_by_column(table,np.array(indices),indices))