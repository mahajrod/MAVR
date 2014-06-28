#!/usr/bin/env python
__author__ = 'mahajrod'
import argparse
import MySQLdb as mdb

#example of usage
#./testMySQL.py --host bioinfinst07.db.8193575.hostedresource.com -u bioinfinst07 -p Bioinfinst123! -d bioinfinst07


def load_data_to_tbl(cursor, table_name, data_file, table_field_list):
    with open(data_file, "r") as in_fd:
        header = in_fd.readline().strip().split(",")

        indexes = {}
        for field in table_field_list:
            indexes[field] = header.index(field)
        for line in in_fd:
            tmp = line.strip().split(",")
            data_string_list = []
            for field in table_field_list:
                data_string_list.append("'" + tmp[indexes[field]] + "'")
                #print (line)
            cursor.execute("INSERT INTO %s (%s) VALUES(%s);" % (table_name, ",".join(table_field_list), ",".join(data_string_list)))


parser = argparse.ArgumentParser()
parser.add_argument('--host', action='store', dest='host', help='Host of database')
parser.add_argument('-u', action='store', dest='user', help='Username')
parser.add_argument('-p', action='store', dest='pwd', help='Password')
parser.add_argument('-d', action='store', dest='db', help='Database')
arg = parser.parse_args()

for argument in [arg.host, arg.user, arg.pwd, arg.db]:
    if not argument:
        raise ValueError("One or more mandatory argument was/were not set")

gene_list_file = "gene_list.csv"
tuberculist_file = "tuberculist.csv"
gene_func_cat_file = "gene_func_cat.csv"

table_set = set(["gene_list", "gene_func_cat", "tuberculist"])
gene_list_fields = ["gene_name",
                    "locus_tag",
                    "gene_id",
                    "start",
                    "end",
                    "strand",
                    "feature_id",
                    "uniprotdb"]

tuberculist_fields = ["locus_tag",
                      "protein_function",
                      "location",
                      "mass",
                      "func_cat"]

gene_func_cat_fields = ["function_id",
                        "function_name"]

cnx = mdb.connect(host=arg.host, user=arg.user, passwd=arg.pwd)
cursor = cnx.cursor()

#Check for presence of database
cursor.execute("SHOW DATABASES;")
existing_db = cursor.fetchall()
existing_db = [line[0] for line in existing_db]
print("Checking presense of database: %s" % arg.db)
if arg.db not in existing_db:
    raise IOError("Database %s was not found on %s host (username: %s)." % (arg.db, arg.host, arg.user))
else:
    print("\tFound.")

cursor.execute("use %s;" % arg.db)
cnx.commit()

#Check for presence of tables
cursor.execute("SHOW TABLES;")
existing_tbl = cursor.fetchall()
existing_tbl = set([line[0] for line in existing_tbl])


if existing_tbl & table_set:
    raise IOError("One or more tables are already present")

cursor.execute("CREATE TABLE gene_list (gene_name varchar(20), locus_tag varchar(20), gene_id varchar(20), start int(255), end int(255), strand varchar(5), feature_id varchar(20), uniprotdb varchar(20) );")
cursor.execute("CREATE TABLE gene_func_cat(function_id int(20),function_name varchar(100));")
cursor.execute("CREATE TABLE tuberculist(locus_tag varchar(20),protein_function varchar(100),location varchar(20),mass float(20),func_cat int(20));")

cursor.execute("ALTER TABLE gene_list ADD CONSTRAINT locus_tag FOREIGN KEY (locus_tag) REFERENCES tuberculist(locus_tag);")
cursor.execute("ALTER TABLE tuberculist ADD CONSTRAINT func_cat FOREIGN KEY (func_cat) REFERENCES gene_func_cat(function_id);")

load_data_to_tbl(cursor, "gene_list", gene_list_file, gene_list_fields)
load_data_to_tbl(cursor, "gene_func_cat", gene_func_cat_file, gene_func_cat_fields)
load_data_to_tbl(cursor, "tuberculist", tuberculist_file, tuberculist_fields)
cursor.close()
cnx.close()


