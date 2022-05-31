import csv
csv_file = 'designs.csv'
txt_file = 'class.txt'
with open(txt_file, "w") as output_file:
    with open(csv_file, "r") as input_file:
        for row in csv.reader(input_file):
            if row[0] == 'group':
                title = row
                continue
            indent_text = "p_data['" + row[0] + "'] = {"
            indent = " "*len(indent_text)
            output_file.write(indent_text)
            output_file.write("'" + title[1] + "':" + row[1])
            i = 2
            while i < len(row): 
                output_file.write(",\n"+indent)
                output_file.write("'" + title[i] + "':" + row[i])
                i = i+1
            output_file.write("\n" + indent + "}\n\n")
    output_file.close()
print('completed')