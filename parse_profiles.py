import csv

multiple_list = list()
with open('profile.txt') as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    for row in reader:
        for k, v in row.items():
            if 'BACT' in k and 'N' == v:
                if k not in multiple_list:
                    #print(row)
                    multiple_list.append(k)

for item in multiple_list:
    print(item)
