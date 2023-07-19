import re
template_root = 'mbank_X27573_2023-07-17-1103-no'
template = open(template_root + '.nex', 'r')
Lines = template.readlines()

output = open(template_root + '-vetu.nex', 'w')

for line in Lines:
      ntax = re.search(r'(NTAXA?\s*=\s*)(\d+)\b', line)
      if ntax:
            output.writelines(re.sub(
                  ntax.group(),
                  ntax.group(1) + str(int(ntax.group(2)) - 3),
                  line))
      elif not re.search(r'Banffia|Vetulicola|Beidazoon|Vetulocystis', line):
          output.writelines(line)

output = open(template_root + '-vetulicolans.nex', 'w')
for line in Lines:
      ntax = re.search(r'(NTAXA?\s*=\s*)(\d+)\b', line)
      if ntax:
            output.writelines(re.sub(
                  ntax.group(),
                  ntax.group(1) + str(int(ntax.group(2)) - 2),
                  line))
      elif not re.search(r'Banffia|Vetulicola|Beidazoon', line):
          output.writelines(line)

                      
