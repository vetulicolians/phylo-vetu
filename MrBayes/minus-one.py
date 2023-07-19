import re
template_root = 'mbank_X27573_2023-07-17-1103-no'
template = open(template_root + '.nex', 'r')
Lines = template.readlines()

i = 0
for line in Lines:
      i += 1
      if re.search(r'^\s*MATRIX\s*$', line):
          matrixStart = i + 1
      elif re.search(r'^\s*;\s*$', line):
          matrixEnd = i - 1
          break

for i in range(matrixStart, matrixEnd):
    taxon = re.search(r'^\s+(\w+)', Lines[i]).group(1)
    newFile = open(template_root + '-' + taxon + '.nex', 'w')
    newFile.writelines(Lines[:i] + Lines[i+1:])
                      
