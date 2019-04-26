import sys

for line in sys.stdin:
    if ',' in line:
        fileSection, moduleSection, funcSection = line.split(',')
        funcSection = funcSection.strip()
        funcName = funcSection.split()[1]
        print(funcName)