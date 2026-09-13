import re
import os

def extract_subroutines(file_path):
    subroutines = set()
    with open(file_path, 'r', encoding='latin-1', errors='ignore') as f:
        for line in f:
            match = re.match(r'^\s*subroutine\s+(\w+)', line, re.IGNORECASE)
            if match:
                subroutines.add(match.group(1).lower())
    return sorted(subroutines)

datcom_f = 'datcom-legacy/datcom.f'
combined_f = 'datcom-legacy/datcom_2000/combined_fortran.f'

if os.path.exists(datcom_f) and os.path.exists(combined_f):
    subs_datcom = extract_subroutines(datcom_f)
    subs_combined = extract_subroutines(combined_f)
    
    print('Subroutines in datcom.f:', len(subs_datcom))
    print('Subroutines in combined_fortran.f:', len(subs_combined))
    
    # Find common, only in datcom, only in combined
    common = set(subs_datcom) & set(subs_combined)
    only_datcom = set(subs_datcom) - set(subs_combined)
    only_combined = set(subs_combined) - set(subs_datcom)
    
    print('Common subroutines:', len(common))
    print('Only in datcom.f:', len(only_datcom))
    print('Only in combined_fortran.f:', len(only_combined))
    
    # Save to file
    with open('datcom-legacy/exist_in_both.txt', 'w') as f:
        f.write('# Subroutine Comparison: datcom.f vs combined_fortran.f\n\n')
        f.write(f'Total in datcom.f: {len(subs_datcom)}\n')
        f.write(f'Total in combined_fortran.f: {len(subs_combined)}\n')
        f.write(f'Common: {len(common)}\n')
        f.write(f'Only in datcom.f: {len(only_datcom)}\n')
        f.write(f'Only in combined_fortran.f: {len(only_combined)}\n\n')
        
        f.write('| Subroutine Name | In datcom.f | In combined_fortran.f |\n')
        f.write('|-----------------|-------------|-------------------------|\n')
        
        all_subs = sorted(set(subs_datcom) | set(subs_combined))
        for sub in all_subs:
            in_datcom = 'Yes' if sub in subs_datcom else 'No'
            in_combined = 'Yes' if sub in subs_combined else 'No'
            f.write(f'| {sub} | {in_datcom} | {in_combined} |\n')
else:
    print('Files not found')