import time

def get_mass_table():
    mass_table = {}
    with open('/Users/richard/Downloads/integer_mass_table.txt') as file:
        for _ in range(20):
            amino,n = file.readline().split()
            mass_table[amino]=int(n)
    return mass_table
    
mass_table = get_mass_table()

def get_prefixMass(peptide):
    res = [0]
    base = 0
    for i in range(len(peptide)):
        base += mass_table[peptide[i]]
        res.append(base)
    return res

def cyclic_spectrum(PrefixMass):
    res = [0]
    for i in range(len(PrefixMass)-1):
        for j in range(i+1,len(PrefixMass)):
            temp = PrefixMass[j]-PrefixMass[i]
            temp1 = PrefixMass[-1]-temp
            res.append(temp)
            ### eliminate double count
            if i>0 and j<len(PrefixMass)-1: res.append(temp1)
    return sorted(res)
    
def get_cyclospectrum(peptide): 
    return cyclic_spectrum(get_prefixMass(peptide))

def convert_to_num_string(aa_string): 
    return '-'.join([str(mass_table[letter]) for letter in aa_string])

def CYCLOPEPTIDESEQUENCING(spectrum,mass_table):
    res = []
    peptides = {aa:mass_table[aa] for aa in mass_table if mass_table[aa] in spectrum}
    while peptides:
        peptides = expand(peptides,mass_table)
        peptides_copy = {pep:peptides[pep] for pep in peptides}
        for peptide in peptides_copy:
            if calc_mass(peptide) in spectrum:
                if set(get_cyclospectrum(peptide)) == spectrum:
                    res.append(peptide)
                    peptides.pop(peptide)
            else:
                peptides.pop(peptide)
    return set(map(convert_to_num_string,res))

def expand(peptides,mass_table):
    new_peptides = {}
    for aa in peptides:
        for aa0 in mass_table:
            new_peptides[aa+aa0] = peptides[aa]+mass_table[aa0]
    return new_peptides

def calc_mass(peptide):
    total = 0
    for letter in peptide:
        total += mass_table[letter]
    return total
       
def read_spectrum(filename):
    with open(filename) as file:
        return map(int,file.readline().split())
        
spectrum = set(read_spectrum("/Users/richard/Downloads/test.txt"))
start = time.time()
res = CYCLOPEPTIDESEQUENCING(spectrum,mass_table)
end = time.time()
print ' '.join(list(res))
print end-start