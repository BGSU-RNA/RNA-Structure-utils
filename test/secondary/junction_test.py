# test cases for junctions

from rnastructure.secondary import dot_bracket as Dot

def show_indices_loops(indices,loops):
    for loop_type in indices.keys():
        for i in range(0,len(indices[loop_type])):
            print(loop_type,loops[loop_type][i],indices[loop_type][i])
    print("")

print("Internal loop and hairpin and external:")
print('01234567890123456')
dot_string = '.((...(..)....)).'
seq        = 'CAGAAACAAGAAAACUC'
print(dot_string)
print(seq)
parser = Dot.Parser(dot_string)
indices = parser.indices(flanking=True)
loops = parser.loops(seq, flanking=True)
show_indices_loops(indices,loops)

print("Bulged loop and exterior:")
print('0123456789012')
dot_string = '((...(..))).'
seq        = 'AGAAACAAGCUC'
print(dot_string)
print(seq)
parser = Dot.Parser(dot_string)
indices = parser.indices(flanking=True)
loops = parser.loops(seq, flanking=True)
show_indices_loops(indices,loops)

print("J3 loop with one empty strand interior and empty hairpin interior:")
print('0123456789012345')
dot_string = '.(...()(..)...)'
seq        = 'CGAAACGCAAGAAAC'
print(dot_string)
print(seq)
parser = Dot.Parser(dot_string)
indices = parser.indices(flanking=True)
loops = parser.loops(seq, flanking=True)
show_indices_loops(indices,loops)

print("J3 loop with longer first stem:")
print(       '0123456789012345678')
dot_string = '((..(..)...(..)..))'
seq        = 'GGAACAAGAAAGAACAACC'
print(dot_string)
print(seq)
parser = Dot.Parser(dot_string)
indices = parser.indices(flanking=True)
loops = parser.loops(seq, flanking=True)
show_indices_loops(indices,loops)

print("J3 with no closing stem:")
print(       '012345678901234')
dot_string = '..(..)...(..)..'
seq        = 'AACAAGAAAGAACAA'
print(dot_string)
print(seq)
parser = Dot.Parser(dot_string)
indices = parser.indices(flanking=True)
loops = parser.loops(seq, flanking=True)
show_indices_loops(indices,loops)

print("J4 loop:")
print(       '012345678901234567890')
dot_string = '(...(..)...(..)(..).)'
seq        = 'CAAACAAGAAACAAGGAACAG'
print(dot_string)
print(seq)
parser = Dot.Parser(dot_string)
indices = parser.indices(flanking=True)
loops = parser.loops(seq, flanking=True)
show_indices_loops(indices,loops)

print("J4 loop with all empty interiors:")
print(       '012345678901234567890')
dot_string = '((..)(..)(..))'
seq        = 'CCAAGCAAGGAACG'
print(dot_string)
print(seq)
parser = Dot.Parser(dot_string)
indices = parser.indices(flanking=True)
loops = parser.loops(seq, flanking=True)
show_indices_loops(indices,loops)

print("J4 loop from RF00162")
dot_string = "(-(((((((----(-(((-(---------------------(((--------)))---------)))-))((((-((-(-(((---------------------------------------------------------------)----)))-----))))))----(-------((((((--------------------------------))))-))----)-)))))))-)--"
print(dot_string)
parser = Dot.Parser(dot_string)
indices = parser.indices(flanking=True)
for key,value in sorted(indices.items()):
    print(key,value)

seq        = "U-UUCUAUCCAGAG-AGG-U-GG-----------------AGGGA--CUGG-CCCUA--UGAA-ACC-UCGGCAACA-------UU-----------------------------------------------------------AU------------UGUGCCA-AUUC--CAG-CAAGC-------GCUA-----------------------GCU-UG-A-AA-GAUAGGA-A--"
loops = parser.loops(seq, flanking=True)
print(dot_string)
show_indices_loops(indices,loops)
seq        = "A-CCUUAUUUUGAG-AAG-C-UG-----------------AGGGA-UUUGG-CCCAU--AGAA-GCU-UCAGCAACC-G-ACU-UUA---------------------------------------------------------AAU----AGC-AC--GGUGCUA-AUAC--CAA-CGAG--------CAA-------------------------CU-CG-A-AU-GAUAAGU-A--"
loops = parser.loops(seq, flanking=True)
print(dot_string)
show_indices_loops(indices,loops)
seq        = "A-ACUUAUCAAGAG-CGG-C-UG-----------------AGGGA--CUGG-ACCUA--UGAA-GCC--CGGCAACC-U-GCA-UAG---------------------------------------------------------UUU----GUA-A---GGUGCUA-CUUC--CAG-CAAAAUG-----AAUUC--------------------CAUUU-UG-A-AA-GAUAAGG-G--"
loops = parser.loops(seq, flanking=True)
print(dot_string)
show_indices_loops(indices,loops)
seq        = "C-UCUUAUCGAGAG-CGG-C-AG-----------------AGGGA--CUGG-CCCGA--UGAA-GCC--CGGCAACC-U-AAC-UUUAUUUAA-----------------------------------------------GCGUAAA----GUG-AA--GGUGCUA-AUUC--CAG-CAAAAUGG---UGUAUU-------------------CCGUUU-UG-G-UA-GAUAAGA-G--"
loops = parser.loops(seq, flanking=True)
print(dot_string)
show_indices_loops(indices,loops)

