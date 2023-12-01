from samfix import is_valid_cigar

valid_cigar_examples = '''
3M2S3H
3H2S3H
3M2I3M
'''

invalid_cigar_examples = '''
3S2P3M
3H2H3M
3S2H3M
'''


for cigar in valid_cigar_examples.strip().split():
    print(cigar)
    c = cigar.strip()
    assert is_valid_cigar(c), str(c)
    
for cigar in invalid_cigar_examples.strip().split():
    print(cigar)
    c = cigar.strip()
    assert not is_valid_cigar(c), str(c)
