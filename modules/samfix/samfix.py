import fileinput
import re
import sys

def printe(x):
    print(x, file=sys.stderr)

qname_regex = re.compile('[!-?A-~]{1,254}')
cigar_regex = re.compile('\*|([0-9]+[MIDNSHPX=])+')

# also check H op can only be first or last op
cigar_h_op_regex = re.compile('\*|[0-9]+[MIDNSHPX=][0-9MIDNSPX=]*[0-9MIDNSHPX=]?')
# s can only have h ops between them and the end of the cigar
cigar_s_op_regex = re.compile('\*|[0-9MIDNHPX=]*[0-9S]*[0-9H]?')


seq_regex = re.compile('\*|[A-Za-z=.]+')

tab = '\t'

cigar_seq_op_regex = re.compile('[0-9]+[MISX=]')

def is_valid_cigar(cigar):
    if cigar_regex.fullmatch(cigar) is None:
        return False
    if cigar_h_op_regex.fullmatch(cigar) is None:
        return False
    if cigar_s_op_regex.fullmatch(cigar) is None:
        return False

    return True

def seq_len_matches_cigar(seq, cigar):
    if seq_regex.fullmatch(seq) is None:
        return False

    if cigar == '*' or seq == '*':
        return True

    cigar_ops = cigar_seq_op_regex.findall(cigar)
    
    op_sum = 0
    for op in cigar_ops:
        op_sum += int(op[:-1])

    if len(seq) != op_sum:
        return False

    return True



def fix_line(line, idx):
    # use qname to detect alignment section
    # headers start with an '@' which cant be in qname
    fields = line.split(tab)

    if qname_regex.fullmatch(fields[0]) is not None:
        # process alignment section

        # sam has 11 or more tab separated fields
        new_fields = fields

        assert len(fields) >= 11
        cigar = fields[5]
        seq = fields[9]

        if not is_valid_cigar(cigar):
            cigar = '*'
            new_fields[5] = cigar

            return ''
        
        if not seq_len_matches_cigar(seq, cigar):
            raise RuntimeError(f'invalid_seq_len in line-no: {idx} with columns: {line}')
            #cigar = '*'
            #new_fields[5] = cigar


        new_line = tab.join(fields)

    else:
        # just return unmodified header
        new_line = line

    return new_line


def main():
    cigar_placeholder_symbol='*'

    idx = 0

    for line in fileinput.input():
        newline = fix_line(line, idx)
        print(newline, end='')
        idx += 1

if __name__ == '__main__':
    main()
