"""
Format is like
{"errors":[],"info":{"apiVersion":1,"dataVersion":1678543586,"deprecationDate":null,"deprecationInfo":null,"acknowledgement":null},"data":[{"insertion":"ins_19398:N","count":1},{"insertion":"ins_22204:GRGCHNRAN","count":1},{"insertion":"ins_28154:A","count":1},{"insertion":"ins_28811:AATTT","count":1},{"insertion":"ins_29557:C","count":1},{"insertion":"ins_29557:A","count":1},{"insertion":"ins_22204:GAGCNRGAR","count":1},{"insertion":"ins_27315:GTGGCCGA","count":1},{"insertion":"ins_19920:G","count":1},{"insertion":"ins_27405:TT","count":5},{"insertion":"ins_22204:GAGMSAGAA","count":1},{"insertion":"ins_22778:GTGTTGATGTT","count":
"""

import json

import argparse

parser = argparse.ArgumentParser(description='Process insertions')

parser.add_argument('--input', dest='input', required=True, help='Input file')

parser.add_argument('--output_csv', dest='output_csv', required=True, help='Output file')

parser.add_argument('--output_fasta', dest='output_fasta', required=True, help='Output file')

args = parser.parse_args()

with open(args.input, 'rt') as f:

    data = json.load(f)

# get the data

data = data['data']



# output as CSV
import requests

def get_epi_isl(insertion):
    url = f"https://lapis.cov-spectrum.org/gisaid/v1/sample/gisaid-epi-isl?country=United+Kingdom&dateFrom=2020-01-06&dateTo=2023-03-07&nucInsertions={insertion}&host=Human&accessKey=9Cb3CqmrFnVjO3XCxQLO6gUnKPd&orderBy=random"
    req = requests.get(url)
    print(req.text)
    data_text = req.text
    # strip
    data_text = data_text.strip()
    return data_text

counter = 0
def process_one_insertion(in2):
    global counter
    counter = counter + 1
   

    insertion = in2['insertion']
    count = in2['count']
    # get the position and the insertion sequence
    pos, seq = insertion.split(':')
    # get the position
    pos = int(pos.split('_')[1])
   
        
    return f"i_{counter}",pos, seq, count, insertion



def filter_fn(insertion):
    # if seq contains an ambiguous base, filter it out
    ambiguous_bases = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
    seq = insertion[2]
    for base in ambiguous_bases:
        if base in seq:
            return False
        
    # if length < 5, filter it out
    if len(seq) < 9:

        return False
    return True

processed = map(process_one_insertion, data)

filtered = list(filter(filter_fn, processed))

def add_episeq(my_tuple):
    epi_isl = ""
    # if count is 1
    if my_tuple[3] == 1:
        epi_isl = get_epi_isl(my_tuple[4])
    return my_tuple + (epi_isl,)


filtered = list(map(add_episeq, filtered))


with open(args.output_csv, 'wt') as f:

    f.write('insertion_id,pos,seq,count,epi_isl\n')

    for insertion in filtered:

        f.write(','.join(map(str, insertion)) + '\n')

# output as FASTA

with open(args.output_fasta, 'wt') as f:
    
        for insertion in filtered:

    
            f.write(f">{insertion[0]}" + '\n')

            f.write(insertion[2] + '\n')

            
    
            


            