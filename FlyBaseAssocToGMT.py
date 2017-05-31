from __future__ import print_function
from sys import argv, stdout
from collections import defaultdict

def parse_obo(obo_file):
    children_of = defaultdict(list)
    current_id = ''
    for line in obo_file:
        if line.strip() == '':
            current_id = ''
        elif line.startswith('id'):
            current_id = line.strip().split()[-1]
        elif current_id and line.startswith('is_a'):
            parent = line.split()[1]
            children_of[parent].append(current_id)
        elif current_id and line.startswith('relationship'):
            parent = line.split()[2]
            children_of[parent].append(current_id)

    return children_of

def all_children_of(children_of, term):
    retval =  set(children_of[term])
    for child_term in children_of[term]:
        retval |= all_children_of(children_of, child_term)
    return retval

def all_terms_of(goterms, children_of, term):
    genes = set(goterms[term])
    for term in all_children_of(children_of, term):
        genes.update(set(goterms[term]))
    return genes


if __name__ == "__main__":
    if len(argv) > 2:
        out = open(argv[2], 'w')
    else:
        out = stdout
    children_of = parse_obo(open('prereqs/go-basic.obo'))
    goterms = defaultdict(list)
    pfcs = {}
    for line in open(argv[1]):
        if line.startswith('!'): continue
        data = line.split('\t')
        FBgn = data[1]
        symbol = data[2]
        GOterm = data[4]
        PFC = data[8] # Process/Function/Component
        goterms[GOterm].append(symbol)
        pfcs[GOterm] = PFC

    for GOterm in goterms.copy():
        goterm_genes = set(goterms[GOterm])
        for child_GOterm in all_children_of(children_of, GOterm):
            goterm_genes |= set(goterms[child_GOterm])
        print(GOterm, pfcs[GOterm], *goterm_genes, sep='\t', file=out)
    if out is not stdout:
        out.close()

