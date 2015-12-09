from __future__ import print_function
from Bio.SeqRecord import SeqRecord
from Bio.Seq import reverse_complement

class OrderedSeqRecord(SeqRecord):
    def __init__(self, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], SeqRecord):
            other = args[0]
            self._per_letter_annotations = other._per_letter_annotations
            self.seq = other.seq
            self.id = other.id
            self.name = other.name
            self.description = other.description
            self.dbxrefs = other.dbxrefs
            self.features = other.features
            self.annotations = other.annotations
            self.letter_annotations = other.annotations

        else:
            super().__init__(*args, **kwargs)
    def __getitem__(self, index):
        if isinstance(index, slice) and index.start is not None and index.stop is not None and index.start > index.stop:
            index = slice(index.stop, index.start, index.step)
            retval = super().__getitem__(index)
            retval.seq = reverse_complement(retval.seq)
            return retval
        else:
            return super().__getitem__(index)
