# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 22:30:29 2020

@author: liusm
"""


from Bio import SeqIO
#Seq object tiene string y alphabet attribute. Muchos string methods apply. (Translation es biological translation)
#Alphabet attribute te dice que secuencia estas viendo. Bio.Alphabet module para ver mas alphabets
#Ejemplo de input local file y print propiedades de seq file. 
for seq_record in SeqIO.parse("C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/S288C_reference_sequence.fsa", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

#Adding records in fasta file to list
list_1 = []
for seq_record in SeqIO.parse("C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/S288C_reference_sequence.fsa", "fasta"):
    list_1.append(seq_record)
#Accessing sequence information in records appended to list and turning to string give you sequence info. Remember
#Slicing no incluye final index    
str(list_1[0].seq[1000:1100])
#Rapid fasta_format......
# fasta_format_string = ">Name\n%s\n" % my_seq
# print(fasta_format_string)
  
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
my_prot = Seq("AGTACACTGGT", IUPAC.protein)
my_prot 
my_prot.alphabet


# for index, letter in enumerate(seq_record):
#     print("%i %s" % (index, letter))
    
    
seq_record.seq[0]
seq_record.seq.count("AAAAAA") #Non overlapping Count

#Seq objects no son mutable. Usa mutable object
#my_seq.tomutable()

#Adding sequences
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
list_of_seqs = [Seq("ACGT", generic_dna), Seq("AACC", generic_dna), Seq("GGTT", generic_dna)]
concatenated = Seq("", generic_dna)
for s in list_of_seqs:
    concatenated += s
sum(list_of_seqs, Seq("", generic_dna))


#Algunas funciones acceptan seq objects



# The SeqRecord (Sequence Record) class is defined in the Bio.SeqRecord module. This class allows higher
# level features such as identifiers and features to be associated with a sequence (see Chapter 3), and is the
# basic data type for the Bio.SeqIO sequence input/output interface (see Chapter 5).
# The SeqRecord class itself is quite simple, and offers the following information as attributes:
# .seq – The sequence itself, typically a Seq object.
# .id – The primary ID used to identify the sequence – a string. In most cases this is something like an
# accession number.
# .name – A “common” name/id for the sequence – a string. In some cases this will be the same as the
# accession number, but it could also be a clone name. I think of this as being analogous to the LOCUS
# id in a GenBank record.
# .description – A human readable description or expressive name for the sequence – a string.
# .letter annotations – Holds per-letter-annotations using a (restricted) dictionary of additional information
# about the letters in the sequence. The keys are the name of the information, and the information is
# contained in the value as a Python sequence (i.e. a list, tuple or string) with the same length as
# the sequence itself. This is often used for quality scores (e.g. Section 20.1.6) or secondary structure
# information (e.g. from Stockholm/PFAM alignment files).
# 34
# .annotations – A dictionary of additional information about the sequence. The keys are the name of
# the information, and the information is contained in the value. This allows the addition of more
# “unstructured” information to the sequence.
# .features – A list of SeqFeature objects with more structured information about the features on a sequence
# (e.g. position of genes on a genome, or domains on a protein sequence). The structure of sequence
# features is described below in Section 4.3.
# .dbxrefs - A list of database cross-references as strings.




# record_iterator = SeqIO.parse("ls_orchid.fasta", "fasta")
# first_record = next(record_iterator)
# print(first_record.id)
# print(first_record.description)
# records = list(SeqIO.parse("ls_orchid.gbk", "genbank"))
# print("Found %i records" % len(records))


# from Bio import SeqIO
#record_iterator = SeqIO.parse("ls_orchid.fasta", "fasta")
#first_record = next(record_iterator)
#first_record.id = "new_id"
#first_record.description = first_record.id + " " + "desired new description"
#print(first_record.format("fasta")[:200])


# from Bio import SeqIO
# with open("ls_orchid.gbk") as handle:
#     print(sum(len(r) for r in SeqIO.parse(handle, "gb")))


# from Bio import Entrez
# from Bio import SeqIO
# Entrez.email = "A.N.Other@example.com"
# with Entrez.efetch(
# db="nucleotide", rettype="gb", retmode="text", id="6273291"
# ) as handle:
#     seq_record = SeqIO.read(handle, "gb") # using "gb" as an alias for "genbank"
#     print("%s with %i features" % (seq_record.id, len(seq_record.features)))

# orchid_dict = SeqIO.to_dict(SeqIO.parse("ls_orchid.gbk", "genbank"))
# orchid_dict = SeqIO.index("ls_orchid.gbk", "genbank")

# import glob
# >>> from Bio import SeqIO
# >>> files = glob.glob("gbvrl*.seq")
# >>> print("%i files to index" % len(files))
# 4
# >>> gb_vrl = SeqIO.index_db("gbvrl.idx", files, "genbank")
# >>> print("%i sequences indexed" % len(gb_vrl))

####Writing to fasta!!!!!!
# rec3 = SeqRecord(
# Seq(
# "MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKEKFKRMC"
# "EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP"
# "KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVKRFMMYQQGCFAGGTVLRMAKDLAEN"
# 66
# "NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDGAAAVIIGSDPIPEVERPLFELV"
# "SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAFQPLGISDWNSLFW"
# "IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFILDEMRKASAKEGLGT"
# "TGEGLEWGVLFGFGPGLTVETVVLHSVAT",
# generic_protein,
# ),
# id="gi|13925890|gb|AAK49457.1|",
# description="chalcone synthase [Nicotiana tabacum]",
# )
# my_records = [rec1, rec2, rec3]

# SeqIO.write(my_records, "my_example.faa", "fasta")

# Bio.SeqIO.convert()

# records = [rec.reverse_complement(id="rc_"+rec.id, description = "reverse complement") \
# ... for rec in SeqIO.parse("ls_orchid.fasta", "fasta") if len(rec)<700]



# from Bio.Alphabet import generic_dna
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.Align import MultipleSeqAlignment
# align1 = MultipleSeqAlignment(
# [
# SeqRecord(Seq("ACTGCTAGCTAG", generic_dna), id="Alpha"),
# SeqRecord(Seq("ACT-CTAGCTAG", generic_dna), id="Beta"),
# SeqRecord(Seq("ACTGCTAGDTAG", generic_dna), id="Gamma"),
# ]
# )
# align2 = MultipleSeqAlignment(
# [
# SeqRecord(Seq("GTCAGC-AG", generic_dna), id="Delta"),
# SeqRecord(Seq("GACAGCTAG", generic_dna), id="Epsilon"),
# SeqRecord(Seq("GTCAGCTAG", generic_dna), id="Zeta"),
# ]
# )
# align3 = MultipleSeqAlignment(
# [
# SeqRecord(Seq("ACTAGTACAGCTG", generic_dna), id="Eta"),
# SeqRecord(Seq("ACTAGTACAGCT-", generic_dna), id="Theta"),
# SeqRecord(Seq("-CTACTACAGGTG", generic_dna), id="Iota"),
# ]
# )
# my_alignments = [align1, align2, align3]