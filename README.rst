=========
RNAlign2D
=========
RNAlign2D is a Python wrapper allowing one to align multiple RNA molecules using
information about secondary structure and sequence.
For this purpose MUSCLE program is used, but one can adjust it
to use any other multiple sequence alignment tool.

This tool firstly remove common RNA modifications, then convert RNA sequences
and secondary structures to pseudo-amino acid sequences (please not mistake for
translation process - that is only technical solution), and use MUSCLE to
align pseudo-amino acid sequences. Finally process is reverted and alignment
with original sequence and secondary structure is restored.

==================
Dot bracket format
==================
Format, where paired residues as marked as brackets () and unpaired residues
are marked as dots.
For example simple structure containing hairpin
(stem of 5 residues length, loop of 6 residues length)
may look like:
........(((((......))))).........

For pseudoknots [] brackets are used. For higher level pseudoknots {}, <>
or even Aa, Bb, Cc ... pairs are used. For the RNAlign2D maximum pseudoknots
level taken into account is Bb.

======================
Input and output files
======================
Input and output files are modified FASTA format.
For each sequence in the file:
* First line contain > character and then optional name
* Second line (or multiple lines) contain sequence
* Third line (or multiple lines) contain structure
Sequence and structure are distinguished on the basis of dot and bracket
content. If there is no secondary structure, linear - unpaired structure
is assumed.
If there is Vienna RNA package installed, this package would automatically
generate missing secondary structures.

Sequence may contain modified nucleotides - these modifications
of nucleotide residues are removed (for example B - 2′-O-methylcytidine
is converted to C - cytidine), and after the alignment they are restored back
to the original state. RNAlign2D uses modifications abbreviations from
MODOMICS database.
Secondary structure may contain higher level pseudoknots, but for the *simple*
mode only fist level - [] pseudoknots are treated differently, and for the
*pseudo* mode only pseudoknots up to level Bb.

=====
Modes
=====

- Mode *simple*
- Mode *pseudo*

Mode *simple* is used in most cases, when there are no pseudoknots or
there are only first level pseudoknots. It uses both sequence and secondary
structure into account. Higher level pseudoknots are treated as
unpaired regions.
Mode *pseudo* on the other hand uses only secondary structure, as it allows to
align higher level pseudoknots. In some cases mode simple may give better
results than mode pseudo even for sequences.
In case of the mode *pseudo* there is a need to explicitly show matrix dedicated
for *pseudo* alignment, otherwise results may not reflect anything.

=======
Options
=======

- *-i* option - input file
- *-o* option - output file
- *-matrix* option - file with custom matrix or default matrix for pseudoknots (data/pseudo_matrix in the project)
- *-mode* option - allow to select one of the modes: *simple* or *pseudo*, default *simple*, please remember to use proper matrix for given mode
- *-gapopen* option - allow to select gap opening penalty for sequence alignment
- *-gapextend* option - allow to select gap extending penalty for sequence alignment

example usage:

``rnalign2d -i my_input -o my_output``

=============
Custom matrix
=============
One can define matrix for himself. As it is simply text file, there is nothing
that cannot be done manually. Hovewer there is script create_matrix, that can
facilitate that process.

* Most probably the strongest score should be given to the same type of
  brackets like ( and ( or [ and [, default value is 7, option: *-same*
* Probably the smallest score should be given to the oposite brackets like
  ( and ), default value is -10, option *-reverse*
* Sometimes there is a bulge in the RNA structure or the loop is bigger by one
  pair of residues, therefore situation like ( and . should be given small
  negative value, default value is -1, option *-bracket-dot*
* In the situation that in one structure there is pseudoknot and in the other
  there is classical Watson-Crick base pair for example ( and [, there should
  be given small positive or negative score, default value is 2, option *-other*
* In case of two unpaired regions . and . there should be positive score, but
  we suggest it to be small, as loops are probably the best place to introduce
  gaps, default value is 3, option *-dot_dot*
* For the *simple* mode there is also possibility to add score for matching
  nucleotide residues in the sequence. Depending on how important the
  sequence is, this value may vary, default value is 1, option *-seq_match_add*
* In case of *pseudo* mode you need to define mode (options are *simple* and
  *pseudo*), default value is *simple*, option *-mode*
* Output name could be defined using *-o* option

exaple usage:

``create_matrix -o my_result -seq_match_add 8``

====================
Remove modifications
====================
One can remove modifications from file for this there is script rm_mod
There are options *-i* to define input file and *-o* to define output file

exaple usage:

``rm_mod -o my_result -i my_input``

============
REQUIREMENTS
============
* Python 3 (tested on Python 3.5)
* MUSCLE (tested on version 3.8.31)
* pytest (tested on version 5.1.3)
* Vienna RNA (optional, tested on version 2.4.14)

=============
EXAMPLE FILES
=============
In the directory *example* there are 3 files:

- dot_bracket.txt - file containing example input sequence
- dot_bracket2.txt - file containing example output
- dot_bracket_cls.txt - file containing only sequences

Files are created using data from T-psi-C-database_tpsic.igcz.poznan.pl

=====
OTHER
=====
In case of WARNING:
``*** WARNING *** Matrix is not symmetrical, �->�=-10, �->�=0``
most probably there is no problem at all,
but if it is your custom matrix, you can check it to ensure it contain desired
data

To test this software Vienna RNA is required, otherwise one test would fail.

========
CITATION
========
If you are using our software in your research - cite us:
Tomasz Woźniak, Małgorzata Sajek, Jadwiga Jaruzelska, Marcin Piotr Sajek; RNAlign2D: a rapid method for combined RNA structure and sequence-based alignment using a pseudo-amino acid substitution matrix; BMC Bioinformatics. 2021 Oct 16;22(1):504. doi: 10.1186/s12859-021-04426-8.
