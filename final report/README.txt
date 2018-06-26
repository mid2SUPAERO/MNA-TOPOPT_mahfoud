ISAE-SUPAERO template for student research projects


INTRODUCTION
------------
This LaTeX report template has been designed in order to help students 
writing their research project reports.


REQUIREMENTS
------------
 A full LaTeX distribution such as TexLive (>2012) or MikTex (>2012).
 An UTF-8 aware text editor such as TexMaker, TexStudio or Emacs.
 Optionally: a bibliography manager such as JabRef can be used in order
to edit the .bib file.


USAGE
-----
 Edit the {.tex,.bib} files using UTF-8 encoding.
 Put your figures in the "images" subdirectory.

The compilation process is the following:
1) Compile the main file ".tex" using PDFLatex.
2) Compile the ".aux" using BibTeX.
3) Compile twice the main file ".tex" using PDFLatex.


FILES MANIFEST
--------------
 "bibliography-report.tex": main LaTeX file of the report.
 "bibliography-report.pdf": main output file of the report. 
 "IEEEtran.cls": document class of the IEEE
 "images": subdirectory used to store the figures.
 "README.txt": the current file.


EXTERNAL DOCUMENTATION
----------------------
Discovering the language with "A not so short introduction to LaTeX":
<http://ctan.mines-albi.fr/info/lshort/english/lshort.pdf>

Doing you own figures in LaTeX using "TikZ for the impatient":
<http://math.et.info.free.fr/TikZ/bdd/TikZ-Impatient.pdf>

Getting inspired by already on-the-shelf figures via "TikZ examples":
<http://www.texample.net/tikz/examples/>

Doing math plotting using "PGFPlots":
<http://pgfplots.sourceforge.net/pgfplots.pdf>

Exporting you Matlab plots thanks to "Matlab2TikZ":
<https://github.com/nschloe/matlab2tikz>


CREDITS
-------
This template has been created by Damien Roque (ISAE-SUPAERO)
based on Michael Shell initial template.
See http://personnel.isae.fr/damien-roque


CHANGELOG
---------
2016/04/14 v1.0 Initial version.
 

