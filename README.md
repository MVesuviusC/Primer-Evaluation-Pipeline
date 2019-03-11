# bsPrimerBlast
	and
		# bsPrimerTree

The point of this code is to make a locally available version of the primer specificity assessment portion of primerBlast. I've copied the algorithm, with a few minor changes, from the primerBlast paper and code. 


Dependencies:
 - Dendroscope
 - getTaxa.pl
 - Megacc - V10+
 - Blastn - not version 2.8.1
   - blast nt database
 - Mafft
 - Imagemagick - optional - to make jpg from svg