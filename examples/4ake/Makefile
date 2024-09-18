test:
	( time ../../scripts/makebloc.pl 4ake.pdb ) > pdb
	( time ../../src/build/rtb2 ) >& rtb.out

clean:
	-rm mov*.mod* matrice.sdij RT.mat resi_lg pdb *.out eigenvec.dat
