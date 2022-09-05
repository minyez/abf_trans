targets = bccHe_k444 \
	bccHe_k333_abacus \
	hbn_k331 hbn_k662 \
	lif_k333 lif_k444_abacus \
	nacl_k333 nacl_k444 \
	scHe_k333_abacus scHe_k444


.PHONY: default $(targets)

default: $(targets)

bccHe_k444:
	cd $@; ./run_coulomb.sh

bccHe_k333_abacus:
	cd $@; ./run_coulomb.sh

hbn_k331:
	cd $@; ./run_coulomb.sh

hbn_k662:
	cd $@; ./run_coulomb.sh

nacl_k333:
	cd $@; ./run_coulomb.sh

nacl_k444:
	cd $@; ./run_coulomb.sh

scHe_k333_abacus:
	cd $@; ./run_coulomb.sh

scHe_k444:
	cd $@; ./run_coulomb.sh

lif_k333:
	cd $@; ./run_coulomb.sh

lif_k444_abacus:
	cd $@; ./run_coulomb.sh
	cd $@; ./run_hamiltonian.sh
