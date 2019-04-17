default:
	@echo Please specify target
	@for t in $$(cat Makefile | grep '^\w' | cut -d: -f1) ; do echo '-' $$t ; done

test:
	python test.py sample/1i7g_ligand.sdf sample/1i7g_protein.pdb.gz --verbose

aconv:
	python atomic_conv.py sample/1i7g_ligand.sdf sample/1i7g_protein.pdb.gz 1i7g_out.npz
