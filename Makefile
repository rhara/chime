default:
	@echo Please specify target
	@for t in $$(cat Makefile | grep '^\w' | cut -d: -f1) ; do echo '-' $$t ; done

# test:
# 	python test.py sample/1i7g_ligand.sdf sample/1i7g_protein.pdb.gz --verbose

aconv_bulk:
	tools/atomic_conv_bulk.py v2015
	tools/atomic_conv_collect.py
