
# ACLOCAL_AMFLAGS=-I m4

# sub directories to check for Makefiles
SUBDIRS = src perl tests .

################################################################
# doxygen
################################################################

if enable_doxygen
## doxygen documentation make target
include $(top_srcdir)/doc/doxygen.am
else
DX_PHONY =
endif

################################################################
# tests
################################################################

# global test directive # ensure all is compiled (need object files for tests)
.PHONY: DX_PHONY tests test
tests: all
	@echo
	@echo " => compiling and running tests ... will take a while ..."
	@echo
	@$(MAKE) check -s -C tests

test: tests

