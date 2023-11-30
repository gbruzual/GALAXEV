install:
#	Use to install public files in GALAXEV subdirectories
	cd $(GALAXEV) ; \
	tar -xvzf ./pygalaxev.tgz ; \
        tar -xvzf ./pylibtar.tgz ; \
        tar -xvzf ./pyauxtar.tgz ; \
        \rm -f ./pylibtar.tgz ./pyauxtar.tgz

update:
#	Use to update public files in GALAXEV subdirectories
	tar -xvzf ./pyglxupd.tgz
