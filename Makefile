#########################################################
#                                                       #
# Ladder Operator Vibrational Configuration Interaction #
#                                                       #
#########################################################

### Compiler settings
CXX=g++ -std=c++14
CXXFLAGS=-static-libgcc -static-libstdc++ -fopenmp -O3
LDFLAGS=-I./src/ -I/usr/include/eigen3/ -I/usr/include/boost/ -I/usr/include/Spectra

#########################################################

### Compile rules for users and devs

install:	title vhcibin manual compdone

clean:	title delbin compdone

#########################################################

### Rules for building various parts of the code

vhcibin:	
	@echo ""; \
	echo "### Compiling the VHCI binary ###"
	$(CXX) $(CXXFLAGS) ./src/VHCI.cpp -o vhci $(LDFLAGS)

manual:	
	@echo ""; \
	echo "### Creating the manual ###"; \
	cp README.md doc/VHCI_manual.txt; \
	echo " [Complete]"

compdone:	
	@echo ""; \
	echo "Done."; \
	echo ""

title:	
	@echo ""; \
	echo "#########################################################"; \
	echo "#                                                       #"; \
	echo "#    Vibrational Heat-Bath Configuration Interaction    #"; \
	echo "#                                                       #"; \
	echo "#########################################################"; \
	echo ""

delbin:	
	@echo ""; \
	echo '     ___'; \
	echo '    |_  |'; \
	echo '      \ \'; \
	echo '      |\ \'; \
	echo '      | \ \'; \
	echo '      \  \ \'; \
	echo '       \  \ \'; \
	echo '        \  \ \       <wrrr vroooom wrrr> '; \
	echo '         \__\ \________'; \
	echo '             |_________\'; \
	echo '             |__________|  ..,  ,.,. .,.,, ,..'; \
	echo ""; \
	echo ""; \
	echo "Removing binaries and manual..."; \
	rm -f ./vhci doc/VHCI_manual.txt
