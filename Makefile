myutil=/home/james.thorpe/bin/myUtils/lib
FC=gfortran
#FC=ifort
FFALGS= -fcheck=bounds -O3 -frecursive

BIN_DIR=./bin
SRC_DIR=./src

.PHONY: all

all: 
	@set -e; \
	for i in $(SRC_DIR); do \
		if [ -d $$i ]; then \
			if [ -f $$i/Makefile ]; then \
				$(MAKE) -C $$i all ;\
			fi; \
		fi; \
	done;
	if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi 
	mv $(SRC_DIR)/myCT $(BIN_DIR)

clean:
	rm -f $(SRC_DIR)/*.o myCT
