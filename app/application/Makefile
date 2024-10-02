CXX=g++
LD=g++

CXXFLAGS += -Wall -std=c++11 -g $(shell root-config --cflags)\
			-I$(WCSIM_BUILD_DIR)/include/WCSim\
			-I${MDTROOT}/cpp/include\
			-I${WCRDROOT}/include

LDFLAGS += $(shell root-config --ldflags) $(shell root-config --libs) -lTreePlayer\
		   -L$(WCSIM_BUILD_DIR)/lib -lWCSimRoot\
		   -L${MDTROOT}/cpp -lMDT\
		   -L${WCRDROOT} -lWCRData

.PHONY: clean Execs

appWCTESingleEvent: appWCTESingleEvent.o PMTResponse3inchR12199_02.o
	$(RM) $@
	$(LD) $^ $(LDFLAGS) -o $@
	@echo "$@ done"

appIWCDSingleEvent: appIWCDSingleEvent.o PMTResponse3inchR12199_02.o
	$(RM) $@
	$(LD) $^ $(LDFLAGS) -o $@
	@echo "$@ done"

appGenPileUpSpill: appGenPileUpSpill.o PMTResponse3inchR12199_02.o
	$(RM) $@
	$(LD) $^ $(LDFLAGS) -o $@
	@echo "$@ done"

%o::%cc
	$(CXX) $(CXXFLAGS) -c $? -o $@

all: appWCTESingleEvent appIWCDSingleEvent appGenPileUpSpill

clean:
	$(RM) *.o
	$(RM) -f appWCTESingleEvent
	$(RM) -f appIWCDSingleEvent
	$(RM) -f appGenPileUpSpill
