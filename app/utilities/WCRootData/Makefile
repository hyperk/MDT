CXX=g++
LD=g++

CXXFLAGS += -Wall -O2 -std=c++11 -g -fPIC -I./include $(shell root-config --cflags) -I$(WCSIM_BUILD_DIR)/include/WCSim -I${MDTROOT}/cpp/include
#LDFLAGS += -shared $(shell root-config --ldflags) $(shell root-config --libs) -lTreePlayer -lMinuit2 -L$(WCSIMDIR) -lWCSimRoot -L${MDTROOT}/cpp -lMDT
LDFLAGS += -shared $(shell root-config --ldflags) $(shell root-config --libs) -L$(WCSIM_BUILD_DIR)/lib -lWCSimRoot -L${MDTROOT}/cpp -lMDT

TARGET=libWCRData.so

SRCDIR = ./src
OBJDIR = ./obj

SRCS = $(wildcard $(SRCDIR)/*.cc)
OBJCTS = $(addprefix $(OBJDIR)/, $(notdir $(SRCS:.cc=.o)))


.PHONY: clean all

all: $(OBJCTS)
	$(LD) $(LDFLAGS) $(OBJCTS) -o $(TARGET)
	@echo "$(TARGET) done"

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	@[ -d $(OBJDIR) ]
	$(LD) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf $(OBJDIR)/*.o $(TARGET)
