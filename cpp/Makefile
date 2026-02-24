CXX=g++
LD=g++

CXXFLAGS += -Wall -O2 -std=c++11 -g -fPIC -I./include $(shell root-config --cflags)
LDFLAGS += -shared $(shell root-config --ldflags) $(shell root-config --libs)

TARGET=libMDT.so
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
