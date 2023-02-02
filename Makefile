CXX ?= g++
BOOST_INSTALL_PATH ?= $(BOOST_ROOT)
OBABEL_INSTALL_PATH ?= $(OPEN_BABEL_ROOT)
CXXFLAGS = -std=c++11
ifeq ($(Debug),Y)
CXXFLAGS += -g -O0
else
CXXFLAGS += -O2
endif

ifeq ($(Symbols),Y)
CXXFLAGS += -g
endif

ifeq ($(OpenMP),Y)
OMPFLAG = -fopenmp
else
OMPFLAG =
endif

CXXFLAGS += $(OMPFLAG)

ifeq ($(BOOST_INSTALL_PATH),)
BOOSTIP =
BOOSTLP =
else
BOOSTIP = -I$(BOOST_INSTALL_PATH)/include
BOOSTLP = -L$(BOOST_INSTALL_PATH)/lib
endif
OBABELIP = -I$(OBABEL_INSTALL_PATH)/include/openbabel-2.0
OBABELLP = -L$(OBABEL_INSTALL_PATH)/lib

.SUFFIXES: .cc .o

_GRID_OBJS = grid_main.o utils.o infile_reader.o Molecule.o Vector3d.o Atom.o AtomEnergyGrid.o EnergyGrid.o EnergyCalculator.o log_writer_stream.o OBMol.o
_CONFDOCK_OBJS = fraggrid_main.o Vector3d.o EnergyGrid.o Molecule.o Fragment.o EnergyGrid.o EnergyCalculator.o infile_reader.o utils.o MoleculeToFragments.o AtomEnergyGrid.o FragmentEnergyGrid.o Atom.o log_writer_stream.o UnionFindTree.o OBMol.o CalcMCFP.o Optimizer.o
_ATOMDOCK_OBJS = nofrag_main.o Vector3d.o EnergyGrid.o Molecule.o EnergyGrid.o EnergyCalculator.o infile_reader.o utils.o AtomEnergyGrid.o Atom.o log_writer_stream.o OBMol.o Optimizer.o
_TEST_OBJS = Atom.o EnergyCalculator.o Molecule.o utils.o Vector3d.o TestEnergyCalculator.o TestMain.o

_EASYTEST_OBJS = easytest_main.o Vector3d.o Molecule.o EnergyGrid.o EnergyCalculator.o infile_reader.o utils.o Atom.o log_writer_stream.o OBMol.o

_SCOREONLY_OBJS = score_only_main.o Vector3d.o Molecule.o EnergyGrid.o EnergyCalculator.o infile_reader.o utils.o Atom.o log_writer_stream.o OBMol.o
_INTRAENERGY_OBJS = intraenergy_main.o Vector3d.o Molecule.o EnergyGrid.o EnergyCalculator.o infile_reader.o utils.o Atom.o log_writer_stream.o OBMol.o

# ALL = objs atomgrid-gen conformer-docking atom-docking unittest easytest-docking score-only intraenergy-only
ALL = atomgrid-gen conformer-docking

GRID_OBJS = $(patsubst %,objs/%,$(_GRID_OBJS))
CONFDOCK_OBJS = $(patsubst %,objs/%,$(_CONFDOCK_OBJS))
ATOMDOCK_OBJS = $(patsubst %,objs/%,$(_ATOMDOCK_OBJS))
TEST_OBJS = $(patsubst %,objs/%,$(_TEST_OBJS))
EASYTEST_OBJS = $(patsubst %,objs/%,$(_EASYTEST_OBJS))
SCOREONLY_OBJS = $(patsubst %,objs/%,$(_SCOREONLY_OBJS))
INTRAENERGY_OBJS = $(patsubst %,objs/%,$(_INTRAENERGY_OBJS))


all: objs $(ALL)
	rm -r objs

objs:
	mkdir -p objs
	mkdir -p objs/test

atomgrid-gen: $(GRID_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(BOOSTLP) $(OBABELLP) -lboost_regex -lboost_program_options -lopenbabel

conformer-docking: $(CONFDOCK_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(BOOSTLP) $(OBABELLP) -lboost_regex -lboost_program_options -lopenbabel

atom-docking: $(ATOMDOCK_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(BOOSTLP) $(OBABELLP) -lboost_regex -lboost_program_options -lopenbabel

unittest: $(TEST_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(BOOSTLP) $(OBABELLP) -lboost_regex -lboost_unit_test_framework -lopenbabel

easytest-docking: $(EASYTEST_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(BOOSTLP) $(OBABELLP) -lboost_regex -lboost_program_options -lopenbabel

score-only: $(SCOREONLY_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(BOOSTLP) $(OBABELLP) -lboost_regex -lboost_program_options -lopenbabel

intraenergy-only: $(INTRAENERGY_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(BOOSTLP) $(OBABELLP) -lboost_regex -lboost_program_options -lopenbabel


objs/%.o: src/%.cc src/%.hpp src/common.hpp
	$(CXX) $(CXXFLAGS) $(BOOSTIP) $(OBABELIP) -o $@ -c $<

objs/%.o: src/%.cc src/common.hpp
	$(CXX) $(CXXFLAGS) $(BOOSTIP) $(OBABELIP) -o $@ -c $<

clean:
	rm -rf objs $(ALL)

