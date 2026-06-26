# ============================================
# KG-Y
# ============================================

FC = h5pfc

TARGET = kg_y.exe
OBJDIR = build
MODDIR = build

FFLAGS = -O3 -march=native -mtune=native -ffast-math -I$(MODDIR)
LIBS   = -llapack -lblas

NPROC  ?= 4
PARAMS ?= params.nml
RUN_DIR = $(basename $(notdir $(PARAMS)))

SRC = \
  main/mpi_lib.f90 \
  main/run_control.f90 \
  MoL/mesh.f90 \
  equations/geometry.f90 \
  numerical/utils.f90 \
  equations/state.f90 \
  numerical/finite_differences.f90 \
  equations/sl_spectrum.f90 \
  equations/initial_data.f90 \
  equations/rhs.f90 \
  MoL/cfl.f90 \
  MoL/time_integrators.f90 \
  numerical/plane_projection.f90 \
  main/hdf5_lib.f90 \
  equations/sl_projector.f90 \
  main/main.f90

OBJ = $(SRC:%.f90=$(OBJDIR)/%.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

$(OBJDIR)/%.o: %.f90
	@mkdir -p $(dir $@) $(MODDIR)
	$(FC) $(FFLAGS) -c $< -o $@

run: $(TARGET)
	@mkdir -p $(RUN_DIR)
	@cp $(PARAMS) $(RUN_DIR)/params.nml
	mpirun --mca btl ^tcp -np $(NPROC) --wdir $(RUN_DIR) $(abspath $(TARGET))

clean:
	rm -rf $(OBJDIR) $(TARGET)

.PHONY: all run clean
