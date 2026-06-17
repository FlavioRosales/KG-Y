# ============================================
#  KG-Y Makefile (AOCC + OpenMPI + HDF5)
# ============================================

SHELL := /bin/bash

FC = h5pfc

MODDIR  = build
OBJDIR  = build
TARGET  = kg_y.exe

FFLAGS_COMMON = -O3 -march=native -mtune=native -ffast-math -J$(MODDIR) -I$(MODDIR)

LIBS = -llapack -lblas

NPROC  ?= 4
PARAMS ?= params.nml
RUN_DIR := $(basename $(notdir $(PARAMS)))

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
  main/hdf5_lib.f90 \
  main/main.f90

OBJ = $(SRC:%.f90=$(OBJDIR)/%.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS_COMMON) -o $@ $^ $(LIBS)

$(OBJDIR)/%.o: %.f90
	@mkdir -p $(dir $@) $(MODDIR)
	$(FC) $(FFLAGS_COMMON) -c $< -o $@

run: $(TARGET)
	@echo ">> Ejecutando con $(NPROC) ranks MPI"
	@echo ">> Archivo de parámetros : $(PARAMS)"
	@echo ">> Directorio de corrida : $(RUN_DIR)"
	@mkdir -p $(RUN_DIR)
	@if [ -f "$(PARAMS)" ]; then \
		echo ">> Copiando $(PARAMS) -> $(RUN_DIR)/params.nml"; \
		cp "$(PARAMS)" "$(RUN_DIR)/params.nml"; \
	else \
		echo ">> [AVISO] No se encontró $(PARAMS). Se usarán defaults."; \
	fi
	HDF5_USE_FILE_LOCKING=FALSE \
	mpirun --bind-to core --map-by slot -np $(NPROC) \
	bash -c "cd $(RUN_DIR) && ../$(TARGET)"

clean:
	rm -rf $(OBJDIR) $(TARGET) *.mod

print-hdf5:
	@echo "which h5pfc ="
	@which h5pfc
	@echo
	@echo "h5pfc -show ="
	@h5pfc -show

print-mpi:
	@echo "which mpif90 ="
	@which mpif90
	@echo
	@echo "mpif90 --showme:command ="
	@mpif90 --showme:command
	@echo
	@echo "mpif90 --showme ="
	@mpif90 --showme

.PHONY: all run clean print-hdf5 print-mpi
