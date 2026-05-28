# ============================================
#  KG-Y Makefile (Ubuntu + OpenMPI + HDF5)
# ============================================

FC = h5pfc

# -------- Flags --------
FFLAGS_COMMON = -O3 -march=native # -Jbuild -fcheck=all 
MODDIR  = build
OBJDIR  = build
TARGET  = kg_y.exe

# -------- MPI run defaults --------
NPROC  ?= 4
PARAMS ?= params.nml
RUN_DIR := $(basename $(notdir $(PARAMS)))

# ============================================
#  Fuentes
# ============================================
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

# ============================================
#  Regla principal
# ============================================
all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS_COMMON) -o $@ $^

# ============================================
#  Compilación
# ============================================
$(OBJDIR)/%.o: %.f90
	@mkdir -p $(dir $@) $(MODDIR)
	$(FC) $(FFLAGS_COMMON) -c $< -o $@

# ============================================
#  Ejecución MPI
# ============================================
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
	mpirun --bind-to core --map-by slot --mca btl ^tcp -np $(NPROC) \
	bash -c "cd $(RUN_DIR) && ../$(TARGET)"

# ============================================
#  Limpieza
# ============================================
clean:
	rm -rf $(OBJDIR) $(MODDIR) $(TARGET)

# ============================================
#  Diagnóstico
# ============================================
print-hdf5:
	@echo "FC = $(FC)"
	@echo "h5pfc -show ="
	@h5pfc -show

print-run:
	@echo "NPROC   = $(NPROC)"
	@echo "PARAMS  = $(PARAMS)"
	@echo "RUN_DIR = $(RUN_DIR)"

.PHONY: all run clean print-hdf5 print-run
