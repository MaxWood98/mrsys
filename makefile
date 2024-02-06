#build arguments  #-Wall -fbounds-check -ffpe-trap=underflow,zero -fopt-info-optimized=$@_opt.dat
# buildargs = -O0 -Wall -fbounds-check -ffpe-trap=underflow,zero,invalid 
# buildargs = -O0 -Wall -fbounds-check
buildargs = -O2 -fopt-info-optimized=$@_opt.dat

#build settings
buildsettings = -J obj 

#object file directory
OBJDIR = obj/

#list object file names -> 1 for each .f90 in the correct compilation order
OBJS = $(addprefix $(OBJDIR), \
		mrsys_data_module.o\
		sparse_mtx_idx_hash_table_module.o\
		sparse_matrix_module.o\
		mrsys_connectivity_module.o\
		io_utilities_module.o\
		mrsys_io_module.o\
		mrsys_ordered_quad_coarsen_module.o\
		mrsys_catmull_clark_module.o\
		mrsys_smooth_reverse_catmull_clark_module.o\
		mrsys_energy_minimised_coarsening_multiresolution_module.o\
		mrsys_main.o\
		)

#object patturn rule -> for every file in $(OBJDIR) of the form var.o make it from src/var.f90
$(OBJDIR)%.o : src/%.f90
	gfortran $(buildsettings) $(buildargs) -c $< -o $@

#main build procedure
build: mrsys_link

#linking procedure
mrsys_link: $(OBJS) $(addprefix $(OBJDIR), mrsys_main.o)
	gfortran -o mrsys $^ $(buildsettings) -I obj $(buildargs) 

#clean procedure 
clean: 
	rm obj/*.mod
	rm obj/*.o 
	rm obj/*.dat