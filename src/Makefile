PROG = Poisson2D
DIR  = .
BIN  = ../bin
FBIN = ../bin


#MOD_SFX= d
MOD_SFX=mod
.SUFFIXES:
.SUFFIXES: .f90 .o .$(MOD_SFX)

LINK    = gfortran -O3 -e
COMP90  = gfortran -O3 
COMP77  = gfortran -O3 
LINK    = gfortran -O3
FFLAGS  = -unroll -fdefault-real-8 -fimplicit-none  
FFLAGS2 = -unroll -ffixed-line-length-132    
#FFLAGS  = -unroll8 -autodouble -W0 -tpp7  
#FFLAGS2 = -unroll8 -autodouble -W0 -tpp7  
#LIBDIR  = -L/usr/local/Matlab6p1/extern/lib/glnx86 
#LIB     = -leng -lmat -lmx -lut
#LIBDIR2 = -L/opt/intel/compiler70/ia32/lib
LIB2    = -lintrins

#LINK    = f77 -O3 -e
#COMP90  = ifc -O3 
#COMP77  = ifc -O3 -132 
#LINK    = ifc -O3
#FFLAGS  = -unroll -autodouble -W0   
#FFLAGS2 = -unroll -autodouble -W0   
##FFLAGS  = -unroll8 -autodouble -W0 -tpp7  
##FFLAGS2 = -unroll8 -autodouble -W0 -tpp7  
#LIBDIR  = -L/usr/local/Matlab6p1/extern/lib/glnx86 
#LIB     = -leng -lmat -lmx -lut
#LIBDIR2 = -L/opt/intel/compiler70/ia32/lib
#LIB2    = -lintrins

# ABSOFT fortran
#COMP90  = f90 -O2
#COMP77  = f77 -O2 -W
#LINK    = f90 -O2  
#FFLAGS  =  -cpu:p7  
#FFLAGS2 =  -cpu:p7
#LIBDIR  = -L/usr/local/Matlab6p1/extern/lib/glnx86
#LIB     = -leng -lmat -lmx -lut 
#-llibmat.so -llibmx.so -llibut.so  


##LIBDIR  = -L/usr/local/matlab6p1/extern/lib/glnx86/
#LIB     =  libeng.so libmat.so libmx.so libut.so  

#LINK    = f77 -fast -xO4 -xsafe=mem -xarch=v8plusa -xchip=ultra
#COMP    = f77 -e -fast -xO4 -xsafe=mem -xarch=v8plusa -xchip=ultra

OBJS = $(FBIN)/Poisson2D.o\
	$(FBIN)/constants.o\
	$(FBIN)/MD2D_Grid.o\
	$(FBIN)/Legendre.o\
	$(FBIN)/splibsrc1.o\
	$(FBIN)/Grid2D_Pack.o\
    $(FBIN)/ERR_Var.o\
	$(FBIN)/State_Var.o\
    $(FBIN)/CG_Pack.o\
	$(FBIN)/State_Pack.o\
    $(FBIN)/Error_Pack.o\
    $(FBIN)/Mat.o\
    $(FBIN)/Source.o\
    $(FBIN)/OutPut_Tec.o\
    $(FBIN)/CG_Var.o\
    $(FBIN)/Multigrid_Var.o\
    $(FBIN)/dsygv.o\
    $(FBIN)/ssygv.o\
    $(FBIN)/blas.o\
    $(FBIN)/Precondition.o\
    $(FBIN)/GMRES_CG_Wrapper.o\
    $(FBIN)/GMRES_Wrapper.o\
    $(FBIN)/HMH_GMRES.o\
    $(FBIN)/GMRES.o\
    $(FBIN)/INPUT.o\
    $(FBIN)/connect.o\
    $(FBIN)/dump.o
#    $(FBIN)/pspm_mg.o






$(PROG):  $(OBJS)
	echo $(BIN)
	@echo ""
	@echo "Linking object files $^ into $@"
	$(LINK) $(OBJS) -o $(PROG) $(LIBDIR) $(LIB)  
	mv $(PROG) ../.

#$(LIBDIR) $(LIB)

#Source file
$(FBIN)/%.o: %.f90 
	@echo ""
	@echo "Compiling source file $< into $@"
	$(COMP90) $(FFLAGS) -c $<
	mv *.o $(FBIN)
#	mv *.mod $(MDIR)


$(FBIN)/%.o: %.for 
	@echo ""
	@echo "Compiling source file $< into $@"
	$(COMP77) $(FFLAGS2) -c $<
	mv *.o $(FBIN)

$(FBIN)/%.o: %.f 
	@echo ""
	@echo "Compiling source file $< into $@"
	$(COMP77) $(FFLAGS2) -c $<
	mv *.o $(FBIN)



#Module 
%.$(MOD_SFX):
	@echo ""
	@if [ ! -s $@ ]; then \
		echo "Error, module file deleted?"; exit 1; fi
	@echo "Module file $@ considered updated through dependence on: $^"
	touch $@

#Dependencies for creating objectfile (USE Statement)
$(FBIN)/Poisson2D.o:\
	$(DIR)/constants.$(MOD_SFX)\
	$(DIR)/md2d_grid.$(MOD_SFX)\
	$(DIR)/legendre.$(MOD_SFX)\
    $(DIR)/err_var.$(MOD_SFX)\
	$(DIR)/state_var.$(MOD_SFX)\
    $(DIR)/cg_var.$(MOD_SFX)\
    $(DIR)/multigrid_var.$(MOD_SFX)\
    $(DIR)/gmres.$(MOD_SFX)\
    $(DIR)/input.$(MOD_SFX)

#$(FBIN)/Test_Pack.o:\
	$(DIR)/constants.$(MOD_SFX)\
	$(DIR)/legendre.$(MOD_SFX)\
	$(DIR)/md2d_grid.$(MOD_SFX)\
	$(DIR)/state_var.$(MOD_SFX)

$(FBIN)/Grid2D_Pack.o:\
	$(DIR)/constants.$(MOD_SFX)\
	$(DIR)/md2d_grid.$(MOD_SFX)\
	$(DIR)/legendre.$(MOD_SFX)
$(FBIN)/Source.o:\
    $(DIR)/constants.$(MOD_SFX)

$(FBIN)/Diff_Pack.o:\
	$(DIR)/legendre.$(MOD_SFX)
#	$(DIR)/metric_var.$(MOD_SFX)


$(FBIN)/Precondition.o:\
    $(DIR)/constants.$(MOD_SFX)\
    $(DIR)/legendre.$(MOD_SFX)\
    $(DIR)/md2d_grid.$(MOD_SFX)\
    $(DIR)/state_var.$(MOD_SFX)\
    $(DIR)/multigrid_var.$(MOD_SFX)\
    $(DIR)/cg_var.$(MOD_SFX)

$(FBIN)/State_Pack.o:\
	$(DIR)/constants.$(MOD_SFX)\
	$(DIR)/legendre.$(MOD_SFX)\
	$(DIR)/md2d_grid.$(MOD_SFX)\
	$(DIR)/state_var.$(MOD_SFX)\
        $(DIR)/input.$(MOD_SFX)

$(FBIN)/Material_Pack.o:\
	$(DIR)/md2d_grid.$(MOD_SFX)\
	$(DIR)/material_var.$(MOD_SFX)

#$(FBIN)/Source_Pack.o:\
	$(DIR)/md2d_grid.$(MOD_SFX)\
	$(DIR)/state_var.$(MOD_SFX)\
	$(DIR)/material_var.$(MOD_SFX)

$(FBIN)/CG_Pack.o:\
	$(DIR)/md2d_grid.$(MOD_SFX)\
    $(DIR)/legendre.$(MOD_SFX)\
	$(DIR)/state_var.$(MOD_SFX)\
    $(DIR)/err_var.$(MOD_SFX)\
    $(DIR)/constants.$(MOD_SFX)\
    $(DIR)/cg_var.$(MOD_SFX)

$(FBIN)/GMRES_CG_Wrapper.o:\
    $(DIR)/legendre.$(MOD_SFX)\
    $(DIR)/md2d_grid.$(MOD_SFX)\
    $(DIR)/state_var.$(MOD_SFX)\
    $(DIR)/multigrid_var.$(MOD_SFX)\
    $(DIR)/cg_var.$(MOD_SFX)\
    $(DIR)/input.$(MOD_SFX)

$(FBIN)/GMRES_Wrapper.o:\
    $(DIR)/legendre.$(MOD_SFX)\
    $(DIR)/md2d_grid.$(MOD_SFX)\
    $(DIR)/state_var.$(MOD_SFX)\
    $(DIR)/multigrid_var.$(MOD_SFX)\
    $(DIR)/cg_var.$(MOD_SFX)

$(FBIN)/HMH_GMRES.o:\
    $(DIR)/gmres.$(MOD_SFX)\
    $(DIR)/md2d_grid.$(MOD_SFX)


$(FBIN)/Error_Pack.o:\
	$(DIR)/md2d_grid.$(MOD_SFX)\
	$(DIR)/state_var.$(MOD_SFX)\
    $(DIR)/err_var.$(MOD_SFX)\
	$(DIR)/legendre.$(MOD_SFX)


$(FBIN)/OutPut_Tec.o:\
	$(DIR)/md2d_grid.$(MOD_SFX)\
	$(DIR)/cg_var.$(MOD_SFX)

$(FBIN)/connect.o:\
	$(DIR)/input.$(MOD_SFX)

$(FBIN)/dump.o:\
	$(DIR)/md2d_grid.$(MOD_SFX)

#$(FBIN)/pspm_mg.o:\
#	$(DIR)/md2d_grid.$(MOD_SFX)\
#	$(DIR)/state_var.$(MOD_SFX)\
#	$(DIR)/multigrid_var.$(MOD_SFX)\
#	$(DIR)/legendre.$(MOD_SFX)

#Location of modules in source file
constants.$(MOD_SFX):$(FBIN)/constants.o
md2d_grid.$(MOD_SFX):$(FBIN)/MD2D_Grid.o 
legendre.$(MOD_SFX):$(FBIN)/Legendre.o
state_var.$(MOD_SFX):$(FBIN)/State_Var.o
err_var.$(MOD_SFX):$(FBIN)/ERR_Var.o
cg_var.$(MOD_SFX):$(FBIN)/CG_Var.o
multigrid_var.$(MOD_SFX):$(FBIN)/Multigrid_Var.o
gmres.$(MOD_SFX):$(FBIN)/GMRES.o
input.$(MOD_SFX):$(FBIN)/INPUT.o

clean:
	rm $(FBIN)/*.o
	rm $(DIR)/*.$(MOD_SFX)
