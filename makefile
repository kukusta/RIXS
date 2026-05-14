MYNAME=RIXS

RIXS_DUM:
	(cd ../.. ; $(MAKE) $(MYNAME))

F90OBJ=createinput.o main.o readrixsdata.o \
readbnd.o readinput.o readsdt.o writerix.o set_nbi_nbf.o set_ncf_ncl.o bz.o \
memory_usage.o geometry.o printdata.o files.o librixs.o time.o readrix.o

OBJ=bzopt.o

OMPOBJ=tetdoss.o

F90OMPOBJ=rixsloss.o set_a.o readrixsmme.o absorption.o set_absorp_rbz.o \
set_absorp_ibz.o kpmesh.o dealloc.o

INCLUDES=PREC.FI

COMMONMODULES=m_rixs.mod

# COBJ= kmesh_wrap.o

programs: m_files.mod $(SYSBIN)/rixs$(VERSION) $(SYSBIN)/jbns$(VERSION)

# $(OMPFLAGS)
m_files.mod:
	$(FC) -c $(FFLAGS) $(OPTFLAGS) $(SRCDIR)/modules.f90

# m_group.mod:
# 	$(FC) -c $(FFLAGS) $(OPTFLAGS) $(OMPFLAGS) $(SRCDIR)/group_operation.f90

# group_operation.o
$(SYSBIN)/rixs$(VERSION): modules.o $(OBJ) $(F90OBJ) $(F90OMPOBJ) $(COBJ) $(OMPOBJ) $(SYSLIB)/libLMTO.a
	$(LINK) -o $@ $(OBJ) $(F90OBJ) $(F90OMPOBJ) $(COBJ) $(OMPOBJ) \
	modules.o $(LFLAGS) -L$(SYSLIB) -lLMTO $(ADDLIB) $(RIXSFLAGS) \
	$(OMPFLAGS) $(MODFLAGS)

modules.o: $(SRCDIR)/modules.f90
	$(FC) -c $(FFLAGS) $(OPTFLAGS) $(RIXSFLAGS) $(OMPFLAGS) $(SRCDIR)/modules.f90

# group_operation.o: $(SRCDIR)/group_operation.f90
# 	$(FC) -c $(FFLAGS) $(OPTFLAGS) $(SRCDIR)/group_operation.f90

MODULES=m_aux.mod m_bz.mod m_results.mod m_sdt.mod m_bnd.mod m_files.mod \
m_params.mod m_time.mod

jointbns.o: $(SRCDIR)/jointbns.f90
	$(FC) -c $(FFLAGS) $(SRCDIR)/jointbns.f90

$(SYSBIN)/jbns$(VERSION): jointbns.o
	$(LINK) -o $@ jointbns.o $(LFLAGS)

objf:
	@echo OPT_FILES=\"$(OPTOBJ)\" > $(FILEOBJ)
	@echo F90OPT_FILES=\"$(F90OPTOBJ)\" >> $(FILEOBJ)
	@echo OPT2_FILES=\"$(OPT2OBJ)\" >> $(FILEOBJ)
	@echo F90OPT2_FILES=\"$(F90OPT2OBJ)\" >> $(FILEOBJ)
	@echo RIXS_FILES=\"$(OBJ)\" >> $(FILEOBJ)
	@echo F90RIXS_FILES=\"$(F90OBJ)\" >> $(FILEOBJ)
	@echo OMP_FILES=\"$(OMPOBJ)\" >> $(FILEOBJ)
	@echo F90RIXSOMP_FILES=\"$(F90OMPOBJ)\" >> $(FILEOBJ)
	@echo NOPT_FILES=\"$(NOPTOBJ)\" >> $(FILEOBJ)
	@echo C_FILES=\"$(COBJ)\" >> $(FILEOBJ)
	@echo INCLUDE_FILES=\"$(INCLUDES)\" >> $(FILEOBJ)
	@echo MODULE_FILES=\"$(MODULES)\" >> $(FILEOBJ)
	@echo COMMON_MODULES=\"$(COMMONMODULES)\" >> $(FILEOBJ)

