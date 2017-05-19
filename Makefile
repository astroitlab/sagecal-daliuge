OUTPUT=
CXX=g++
CXXFLAGS=-O3 -Wall -g
CASA_LIBDIR=/usr/local/lib
CASA_INCDIR=-I/usr/local/include/casacore
CASA_LIBS=-lcasa_casa -lcasa_tables -lcasa_measures -lcasa_ms
LAPACK=-lopenblas -lgfortran -lpthread
LAPACK_DIR=/usr/lib64

LDFLAGS=-Wl,--rpath,/usr/lib64,--rpath,${CASA_LIBDIR}

MY_LIBS=-lm -lsagecal
INCLUDES=-I. -I./lib $(CASA_INCDIR) -I/usr/include
LIBPATH=-L$(LAPACK_DIR) -L$(CASA_LIBDIR)  -L./lib

#### glib
GLIBI=-I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include -I/usr/lib/x86_64-linux-gnu/glib-2.0/include/
GLIBL=-lglib-2.0

default: all

data.o: data.cpp data.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GLIBI) -c $<

utils.o: utils.cpp utils.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GLIBI) -c $<

aux_reader.o: aux_reader.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GLIBI) -c $<

admm_master.o: admm_master.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GLIBI) -c $<

admm_slave.o: admm_slave.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GLIBI) -c $<

fratio_master.o: fratio_master.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GLIBI) -c $<

fratio_slave.o: fratio_slave.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GLIBI) -c $<

coh_slave.o: coh_slave.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GLIBI) -c $<

sagefit_slave.o: sagefit_slave.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GLIBI) -c $<

update_z_master.o: update_z_master.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GLIBI) -c $<

update_y_slave.o: update_y_slave.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GLIBI) -c $<

write_z_master.o: write_z_master.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GLIBI) -c $<

write_residual_slave.o: write_residual_slave.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(GLIBI) -c $<

aux_reader: aux_reader.o data.o utils.o ./lib/libsagecal.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $(GLIBI) $(LIBPATH)  -o $@ aux_reader.o data.o utils.o $(MY_LIBS) $(LAPACK) $(CASA_LIBS)  $(GLIBL)

admm_master: admm_master.o data.o utils.o ./lib/libsagecal.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $(GLIBI) $(LIBPATH)  -o $@ admm_master.o data.o utils.o $(MY_LIBS) $(LAPACK) $(CASA_LIBS)  $(GLIBL)

admm_slave: admm_slave.o data.o utils.o ./lib/libsagecal.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $(GLIBI) $(LIBPATH)  -o $@ admm_slave.o data.o utils.o $(MY_LIBS) $(LAPACK) $(CASA_LIBS)  $(GLIBL)

fratio_master: fratio_master.o data.o utils.o ./lib/libsagecal.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $(GLIBI) $(LIBPATH)  -o $@ fratio_master.o data.o utils.o $(MY_LIBS) $(LAPACK) $(CASA_LIBS)  $(GLIBL)

fratio_slave: fratio_slave.o data.o utils.o ./lib/libsagecal.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $(GLIBI) $(LIBPATH)  -o $@ fratio_slave.o data.o utils.o $(MY_LIBS) $(LAPACK) $(CASA_LIBS)  $(GLIBL)

sagefit_slave: sagefit_slave.o data.o utils.o ./lib/libsagecal.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $(GLIBI) $(LIBPATH)  -o $@ sagefit_slave.o data.o utils.o $(MY_LIBS) $(LAPACK) $(CASA_LIBS)  $(GLIBL)

update_z_master: update_z_master.o data.o utils.o ./lib/libsagecal.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $(GLIBI) $(LIBPATH)  -o $@ update_z_master.o data.o utils.o $(MY_LIBS) $(LAPACK) $(CASA_LIBS)  $(GLIBL)

update_y_slave: update_y_slave.o data.o utils.o ./lib/libsagecal.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $(GLIBI) $(LIBPATH)  -o $@ update_y_slave.o data.o utils.o $(MY_LIBS) $(LAPACK) $(CASA_LIBS)  $(GLIBL)

write_z_master: write_z_master.o data.o utils.o ./lib/libsagecal.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $(GLIBI) $(LIBPATH)  -o $@ write_z_master.o data.o utils.o $(MY_LIBS) $(LAPACK) $(CASA_LIBS)  $(GLIBL)

write_residual_slave: write_residual_slave.o data.o utils.o ./lib/libsagecal.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $(GLIBI) $(LIBPATH)  -o $@ write_residual_slave.o data.o utils.o $(MY_LIBS) $(LAPACK) $(CASA_LIBS)  $(GLIBL)

coh_slave: coh_slave.o data.o utils.o ./lib/libsagecal.a
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $(GLIBI) $(LIBPATH)  -o $@ coh_slave.o data.o utils.o $(MY_LIBS) $(LAPACK) $(CASA_LIBS)  $(GLIBL)

all: aux_reader admm_master admm_slave fratio_master fratio_slave sagefit_slave update_z_master write_z_master write_residual_slave update_y_slave coh_slave

clean:
	rm *.o aux_reader admm_master admm_slave fratio_master fratio_slave sagefit_slave update_z_master write_z_master write_residual_slave update_y_slave coh_slave
