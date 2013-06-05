
all: gns_velocity_linz gns_velocity

FOPTS=-g -ffixed-line-length-none -fno-align-commons

gns_velocity_linz: gns_velocity_linz.f makefile
	f95 $(FOPTS) gns_velocity_linz.f -o gns_velocity_linz


gns_velocity: gns_velocity.f makefile
	f95 $(FOPTS) gns_velocity.f -o gns_velocity

