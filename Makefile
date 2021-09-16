include ./make.inc

SOURCES=Definition.f90 \
	Riemann_Module.f90 \
	Nuclear_module.f90 \
	FlameTable_module.f90 \
	EcapTable_module.f90 \
	Flame_module.f90 \
	NuSpec_module.f90 \
	Neutrino_module.f90 \
	WENO_Module.f90 \
	PPM_Module.f90 \
	MP5_Module.f90 \
	PPT_module.f90 \
	Main.f90 \
	AlphaSplit.f90 \
	Boundary1D.f90 \
	Boundary2D.f90 \
	Boundary2D_X.f90 \
	BuildWENO.f90 \
	BuildHydro.f90 \
	CheckNan.f90 \
	CheckRho.f90 \
	EOSTABLE.f90 \
	FermiMo.f90 \
	FindAsh.f90 \
	FindCentralDensity.f90 \
	FindDt.f90 \
	FindEnergy.f90 \
	FindMass.f90 \
	FindPotential.f90 \
	FindPressure.f90 \
	FromRVEToU.f90 \
	GetEpsilon.f90 \
	GetGrid.f90 \
	GetRho.f90 \
	GetRho_EOSPToR.f90 \
	GetRho_EOSRToP.f90 \
	GetRho_simple.f90 \
	GetSponge.f90 \
	GetVel.f90 \
	Helmeos.f90 \
	Ini_Der.f90 \
	Initial.f90 \
 	Interpolation.f90 \
	Multipole.f90 \
	Openfile.f90 \
	Output.f90 \
	Public_NSE.f90 \
	RungeKutta.f90 \
	SedovExplosion.f90 \
	Spatial.f90 \
	TVD.f90 \
	Update.f90 
	
OBJECTS=$(SOURCES:.f90=.o )

CODE_1D: Definition.o $(OBJECTS)  
	$(F90) $(LDFLAGS) -o CUCODE1D $(OBJECTS) 

$(OBJECTS): %.o: %.f90 
	$(F90) $(F90FLAGS) -c $< -o $@

Definition.mod: Definition.f90        
	$(F90) -C Definition.f90

Riemann.mod: Riemann_Module.f90        
	$(F90) -C Riemann_Module.f90

weno.mod: WENO_Module.f90
	$(F90) -C WENO_Module.f90

mp5.mod: MP5_Module.f90
	$(F90) -C MP5_Module.f90

ppm.mod: PPM_Module.f90
	$(F90) -C PPM_Module.f90
 
nuclear.mod: Nuclear_module.f90
	$(F90) -C Nuclear_module.f90
 
ecaptable.mod: EcapTable_module.f90
	$(F90) -C EcapTable_module.f90
 
flametable.mod: FlameTable_module.f90
	$(F90) -C FlameTable_module.f90

flame.mod: Flame_module.f90
	$(F90) -C Flame_module.f90

Neutrino.mod: Neutrino_module.f90
	$(F90) -C Neutrino_module.f90

Nuspec.mod: NuSpec_module.f90
	$(F90) -C NuSpec_module.f90

PPT.mod: PPT_module.f90
	$(F90) -C PPT_module.f90

clean:
	rm -rf Definition
	rm -rf *.o 	
	rm -rf *.mod
	rm tmp.txt
cleanpic:
	rm -rf ./Plot/Flame/*.png
	rm -rf ./Plot/Hydro/*.png
	rm -rf ./Plot/Isotope/*.png
	rm -rf ./Plot/Level-Set/*.png
	rm -rf ./Plot/Neutrino/*.png
	rm -rf ./Plot/Tracer/*.png
cleanfile:
	rm -rf ./Outfile/Flame/*.dat
	rm -rf ./Outfile/Hydro/*.dat
	rm -rf ./Outfile/Isotope/*.dat
	rm -rf ./Outfile/Level-Set/*.dat
	rm -rf ./Outfile/Neutrino/*.dat
	rm -rf ./Outfile/Tracer/*.dat