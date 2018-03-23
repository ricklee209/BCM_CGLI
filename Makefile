OBJS = main.o BCM_Grid.o BCM_Immersed_boundary.o BCM_Initial_condition.o BCM_Interface_table.o BCM_Interface.o \
BCM_FWS_Interface.o BCM_Reading_IBM.o \
BCM_Ghostcell_minus.o BCM_Ghostcell_plus.o BCM_Ghostcell_minus_Tem.o BCM_Ghostcell_plus_Tem.o \
BCM_X_boundary_condition.o BCM_Y_boundary_condition.o BCM_Z_boundary_condition.o \
BCM_Abs_XYZ_boundary_condition.o BCM_Abs_XY_boundary_condition.o BCM_Abs_X_boundary_condition.o BCM_Abs_Y_boundary_condition.o BCM_Abs_Z_boundary_condition.o \
BCM_Flux_XYZ_Viscous_Runge_kutta.o BCM_Flux_XYZ_Viscous_DPLUSGS.o BCM_Flux_XYZ_Viscous_LUSGS.o\
BCM_Output.o BCM_Output_coarse.o BCM_Statistic.o BCM_ADM_filter.o \
BCM_Mean_pressure_coefficient_Sphere.o BCM_Point_probe.o BCM_Slice_output.o\

HEAD = BCM.h MPI_prm.h prm.h Resolution.h

.SUFFIXES: .o .cpp

.cpp.o:
	mpicxx -O3 -c  $< 


a.out : $(OBJS)

	mpicxx -O3 -o a.out  $(OBJS) 


main.o: main.h Array.h $(HEAD) Pre_selection.h

BCM_Grid.o: $(HEAD) 

BCM_Immersed_boundary.o: $(HEAD) BCM_Immersed_boundary.h 

BCM_Initial_condition.o: $(HEAD) Pre_selection.h

BCM_Interface_table.o: $(HEAD)

BCM_Interface.o: $(HEAD) 

BCM_FWS_Interface.o: $(HEAD) 

BCM_Reading_IBM.o: $(HEAD) 

BCM_Ghostcell_minus.o: $(HEAD)  

BCM_Ghostcell_plus.o: $(HEAD) 

BCM_Ghostcell_minus_Tem.o: $(HEAD)  

BCM_Ghostcell_plus_Tem.o: $(HEAD) 

BCM_X_boundary_condition.o: $(HEAD) 

BCM_Y_boundary_condition.o: $(HEAD) 

BCM_Z_boundary_condition.o: $(HEAD) 

BCM_Abs_XYZ_boundary_condition.o: $(HEAD) Pre_selection.h

BCM_Abs_XY_boundary_condition.o: $(HEAD)

BCM_Abs_X_boundary_condition.o: $(HEAD) 

BCM_Abs_Y_boundary_condition.o: $(HEAD) 

BCM_Abs_Z_boundary_condition.o: $(HEAD) 

BCM_Flux_XYZ_Viscous_Runge_kutta.o: $(HEAD) Pre_selection.h

BCM_Flux_XYZ_Viscous_DPLUSGS.o: $(HEAD) Pre_selection.h implicit.h

BCM_Flux_XYZ_Viscous_LUSGS.o: $(HEAD) Pre_selection.h implicit.h

BCM_Output.o: $(HEAD)

BCM_Output_coarse.o: $(HEAD) Pre_selection.h

BCM_Statistic.o: $(HEAD)

BCM_ADM_filter.o: $(HEAD)

BCM_Point_probe.o: $(HEAD)

BCM_Slice_output.o: $(HEAD)


.PHONY: clean
clean:
	rm $(OBJS)



