# Load a structure file and a trajectory file
mol new x.prmtop
mol addfile centered.nc waitfor all

# Get the total number of frames
set num_frames [molinfo top get numframes]

# Get the total number of frames
set num_frames [molinfo top get numframes]
set pi [expr {atan(1) * 4}]

# Pg is selected by index 7064
set selection_list {"(atomicnumber 8 and resname WAT) and (same resid as within 2.7 of index 7064)"}

# Angle 0: all angles, Angle 160: between 160-180
set angle_list {0 160}

# This example is for RNA-7A
puts "rna7"
foreach selection $selection_list {
	foreach angle_cutoff $angle_list {
	   set final_frames {}
	   set total_resid {}
	   set control_frames_sin_aguas 0
	   set control_frames_mal_angulo {}

	   # Loop over all the frames
	   for {set i 0} {$i < $num_frames} {incr i} {
	       # Go to the i-th frame
	       animate goto $i
	   
	       # Selection for PG atoms
	       set sel_pg [atomselect top "index 7064"]

	       # Selection for PB atoms
	       set sel_ob [atomselect top "index 7067"]

	       # Get positions of PG and PB atoms
	       set pg_pos [$sel_pg get {x y z}]
	       set ob_pos [$sel_ob get {x y z}]
	    
	       # Select WAT
	       set sel_wat [atomselect top $selection]
	       set wat_resid_lista [$sel_wat get resid]
	       set number_waters [llength $wat_resid_lista]
	       set ow_pos [$sel_wat get {x y z}]
	       
	       if {$number_waters > 0} {
		  
		  # Loop over all OW atoms in the selected WAT
		  for { set j 0 } { $j < $number_waters } { incr j } {
		     
		     set agua_resid [lindex $wat_resid_lista $j]
		     
		     # Calculate the vectors PG-OW and PG-PB
		     set ow_pg_vec [vecsub [lindex $ow_pos $j] [lindex $pg_pos 0]]
		     set ob_pg_vec [vecsub [lindex $ob_pos 0] [lindex $pg_pos 0]]		    
		     
		     # Calculate the angle between the vectors OW-PG and PB-PG
		     set nvec3_2 [vecnorm $ow_pg_vec]
		     set nvec3_4 [vecnorm $ob_pg_vec]
		     set angle_top [vecdot $nvec3_2 $nvec3_4]
		     set angle_rad [expr acos($angle_top)]
		     set angle [expr $angle_rad*(180/$pi)]
		     
		     # Check the angle
		     if { ($angle >= $angle_cutoff) && ($agua_resid ni $total_resid)} {
			lappend total_resid $agua_resid
		     }   
		     if {($angle >= $angle_cutoff) && ($i ni $final_frames)} {
			lappend final_frames $i
		   
		     } 
		  }
		  # Si no pudiste agregar el frame a final_frames, quiere decir que no le daba el angulo 	  
		  if {$i ni $final_frames} {
		     lappend control_frames_mal_angulo $i
		  }	  
		    
		  
	       } else {
		    incr control_frames_sin_aguas
	       }
	       # Delete selections to free memory
	       $sel_wat delete
	       $sel_pg delete
	       $sel_ob delete
	   }		
	#   puts $total_resid
	#   puts $final_frames
	#   set r [llength $total_resid]
	#   puts "This is the length of total_resid: $r"
	   set l [llength $final_frames]
	   puts "This is the length of final_frames: $l"
	#   puts "Hay $control_frames_sin_aguas frames sin agua"
	   set m [llength $control_frames_mal_angulo]
	   set suma [expr $l + $control_frames_sin_aguas]
	   set suma_total [expr $suma + $m]
	#   puts "Hay $m frames con mal angulo"
	   puts "Control: Este numero tiene que ser $num_frames: $suma_total"
	   puts "El angulo cut off fue $angle_cutoff"
	   puts "El filtro usado fue $selection"
	   puts "------------------------------------------------> Proximo angulo"
	}
}

exit

