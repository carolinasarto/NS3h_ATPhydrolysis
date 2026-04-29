# vmd -dispdev text -e eta_clust.tcl > eta.out

# date of last modification
set DATE 13/07/2019

###########################################################################################################################
#########################################   ADJUSTABLE PARAMETERS   #######################################################
###########################################################################################################################

# topology file
set TOP ../rep1/x.prmtop
# trajectory file
set TRAJ cent.nc
# trajectory file type
set TYPE netcdf
# reference PDB structure (e.g. a protein-ligand complex)
set REF ref.pdb
# residue name of the cosolvent
set COSOLV ETA
# atom names to clusterize
set solv_atom_list {C1 O1}
# number of MD snapshots
set first 0
set last 1200
set step 1
# residues considered in the analysis (alignment + clustering region) -e.g. binding site-
set resid_list {1 to 439}
#set resid_list {4 5 6 18 19 26 27 28 29 30 31 45 48 49 50 52 54 55 57 97 103 116}
# maximum distance between consecutive probes #0.28 default
set dist 0.28
# minimum number of molecule probes to consider a cluster (reference value = 10% of snapshots)
set watnumbermin 120
# delta distance for plots: WFR(r) and population(r)
set dr 0.20
# radius for WFRr
set WFRr 0.60
# porcentaje of population that defines the radius that contains it (between 0 and 1)
set pop 0.90

#############################################################################################################################
##########################################   Useful subrutines   ############################################################
#############################################################################################################################

# VMD directory
set defaultvmddir /usr/local/lib/vmd

# add topology file
mol new ${TOP} type parm7 first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
# load trajectory
mol addfile ${TRAJ} type ${TYPE} first $first last $last step $step filebonds 1 autobonds 1 waitfor all
# load reference
mol new ${REF} type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

# trajectory working ID in vmd
set workingid 0
# reference working ID in vmd
set workingid2 1
# select residues
set sel [atomselect top "resid $resid_list"]

#######################################################################################################

# Prints the RMSD of the protein atoms between each timestep
# and the first timestep for the given molecule id (default: top)
proc print_rmsd_through_time {mol mol_ref resid_list} {
		# use frame 0 for the reference
		set reference [atomselect $mol_ref "resid $resid_list and noh" frame 0]
		# the frame being compared
		set compare [atomselect $mol "resid $resid_list and noh"]

		set num_steps [molinfo $mol get numframes]
		for {set frame 0} {$frame < $num_steps} {incr frame} {
				# get the correct frame
				$compare frame $frame

				# compute the transformation
				set trans_mat [measure fit $compare $reference]
				# do the alignment
				set compare1 [atomselect $mol "all" frame $frame]
				$compare1 move $trans_mat
				# compute the RMSD
				set rmsd [measure rmsd $compare $reference]
				# print the RMSD
				puts "RMSD of $frame is $rmsd"
		}
}

print_rmsd_through_time $workingid $workingid2 $resid_list

# geometric center of selected residues (useful for entropy calculations)
proc geom_center {selection} {
        # set the geometrical center to 0
        set gc [veczero]
        # [$selection get {x y z}] returns a list of {x y z} 
        #    values (one per atoms) so get each term one by one
        foreach coord [$selection get {x y z}] {
           # sum up the coordinates
           set gc [vecadd $gc $coord]
        }
        # and scale by the inverse of the number of atoms
        return [vecscale [expr 1.0 /[$selection num]] $gc]
}

set vector_resid [geom_center $sel]

foreach atom_solv $solv_atom_list {

set outfolder ${COSOLV}_${atom_solv}_DATA ; # directory to put output files

#############################################################################################################################
############################################# ANALYSIS ######################################################################
#############################################################################################################################

#Borra la carpeta $outfolder
     file delete -force $outfolder
#Genera la carpeta $outfolder
     file mkdir $outfolder

#hace que sea top donde está la trayectoria a anlizar
     mol top $workingid

#Saca el número de frames totales
     set n [molinfo top get numframes]

#Genera una lista vacía para su uso posterior
     set outwatsel {}
#Itera sobre cada frame
     for {set f 0} {$f < $n} {incr f} {
     #Selecciona los O de las aguas cercanas al residuo del frame "f"
        # set watsel [atomselect top "noh and water and same residue as within 4 of resid ${resid_list}" frame $f]
          set watsel [atomselect top "resname ${COSOLV} and name ${atom_solv} and same residue as within 4 of resid ${resid_list}" frame $f] 
     #Agregar a beta factor la info del frame
        $watsel set beta [expr {$f * 0.01}]
     #Guarda en un pdb la selección
        animate write pdb watsel.pdb beg $f end $f skip 1 waitfor all sel $watsel top
     #Se incorpora en una lista el pdb generado si fueron seleccionadas aguas
        if {[llength [$watsel get serial]] < 1} {} else {
           append outwatsel "[exec grep ATOM watsel.pdb]\n"
        }
     #Borra la selección
        $watsel delete
     #Escribe en la consola el frame hecho
       # puts "the waters of ${f}/${n} frame was done"
     }
#borra el pdb del último frame generado
     exec rm watsel.pdb
#Guardo en el archivo watall_X.pdb la lista de aguas 
     set filewatsel [open watall.pdb w]
     puts $filewatsel "$outwatsel"
     close $filewatsel
#Carga el pdb que contiene los atomos de oxigeno de las aguas previemante elegidas foto a foto
     mol load pdb watall.pdb
     set cluster 1
#Selecciona todos los átomos de O
     set sel [atomselect top "all"]
#Extrae el serial de todos los átomos de O
     set serialsel [$sel get serial]
#Borra la selección
     $sel delete
#COMIENZA LA CLUSTERIZACION
#Genera una lista
     set outwatclust {}
#Genera una lista
     set outwatcent {}
#extrae el número de aguas de la lista serialsel
     set lserialsel [llength $serialsel]
#Itera sobre la lista lserialsel hasta que no quede ningún elemento
     while {$lserialsel > 0} {
     #Extrae el primer elemento de la lista serialsel y lo llama i
        set i [lindex $serialsel 0]
     #Genera una lista que inicia con el elemento i
        set serialisel $i
     #Cuenta cuantos elementos posee la lista serialisel
        set lserialisel [llength $serialisel]
     #Genera una lista que inicia con el elemento i
        set serialiselall $i
     #Itera sobre la lista lserialisel hasta que no quede ningún elemento
        while {$lserialisel > 0} {
        #Extrae el primer elemento de la lista serialisel y lo llama j
           set j [lindex $serialisel 0]
        #Borra de la lista serialisel el primer elemento
           set serialisel [lreplace $serialisel 0 0]
        #Busca en la lista serialsel el elemento $j
           set serialsellocation [lsearch $serialsel $j]
        #Borra de la lista serialsel el elemento $j
           set serialsel [lreplace $serialsel $serialsellocation $serialsellocation]
        #Selecciona las aguas que están próximas al elemento $j que no fueron seleccionadas previamente
           set jsel [atomselect top "not serial $serialiselall and within $dist of serial $j"]
        #Extrae el serial de la selección
           set serialjsel [$jsel get serial]
        #Borra la selección
           $jsel delete

        #Itera sobre las aguas seleccionadas
           foreach k $serialjsel {
           #Inserta en la lista serialisel las aguas seleccionadas
              set serialisel [linsert $serialisel end $k]
           #Inserta en la lista serialiselall las aguas seleccionadas
              set serialiselall [linsert $serialiselall end $k]
           }

        #Cuenta nuevamente el número de aguas en la lista serialisel
           set lserialisel [llength $serialisel]

        }

     #Obtiene el número de aguas del cluster
        set watnumber [llength $serialiselall]
     #Considera el cluster si el número de aguas es >= a watnumbermin previamente definido
        if {$watnumber >= $watnumbermin} {
        #GUARDA .PDB CON LAS AGUAS DEL CLUSTER Y CALCULA EL CENTRO DE MASA
        #Selecciona las aguas del cluster
           set selcluster [atomselect top "serial $serialiselall"]
        #Guarda en un archivo .pdb los átomos de O del cluster
           animate write pdb $outfolder/watclust${cluster}.pdb skip 1 waitfor all sel $selcluster top
        #Calcula el centro de masa del cluster
           set clustcenter [measure center $selcluster]
        #Extrae las coordenas del centro de masa
           set x [lindex $clustcenter 0]
           set y [lindex $clustcenter 1]
           set z [lindex $clustcenter 2]
        #GENERA EL PDB DEL CENTRO DE MASA
        #Selecctiona la primer agua del cluster
           set selfirst [atomselect top "serial [lindex $serialiselall 0]"]
        #Cambia las coordenadas xyz del agua seleccionada por las del centro de masa
           $selfirst set x $x
           $selfirst set y $y
           $selfirst set z $z
           $selfirst set resid $cluster
        #Guarda en un archivo .pdb el átomos de O con las coordenadas xyz del centro de masa
           animate write pdb center.pdb skip 1 waitfor all sel $selfirst top
        #Borra la selección
           $selfirst delete
        #Incorpora en la lista outwatclust los datos a guardar
           append outwatcent "[exec grep ATOM center.pdb]\n"
        #Borra el file center.pdb
           exec rm -f center.pdb
        #GENERA EL WFR
        #Selecciona las aguas que están a un radio de 0.6A del centro del cluster
           set selWFR [atomselect top "serial $serialiselall and ((x-$x)*(x-$x) + (y-$y)*(y-$y) + (z-$z)*(z-$z) < (${WFRr}*${WFRr}))"]
        #Extrae el serial de la selección
           set serialWFR [$selWFR get serial]
        #Borra la selección
           $selWFR delete
        #Cuenta las aguas de la selección
           set lserialWFR [llength $serialWFR]
        #Ecuación de WFR para r=0.6A
        #Wat(r) / (#fotos * 4/3 * phi * r^3 * densidad)
           set WFR [expr {$lserialWFR / ($n * 4/3 * ${M_PI} * $WFRr * $WFRr * $WFRr * 0.00215)}]
        #GRAFICA LA POBLACIÓN DEL CLUSTER EN FUNSIÓN DEL RADIO
        #Genera una lista vacía para su uso posterior
           set outwatr {}
        #Inicia con la pob a radio 0.00A
           set popr {0.00}
        #Inicia con el radio en 0.00A
           set r {0.00}
        #Itera sobre el radio hasta hallar el 100 % de la población
           while {$popr < 1.0} {
       ####    #puts $popr
           #Incrementa de a "dr" A el radio
              set r [expr {$r + $dr}]
           #Selecciona las aguas del cluster que se encuentran en la capa entre los radios r-dr y r para el cálculo del g(r)
              if {$r == $dr} {
                 set selgr [atomselect top "serial $serialiselall and ((x-$x)*(x-$x) + (y-$y)*(y-$y) + (z-$z)*(z-$z) < (${r}*${r}))"]
                 } else {
                set selgr [atomselect top "serial $serialiselall and ((x-$x)*(x-$x) + (y-$y)*(y-$y) + (z-$z)*(z-$z) < (${r}*${r})) and not ((x-$x)*(x-$x) + (y-$y)*(y-$y) + (z-$z)*(z-$z) <= ([expr {$r-$dr}]*[expr {$r-$dr}]))"]
                 }
           #Selecciona las aguas del cluster que se encuentran en un radio r para el cálculo del wrf(r)
              set selwfrr [atomselect top "serial $serialiselall and ((x-$x)*(x-$x) + (y-$y)*(y-$y) + (z-$z)*(z-$z) < (${r}*${r}))"]
           #Extrae el serial de la selección
              set serialgr [$selgr get serial]
              set serialwfrr [$selwfrr get serial]
           #Borra la selección
              $selgr delete
              $selwfrr delete
           #Ecuación de g(r)
           #Wat(dr) / (4/3 * phi * (r^3 - (r-dr)^3) * densidad del solvente bulk * número de fotos)
              set gr [expr {[llength $serialgr] / (4/3 * ${M_PI} * ( ($r * $r * $r) - (($r-$dr) * ($r-$dr) * ($r-$dr)) ) * 0.00215 * $n)}]
           #Ecuación de WFR(r)
           #Wat(r) / (#fotos * 4/3 * phi * r^3 * densidad)
              set wfr [expr {[llength $serialwfrr] / ($n * 4/3 * ${M_PI} * ($r * $r * $r) * 0.00215)}]
           #Incorpora en la lista outwatr los datos a guardar
              append outwatr [format "%.3f\t%.3f\t%.3f\n" [expr {$r - $dr/2}] $gr $wfr]
           #Calcula el porcentaje de pop de WAT utilizado dentro del cluster para el radio r
              set popr [expr {[llength $serialwfrr] * 1.000 / $watnumber}]
           #ENCUENTRA EL RADIO AL CUAL SE HALLA DETERMADA CANTIDAD DE LA POBLACION DE AGUAS DEL CLUSTER
           #Entra al loop siempre que la población a dado r sea <= que la pop determinada los "parámetros ajustables"
              if {$popr <= $pop} {set rpop $r} else {}
           }

        #Guarda en un file .dat la info f(r) de cada cluster
           set filewatr [open $outfolder/watr${cluster}.dat w]
           puts $filewatr "$outwatr"
           close $filewatr
         #Incorpora en la lista outwatclust los datos a guardar
           if {$cluster == 1} {
              append outwatclust [format "WS\t#TotWat\tWFP(r=${WFRr})\tr(pop=${pop})\n"]
              append outwatclust [format "WS%.0f\t%.0f\t%.2f\t\t%.2f\n" ${cluster} $watnumber $WFR $rpop]
           } else {
             append outwatclust [format "WS%.0f\t%.0f\t%.2f\t\t%.2f\n" ${cluster} $watnumber $WFR $rpop]
           }
               

##################
        # writes on terminal...
#           puts "cluster $cluster done\nanalysing if there is another cluster of waters...\n"
        #Incrementa el número de cluster
           set cluster [expr {$cluster + 1}]
           } else {}
        #Cuenta nuevamente el número de aguas en la lista serialsel
           set lserialsel [llength $serialsel]
       }
     #Borra la última representación cargada
        mol delete top
     #Borra el file watall.pdb
        exec rm watall.pdb
     #Guarda en un archivo .pdb las coordenas del centro de masa de todos los clusters de cada resid
        set filewatcent [open $outfolder/watcent.pdb w]
        puts $filewatcent "$outwatcent"
        close $filewatcent
     #Carga los centros de masa de todos los watersites y los resperesnta por VDW y colores
        mol load pdb $outfolder/watcent.pdb
        mol delrep 0 top
        mol color type
        mol representation VDW
        mol material Transparent
        mol addrep top
     #Guarda en un file .dat toda la info de cada WS
        set filewatclust [open $outfolder/watclust.dat w]
        puts $filewatclust "$outwatclust"
        close $filewatclust

#################################################################################################################
######################################### SCRIPT END ############################################################
#################################################################################################################
}

quit

