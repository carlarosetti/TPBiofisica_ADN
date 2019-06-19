# Auth: RCarr 4/27/2009
#


# ---  Load all the molecules --- #

set molID [mol new pqrm1_editado.pqr type pqr first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all ]
mol top $molID

# --- --- --- --- --- --- --- --- #
# --  Load the representations -- #

# Do the Front molecule

mol delrep 0 $molID
set repno 0
mol representation Licorice 0.1 33 33
mol color "ResName"
mol selection {chain A B}
mol material Opaque
mol addrep $molID
mol smoothrep $molID $repno 2
incr repno

mol representation Licorice 0.2 33.0 33.0
mol color "Name"
mol selection {chain A B and backbone}
mol material Opaque
mol addrep $molID
mol smoothrep $molID $repno 2
mol showrep $molID $repno off
incr repno

mol representation Licorice 0.1 33.0 33.0
mol color "Charge"
mol selection {chain A B}
mol material Opaque
mol addrep $molID
mol smoothrep $molID $repno 2
mol showrep $molID $repno off
incr repno

mol representation Hbonds 3.3 31 4
mol color "Name"
mol selection {chain A B}
mol addrep $molID
mol smoothrep $molID $repno 2
mol showrep $molID $repno off
incr repno

mol representation Tube 5 1
mol color "ColorID" 1
mol selection {chain A B and backbone}
mol material Opaque
mol addrep $molID
mol smoothrep $molID $repno 2
mol showrep $molID $repno off
incr repno

mol representation NewCartoon 0.3 10 4.1 0 
mol color "Chain"
mol selection {chain E F}
mol material Opaque
mol addrep $molID
mol smoothrep $molID $repno 2
mol showrep $molID $repno off
incr repno

mol representation VDW 0.3 41
mol color "ColorID" 0
mol selection {chain E F and basic}
mol material Opaque
mol addrep $molID
mol smoothrep $molID $repno 2
mol showrep $molID $repno off
incr repno

mol representation VDW 0.3 41
mol color "Charge"
mol selection {resname DG and name N7 O6 or resname DA and name N7 or resname DT and name C7 H71 H72 H73 O4 or resname DC and name H5}
mol material Opaque
mol addrep $molID
mol smoothrep $molID $repno 2
mol showrep $molID $repno off
incr repno

mol representation VDW 0.8 41
mol color "Charge"
mol selection {resname DA and name H62 H61 or resname DC and name H42}
mol material Opaque
mol addrep $molID
mol smoothrep $molID $repno 2
mol showrep $molID $repno off
incr repno

mol representation VDW 0.3 41
mol color "Charge"
mol selection {resname DG and name N3 or resname DA and name N3 H2 or resname DC and name O2 or resname DT and name O2}
mol material Opaque
mol addrep $molID
mol smoothrep $molID $repno 2
mol showrep $molID $repno off
incr repno


mol representation VDW 0.8 41
mol color "Charge"
mol selection {resname DG and name H22 or resname DC and name N2}
mol material Opaque
mol addrep $molID
mol smoothrep $molID $repno 2
mol showrep $molID $repno off
incr repno

mol representation Licorice 0.3 33 33 
mol color "ResName"
mol selection {resname PHE and residue 203 or residue 96 97}
mol material Opaque
mol addrep $molID
mol smoothrep $molID $repno 2
mol showrep $molID $repno off
incr repno
# --- --- -- Make it nice --- --- #


mol rename $molID "TP ADN"

display resetview
rotate x by -90



# --- --- --- --- --- --- --- --- #
