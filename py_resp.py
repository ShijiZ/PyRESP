#!/usr/bin/env python3

#-----------------------------------------------------------------------
#    History Versions and Authors
#-----------------------------------------------------------------------
#
#      PyRESP version 1.0     Feb 2022 - Shiji Zhao
#
#      RESP   version 2.4     Nov 2013 - q4md-forcefieldtools.org
#      RESP   version 2.2     Jan 2011 - q4md-forcefieldtools.org
#      RESP   version 2.1     Oct 1994 - Jim Caldwell
#      RESP   version 2.0     Sep 1992 - Christopher Bayly
#
#      ESPFIT version 1.0 (modified)   - Ian Gould
#      ESPFIT version 1.0              - U.Chandra Singh and P.A.Kollman
#
#-----------------------------------------------------------------------
#    Affiliations
#-----------------------------------------------------------------------
#
#      Shiji Zhao:
#                 Center for Complex Biological Systems
#                 University of California, Irvine
#                 Irvine, CA 92697-2280
#
#      All other authors:
#                 Department of Pharmaceutical Chemistry
#                 School of Pharmacy
#                 University of California, San Francisco
#                 San Francisco, CA 94143
#
#---------------------------------------------------------------------------------------------
#
#    This program fits the quantum mechanically calculated electrostatic potential
#    at molecular surfaces using electrostatic models with atom-centered (1) permanent
#    charges and (2) induced dipoles and (3) permanent dipoles. The molecular surfaces
#    should be generated beyond Van der Waal surface in order to minimize other 
#    contributions such as exchange repulsion and charge transfer. 
#
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#
#    Unit -i (input) input of general information
#
#---------------------------------------------------------------------------------------------
#
#    -1st line-  TITLE          a character string
#
#---------------------------------------------------------------------------------------------
#
#    -2nd section- control parameters
#       begin with " &cntrl"; end with " &end"
#
#       nmol     =  the number of structure(s) in a multiple structure fit (default 1)
#                   structure(s): orientation(s), conformation(s) or molecule(s)
#  
#       iqopt    =  1  ... reset all initial charges (and permanent dipoles) to zero (default)
#                =  2  ... read in new initial charges (and permanent dipoles) from -q unit
#
#       ihfree   =  0  ... all atoms are restrained
#                =  1  ... hydrogens not restrained (default)
#
#       irstrnt  =  0  ... harmonic restraints (old style)
#                =  1  ... hyperbolic restraint to parameters of zero (default)
#                =  2  ... only analysis of input parameters; no parameterization is carried out
#                          
#       qwt      =  restraint weight for charges; default is 0.0005
#
#       ioutopt  =  0  ... normal run
#                =  1  ... write restart info of new esp to -s unit (default)
#
#       ireornt  =  0  ... normal run (default)
#                =  1  ... reorient molecule to standard reorientation in Gaussian definition before 
#                          calculating molecular dipole and quadrupole moments
#
#       iquad    =  0  ... report molecular quadrupole moment in Buckingham definition
#                =  1  ... report molecular quadrupole moment in Gaussian definition (default)
#
#       ipol     =  0  ... additive RESP model; no atomic dipole calculations
#                   1  ... Applequist scheme without damping
#                   2  ... Tinker-exponential damping scheme
#                   3  ... exponential damping scheme
#                   4  ... linear damping scheme
#                   5  ... pGM damping scheme(default)
#               
#       igdm     =  0  ... normal run
#                =  1  ... use distributed pGM charges and dipoles in ESP fitting (default)
#                          only use with ipol = 5
#
#       exc12    =  0  ... include 1-2 interactions for electric field calculations (default)
#                =  1  ... exclude 1-2 interactions
#
#       exc13    =  0  ... include 1-3 interactions for electric field calculations (default)
#                =  1  ... exclude 1-3 interactions
#
#       ipermdip =  0  ... RESP-ind model; do not calculate permanent dipole
#                =  1  ... RESP-perm model; calculate permanent dipoles (default)
#
#       pwt      =  restraint weight for permanent dipoles; default is 0.0005
#
#       virtual  =  0  ... normal run (default)
#                =  1  ... enable permanent dipoles for 1-3 virtual bonds
#
#---------------------------------------------------------------------------------------------
#
#     -3rd line- wtmol ... relative weight for the structure if multiple structure fit (1.0 otherwise) 
#                    
#---------------------------------------------------------------------------------------------
#
#     -4th line- subtitle for the structure (a character string)
#
#---------------------------------------------------------------------------------------------
#
#     -5th line- charge   iuniq   (iuniq_p)
#        charge  = total charge value for this structure (-99 if no total charge constraint)
#        iuniq   = total number of atoms for this structure
#        iuniq_p = total number of permanent dipoles for this structure
#
#---------------------------------------------------------------------------------------------
#
#     -6th section- one line for each atom
#        element number = element number in periodic table
#        ivary          = control charge variations of each center
#        (ivary_p       = control permanent dipole variations of each center)
#        Note: The permanent dipoles of each atom are ordered with the atom number of reference atoms.
#              If virtual = 1, real permanent dipoles come before all virtual dipoles for each atom.
#        ivary & ivary_p 
#                       =  0 current parameters fitted independently of other centers
#                       = -1 current parameters frozen at "initial stage" value typically read in from -q unit
#                       =  n current parameters fitted and equivalenced to that of center "n"
#
#---------------------------------------------------------------------------------------------
#
#     -7th section- intra-molecular charge constraint(s); blank line if no constriant
#        ngrp = the number of charge centers in the group of atoms associated with this constraint
#        grpchg = charge value to which the associated group of atoms (given on the next section) is to be constrained
#
#     -7.1th section-
#        imol  iatom (repeat if more than 8 centers)
#        the list ("ngrp" long) of the atoms to be constrained to the charge specified on the previous line.
#
#        blank to end
#
#---------------------------------------------------------------------------------------------
#
#     -8th section- inter-molecular charge constraint(s)
#        same format as intra-molecular charge constraint(s) - see the 7th & 7.1th sections
#
#        blank to end
#
#---------------------------------------------------------------------------------------------
#
#     -9th section- multiple structure atom charge equivalencing
#        ngrp = the number of charge enters in the group of atoms for equivalencing
#
#     -9.1th section-
#        imol  iatom (repeat if more than 8 centers)
#        the list ("ngrp" long) of the atoms to be equivalenced
#
#        blank to end
#
#---------------------------------------------------------------------------------------------
#
#     -10th section- multiple structure permanent dipole equivalencing
#        ngrp = the number of permanent dipole centers for the group of permanent dipoles for equivalencing
#
#     -10.1th section-
#        imol  idip (repeat if more than 8 centers)
#        the list ("ngrp" long) of the permanent dipoles to be equivalenced
#
#        blank to end
#
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#
#    Unit -e (espot) input of ESP and coordinates
#
#
#     -1st line- natom   nesp (total number of atoms & ESP points)
#
#
#     -2nd line up to natom+1 line-
#        atom coordinates X Y Z (in Bohrs) & element number & atom type
#        Note: atom type can be generated with the espgen program by setting -p 1 
#
#
#     -natom+2 line up to natom+2+nesp line- ESP & coordinates
#        espot X Y Z (in a.u. & Bohrs)
#
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#
#    Unit -q (qin) input of replacement parameters if requested
#    (note same format as that produced by -t unit)
#
#
#     %FLAG TITLE: a character string
#     subtitle for the structure 
#
#
#     %FLAG ATOM CRD: I4,3E16.7
#     atm.no: atom number 
#     X, Y, Z: atom coordinates X Y Z 
#
#
#     %FLAG ATOM CHRG: 2(I4,X7),I4,X2,E16.7
#     atm,no: atom number
#     element.no: element number in periodic table
#     ivary: charge variations
#     q(opt): optimized atomic charge
#
#
#     %FLAG PERM DIP LOCAL: 3(I4,X5),I4,X2,E16.7 (for ipermdip = 1)
#     dip.no: dipole number
#     atm.no: atom number
#     ref.no: reference atom number
#     ivary: permanent dipole variations
#     p(opt): optimized permanent dipole in local frame 
#
#
#     %FLAG PERM DIP GLOBAL: I4,3E16.7 (for ipermdip = 1)
#     atm.no: atom number 
#     X, Y, Z: optimized permanent dipole in X Y Z directions of global frame
#
#
#     %FLAG IND DIP GLOBAL: I4,3E16.7 (for ipol > 0)
#     atm.no: atom number 
#     X, Y, Z: induced dipole in X Y Z directions of global frame
#
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#
#    Unit -ip (polariz) input of atomic polarizabilities if requested
#
#
#     -1st line- Some comments
#
#
#     -2st section- atom type & atomic polarizability (a.u.) & damping factor (a.u.)
#        Note: For each of the provided damping schemes (Tinker-exponential. exponential, linear, and
#             pGM), the damping factor has different meanings. See PyRESP publication for reference.        
#
#        end with a line starting with 'a' (unused)
#
#
#     -3rd section- atom type equivalence information
#        EQ & a list of atome types sharing identical polarizability and damping factors
#

import argparse as ap
import numpy as np
import sys
import f90nml
import math
from scipy.linalg import lu_factor, lu_solve

###### parameters ######
maxq = 8000

###### constants ######
au2A = 0.52917725   # constant converting from a.u. (bohr) to angstrom
au2D = 2.541765     # constant converting from a.u. (e*bohr) to debye

###### iostuf ######
ioutopt,iqopt = 0,1

###### files ######
Input,output,qin,qout,espot,esout,polariz = [None]*7
input_file,output_file,qin_file,qout_file,espot_file,esout_file,polariz_file = [None]*7

###### infoa ######
iuniq,iuniq_p,nqpl,ihfree,irstrnt,ireornt,iquad = 0,0,0,1,1,0,1

###### runlab ######
title = None

###### espcom ######
a_pot,b_pot = [None]*2
ssvpot_list,tot_nesp = None,0

###### calcul ###### 
qcal,a,b,qwtval,pwtval,iqpcntr = [None]*6

###### polarization ######
ipol,ipermdip,igdm = 5,1,1
pol_dict,AinvQ_list,AinvP_L2G_list,L2G_list,mat_pot = [None]*5

###### lagrng ###### 
grpchg = np.zeros((1))
lgrcnt = np.zeros((1,maxq))
nlgrng = 0

###### orig ###### 
q0,p0,crd,izan,izan_p,ivary,ivary_p,qwt,pwt = [None]*9

###### worker ###### 
awork,bwork = [None]*2

###### propty ######
cmas_mol,co,dipol_mol,dipol_mol_com,dipmom_mol,quad_mol,quad_mol_com = [None]*7
dipind,dipperm,dipindperm,dipindperm_com = [None]*4

###### mltmol ###### 
wtmol,ibeg,iend,ibeg_p,iend_p = [None]*5
nmol = 1

##### geometry #####
n12_list,dict12_list,n13_list,dict13_list,exc12,exc13,virtual = [None]*7

def file_in():
	parser = ap.ArgumentParser(usage='py_resp.py [-O] -i input -o output [-q qin] [-ip polariz] -t qout -e espot -s esout')
	parser.add_argument("-O", action='store_true', help="Overwrite output files if they exist")
	parser.add_argument("-i", "--input", required=True, help="type: input, required; description: input of general information")
	parser.add_argument("-e", "--espot", required=True, help="type: input, required; description: input of ESP and coordinates")
	parser.add_argument("-q", "--qin", help="type: input, optional; description: replacement parameters")
	parser.add_argument("-o", "--output", required=True, help="type: output, always produced; description: output of results")
	parser.add_argument("-t", "--qout", required=True, help="type: output, always produced; description: output of parameters")
	parser.add_argument("-s", "--esout", help="type: output, optional; description: generated ESP values for new parameters")
	parser.add_argument("-ip", "--polariz", help="type: input, optional; description: atomic polarizabilities")
	
	args = parser.parse_args()
	Input,espot,qin,polariz = args.input,args.espot,args.qin,args.polariz
	output,qout,esout = args.output,args.qout,args.esout

	return Input,output,qin,polariz,qout,espot,esout

def read_in():
	#############################################
	# This function reads in control parameters.
 	#
 	# called from Main
 	#############################################
	global nmol,iqopt,ihfree,irstrnt,qwt,ioutopt,ireornt,iquad
	global ipol,igdm,exc12,exc13,ipermdip,pwt,virtual

	output_file.write('\n -----------------------------------------------')
	output_file.write('\n              Py_RESP Version 1.0')
	output_file.write('\n -----------------------------------------------')
	output_file.write('\n '+input_file.readline())
	output_file.write(' -----------------------------------------------\n')
	output_file.write('\n ------------------')
	output_file.write('\n Control Parameters')
	output_file.write('\n ------------------')

	# Read in control parameters
	if 'cntrl' in f90nml.read(Input):
		cntrl_nml = f90nml.read(Input)['cntrl']
		# Common control parameters
		nmol = cntrl_nml['nmol'] if 'nmol' in cntrl_nml else 1
		iqopt = cntrl_nml['iqopt'] if 'iqopt' in cntrl_nml else 1
		ihfree = cntrl_nml['ihfree'] if 'ihree' in cntrl_nml else 1
		irstrnt = cntrl_nml['irstrnt'] if 'irstrnt' in cntrl_nml else 1
		qwt = cntrl_nml['qwt'] if 'qwt' in cntrl_nml else 0.0005
		ioutopt = cntrl_nml['ioutopt'] if 'ioutopt' in cntrl_nml else 0
		ireornt = cntrl_nml['ireornt'] if 'ireornt' in cntrl_nml else 0
		iquad = cntrl_nml['iquad'] if 'iquad' in cntrl_nml else 1

		# With dipole
		ipol = cntrl_nml['ipol'] if 'ipol' in cntrl_nml else 5
		igdm = cntrl_nml['igdm'] if 'igdm' in cntrl_nml else 1
		exc12 = cntrl_nml['exc12'] if 'exc12' in cntrl_nml else 0
		exc13 = cntrl_nml['exc13'] if 'exc13' in cntrl_nml else 0

		# With permanent dipole
		ipermdip = cntrl_nml['ipermdip'] if 'ipermdip' in cntrl_nml else 1
		pwt = cntrl_nml['pwt'] if 'pwt' in cntrl_nml else 0.0005
		virtual = cntrl_nml['virtual'] if 'virtual' in cntrl_nml else 0
	else:
		output_file.write('\n Error: Must use namelist input\n')
		sys.exit()

	output_file.write('\n nmol        = %d   iqopt       = %d'%(nmol, iqopt))
	output_file.write('\n ihfree      = %d   irstrnt     = %d'%(ihfree, irstrnt))
	output_file.write('\n ioutopt     = %d   ireornt     = %d'%(ioutopt, ireornt))
	output_file.write('\n iquad       = %d   qwt         = %.8f'%(iquad, qwt))
	if ipol > 0:
		output_file.write('\n ipol        = %d   igdm        = %d'%(ipol, igdm))
		output_file.write('\n exc12       = %d   exc13       = %d'%(exc12, exc13))
		if ipermdip > 0:
			output_file.write('\n ipermdip    = %d   pwt         = %.8f'%(ipermdip, pwt))
			output_file.write('\n virtual     = %d'%(virtual))

	for line in input_file:
		if "&end" in line:
			break

	if ipermdip > 0 and ipol <= 0:
		output_file.write('\n Error: permanent dipole can be enabled only if ipol > 0\n')
		sys.exit()

def sing_mol():
	##############################################################################################
	# This function reads in single molecule input. nmol = 1 caused this function to be called.
 	#
 	# The molecule control decks contains
 	#  - ich, iuniq (and icntrs_p)
 	#  - izan, ivary (and ivary_p)
 	#  - (Largrange information for individual-molecule)
 	#
 	# called from Main
 	###############################################################################################
	global nlgrng,iuniq,title,wtmol,ibeg,iend,izan,ivary
	if ipermdip > 0:
		global iuniq_p,ibeg_p,iend_p,izan_p,ivary_p

	output_file.write('\n\n --------------------')
	output_file.write('\n Molecule Information')
	output_file.write('\n --------------------')
	output_file.write("\n Single-molecule run")

	title = []
	wtmol = np.ndarray((1))

	# Read in fitting weight for q0 and esp point weighting
	wtmol[0] = float(input_file.readline())
	output_file.write("\n\n Molecule 1 weight:%10.3f \n"%wtmol[0])
	title.append(input_file.readline())
	output_file.write(" Molecule 1 name: %s"%title[0])

	# Read in charge, number of charge centers
	line = input_file.readline().split()
	ich, iuniq = int(line[0]), int(line[1])
	output_file.write(' Total charge (ich):%3d'%ich)
	output_file.write('\n Number of centers:%3d'%iuniq)

	# Read in number of permanent dipoles
	if ipermdip > 0:
		iuniq_p = int(line[2])
		output_file.write('\n Number of permanent dipoles:%3d'%iuniq_p)

	# Initialize izan, ivary (, izan_p and ivary_p)
	# Build ibeg, iend (, ibeg_p and iend_p)
	ibeg = np.ndarray((1), dtype=int)
	iend = np.ndarray((1), dtype=int)
	izan = np.ndarray((iuniq), dtype=int)
	ivary = np.ndarray((iuniq), dtype=int)
	ibeg[0] = 0
	iend[0] = iuniq-1
	wtmol[0] = 1.0
	dip_num = []
	if ipermdip > 0:
		ibeg_p = np.ndarray((1), dtype=int)
		iend_p = np.ndarray((1), dtype=int)
		izan_p = np.ndarray((iuniq_p), dtype=int)
		ivary_p = np.ndarray((iuniq_p), dtype=int)
		ibeg_p[0] = 0
		iend_p[0] = iuniq_p-1
		p_cnt = 0

	# Build atomic number izan[i] and atom equivalencing info ivary[i] for charges
	for i in range(iuniq):
		line = input_file.readline().split()
		izan[i], ivary[i] = int(line[0]), int(line[1])
		output_file.write("\n%5d%5d%5d"%(i+1, izan[i], ivary[i]))

		# Build atomic number izan_p[i] and equivalencing info ivary_p[i] for permanent dipoles
		if ipermdip > 0:
			dip_num.append(len(line)-2)
			for j in range(len(line)-2):
				izan_p[p_cnt] = int(line[0])
				ivary_p[p_cnt] = int(line[j+2])
				output_file.write("%5d"%(ivary_p[p_cnt]))
				p_cnt += 1

	# Read in lagrange constraints (including total and intra-molecular charge constraints)
	lagrange(ich, 0)

	input_summary(dip_num)

def mult_mol():
	##############################################################################################
	# This function reads in multiple molecule input. nmol > 1 caused this function to be called.
	#
 	# The input form for molecules is that their entire control decks are appended to the control 
 	# parameters just read in. Each control deck is separated by a blank line. Then comes the
 	# multiple-molecule specific input, which is
 	#  - the lagrange constraints to be applied between molecules
 	#  - equivalencing of centers between molecules in the series
 	#
 	# Each molecule control decks contains
 	#  - ich, icntrs (and icntrs_p)
 	#  - izan, ivary (and ivary_p)
 	#  - (Largrange information for individual-molecule)
 	#
 	# called from Main
 	###############################################################################################
	global nlgrng,iuniq,title,wtmol,ibeg,iend,izan,ivary
	if ipermdip > 0:
		global iuniq_p,ibeg_p,iend_p,izan_p,ivary_p

	output_file.write('\n\n --------------------')
	output_file.write('\n Molecule Information')
	output_file.write('\n --------------------')
	output_file.write("\n Multiple-molecule run of %d molecules"%nmol)

	title = []
	wtmol = np.ndarray((nmol))

	# Initialize izan, ivary, ibeg, iend (, izan_p, ivary_p, ibeg_p and iend_p)
	ibeg = np.ndarray((nmol), dtype=int)
	iend = np.ndarray((nmol), dtype=int)
	izan = np.ndarray((1), dtype=int)
	ivary = np.ndarray((1), dtype=int)
	dip_num = []
	if ipermdip > 0:
		ibeg_p = np.ndarray((nmol), dtype=int)
		iend_p = np.ndarray((nmol), dtype=int)
		izan_p = np.ndarray((1), dtype=int)
		ivary_p = np.ndarray((1), dtype=int)
		p_cnt = 0

	for imol in range(nmol):
		# Read in fitting weight for q0 and esp point weighting
		wtmol[imol] = float(input_file.readline())
		output_file.write("\n\n Molecule %d weight:%10.3f\n"%(imol+1,wtmol[imol]))
		title.append(input_file.readline())
		output_file.write(" Molecule %d Name: %s"%(imol+1,title[imol]))

		# Read in charge, number of charge centers
		line = input_file.readline().split()
		ich, icntrs = int(line[0]), int(line[1])
		output_file.write(' Total charge (ich):%3d'%ich)
		output_file.write('\n Number of centers:%3d'%icntrs)

		# Read in number of permanent dipoles
		if ipermdip > 0:
			icntrs_p = int(line[2])
			output_file.write('\n Number of permanent dipoles:%3d'%icntrs_p)

		# Now some book-keeping: iuniq is the global variable for the total number of centers over all molecules.
		# The first center of this mol therefore starts in iuniq and goes to iuniq+icntrs-1.
		#
		# Same for permanent dipoles.
		ibeg[imol] = iuniq
		iend[imol] = iuniq+icntrs-1
		if ipermdip > 0:
			ibeg_p[imol] = iuniq_p
			iend_p[imol] = iuniq_p+icntrs_p-1

		# Trap for having too many centers
		if iend[imol]+1 > maxq:
			output_file.write('\n ERROR: more than %5d centers'%maxq)
			sys.exit()

		izan.resize((iuniq+icntrs), refcheck=False)
		ivary.resize((iuniq+icntrs), refcheck=False)
		iuniq += icntrs
		if ipermdip > 0:
			izan_p.resize((iuniq_p+icntrs_p), refcheck=False)
			ivary_p.resize((iuniq_p+icntrs_p), refcheck=False)
			iuniq_p += icntrs_p

		# Read in atomic number izan[i] and ivary[i]
		# Since ivary[i] is supposed to correspond to a center-number in the same molecule, this has to be adjusted 
		# to ivary[i]+ibeg[imol]
		for i in range(ibeg[imol], iend[imol]+1):
			line = input_file.readline().split()
			izan[i], ivary[i] = int(line[0]), int(line[1])
			output_file.write("\n%5d%5d%5d"%(i+1, izan[i], ivary[i]))
			if ivary[i] > 0:
				ivary[i] += ibeg[imol]

			# Read in atomic number izan_p[i] and equivalencing info ivary_p[i] for permanent dipoles
			if ipermdip > 0:
				dip_num.append(len(line)-2)
				for j in range(len(line)-2):
					izan_p[p_cnt] = int(line[0])
					ivary_p[p_cnt] = int(line[j+2])
					output_file.write("%5d"%(ivary_p[p_cnt]))
					if ivary_p[p_cnt] > 0:
						ivary_p[p_cnt] += ibeg_p[imol]
					p_cnt += 1		

		# Now read in the lagrange constraints for this molecule (including total and intra-molecular charge constraints)
		lagrange(ich, imol)

	# End of molecule input, now do other preparation stuff

	# Read past a blank line after the final molecule control deck and then read in inter-molecule lagrange constraints.
	# The "-99" for the total charge tells lagrange to drop the total charge constraint. The "-1" for the molecule number
	# tells lagrange that this is inter-molecular charge constraints.
	lagrange(-99, -1)

	# Perform inter-molecule equivalencing
	mol_equiv()

	input_summary(dip_num)

def lagrange(ncharge, imol):
	#################################################################
	# This function read in and assign lagrange constraint pointers.
	#
	# called from sing_mol() and mult_mol() 
	#################################################################
	global nlgrng

	line = input_file.readline().split()
	# If line is not empty, implement intra- or inter- molecular charge constraint(s)
	if line:
		output_file.write("\n ---------------------------------------------------------")
		if imol > -1:
			output_file.write("\n Intra-Molecular Charge Constraint(s) Info for Molecule %d"%(imol+1))
		else:
			output_file.write("\n Inter-Molecular Charge Constraint(s) Info")
		output_file.write("\n ---------------------------------------------------------")
		while line:  # line is not empty
			nlgrng += 1
			grpchg.resize((nlgrng), refcheck=False)
			lgrcnt.resize((nlgrng,maxq), refcheck=False)

			ngrp, grpchg[nlgrng-1] = int(line[0]), float(line[1])
			imoll = np.ndarray((ngrp), dtype=int)
			iatm = np.ndarray((ngrp), dtype=int)
			cnt = ngrp
			for row in range(math.ceil(ngrp/8)):
				line = input_file.readline().split()
				output_file.write("\n")
				for j in range(8):
					idx = row*8+j
					imoll[idx] = int(line[2*j])
					iatm[idx] = int(line[2*j+1])
					output_file.write("%5d%5d"%(imoll[idx], iatm[idx]))
					cnt -= 1
					if cnt == 0:
						break

			for i in range(ngrp):
				iatm[i] += ibeg[imoll[i]-1]     # -1 since python is 0 indexed
				lgrcnt[nlgrng-1][iatm[i]-1] = 1 # -1 since python is 0 indexed

			line = input_file.readline().split()

	# As long as ncharge is not -99, implement the "total charge" constraint
	if ncharge > -99:
		nlgrng += 1
		grpchg.resize((nlgrng), refcheck=False)
		lgrcnt.resize((nlgrng,maxq), refcheck=False)

		grpchg[nlgrng-1] = float(ncharge)
		for j in range(ibeg[imol], iend[imol]+1):
			lgrcnt[nlgrng-1][j] = 1

def mol_equiv():
	###############################################################################################
	# This function carries out the inter-molecule charge (and permanent dipole) equivalencing by
	#
	# First : read the cards saying how many centers will be read in in the next card.
	# Second: read the first-occurrence-in-each-molecule of the centers to be equivalenced.
	#
	#        - The specifcations MUST be in ascending order.
	#        - The expanding of the centers within each molecule is based on the ivary values for 
	#          the individual mol.
	#
	# called from mult_mol()
	###############################################################################################
	output_file.write("\n\n -------------------------------")
	output_file.write("\n Charges Equivalence Information ")
	output_file.write("\n -------------------------------")

	line = input_file.readline().split()
	while line:  # line is not empty
		ngrp = int(line[0])
		
		imoll = np.ndarray((ngrp), dtype=int)
		iatm = np.ndarray((ngrp), dtype=int)
		cnt = ngrp
		for row in range(math.ceil(ngrp/8)):
			line = input_file.readline().split()
			output_file.write("\n")
			for j in range(8):
				idx = row*8+j
				imoll[idx] = int(line[2*j])
				iatm[idx] = int(line[2*j+1])
				output_file.write("%5d%5d"%(imoll[idx], iatm[idx]))
				cnt -= 1
				if cnt == 0:
					break

		for i in range(ngrp):
			iatm[i] += ibeg[imoll[i]-1]   # -1 since python is 0 indexed

		for i in range(1, ngrp):
			ivary[iatm[i]-1] = iatm[0]   # -1 since python is 0 indexed

		line = input_file.readline().split()

	# If ivary for mol 2+ is greater than zero it is replaced with the atom number of mol 1.
	for i in range(iuniq):
		if ivary[i] > 0:
			vary = ivary[ivary[i]-1]
			if vary > 0:
				ivary[i] = vary

	if ipermdip > 0:
		output_file.write("\n\n -----------------------------------------")
		output_file.write("\n Permanent Dipoles Equivalence Information")
		output_file.write("\n -----------------------------------------")

		line = input_file.readline().split()
		while line:  # line is not empty
			ngrp = int(line[0])
			
			imoll = np.ndarray((ngrp), dtype=int)
			idip = np.ndarray((ngrp), dtype=int)
			cnt = ngrp
			for row in range(math.ceil(ngrp/8)):
				line = input_file.readline().split()
				output_file.write("\n")
				for j in range(8):
					idx = row*8+j
					imoll[idx] = int(line[2*j])
					idip[idx] = int(line[2*j+1])
					output_file.write("%5d%5d"%(imoll[idx], idip[idx]))
					cnt -= 1
					if cnt == 0:
						break
	
			for i in range(ngrp):
				idip[i] += ibeg_p[imoll[i]-1]   # -1 since python is 0 indexed
	
			for i in range(1, ngrp):
				ivary_p[idip[i]-1] = idip[0]   # -1 since python is 0 indexed
	
			line = input_file.readline().split()

		# If ivary_p for mol 2+ is greater than zero it is replaced with the atom number of mol 1.
		for i in range(iuniq_p):
			if ivary_p[i] > 0:
				vary = ivary_p[ivary_p[i]-1]
				if vary > 0:
					ivary_p[i] = vary

def input_summary(dip_num):
	#####################################################################################
	# This function outputs summary info for ivary, iuniq, iuniq_p, qwt, pwt and nlgrng.
	#
	# called from sing_mol() and mult_mol()
	#####################################################################################
	output_file.write("\n\n -----------------------------")
	if ipermdip > 0:
		output_file.write("\n Atom/Dipole Ivary Information")
		pcnt = 0
	else:
		output_file.write("\n Atom Ivary Information")
	output_file.write("\n -----------------------------")
	icnt = 0
	jcnt = 0
	for iat in range(iuniq):
		output_file.write("\n%5d%5d"%(izan[iat], ivary[iat]))
		if ipermdip > 0:
			for j in range(dip_num[iat]):
				output_file.write("%5d"%(ivary_p[pcnt]))
				pcnt += 1
		jcnt += 1
		if (jcnt > iend[icnt]):
			output_file.write("\n")
			icnt += 1
	
	output_file.write("\n -------------------------")
	output_file.write("\n Input Information Summary")
	output_file.write("\n -------------------------")
	output_file.write("\n Total number of atoms      =%5d"%iuniq)
	if ipermdip > 0:
		output_file.write("\n Total number of permanent dipoles      =%5d"%iuniq_p)
	output_file.write("\n Weight factor on initial charge restraints=%10.6f"%qwt)
	if ipermdip > 0:
		output_file.write("\n Weight factor on initial dipole restraints=%10.6f"%pwt)
	output_file.write("\n There are%3d charge constraints"%nlgrng)

def read_pol_dict():
	####################################################################################
	# This function builds polarizability dictionary pol_dict for each atom type by
	# reading the polarizability file.
	#
	# called from Main
	####################################################################################
	global pol_dict

	output_file.write("\n\n Read polarizability and radii information from:")
	output_file.write("\n %s"%polariz)
	output_file.write("\n Assume the units are in a.u.")

	pol_dict = {}
	polariz_file = open(polariz, 'r')
	polariz_file.readline()

	line = polariz_file.readline().split()
	while (line[0] != 'a'):
		pol_dict[line[0].lower()] = [float(line[1]), float(line[2])]
		line = polariz_file.readline().split()

	line = polariz_file.readline().split()
	while (line and line[0] == 'EQ'):
		for i in range(2, len(line)):
			pol_dict[line[i].lower()] = pol_dict[line[1].lower()]
		line = polariz_file.readline().split()

def neighbors(imol):
	####################################################################################
	# This function builds the following for molecule imol:
	#  - n12: see build12()
	#  - dict12: see build12()
	#  - n13: see build13()
	#  - dict13: see build13()
	#
	# called from matpot()
	####################################################################################
	natm = iend[imol]-ibeg[imol]+1

	n12, dict12, ibond, nbonds = build12(imol, natm)
	n13, dict13 = build13(ibond, nbonds, natm)
	#build14()
	return n12, dict12, n13, dict13

def build12(imol, natm):
	####################################################################################
	# This function builds the following for molecule imol:
	#  - n12: list for # of 12-connecting atoms of each atom
	#  - dict12: dictionary for 12-connecting atoms of each atom
	#  - ibond: see build_bonds()
	#  - nbond: see build_bonds()
	#
	# called from neighbors()
	####################################################################################
	nbonds, ibond = build_bonds(imol, natm)

	n12 = np.zeros(natm, dtype=int)
	dict12 = {}
	for l in range(nbonds):
		i = ibond[0][l]
		j = ibond[1][l]
		n12[i] += 1
		n12[j] += 1
		if i in dict12:
			dict12[i].append(j)
		else:
			dict12[i] = [j]

		if j in dict12:
			dict12[j].append(i)
		else:
			dict12[j] = [i]
	return n12, dict12, ibond, nbonds

def build_bonds(imol, natm):
	####################################################################################
	# This function builds the following for molecule imol by scanning distance between
	# atom pairs:
	#  - nbond: total # of 12-connecting pairs
	#  - ibond: list of all 12-connecting pairs
	#
	# called from build12()
	####################################################################################
	ibond = np.ndarray((2, 4*natm), dtype=int)
	cutoff = [0.8,1.95,1.6,1.85,1.58,1.50,1.2]    # in angstrom
	cutoff2 = [(i/au2A)**2 for i in cutoff]

	nbonds = 0
	for i in range(ibeg[imol], iend[imol]):
		izi = izan[i]
		for j in range(i+1, iend[imol]+1):
			dist2 = (crd[0][j]-crd[0][i])**2 + (crd[1][j]-crd[1][i])**2 + (crd[2][j]-crd[2][i])**2
			bonded = False
			if dist2 < cutoff2[0]:
				output_file.write("\n Molecule %d has wrong geometry. Stop!!"%imol)
				sys.exit()
			if dist2 <= max(cutoff2):
				izj = izan[j]
				if izi == 6 or izj == 6:
					izx = 0
					if izi == 6:
						izx = izj
					else:
						izx = izi

					if izx == 35 or izx == 53:
						bonded = True
					elif izx == 16:
						bonded = dist2 < cutoff2[3]
					elif izx < 35:
						bonded = dist2 < cutoff2[4]
				elif (izi == 7 and izj == 8) or (izi == 8 and izj == 7):
					bonded = dist2 < cutoff2[5]
				elif (izi == 15 and izj == 8) or (izi == 8 and izj == 15):
					bonded = dist2 < cutoff2[1]
				elif izi == 1 or izj == 1:
					bonded = dist2 < cutoff2[6]
				if bonded:
					ibond[0][nbonds] = i-ibeg[imol]
					ibond[1][nbonds] = j-ibeg[imol]
					nbonds += 1
	return nbonds, ibond

def build13(ibond, nbonds, natm):
	####################################################################################
	# This function builds the following:
	#  - n13: list for # of 13-connecting atoms of each atom
	#  - dict13: dictionary for 13-connecting atoms of each atom
	#
	# called from neighbors()
	####################################################################################
	nangls, iangle = build_angle(ibond, nbonds)

	n13 = np.zeros(natm, dtype=int)
	dict13 = {}
	for l in range(nangls):
		i = iangle[0][l]
		j = iangle[2][l]
		n13[i] += 1
		n13[j] += 1
		if i in dict13:
			dict13[i].append(j)
		else:
			dict13[i] = [j]

		if j in dict13:
			dict13[j].append(i)
		else:
			dict13[j] = [i]
	return n13, dict13

def build_angle(ibond, nbonds):
	####################################################################################
	# This function builds the following by scanning ibond:
	#  - nangls: total # of 13-connecting pairs
	#  - iangle: list of all 13-connecting pairs
	#
	# called from build13()
	####################################################################################
	iangle = np.ndarray((3, 2*nbonds), dtype=int)

	nangls = 0
	for i in range(nbonds-1):
		for j in range(i+1, nbonds):
			if ibond[0][i] == ibond[0][j]:
				iangle[0][nangls] = ibond[1][i]
				iangle[1][nangls] = ibond[0][i]  # The middle atom
				iangle[2][nangls] = ibond[1][j]
				nangls += 1
			if ibond[0][i] == ibond[1][j]:
				iangle[0][nangls] = ibond[1][i]
				iangle[1][nangls] = ibond[0][i]  # The middle atom
				iangle[2][nangls] = ibond[0][j]
				nangls += 1
			if ibond[1][i] == ibond[0][j]:
				iangle[0][nangls] = ibond[0][i]
				iangle[1][nangls] = ibond[1][i]  # The middle atom
				iangle[2][nangls] = ibond[1][j]
				nangls += 1
			if ibond[1][i] == ibond[1][j]:
				iangle[0][nangls] = ibond[0][i]
				iangle[1][nangls] = ibond[1][i]  # The middle atom
				iangle[2][nangls] = ibond[0][j]
				nangls += 1
	return nangls, iangle

def matpot():
	#####################################################################################################
	# This function read in the electrostatic potential points used in the fitting, building up the 
	# matrices for LU decomposition:
	#  - a_pot: Weighted symmetric atom-esp distance matrix
	#  - b_pot: Weighted potential vector
	#
	# Other matrices include:
	#  - crd: Atom coordinate matrix
	#  - mat_pot: Original atom-esp distance matrix
	#
	# called from Main
	#####################################################################################################
	global a_pot, b_pot, crd, mat_pot, ssvpot_list, tot_nesp
	if ipol > 0:
		global n12_list, dict12_list, n13_list, dict13_list, AinvQ_list
		if ipermdip > 0:
			global AinvP_L2G_list, L2G_list

	###########################################
	# Section 1. Initiaize required variables #
	###########################################
	a_pot = np.zeros((iuniq+iuniq_p,iuniq+iuniq_p))
	b_pot = np.zeros((iuniq+iuniq_p))
	crd = np.zeros((3,iuniq))
	mat_pot = np.ndarray((0, iuniq+iuniq_p))
	ssvpot_list = np.zeros((nmol))

	if ipol > 0:
		n12_list, dict12_list = [], []
		n13_list, dict13_list = [], []
		AinvQ_list = []
		if ipermdip > 0:
			AinvP_L2G_list, L2G_list = [], []

	#######################################################################
	# Section 2. Scan the espot file and build a_pot, b_pot, crd, mat_pot #
	#######################################################################
	output_file.write("\n\n --------------------------------------")
	output_file.write("\n Atomic and ESP Coordinates Information")
	output_file.write("\n --------------------------------------")
	espot_file = open(espot, 'r')
	ioff = 0      # local variable
	max_nesp = 0  # local variable
	for imol in range(nmol):
		###################################################
		# Section 2.1. Read the metainfo of molecule imol #
		###################################################
		if ipol > 0:
			atype = []    # local variable
		line = espot_file.readline()
		natm, nesp = int(line[0:5]), int(line[5:10])

		tot_nesp += nesp
		max_nesp = max(max_nesp, nesp)
		mat_pot.resize((max_nesp, iuniq+iuniq_p), refcheck=False)

		output_file.write("\n Reading esp's for molecule %3d"%(imol+1))
		output_file.write("\n total number of atoms      = %5d"%natm)
		output_file.write("\n total number of esp points = %5d"%nesp)
		output_file.write("\n\n center       X               Y               Z")

		#################################
		# Section 2.2. Build matrix crd #
		#################################
		for i in range(natm):
			# Read the atomic coordinate (and atom type) info of molecule imol
			line = espot_file.readline().split()
			crd[0][ioff], crd[1][ioff], crd[2][ioff] = float(line[0]),float(line[1]),float(line[2])
			if ipol > 0:
				atype.append(line[4].lower())
			output_file.write("\n {:4d}{:16.7E}{:16.7E}{:16.7E}".format(i+1, crd[0][ioff], crd[1][ioff], crd[2][ioff]))
			ioff += 1
		output_file.write("\n\n")

		###########################################################################################
		# Section 2.3. Initialize dismatq, matq, (dismatdip and matp. Build AinvQ, AinvP and L2G) #
		###########################################################################################
		dismatq = np.zeros((natm))
		matq = np.zeros((natm))
		if ipol > 0:
			dismatdip = np.zeros((1, 3*natm))

			# Build n12, dict12, n13, dict13 for imol, and append into corresponding lists
			n12, dict12, n13, dict13 = neighbors(imol)
			n12_list.append(n12)
			dict12_list.append(dict12)
			n13_list.append(n13)
			dict13_list.append(dict13)

			# Considering permanent dipole, we need to build AinvQ, AinvP, L2G; without permanent dipole
			# we only need to build AinvQ. See bld_AinvQP() for the meaning of each matrix.
			if ipermdip > 0:
				npermdip = iend_p[imol]-ibeg_p[imol]+1
				matp = np.zeros((npermdip))

				AinvQ, AinvP, L2G = bld_AinvQP(imol, atype)
				AinvP_I = AinvP + np.identity(3*natm)

				AinvP_L2G_list.append(np.matmul(AinvP,L2G))
				L2G_list.append(L2G)
			else:
				AinvQ = bld_AinvQP(imol, atype)

			AinvQ_list.append(AinvQ)

		###############################################
		# Section 2.4. build a_pot, b_pot and mat_pot #
		###############################################
		wt = wtmol[imol]
		wt2 = wt*wt
		xij = np.ndarray((3)) # local variable for the distance elements between esp point i and atom j
		for i in range(nesp):
			# Read the esp point coordinate info of molecule imol
			line = espot_file.readline().split()
			if not line:
				output_file.write("\n Error: premature end of potential file")
				sys.exit()
			espi, xi, yi, zi = float(line[0]),float(line[1]),float(line[2]),float(line[3])
			ssvpot_list[imol] += wt2*espi*espi

			# Calculate the distance between esp point i and atom j
			for j in range(ibeg[imol], iend[imol]+1):
				xij[0] = xi-crd[0][j]
				xij[1] = yi-crd[1][j]
				xij[2] = zi-crd[2][j]
				rij2 = xij[0]**2 + xij[1]**2 + xij[2]**2
				rij = math.sqrt(rij2)
				rij3 = rij*rij2

				j_idx = j-ibeg[imol]

				# Calculate damping factors fe and f0 for the pGM model
				if ipol == 5 and igdm > 0:
					polj, radj = pol_dict[atype[j_idx]][0], pol_dict[atype[j_idx]][1]
					fe, ft, f0 = damp_facts(rij, 1, polj, 0, radj)
					rij /= f0
					rij3 /= fe

				# Build dismatq (and dismatdip)
				dismatq[j_idx] = 1/rij
				if ipol > 0:
					for k in range(3):
						dismatdip[0][3*j_idx+k] = xij[k]/rij3

			# Diagonal elements of a_pot (charge)
			if ipol > 0:
				dismatdip_AinvQ = np.matmul(dismatdip, AinvQ)
			for j in range(ibeg[imol], iend[imol]+1):
				j_idx = j-ibeg[imol]

				if ipol > 0:
					matq[j_idx] = dismatq[j_idx] + dismatdip_AinvQ[0][j_idx]
				else:
					matq[j_idx] = dismatq[j_idx]

				mat_pot[i][j] = matq[j_idx]
				b_pot[j] += espi*matq[j_idx]
				a_pot[j][j] += matq[j_idx]**2

			# Off-diagonal elements of a_pot (charge)
			for j in range(ibeg[imol], iend[imol]):
				j_idx = j-ibeg[imol]
				for k in range(j+1, iend[imol]+1):
					k_idx = k-ibeg[imol]

					a_pot[j][k] += matq[j_idx]*matq[k_idx]

			if ipermdip > 0:
				# Diagonal elements of a_pot (dipole)
				dismatdip_AinvP_L2G = np.matmul(np.matmul(dismatdip, AinvP_I),L2G)
				for j in range(ibeg_p[imol], iend_p[imol]+1):
					j_idx1 = j-ibeg_p[imol]
					j_idx2 = j+iuniq
					matp[j_idx1] = dismatdip_AinvP_L2G[0][j_idx1]

					mat_pot[i][j_idx2] = matp[j_idx1]
					b_pot[j_idx2] += espi*matp[j_idx1]
					a_pot[j_idx2][j_idx2] += matp[j_idx1]**2

				# Off-diagonal elements of a_pot (dipole)
				for j in range(ibeg_p[imol], iend_p[imol]):
					j_idx1 = j-ibeg_p[imol]
					j_idx2 = j+iuniq
					for k in range(j+1, iend_p[imol]+1):
						k_idx1 = k-ibeg_p[imol]
						k_idx2 = k+iuniq
	
						a_pot[j_idx2][k_idx2] += matp[j_idx1]*matp[k_idx1]

				# Off-diagonal elements of a_pot (charge * dipole)
				for j in range(ibeg[imol], iend[imol]+1):
					j_idx = j-ibeg[imol]
					for k in range(ibeg_p[imol], iend_p[imol]+1):
						k_idx1 = k-ibeg_p[imol]
						k_idx2 = k+iuniq

						a_pot[j][k_idx2] += matq[j_idx]*matp[k_idx1]

		# Apply weights on a_pot and b_pot
		for j in range(ibeg[imol], iend[imol]+1):
			# Diagonal elements of a_pot (charge)
			b_pot[j] *= wt2
			a_pot[j][j] *= wt2
			for k in range(j+1, iend[imol]+1):
				# Off-diagonal elements of a_pot (charge)
				a_pot[j][k] *= wt2

		if ipermdip > 0:
			for j in range(ibeg_p[imol], iend_p[imol]+1):
				# Diagonal elements of a_pot (dipole)
				j_idx = j+iuniq
				b_pot[j_idx] *= wt2
				a_pot[j_idx][j_idx] *= wt2
				for k in range(j+1, iend_p[imol]+1):
					# Off-diagonal elements of a_pot (dipole)
					k_idx = k+iuniq
					a_pot[j_idx][k_idx] *= wt2

			for j in range(ibeg[imol], iend[imol]+1):
				for k in range(ibeg_p[imol], iend_p[imol]+1):
					# Off-diagonal elements of a_pot (charge * dipole)
					k_idx = k+iuniq
					a_pot[j][k_idx] *= wt2

	# Symmetrize the a_pot
	for j in range(iuniq+iuniq_p-1):
		for k in range(j+1, iuniq+iuniq_p):
			a_pot[k][j] = a_pot[j][k]

	espot_file.close()

def bld_AinvQP(imol, atype):
	######################################################################################
	# This function builds the following for molecule imol:
	#  - A: Applequist matrix (Relay matrix)
	#  - Ainv: Inverse of Applequist matrix. This is the molecular polarizability matrix.
	#  - QFld: Charge field matrix
	#  - PFld: Dipole field matrix
	#  - AinvQ: Ainv * QFld
	#  - AinvP: Ainv * PFld
	#  - L2G: Matrix that converts permanent dipoles from local frames to global frames.
	#
	# called from matpot()
	######################################################################################

	##########################################################
	# Section 1. Initiaize matrices A, QFld (, PFld and L2G) #
	##########################################################
	natm = iend[imol]-ibeg[imol]+1
	A = np.zeros((3*natm, 3*natm))
	QFld = np.zeros((3*natm, natm))
	xij = np.ndarray((3)) # local variable for the distance elements between atoms i and j

	if ipermdip > 0:
		npermdip = iend_p[imol]-ibeg_p[imol]+1
		PFld = np.zeros((3*natm, 3*natm))
		L2G = np.zeros((3*natm, npermdip))

	################################################
	# Section 2. Build matrices A, QFld (and PFld) #
	################################################

	# Fill the diagonal matrices of A
	for i in range(natm):
		for j in range(3):
			idx = 3*i+j
			poli = pol_dict[atype[i]][0]
			A[idx][idx] = 1/poli

	for i in range(natm-1):
		i_idx = ibeg[imol]+i
		for j in range(i+1, natm):
			j_idx = ibeg[imol]+j

			# Calculate the distance between atoms i and j
			for k in range(3):
				xij[k] = crd[k][i_idx]-crd[k][j_idx]   # !!! from j to i

			rij2 = xij[0]**2 + xij[1]**2 + xij[2]**2
			rij = math.sqrt(rij2)
			rij3 = rij*rij2
			rij5 = rij2*rij3

			# Calculate damping factors fe and ft for selected Applequist, Thole or pGM models
			poli, polj = pol_dict[atype[i]][0], pol_dict[atype[j]][0]
			radi, radj = pol_dict[atype[i]][1], pol_dict[atype[j]][1]
			fe, ft, f0 = damp_facts(rij, poli, polj, radi, radj)

			for k in range(3):
				idx1 = 3*i+k
				idx2 = 3*j+k

				# Fill the off-diagonal vectors of QFld
				element1 = xij[k]*fe/rij3
				QFld[idx1][j] = element1
				QFld[idx2][i] = -element1

				# Fill the diagonal elements of the off-diagonal matrices of A
				element2 = fe/rij3 - 3*(xij[k]**2)*ft/rij5
				A[idx1][idx2] = element2
				A[idx2][idx1] = element2

				# Fill the diagonal elements of the off-diagonal matrices of PFld
				if ipermdip > 0:
					PFld[idx1][idx2] = -element2
					PFld[idx2][idx1] = -element2

			for k in range(2):
				for l in range(k+1, 3):
					idx1 = i*3+k
					idx2 = j*3+l
					idx3 = i*3+l
					idx4 = j*3+k

					# Fill the off-diagonal elements of the off-diagonal matrices of A
					element = -3*xij[k]*xij[l]*ft/rij5
					A[idx1][idx2] = element
					A[idx3][idx4] = element
					A[idx2][idx1] = element
					A[idx4][idx3] = element

					# Fill the off-diagonal elements of the off-diagonal matrices of PFld
					if ipermdip > 0:
						PFld[idx1][idx2] = -element
						PFld[idx3][idx4] = -element
						PFld[idx2][idx1] = -element
						PFld[idx4][idx3] = -element

	###############################
	# Section 3. Build matrix L2G #
	###############################
	if ipermdip > 0:
		p_cnt = 0
		for i in range(natm):
			# Fill the 12-connecting elements of L2G
			if i in dict12_list[imol]:
				i_idx = ibeg[imol]+i
				for j in dict12_list[imol][i]:
					j_idx = ibeg[imol]+j

					for k in range(3):
						xij[k] = crd[k][j_idx]-crd[k][i_idx]  # !!! from i to j

					rij = math.sqrt(xij[0]**2 + xij[1]**2 + xij[2]**2)

					for k in range(3):
						L2G[3*i+k][p_cnt] = xij[k]/rij
					p_cnt += 1
			if virtual > 0:
				# If virtual bond is enabled, fill the 13-connecting elements of L2G
				if i in dict13_list[imol]:
					i_idx = ibeg[imol]+i
					for j in dict13_list[imol][i]:
						j_idx = ibeg[imol]+j
	
						for k in range(3):
							xij[k] = crd[k][j_idx]-crd[k][i_idx]  # !!! from i to j
	
						rij = math.sqrt(xij[0]**2 + xij[1]**2 + xij[2]**2)
	
						for k in range(3):
							L2G[3*i+k][p_cnt] = xij[k]/rij
						p_cnt += 1

	############################################################
	# Section 4. Modify QFld and PFld based on exc12 and exc13 #
	############################################################
	if exc12 == 1:
		for i in range(natm):
			if i in dict12_list[imol]:
				for j in dict12_list[imol][i]:
					for k in range(3):
						idx1 = 3*i+k

						# Exclude 1-2 interactions in QFld
						QFld[idx1][j] = 0.0

						# Exclude 1-2 interactions in PFld
						if ipermdip > 0:
							for m in range(3):
								idx2 = 3*j+m
								PFld[idx1][idx2] = 0.0

	if exc13 == 1:
		for i in range(natm):
			if i in dict13_list[imol]:
				for j in dict13_list[imol][i]:
					for k in range(3):
						idx1 = 3*i+k

						# Exclude 1-3 interactions in QFld
						QFld[idx1][j] = 0.0

						# Exclude 1-3 interactions in PFld
						if ipermdip > 0:
							for m in range(3):
								idx2 = 3*j+m
								PFld[idx1][idx2] = 0.0

	#####################################################
	# Section 5. Finalize by building AinvQ (and AinvP) #
	#####################################################
	Ainv = np.linalg.inv(A)

	AinvQ = np.matmul(Ainv, QFld)
	if ipermdip > 0:
		AinvP = np.matmul(Ainv, PFld)
		return AinvQ, AinvP, L2G
	else:
		return AinvQ

def damp_facts(r, poli, polj, radi, radj):
	######################################################################################
	# This function calculates damping factors fe, ft (and f0, pGM model only) for Thole
	# and pGM models.
	#
	# called from bld_AinvQP() and matpot()
	######################################################################################
	fe, ft, f0 = 1, 1, 1  # Applequist model

	if ipol == 2:  # Tinker-exponential model
		coef1 = math.sqrt(radi/poli)
		coef2 = math.sqrt(radj/polj)
		v = r**3 * coef1 * coef2
		exp_v = math.exp(-v)
		fe = 1 - exp_v
		ft = fe - v * exp_v

	elif ipol == 3:  # Exponential-Thole model
		coef1 = math.sqrt(radi) / poli**(1/6)
		coef2 = math.sqrt(radj) / polj**(1/6)
		v = r * coef1 * coef2
		exp_v = math.exp(-v)
		v2 = v**2
		v3 = v * v2
		fe = 1 - (v2/2 + v + 1) * exp_v
		ft = fe - (v3/6) * exp_v

	elif ipol == 4:  # linear Thole model
		coef1 = math.sqrt(radi) * poli**(1/6)
		coef2 = math.sqrt(radj) * polj**(1/6)
		v = r / (coef1 * coef2)
		if v < 1:
			v3 = v**3
			v4 = v3*v
			fe = 4*v3 - 3*v4
			ft = v4

	elif ipol == 5:  # pGM model
		beta = 1/math.sqrt((radi**2 + radj**2)*2)  # This doesn't match pGM19 paper eq(11)
		s1 = beta * r
		if s1 < 20:  # why < 20?
			sqrt_pi = 1.7724538509055200
			s2 = s1 * s1
			es1= math.erf(s1)
			es2= s1*math.exp(-s2)/sqrt_pi
			fe = es1 - 2*es2
			ft = fe - 4*s2*es2/3
			f0 = es1

	return fe, ft, f0

def init_q0_p0():
	##################################################################
	# This function initializes charges q0 and permanent dipoles p0.
	#
	# called from Main
	##################################################################
	global q0
	
	q0 = np.ndarray((iuniq))
	if ipermdip > 0:
		global p0
		p0 = np.ndarray((iuniq_p))

	output_file.write("\n ------------------------")
	output_file.write("\n Optimization Information")
	output_file.write("\n ------------------------")
	if iqopt > 1:
		# Replace initial charges q0 and permanent dipoles p0 from qin if iqopt>1
		qin_file = open(qin, 'r')
		if ipermdip > 0:
			output_file.write("\n Since iqopt>1, %4d new q0 values and %4d new p0 values will be"%(iuniq, iuniq_p))
			output_file.write("\n read from file %s\n"%qin)
		else:
			output_file.write("\n Since iqopt>1, %4d new q0 values will be read from file"%iuniq)
			output_file.write("\n %s\n"%qin)

		q0_idx = 0
		if ipermdip > 0:
			p0_idx = 0

		line = qin_file.readline().split() # skip 1st line
		for imol in range(nmol):
			natm = iend[imol]-ibeg[imol]+1

			for i in range(natm+9):
				line = qin_file.readline().split() # skip ATOM CRD

			while line:
				q0[q0_idx] = float(line[3])
				q0_idx += 1
				line = qin_file.readline().split() # Read in replacement charges

			if ipol > 0:
				if ipermdip > 0:
					for i in range(3):
						line = qin_file.readline().split()
	
					while line:
						p0[p0_idx] = float(line[4])
						p0_idx += 1
						line = qin_file.readline().split() # Read in replacement permanent dipoles
	
					for i in range(2*natm+6):
						line = qin_file.readline() # skip PERM DIP GLOBAL and IND DIP GLOBAL
				else:
					for i in range(natm+3):
						line = qin_file.readline() # skip IND DIP GLOBAL

		qin_file.close()
	else:
		# Set initial charges (and permanent dipoles) to 0
		q0.fill(0.0)
		if ipermdip > 0:
			output_file.write("\n iqopt=1, all q0 and p0 values will be set to 0\n")
			p0.fill(0.0)
		else:
			output_file.write("\n iqopt=1, all q0 values will be set to 0\n")

def data_prep():
	##############################################################################################
	# This function setups pointers for groups of charges based on "ivary" info.
	#
	#************************************************************************
	# begin section: set lists for combined and frozen charges
	#
	# ivary[i] = 0, it is a new charge center to be fitted
	# ivary[i] =+n, it is a charge center to be fitted with center n (center n must be a previous 
	#                    center entered with ivary[n] = 0
	# ivary[i] =-n, it is a frozen charge center to be kept at q0[i]
	#
	# Similarly:
	#
	# ivary_p[i] = 0, it is a new permanent dipole center to be fitted center n (center n must be a
	# ivary_p[i] =+n, it is a permanent dipole center to be fitted with 
	#                    previous center entered with ivary_p[n] = 0
	# ivary_p[i] =-n, it is a frozen permanent dipole center to be kept at p0[i]
	#************************************************************************
	#
	# called from Main
	###############################################################################################
	global iqpcntr, nqpl

	nqp = 0
	iqpcntr = np.zeros((iuniq+iuniq_p+nlgrng), dtype=int)

	# Fill in charge equivalence information in iqpcntr
	for i in range(iuniq):
		if ivary[i] == 0:
			nqp += 1
			iqpcntr[i] = nqp
		elif ivary[i] > 0:
			iqpcntr[i] = iqpcntr[ivary[i]-1]

			if iqpcntr[i] > nqp:
				output_file.write("\n Error: data_prep() charge equivalence input is screwy")
				sys.exit()
		else:
			iqpcntr[i] = -1

	if nqp == 0:
		output_file.write("\n Warning: ALL charges are frozen!!!")
	else:
		output_file.write("\n Number of unique UNfrozen charge centers = %5d"%nqp)
	nq = nqp

	# Fill in permanent dipole equivalence information in iqpcntr
	if ipermdip > 0:
		for i in range(iuniq_p):
			if ivary_p[i] == 0:
				nqp += 1
				iqpcntr[iuniq+i] = nqp
			elif ivary_p[i] > 0:
				iqpcntr[iuniq+i] = iqpcntr[iuniq+ivary_p[i]-1]

				if iqpcntr[iuniq+i] > nqp:
					output_file.write("\n Error: data_prep() permanent dipole equivalence input is screwy")
					sys.exit()
			else:
				iqpcntr[iuniq+i] = -1

		if (nqp-nq) == 0:
			output_file.write("\n Warning: ALL permanent dipoles are frozen!!!")
		else:
			output_file.write("\n Number of unique UNfrozen permanent dipoles = %5d"%(nqp-nq))

	# Fill in Lagrange constraints information in iqpcntr
	for i in range(nlgrng):
		iqpcntr[iuniq+iuniq_p+i] = nqp+i+1

	# Set nqpl to the total # of row elements (charges/permanent dipoles to be independantly fit + constraints)
	# in the fitting matrix
	nqpl = nqp + nlgrng

def charge_opt():
	#########################################################################
	# This function is the driver for charge (and permanent dipole) fitting
	#
	# called from Main
	#########################################################################
	global irstrnt, awork, bwork

	qold = np.zeros((iuniq))   # local variable
	irsave = 0                 # local variable
	nitern = 0                 # local variable

	# criteria for convergence & maximum iterations for the non-linear optimizations
	qtol = 0.000001
	maxit = 42

	# Only on first pass through this function (indicated by nitern= 0), if irstrnt > 0, transfer irstrnt to irsave
	# and reset irstrnt to 0, in order to get an initial guess using a harmonic constraint. This is done so the
	# restraint function rstran() will use a harmonic restraint.
	if irstrnt > 0:
		irsave = irstrnt
		irstrnt = 0
		output_file.write("\n\n Non-linear optimization requested.")

	# Now go do a "harmonic restraint" run, restraint= qwt(qcal[i]-q0[i])**2
	# -- loop to convergence
	while nitern < maxit:
		matbld()
	
		# Solve (Ax = b) where A and b are input, x is solution awork x = bwork          
		# 
		# the solution "x" is returned in "b" (bwork)
		# 
		# -- condition the matrix diagonal to avoid DGETRF() detectingNN (wrapped in lu_factor) singularity 
		for jn in range(nqpl):
			if abs(awork[jn][jn]) < 1.0E-10:
				awork[jn][jn] = 1.0E-10
	
		awork, piv = lu_factor(awork, overwrite_a=True)
		bwork = lu_solve((awork, piv), bwork, overwrite_b=True)
	
		# -- copy solution vector "bwork" to 'calculated parameters' qcal and pcal
		for k in range(iuniq):
			icntr = iqpcntr[k]
			if icntr >= 1:
				# -- new charge
				qcal[k] = bwork[icntr-1]
			else:
				# -- frozen charge
				qcal[k] = q0[k]

		for k in range(iuniq_p):
			icntr = iqpcntr[iuniq+k]
			if icntr >= 1:
				# -- new permanent dipole
				pcal[k] = bwork[icntr-1]
			else:
				# -- frozen permanent dipole
				pcal[k] = p0[k]
	
		# -- A quick check from rstrn: if irstrnt is now negative, there are no restraints because no qwtval(i) > 
		#    0.1e-10, so reset irsave= 0
		if irstrnt < 0:
			irsave = 0
			output_file.write("\n\n WARNING: Restraints were requested, but the restraint weights were all zero\n")
	
		# -- We're finished if it's only a "harmonic restraint" run. But if it's a non-linear optimization (irsave>0),
		#    we've only just begun (i.e. we have our initial guess)
		if irsave <= 0:
			return
		else:
			# -- It's a non-linear optimization: reset irstrnt (to now calculate the proper non-linear restraint
			#    derivatives in routine rstran)
			irstrnt = irsave
	
		# -- Begin iterative optimization loop with comparison of old & new charges; calculate the convergence and
		#    replace the old charges with the new
		qchnge = 0.0
		for i in range(iuniq):
			qdiff = qcal[i] - qold[i]
			qchnge += qdiff*qdiff
			qold[i] = qcal[i]
	
		qchnge = math.sqrt(qchnge)/iuniq
		output_file.write("\n  qchnge ={:20.10E}".format(qchnge))
	
		# -- If this is less than qtol then we're done
		if qchnge < qtol and nitern > 1:
			output_file.write("\n\n Convergence in%5d iterations\n"%nitern)
			return

		# Loop again
		nitern += 1

	output_file.write("\n Warning: after %5d iterations, no convergence!\n"%maxit)

def matbld():
	##########################################################################################################
	# This function build up matrices for LU decomposition:
	#
	#   stage 1: copy weighted matrices a_pot and b_pot to work arrays awork and bwork (which are destroyed in
	#            the LU decomp & back subst)
	#
	#   stage 2: if charge (and permanent dipole) are to be restrained, modify awork and bwork appropriately  
	#
	# called from charge_opt()
	##########################################################################################################
	global a, b, awork, bwork

	a = np.zeros((iuniq+iuniq_p+nlgrng,iuniq+iuniq_p+nlgrng))
	b = np.zeros((iuniq+iuniq_p+nlgrng))

	for k in range(iuniq+iuniq_p):
		b[k] = b_pot[k]
		for j in range(iuniq+iuniq_p):
			a[j][k] = a_pot[j][k]

	# Fill in the final columns & rows of A with the Lagrange constraints which keep the charge on groups of 
	# atoms to a constant
	for i in range(nlgrng):
		b[iuniq+iuniq_p+i] = grpchg[i]
		for j in range(iuniq+iuniq_p+nlgrng):
			a[iuniq+iuniq_p+i][j] = lgrcnt[i][j]
			a[j][iuniq+iuniq_p+i] = lgrcnt[i][j]

	# Add restraint to initial charge q0 and permanent dipole p0:
	rstran()

	# Build awork and bwork based on "combined and frozen centers" info:
	#
	# 1) Frozen centers do not appear in the matrix of fitted charges (or permanent dipoles)
	# 2) Combined centers appear as one single charge (or permanent dipole) center for fitting
	#
	# First, since we accumulate values, zero out awork & bwork up to nqpl (the independant + contraint number):
	awork = np.zeros((nqpl,nqpl))
	bwork = np.zeros((nqpl))

	# Loop over all centers, building awork & bwork from A and B based on iqpcntr: for each center, iqpcntr[i]
	# dictates which of the fitted charges (or permanent dipoles) it is and therefore where it goes in the matrices.
	# If iqpcntr[j] < 1, this center is a frozen charge (or permanent dipole) and it is skipped as far as forming a
	# row in awork, and its esp contribution is subtracted from bwork to take care of it's awork jth column-element
	# for each i.
	for i in range(iuniq+iuniq_p+nlgrng):
		icntr = iqpcntr[i]
		if icntr > 0:
			# i is active
			bwork[icntr-1] += b[i]
			for j in range(iuniq+iuniq_p+nlgrng):
				jcntr = iqpcntr[j]
				if jcntr > 0:
					# j is active
					awork[icntr-1][jcntr-1] += a[i][j]
				else:
					if j < iuniq:
						# j is a frozen charge
						bwork[icntr-1] -= q0[j]*a[i][j]
					else:
						# j is a frozen permanent dipole
						bwork[icntr-1] -= p0[j-iuniq]*a[i][j]

def rstran():
	##########################################################################################################
	# This function assigns the retraint weights to the diagonal of A and to B.
	#
	#----------------------------------------------------------------------
	# Two kinds of restraints are available:
	#
	# a) Harmonic restraint to the initial charges (and permanent dipoles). Fine as long as there aren't any
	#  large charges (and permanent dipoles) that SHOULD be large... these really feel a strong force if they
	#  are restrained to a low value.
	#
	# b) Hyperbolic restraint to charges (and permanent dipoles) towards 0. This gets asymptotic at "large"
	#  values, so "large" charges (and permanent dipoles) aren't pulled down any stronger than some (reasonable)
	#  limiting force.  This is a non-linear weighting function, so the fit procedure is iterative.
	#
	# Other options for restraints to initial charges (and permanent dipoles):
	# if requested, restrain the charges (and permanent dipoles) by modifying the sum-of-squares cost function
	# derivative. The scheme for doing this is as follows:
	#
	# If control variable ihfree > 0, let hydrogen charges (and permanent dipoles) free (i.e. reset their qwtval
	# and pwtval to 0.0).
	#-----------------------------------------------------------------------
	#
	# called from matbld()
	##########################################################################################################
	global irstrnt
	qwtval.fill(qwt)
	if ipermdip > 0:
		pwtval.fill(pwt)

	for i in range(iuniq):
		if ihfree > 0 and izan[i] == 1:
			qwtval[i] = 0.0

		# use of harmonic restraints
		if irstrnt == 0:
			a[i][i] += qwtval[i]
			b[i] += qwtval[i]*q0[i] # q0 has the initial and/or frozen charge
		# use of hyperbolic restraints
		elif irstrnt > 0 and qwtval[i] > 0.1E-10:
			qwtval[i] = qwt/math.sqrt(qcal[i]*qcal[i] + 0.01) # qcal has the current (calculated) charge
			a[i][i] += qwtval[i]

	for i in range(iuniq_p):
		i_idx = iuniq+i
		if ihfree > 0 and izan_p[i] == 1:
			pwtval[i] = 0.0

		# use of harmonic restraints
		if irstrnt == 0:
			a[i_idx][i_idx] += pwtval[i]
			b[i_idx] += pwtval[i]*p0[i] # p0 has the initial and/or frozen dipole
		# use of hyperbolic restraints
		elif irstrnt > 0 and pwtval[i] > 0.1E-10:
			pwtval[i] = pwt/math.sqrt(pcal[i]*pcal[i] + 0.01) # pcal has the current (calculated) dipole
			a[i_idx][i_idx] += pwtval[i]

	# If all qwtval[i] and pwtval[i] are 0.0, no restraints so reset irstrnt= -1
	for i in range(iuniq):
		if qwtval[i] > 0.1E-10:
			return
	for i in range(iuniq_p):
		if pwtval[i] > 0.1E-10:
			return

	irstrnt = -1

def calc_dip():
	#############################################################################
	# This function calculates the dipole moments, including
	#  - dipol_mol: Molecular dipole components
	#  - dipmom_mol: Molecular dipole magnitude
	#  - dipind: Atomic induced dipole components (global frame)
	#  - dipperm: Atomic permanent dipole components (global frame)
	#  - dipindperm: Atomic induced + permanent dipole components (global frame)
	#
	# called from Main
	#############################################################################
	global dipol_mol, dipmom_mol, dipind, dipperm, dipindperm

	dipol_mol = np.zeros((3,nmol))
	dipmom_mol = np.zeros((nmol))
	dipind = np.zeros((3,iuniq))
	dipperm = np.zeros((3,iuniq))
	dipindperm = np.zeros((3,iuniq))

	# Calculate atomic induced and permanent dipoles (global frame)
	if ipol > 0:
		for imol in range(nmol):
			natm = iend[imol]-ibeg[imol]+1
			if ipermdip > 0:
				npermdip = iend_p[imol]-ibeg_p[imol]+1

			for i in range(natm):
				i_idx = ibeg[imol]+i

				for j in range(natm):
					j_idx = ibeg[imol]+j
					for k in range(3):
						# The induced dipole due to permanent charges
						dipind[k][i_idx] += qcal[j_idx]*AinvQ_list[imol][3*i+k][j]

				if ipermdip > 0:
					for j in range(npermdip):
						j_idx = ibeg_p[imol]+j
						for k in range(3):
							# The induced dipole due to permanent dipoles
							dipind[k][i_idx] += pcal[j_idx]*AinvP_L2G_list[imol][3*i+k][j]

							# The permanent dipole
							dipperm[k][i_idx] += pcal[j_idx]*L2G_list[imol][3*i+k][j]

	# Calculate atomic induced + permanent dipole and molecular dipoles
	for imol in range(nmol):
		for i in range(ibeg[imol], iend[imol]+1):
			for j in range(3):
				dipindperm[j][i] = dipind[j][i] + dipperm[j][i]
				dipol_mol[j][imol] += qcal[i]*crd[j][i] + dipindperm[j][i]

		dipmom_mol[imol] = math.sqrt(dipol_mol[0][imol]**2 + dipol_mol[1][imol]**2 + dipol_mol[2][imol]**2)

	# Convert molecular dipole from a.u. to debyes
	dipol_mol *= au2D
	dipmom_mol *= au2D

def reornt():
	############################################################################################
	# This function translates molecule to center of mass and reorients it along principal axes
	# of rotation in preparation for dipole and quadrupole moment calculation.
	#
	# called from Main
	############################################################################################

	# ---- ATOMIC WEIGHT ARRAY FOR CENTER OF MASS ----
	#
	#  20 elements were originally handled: H(1) - Ca(20)
	# 103 elements are now considered:      H(1) - Lr(103)
	#
	# Atomic weight from The Merck Index - Thirteeth edition, Merck & Co., INC., Whitehouse Station,
	# NJ, 2001
	#
	# F.-Y. Dupradeau & P. Cieplak, http://q4md-forcefieldtools.org/
	global co, cmas_mol, dipol_mol_com, dipindperm_com

	cmas_mol = np.ndarray((3,nmol))
	dipol_mol_com = np.zeros((3,nmol))
	co = np.ndarray((3,iuniq))
	dipindperm_com = np.zeros((3,iuniq))

	wt = [1.0079,4.0026,6.9410,9.0122,10.8110,12.0107,14.0067,15.9994,18.9984,20.1797,
	22.9898,24.3050,26.9815,28.0855,30.9738,32.0650,35.4530,39.9480,39.0983,40.0780,
	44.9559,47.8670,50.9415,51.9961,54.9380,55.8450,58.9332,58.6934,63.5460,65.3900,
	69.7230,72.6400,74.9216,78.9600,79.9040,83.8000,85.4678,87.6200,88.9058,91.2240,
	92.9064,95.9400,97.9072,101.0700,102.9055,106.4200,107.8682,112.4110,114.8180,118.7100,
	121.7600,127.6000,126.9045,131.2930,132.9054,137.3270,138.9055,140.1160,140.9076,144.2400,
	144.9127,150.3600,151.9640,157.2500,158.9253,162.5000,164.9303,167.2590,168.9342,173.0400,
	174.9670,178.4900,180.9479,183.8400,186.2070,190.2300,192.2170,195.0780,196.9665,200.5900,
	204.3833,207.2000,208.9804,208.9824,209.9871,222.0176,223.0197,226.0254,227.0277,232.0381,
	231.0359,238.0289,237.0482,244.0642,243.0614,247.0704,247.0703,251.0796,252.0830,257.0951,
	258.0984,259.1010,262.1097]

	for imol in range(nmol):
		# ----- CALCULATE THE CENTER OF MASS -----
		xc, yc, zc = cmass(wt, imol)
		cmas_mol[0][imol] = xc
		cmas_mol[1][imol] = yc
		cmas_mol[2][imol] = zc
	
		# ----- MOVE THE ORIGIN TO CENTER OF MASS ----- 
		cmove(xc, yc, zc, imol)
	
		# ----- CALCULATE THE MOMENT OF INERTIA -----
		s = momin(wt, imol)
	
		# ----- ROTATE ALONG THE PRINCIPLE AXES OF MOMENT OF INERTIA -----
		mominrot(s, co, co, imol)
		mominrot(s, dipol_mol, dipol_mol_com, imol)
	
		if ipol > 0:
			mominrot(s, dipindperm, dipindperm_com, imol)
	
	# Convert center of mass from bohr to angstroms
	cmas_mol *= au2A

def cmass(wt, imol):
	#################################################################
	# This function calculates the center of mass of the molecule.
	#
	# called from reornt()
	#################################################################
	sumx, sumy, sumz, Sum = 0.0, 0.0, 0.0, 0.0
	for i in range(ibeg[imol], iend[imol]+1):
		idx = izan[i]-1
		sumx += crd[0][i] * wt[idx]
		sumy += crd[1][i] * wt[idx]
		sumz += crd[2][i] * wt[idx]
		Sum += wt[idx]
	return sumx/Sum, sumy/Sum, sumz/Sum

def cmove(xc, yc, zc, imol):
	#################################################################
	# This function moves the origin to the center of mass.
	#
	# called from reornt()
	#################################################################
	for i in range(ibeg[imol], iend[imol]+1):
		co[0][i] = crd[0][i] - xc
		co[1][i] = crd[1][i] - yc
		co[2][i] = crd[2][i] - zc

def momin(wt, imol):
	################################################################################################
	# This function calculates the moments of inertia tensor and the principal axes of moment of
	# inertia of the molecule. It then reorients the molecule along the principal axes (with the 
	# origin at the center of mass) in preparation for calculation of quadrapole moment components.
	#  - ain: The moments of inertia tensor
	#  - s: The principal axes of moment of inertia matrix
	#
	# called from reornt()
	################################################################################################
	ain = np.ndarray((3,3))
	sxx, syy, szz, sxy, sxz, syz = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	for i in range(ibeg[imol], iend[imol]+1):
		idx = izan[i]-1
		xx = co[0][i]**2
		yy = co[1][i]**2
		zz = co[2][i]**2
		sxx += wt[idx]*(yy+zz)
		syy += wt[idx]*(xx+zz)
		szz += wt[idx]*(xx+yy)
		sxy += wt[idx]*co[0][i]*co[1][i]
		sxz += wt[idx]*co[0][i]*co[2][i]
		syz += wt[idx]*co[1][i]*co[2][i]
	ain[0][0] = sxx
	ain[1][1] = syy
	ain[2][2] = szz
	ain[0][1] = -sxy
	ain[0][2] = -sxz
	ain[1][2] = -syz
	ain[1][0] = ain[0][1]
	ain[2][0] = ain[0][2]
	ain[2][1] = ain[1][2]

	#  ----- CALCULATE PRINCIPAL AXES OF MOMENT OF INERTIA -----
	eigvals, s = np.linalg.eigh(ain)
	return s

def mominrot(s, vecin, vecout, imol):
	########################################################################################
	# This function rotate the following along the principle axes of moment of inertia.
	#  - dipol_mol: Molecular dipole components
	#  - dipindperm: Atomic induced + permanent dipole components
	#  - co: Atomic coordinates
	#
	# called from reornt()
	########################################################################################
	if len(vecin[0]) == nmol:
		# Rotates dipol_mol
		xs = vecin[0][imol]*s[0][0] + vecin[1][imol]*s[1][0] + vecin[2][imol]*s[2][0]
		ys = vecin[0][imol]*s[0][1] + vecin[1][imol]*s[1][1] + vecin[2][imol]*s[2][1]
		zs = vecin[0][imol]*s[0][2] + vecin[1][imol]*s[1][2] + vecin[2][imol]*s[2][2]
		vecout[0][imol] = xs
		vecout[1][imol] = ys
		vecout[2][imol] = zs
	else:
		# Rotates dipindperm or co
		for i in range(ibeg[imol], iend[imol]+1):
			xs = vecin[0][i]*s[0][0] + vecin[1][i]*s[1][0] + vecin[2][i]*s[2][0]
			ys = vecin[0][i]*s[0][1] + vecin[1][i]*s[1][1] + vecin[2][i]*s[2][1]
			zs = vecin[0][i]*s[0][2] + vecin[1][i]*s[1][2] + vecin[2][i]*s[2][2]
			vecout[0][i] = xs
			vecout[1][i] = ys
			vecout[2][i] = zs

def calc_quad(coord, dip):
	############################################################################################
	# This function calculates the molecular quadrupole moments (before or after reorientation).
	#
	# called from Main
	############################################################################################
	quad = np.zeros((6,nmol))
	for imol in range(nmol):
		for i in range(ibeg[imol], iend[imol]+1):
			# Calculates quadrupole contributed by permanent charge
			qqxx = 3.0*qcal[i]*coord[0][i]*coord[0][i]
			qqyy = 3.0*qcal[i]*coord[1][i]*coord[1][i]
			qqzz = 3.0*qcal[i]*coord[2][i]*coord[2][i]
			qqiso = (qqxx + qqyy + qqzz)/3
			quad[0][imol] += qqxx - qqiso # XX
			quad[1][imol] += qqyy - qqiso # YY
			quad[2][imol] += qqzz - qqiso # ZZ
			quad[3][imol] += 3.0*qcal[i]*coord[0][i]*coord[1][i] # XY
			quad[4][imol] += 3.0*qcal[i]*coord[0][i]*coord[2][i] # XZ
			quad[5][imol] += 3.0*qcal[i]*coord[1][i]*coord[2][i] # YZ

			if ipol > 0:
				# Calculates quadrupole contributed by permanent and induced dipoles
				qdxx = 6.0*dip[0][i]*coord[0][i]
				qdyy = 6.0*dip[1][i]*coord[1][i]
				qdzz = 6.0*dip[2][i]*coord[2][i]
				qdiso = (qdxx + qdyy + qdzz)/3
				quad[0][imol] += qdxx - qdiso # XX
				quad[1][imol] += qdyy - qdiso # YY
				quad[2][imol] += qdzz - qdiso # ZZ
				quad[3][imol] += 3.0*(dip[0][i]*coord[1][i] + dip[1][i]*coord[0][i]) # XY
				quad[4][imol] += 3.0*(dip[0][i]*coord[2][i] + dip[2][i]*coord[0][i]) # XZ
				quad[5][imol] += 3.0*(dip[1][i]*coord[2][i] + dip[2][i]*coord[1][i]) # YZ

	# Convert quadrupoles from a.u. to debye*angstroms
	quad *= au2D*au2A
	return quad

def wrt_qout():
	########################################################################################
	# This function writes out the qout file. The follwoing sections appear molecule-wisely.
	#
	# called from Main
	########################################################################################
	qout_file = open(qout, 'w')
	qout_file.write("All values are reported in atomic units\n")

	ioff = 0   # local variable
	for imol in range(nmol):
		# 1. Title
		qout_file.write("%FLAG TITLE\n")
		qout_file.write(" molecule %d: %s"%(imol+1, title[imol]))

		# 2. Atomic coordinates
		qout_file.write("\n%FLAG ATOM CRD: I4,3E16.7\n")
		qout_file.write(" atm.no      X               Y               Z\n")
		for i in range(ibeg[imol], iend[imol]+1):
			qout_file.write("{:4d}{:16.7E}{:16.7E}{:16.7E}\n".format(i+1, crd[0][i], crd[1][i], crd[2][i]))

		# 3. Atomic charges
		qout_file.write("\n%FLAG ATOM CHRG: 2(I4,X7),I4,X2,E16.7\n")
		qout_file.write(" atm.no   element.no   ivary      q(opt)\n")
		for i in range(ibeg[imol], iend[imol]+1):
			qout_file.write("{:4d}       {:4d}       {:4d}  {:16.7E}\n".format(i+1, izan[i], ivary[i], qcal[i]))

		if ipol > 0:
			if ipermdip > 0:
				# 4. Permanent dipoles (local frame)
				qout_file.write("\n%FLAG PERM DIP LOCAL: 3(I4,X5),I4,X2,E16.7\n")
				qout_file.write(" dip.no   atm.no   ref.no   ivary      p(opt)\n")
				natm = iend[imol]-ibeg[imol]+1
				for atm1 in range(natm):
					if atm1 in dict12_list[imol]:
						for atm2 in dict12_list[imol][atm1]:
							qout_file.write("{:4d}     {:4d}     {:4d}     {:4d}  {:16.7E}\n".format(ioff+1, atm1+1, atm2+1, ivary_p[ioff], pcal[ioff]))
							ioff += 1
					if virtual > 0:
						if atm1 in dict13_list[imol]:
							for atm2 in dict13_list[imol][atm1]:
								qout_file.write("{:4d}     {:4d}     {:4d}     {:4d}  {:16.7E}\n".format(ioff+1, atm1+1, atm2+1, ivary_p[ioff], pcal[ioff]))
								ioff += 1

				# 5. Permanent dipoles (global frame)
				qout_file.write("\n%FLAG PERM DIP GLOBAL: I4,3E16.7\n")
				qout_file.write(" atm.no      X               Y               Z\n")
				for i in range(ibeg[imol], iend[imol]+1):
					qout_file.write("{:4d}{:16.7E}{:16.7E}{:16.7E}\n".format(i+1, dipperm[0][i], dipperm[1][i], dipperm[2][i]))

			# 6. Induced dipoles (global frame)
			qout_file.write("\n%FLAG IND DIP GLOBAL: I4,3E16.7\n")
			qout_file.write(" atm.no      X               Y               Z\n")
			for i in range(ibeg[imol], iend[imol]+1):
				qout_file.write("{:4d}{:16.7E}{:16.7E}{:16.7E}\n".format(i+1, dipind[0][i], dipind[1][i], dipind[2][i]))
		qout_file.write("\n")

	qout_file.close()

def evlchi():
	###############################################################################################################
	# This function valuates chi-square for linear function espclci = sum_j(qj*termij), where j = number of terms,
	# qj is the coefficient to termij, and chi-square is the merit function: chi-square = sum_i((espi-espclci)**2),
	# where i is the number of data points for which esp is known.
	#
	# called from wrt_out()
	###############################################################################################################
	espot_file = open(espot, 'r')
	if ioutopt == 1:
		esout_file = open(esout, 'w')

	chipot = 0.0
	for imol in range(nmol):
		# Read the metainfo of molecule imol
		line = espot_file.readline()
		natm, nesp = int(line[0:5]), int(line[5:10])
		if ioutopt == 1:
			esout_file.write("molecule %d: %s"%(imol+1, title[imol]))
			esout_file.write("nesp: %d   natm: %d\n"%(nesp, natm))
			esout_file.write("      X         Y         Z        esp_qm      esp_clc       diff\n")
	
		# Skip the atomic coordinate and atom type info of molecule imol
		for i in range(natm):
			line = espot_file.readline()
	
		chipot_mol = 0.0
		for i in range(nesp):
			# Read in the esp points used in the fitting
			line = espot_file.readline().split()
			if not line:
				output_file.write("\n Error: unexpected eof in %s"%espot)
				sys.exit()
			espqmi, xi, yi, zi = float(line[0]),float(line[1]),float(line[2]),float(line[3])

			espclc = 0.0
			# Calculate esp contributed from fitted charges
			for j in range(ibeg[imol], iend[imol]+1):
				espclc += mat_pot[i][j]*qcal[j]

			# Calculate esp contributed from fitted permanent dipoles
			if ipermdip > 0:
				for j in range(ibeg_p[imol], iend_p[imol]+1):
					espclc += mat_pot[i][iuniq+j]*pcal[j]
	
			# Caculate the esp residual
			vresid = espqmi - espclc
			chipot_mol += vresid**2

			# Write the coords, qm esps, calculated esps and residuals in the esout file
			if ioutopt == 1:
				esout_file.write("%10.5f%10.5f%10.5f%12.5f%12.5f%12.5f\n"%(xi,yi,zi,espqmi,espclc,vresid))

		chipot += chipot_mol

		if ioutopt == 1:
			rmse_mol = math.sqrt(chipot_mol/nesp)
			rrmse_mol = math.sqrt(chipot_mol/ssvpot_list[imol])
			esout_file.write("molecule {:d}  rss={:13.7E} rmse={:13.7E} rrmse={:13.7E}\n\n".format(imol+1,chipot_mol,rmse_mol,rrmse_mol))

	espot_file.close()
	if ioutopt == 1:
		esout_file.close()

	rmse = math.sqrt(chipot/tot_nesp)
	rrmse = math.sqrt(chipot/sum(ssvpot_list))

	return chipot, rmse, rrmse

def wrt_quad(quad_mol, prncpl):
	if prncpl:
		for imol in range(nmol):
			quad_orig = np.ndarray((3,3))
			quad_orig[0][0] = quad_mol[0][imol]
			quad_orig[0][1] = quad_mol[3][imol]
			quad_orig[0][2] = quad_mol[4][imol]
			quad_orig[1][0] = quad_mol[3][imol]
			quad_orig[1][1] = quad_mol[1][imol]
			quad_orig[1][2] = quad_mol[5][imol]
			quad_orig[2][0] = quad_mol[4][imol]
			quad_orig[2][1] = quad_mol[5][imol]
			quad_orig[2][2] = quad_mol[2][imol]
			quad_diag, s = np.linalg.eigh(quad_orig)
			quad_diag = sorted(quad_diag, reverse=True)
			output_file.write("\n #MOL          X          Y          Z")
			if iquad == 0:
				output_file.write("\n %3d    X %10.5f"%(imol+1, quad_diag[0]))
				output_file.write("\n        Y %10.5f %10.5f"%(0.0, quad_diag[1]))
				output_file.write("\n        Z %10.5f %10.5f %10.5f"%(0.0, 0.0, quad_diag[2]))
			else:
				output_file.write("\n %3d    X %10.5f"%(imol+1, quad_diag[0]/3))
				output_file.write("\n        Y %10.5f %10.5f"%(0.0, quad_diag[1]/3))
				output_file.write("\n        Z %10.5f %10.5f %10.5f"%(0.0, 0.0, quad_diag[2]/3))
	else:
		for imol in range(nmol):
			output_file.write("\n #MOL          X          Y          Z")
			if iquad == 0:
				output_file.write("\n %3d    X %10.5f"%(imol+1,quad_mol[0][imol]))
				output_file.write("\n        Y %10.5f %10.5f"%(quad_mol[3][imol],quad_mol[1][imol]))
				output_file.write("\n        Z %10.5f %10.5f %10.5f"%(quad_mol[4][imol], quad_mol[5][imol], quad_mol[2][imol]))
			else:
				output_file.write("\n %3d    X %10.5f"%(imol+1,quad_mol[0][imol]/3))
				output_file.write("\n        Y %10.5f %10.5f"%(quad_mol[3][imol]/3,quad_mol[1][imol]/3))
				output_file.write("\n        Z %10.5f %10.5f %10.5f"%(quad_mol[4][imol]/3, quad_mol[5][imol]/3, quad_mol[2][imol]/3))

def wrt_out():
	########################################################################################
	# This function writes out the fitting results and statistics into the output file.
	#
	# called from Main
	########################################################################################
	global rmse, rrmse

	output_file.write("\n -----------------------")
	output_file.write("\n Fitting Results Summary")
	output_file.write("\n -----------------------\n")

	for imol in range(nmol):
		output_file.write(" molecule %d: %s"%(imol+1, title[imol]))

	# Print the charges
	output_file.write("\n          Point Charges Before & After Optimization")
	output_file.write("\n atm.no   element.no   q(init)        q(opt)     ivary    d(rstr)/dq")

	icnt, jcnt = 0, 0
	chge = 0.0
	for j in range(iuniq):
		output_file.write("\n %4d      %4d      %10.6f    %10.6f%7d%15.6f"%(j+1,izan[j],q0[j],qcal[j],ivary[j],qwtval[j]))
		chge += qcal[j]
		jcnt += 1
		if (jcnt > iend[icnt]):
			output_file.write("\n")
			icnt += 1

	output_file.write("\n Sum over the calculated charges: %10.3f\n"%chge)

	# Print the permanent dipoles
	if ipermdip > 0:
		output_file.write("\n       Permanent Dipoles Before & After Optimization")
		output_file.write("\n dip.no   element.no   p(init)        p(opt)     ivary    d(rstr)/dp")

		icnt, jcnt = 0, 0
		for j in range(iuniq_p):
			output_file.write("\n %4d      %4d      %10.6f    %10.6f%7d%15.6f"%(j+1,izan_p[j],p0[j],pcal[j],ivary_p[j],pwtval[j]))
			jcnt += 1
			if (jcnt > iend_p[icnt]):
				output_file.write("\n")
				icnt += 1

	# Calculate statistics for the esp's
	chipot,rmse,rrmse = evlchi()

	# Write statistics out
	output_file.write("\n --------------------------")
	output_file.write("\n Fitting Statistics Summary")
	output_file.write("\n --------------------------")
	output_file.write("\n The initial sum of squares  (ssvpot = sum_i(espi**2))          {:13.7E}".format(sum(ssvpot_list)))
	output_file.write("\n The residual sum of squares (RSS = sum_i((espi-espclci)**2))   {:13.7E}".format(chipot))
	output_file.write("\n The root-mean-squared error (RMSE = sqrt(RSS/N))               {:13.7E}".format(rmse))
	output_file.write("\n The relative RMSE           (RRMSE = sqrt(RSS/ssvpot))         {:13.7E}".format(rrmse))

	# Print the center of mass, (original and reoriented) molecular dipole and quadrupole
	output_file.write("\n\n ---------------------------")
	output_file.write("\n Molecular Multipole Summary")
	output_file.write("\n ---------------------------")
	if ireornt > 0:
		output_file.write("\n Center of Mass (Angst.):")
		for imol in range(nmol):
			output_file.write("\n #MOL          X          Y          Z")
			output_file.write("\n %3d     %10.5f %10.5f %10.5f"%(imol+1,cmas_mol[0][imol], cmas_mol[1][imol], cmas_mol[2][imol]))
		output_file.write("\n\n Dipole Moments Reoriented (Debye):")
		for imol in range(nmol):
			output_file.write("\n #MOL         D          Dx         Dy         Dz")
			output_file.write("\n %3d     %10.5f %10.5f %10.5f %10.5f"%(imol+1,dipmom_mol[imol],dipol_mol_com[0][imol],dipol_mol_com[1][imol],dipol_mol_com[2][imol]))
		output_file.write("\n\n Traceless Quadrupole Moments Reoriented (Debye*Angst.):")
		wrt_quad(quad_mol_com, False)
		output_file.write("\n\n Traceless Quadrupole Moments in Principal Axes (Debye*Angst.):")
		wrt_quad(quad_mol_com, True)
	else:
		output_file.write("\n Dipole Moments (Debye):")
		for imol in range(nmol):
			output_file.write("\n #MOL         D          Dx         Dy         Dz")
			output_file.write("\n %3d     %10.5f %10.5f %10.5f %10.5f"%(imol+1,dipmom_mol[imol],dipol_mol[0][imol],dipol_mol[1][imol],dipol_mol[2][imol]))
		output_file.write("\n\n Traceless Quadrupole Moments (Debye*Angst.):")
		wrt_quad(quad_mol, False)
		output_file.write("\n\n Traceless Quadrupole Moments in Principal Axes (Debye*Angst.):")
		wrt_quad(quad_mol, True)

#-----------------------------------------------------------------------
# The beginning of the main program
#-----------------------------------------------------------------------

print("Reading input file ...")
###### Get the file names ######
Input,output,qin,polariz,qout,espot,esout = file_in()

input_file = open(Input, 'r')
output_file = open(output, 'w')

###### Read the control parameters ######
read_in()

# If nmol > 1, this is a multiple molecule run;
# otherwise it is a single molecule run
if nmol > 1:
	mult_mol()
else:
	sing_mol()

input_file.close()

if ipol > 0:
	read_pol_dict()

print("Reading espot file and building matrices ...")
###### Read in the qm esp, forming the matrices a_pot and b_pot
matpot()

###### Initialize q0 and p0 according to iqopt
init_q0_p0()

###### Process the input (freezing, equivalencing charges and permanent dipoles)
data_prep()

qcal = np.zeros((iuniq))
qwtval = np.zeros((iuniq))
if ipermdip > 0:
	pcal = np.zeros((iuniq_p))
	pwtval = np.zeros((iuniq_p))

if irstrnt == 2:
	# If irstrnt= 2 then we just want to calculate esp's of q0's (and p0's)
	for k in range(iuniq):
		qcal[k] = q0[k]
	qwt = 0.0

	if ipermdip > 0:
		for k in range(iuniq_p):
			pcal[k] = p0[k]
		pwt = 0.0
else:
	# Do the least-squares fitting
	print("Performing least-squares fitting ...")
	charge_opt()

print("Outputting statistics and molecular moments ...")
###### Calculate dipole moments ######
calc_dip()

###### Center & reorient molecule/dipole and calculate quadrupole moments ######
if ireornt > 0:
	reornt()
	quad_mol_com = calc_quad(co, dipindperm_com)
else:
	np.zeros((6,nmol))
	quad_mol = calc_quad(crd, dipindperm)

###### Write charges & permanent dipoles ######
wrt_qout()

###### Calculate and print sum-of-squares, rmse, and rrmse ######
wrt_out()

output_file.close()

print("-----------------------------------------------------------------")
print("To cite PyRESP use:")
print()
print("  Shiji Zhao, Haixin Wei, Piotr Cieplak, Yong Duan, and Ray Luo, \n"
	"  \"PyRESP: A Program for Electrostatic Parameterizations of Additive \n"
	"  and Induced Dipole Polarizable Force Fields\". J. Chem. Theory \n"
	"  Comput. 2022, 18, 6, 3654-3670.")
#-----------------------------------------------------------------------
# The end of the main program
#-----------------------------------------------------------------------
