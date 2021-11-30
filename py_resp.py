#!/usr/bin/env python3

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
iuniq,iuniq_p,nesp,nqpl,ihfree,irstrnt = 0,0,0,0,1,1

###### runlab ######
title = None

###### espcom ######
apot,awt,bpot,bwt = [None]*4
ssvpot,rmse,rrmse,tot_nesp,max_nesp = 0.0,0.0,0.0,0,0

###### calcul ###### 
qcal,a,b,qwtval,pwtval,iqpcntr = [None]*6

###### polarization ######
pol_dict,ipol,ipermdip,igdm = [None]*4
AinvQ_list,AinvP_L2G_list,L2G_list,mat_out = [None]*4

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
	parser.add_argument("-i", "--input", required=True, help="type: input, required; description: input options")
	parser.add_argument("-e", "--espot", required=True, help="type: input, required; description: input of MEP and coordinates")
	parser.add_argument("-q", "--qin", help="type: input, optional; description: replacement charges")
	parser.add_argument("-o", "--output", required=True, help="type: output, always produced; description: output of results")
	parser.add_argument("-t", "--qout", required=True, help="type: output, always produced; description: output of current charges")
	parser.add_argument("-s", "--esout", help="type: output, optional; description: generated MEP values for new charges")
	parser.add_argument("-ip", "--polariz", help="type: input, optional; description: atomic polarizabilities")
	
	args = parser.parse_args()
	Input,espot,qin,polariz = args.input,args.espot,args.qin,args.polariz
	output,qout,esout = args.output,args.qout,args.esout

	return Input,output,qin,polariz,qout,espot,esout

def read_in():
	global ioutopt,nmol,iqopt,irstrnt,ihfree,qwt,pwt
	global ipol,ipermdip,igdm,exc12,exc13,virtual

	# start of molecule input
	output_file.write('\n -----------------------------------------------')
	output_file.write('\n      Py_RESP Alpha Version  ')
	output_file.write('\n -----------------------------------------------')
	output_file.write('\n '+input_file.readline())
	output_file.write(' -----------------------------------------------\n\n')

	# read in charge, number of charge centers, and control parameters
	if 'cntrl' in f90nml.read(Input):
		cntrl_nml = f90nml.read(Input)['cntrl']
		ioutopt = cntrl_nml['ioutopt'] if 'ioutopt' in cntrl_nml else 0
		nmol = cntrl_nml['nmol'] if 'nmol' in cntrl_nml else 1
		iqopt = cntrl_nml['iqopt'] if 'iqopt' in cntrl_nml else 1
		irstrnt = cntrl_nml['irstrnt'] if 'irstrnt' in cntrl_nml else 1
		ihfree = cntrl_nml['ihfree'] if 'ihree' in cntrl_nml else 1
		qwt = cntrl_nml['qwt'] if 'qwt' in cntrl_nml else 0.0005
		pwt = cntrl_nml['pwt'] if 'pwt' in cntrl_nml else 0.0005
		ipol = cntrl_nml['ipol'] if 'ipol' in cntrl_nml else 0
		ipermdip = cntrl_nml['ipermdip'] if 'ipermdip' in cntrl_nml else 0
		igdm = cntrl_nml['igdm'] if 'igdm' in cntrl_nml else 1
		exc12 = cntrl_nml['exc12'] if 'exc12' in cntrl_nml else 1
		exc13 = cntrl_nml['exc13'] if 'exc13' in cntrl_nml else 1
		virtual = cntrl_nml['virtual'] if 'virtual' in cntrl_nml else 0
	else:
		output_file.write('Sorry, you must use namelist input\n')
		sys.exit()

	output_file.write('\n nmol        = %d   iqopt       = %d'%(nmol, iqopt))
	output_file.write('\n ihfree      = %d   irstrnt     = %d'%(ihfree, irstrnt))
	output_file.write('\n ioutopt     = %d   qwt         = %.8f'%(ioutopt, qwt))
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
		output_file.write('Error: permanent dipole can be enabled only if ipol > 0\n')
		sys.exit()

def sing_mol():
	global nlgrng,iuniq,title,wtmol,ibeg,iend,izan,ivary
	if ipermdip > 0:
		global iuniq_p,ibeg_p,iend_p,izan_p,ivary_p

	output_file.write("\n\n single-molecule run")

	title = []
	wtmol = np.ndarray((1))

	# read in fitting weight for q0 and esp point weighting
	wtmol[0] = float(input_file.readline())
	output_file.write("\n\n Reading input for molecule 1 weight:%10.3f \n"%wtmol[0])
	title.append(input_file.readline())
	output_file.write(title[0])

	# read in charge, number of charge centers (and number of permanent dipoles)
	line = input_file.readline().split()
	ich, iuniq = int(line[0]), int(line[1])
	output_file.write('\n Total charge (ich):%3d'%ich)
	output_file.write('\n Number of centers:%3d'%iuniq)
	if ipermdip > 0:
		iuniq_p = int(line[2])
		output_file.write('\n Number of permanent dipoles:%3d'%iuniq_p)

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

	print("ibeg",ibeg)
	print("ibeg_p",ibeg_p)
	print("iend",iend)
	print("iend_p",iend_p)

	# read in atomic number izan[i] and ivary[i]
	for i in range(iuniq):
		line = input_file.readline().split()
		izan[i], ivary[i] = int(line[0]), int(line[1])
		output_file.write("\n%5d%5d%5d"%(i+1, izan[i], ivary[i]))
		if ipermdip > 0:
			dip_num.append(len(line)-2)
			for j in range(len(line)-2):
				izan_p[p_cnt] = int(line[0])
				ivary_p[p_cnt] = int(line[j+2])
				output_file.write("%5d"%(ivary_p[p_cnt]))
				p_cnt += 1

	# read in lagrange constraints (including total and intra-molecular 
	# charge constraints)
	lagrange(ich, 0)

	input_summary(dip_num)

def mult_mol():
	######################################################################
	# this function reads in multiple molecule input. In function readin
 	# it has already read the control variable input for the run, namely:
 	# ich, iuniq, inopt, iqopt, ihfree, irstrnt, nlgrng, nmol
 	# where nmol > 1 caused this function to be called.
 	# The input form for the other molecules is that their entire control
 	# decks are appended to the initial 2 control lines just read in,
 	# each control deck separated by a blank line, and then comes the
 	# multiple-molecule specific input, which is
 	#  - equivalencing of centers between molecules in the series
 	#  - the lagrange constraints to be applied between molecules
 	#
 	# the control characters read in the individual job decks are ignored
 	# except for ich, and icntrs. The lagrange (charge) constraints
 	# contained in the individual-molecule inputs ARE included.
 	#
 	######################################################################
	global nlgrng,iuniq,title,wtmol,ibeg,iend,izan,ivary
	if ipermdip > 0:
		global iuniq_p,ibeg_p,iend_p,izan_p,ivary_p

	output_file.write("\n\n multiple-molecule run of %d molecules"%nmol)

	title = []
	wtmol = np.ndarray((nmol))
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
		# read in fitting weight for q0 and esp point weighting
		wtmol[imol] = float(input_file.readline())
		output_file.write("\n\n Reading input for molecule %d weight:%10.3f\n"%(imol+1,wtmol[imol]))
		title.append(input_file.readline())
		output_file.write(title[imol])

		# read in charge, number of charge centers (and number of permanent dipoles)
		line = input_file.readline().split()
		ich, icntrs = int(line[0]), int(line[1])
		output_file.write('\n Total charge (ich):%3d'%ich)
		output_file.write('\n Number of centers:%3d'%icntrs)
		if ipermdip > 0:
			icntrs_p = int(line[2])
			output_file.write('\n Number of permanent dipoles:%3d'%icntrs_p)

		# now some book-keeping: iuniq is the global variable for the total
		# number of centers over all molecules. The first center of this
		# mol therefore starts in iuniq and goes to iuniq+icntrs-1.
		#
		# Same for permanent dipoles.
		ibeg[imol] = iuniq
		iend[imol] = iuniq+icntrs-1
		if ipermdip > 0:
			ibeg_p[imol] = iuniq_p
			iend_p[imol] = iuniq_p+icntrs_p-1

		# trap for having too many centers
		if iend[imol]+1 > maxq:
			output_file.write('\n ERROR: more than %5d centers'%maxq)
			sys.exit()

		# Read in atomic number izan[i] and ivary[i]
		# Since ivary[i] is supposed to correspond to a center-number in the
		# same molecule, this has to be adjusted to ivary[i]+ibeg[imol]
		# convert angstroms to bohrs if necessary
		izan.resize((iuniq+icntrs), refcheck=False)
		ivary.resize((iuniq+icntrs), refcheck=False)
		iuniq += icntrs
		if ipermdip > 0:
			izan_p.resize((iuniq_p+icntrs_p), refcheck=False)
			ivary_p.resize((iuniq_p+icntrs_p), refcheck=False)
			iuniq_p += icntrs_p

		for i in range(ibeg[imol], iend[imol]+1):
			line = input_file.readline().split()
			izan[i], ivary[i] = int(line[0]), int(line[1])
			output_file.write("\n%5d%5d%5d"%(i+1, izan[i], ivary[i]))
			if ivary[i] > 0:
				ivary[i] += ibeg[imol]

			if ipermdip > 0:
				dip_num.append(len(line)-2)
				for j in range(len(line)-2):
					izan_p[p_cnt] = int(line[0])
					ivary_p[p_cnt] = int(line[j+2])
					output_file.write("%5d"%(ivary_p[p_cnt]))
					if ivary_p[p_cnt] > 0:
						ivary_p[p_cnt] += ibeg_p[imol]
					p_cnt += 1		

		# now read in the lagrange constraints for this molecule (including total  
		# and intra-molecular charge constraints)
		lagrange(ich, imol)

	# end of molecule input, now do other preparation stuff

	# read past a blank line after the final molecule job deck and then
	# read in inter-molecule lagrange constraints.
	# The "-99" for the total charge tells lagrange to drop the total charge
	# constraint.
	lagrange(-99)

	# inter-molecule equivalencing
	print("before equiv, ivary",ivary)
	print("before equiv, ivary_p",ivary_p)
	mol_equiv()
	print("after equiv, ivary",ivary)
	print("after equiv, ivary_p",ivary_p)

	input_summary(dip_num)

def lagrange(ncharge, imol=-1):
	###################################################
	# read in and assign lagrange constraint pointers #
	# called from "readin" and "mult_mol"             #
	###################################################
	global nlgrng

	line = input_file.readline().split()
	# If line is not empty, implement intra- or inter- molecular charge constraint(s)
	if line:
		output_file.write("\n ------------------------------------------------------------")
		output_file.write("\n reading intra- or inter- molecular charge constraint(s) info")
		output_file.write("\n ------------------------------------------------------------")
		while line:
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

	# as long as ncharge is not -99, implement the "total charge" constraint
	if ncharge > -99:
		nlgrng += 1
		grpchg.resize((nlgrng), refcheck=False)
		lgrcnt.resize((nlgrng,maxq), refcheck=False)

		grpchg[nlgrng-1] = float(ncharge)
		for j in range(ibeg[imol], iend[imol]+1):
			lgrcnt[nlgrng-1][j] = 1

def mol_equiv():
	# This function carries out the inter-molecule equivalencing by
	#
	# First : read the cards saying how many centers will be read in in the
	#         next card. a zero means we have finished input
	#
	# Second: read the first-occurrence-in-each-molecule of the centers to
	#        be equivalenced.
	#
	#        The specifcations MUST be in ascending order.
	#
	#        The expanding of the centers within each
	#        molecule is based on the ivary values for the individual mol.
	#
	#        if ivary for mol 2+ is zero it is replaced with the atom number
	#        of mol 1.
	output_file.write("\n ------------------------------------")
	output_file.write("\n reading equivalence info for charges")
	output_file.write("\n ------------------------------------")

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
			ivary[iatm[i]-1] = iatm[0]

		line = input_file.readline().split()

	for i in range(iuniq):
		if ivary[i] > 0:
			vary = ivary[ivary[i]-1]
			if vary > 0:
				ivary[i] = vary

	if ipermdip > 0:
		output_file.write("\n ----------------------------------------------")
		output_file.write("\n reading equivalence info for permanent dipoles")
		output_file.write("\n ----------------------------------------------")

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
				ivary_p[idip[i]-1] = idip[0]
	
			line = input_file.readline().split()

		for i in range(iuniq_p):
			if ivary_p[i] > 0:
				vary = ivary_p[ivary_p[i]-1]
				if vary > 0:
					ivary_p[i] = vary

def input_summary(dip_num):
	output_file.write("\n\n  --------------------")
	if ipermdip > 0:
		output_file.write("\n   Atom/Dipole Ivary")
		pcnt = 0
	else:
		output_file.write("\n     Atom   Ivary")
	output_file.write("\n  --------------------")
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
	
	output_file.write("\n\n Total number of atoms      =%5d"%iuniq)
	if ipermdip > 0:
		output_file.write("\n Total number of permanent dipoles      =%5d"%iuniq_p)
	output_file.write("\n Weight factor on initial charge restraints=%10.6f\n"%qwt)
	if ipermdip > 0:
		output_file.write("\n Weight factor on initial dipole restraints=%10.6f\n"%pwt)
	output_file.write("\n\n There are%3d charge constraints"%nlgrng)

def read_pol_dict():
	global pol_dict

	output_file.write("\n\n\n Read polarizability and radii information from %s"%polariz)
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
	print("pol_dict",pol_dict)

def neighbors(imol):
	natm = iend[imol]-ibeg[imol]+1

	n12, dict12, ibond, nbonds = build12(imol, natm)
	n13, dict13 = build13(ibond, nbonds, natm)
	#build14()
	return n12, dict12, n13, dict13

def build12(imol, natm):
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
				elif izi == 1 or izj == 1:
					bonded = dist2 < cutoff2[6]
				if bonded:
					ibond[0][nbonds] = i-ibeg[imol]
					ibond[1][nbonds] = j-ibeg[imol]
					nbonds += 1
	return nbonds, ibond

def build13(ibond, nbonds, natm):
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
	# read in the electrostatic potential points used in the fitting,
	# building up as we go the matrices for LU decomposition
	#
	# called from Main
	global apot, awt, bpot, bwt, crd, ssvpot, nesp, mat_out, tot_nesp, max_nesp
	if ipol > 0:
		global n12_list, dict12_list, n13_list, dict13_list, AinvQ_list
		if ipermdip > 0:
			global AinvP_L2G_list, L2G_list
	
	apot = np.zeros((iuniq+iuniq_p,iuniq+iuniq_p))
	awt = np.zeros((iuniq+iuniq_p,iuniq+iuniq_p))
	bpot = np.zeros((iuniq+iuniq_p))
	bwt = np.zeros((iuniq+iuniq_p))
	crd = np.zeros((3,iuniq))
	mat_out = np.ndarray((0, iuniq+iuniq_p))

	espot_file = open(espot, 'r')

	if ipol > 0:
		n12_list, dict12_list = [], []
		n13_list, dict13_list = [], []
		AinvQ_list = []
		if ipermdip > 0:
			AinvP_L2G_list, L2G_list = [], []

	ioff = 0     # local variable
	for imol in range(nmol):
		if ipol > 0:
			atype = []        # local variable
		line = espot_file.readline().split()
		natm, nesp = int(line[0]), int(line[1])
		tot_nesp += nesp
		max_nesp = max(max_nesp, nesp)
		output_file.write("\n\n\n Reading esp's for molecule %3d"%(imol+1))
		output_file.write("\n total number of atoms      = %5d"%natm)
		output_file.write("\n total number of esp points = %5d"%nesp)

		output_file.write("\n\n center       X               Y               Z")
		for i in range(natm):
			line = espot_file.readline().split()
			crd[0][ioff], crd[1][ioff], crd[2][ioff] = float(line[0]),float(line[1]),float(line[2])
			if ipol > 0:
				atype.append(line[3].lower())
			output_file.write("\n {:4d}{:16.7E}{:16.7E}{:16.7E}".format(i+1, crd[0][ioff], crd[1][ioff], crd[2][ioff]))
			ioff += 1

		wt = wtmol[imol]
		wt2 = wt*wt

		dismatq = np.zeros((natm))
		matq = np.zeros((natm))
		if ipol > 0:
			dismatdip = np.zeros((1, 3*natm))

			n12, dict12, n13, dict13 = neighbors(imol)
			n12_list.append(n12)
			dict12_list.append(dict12)
			n13_list.append(n13)
			dict13_list.append(dict13)

			print()
			print("n12_list, dict12_list")
			print(n12_list, dict12_list)
			print("n13_list, dict13_list")
			print(n13_list, dict13_list)

			if ipermdip > 0:
				npermdip = iend_p[imol]-ibeg_p[imol]+1
				matp = np.zeros((npermdip))

				AinvQ, AinvP, L2G = bld_AinvQP(imol, atype)
				AinvP_I = AinvP + np.identity(3*natm)
				print()
				print("AinvP_I")
				print(AinvP_I)

				AinvP_L2G_list.append(np.matmul(AinvP,L2G))
				L2G_list.append(L2G)
			else:
				AinvQ = bld_AinvQP(imol, atype)

			AinvQ_list.append(AinvQ)
			print()
			print("AinvQ",AinvQ)
			print()
			print("AinvP",AinvP)

		mat_out.resize((max_nesp, iuniq+iuniq_p), refcheck=False)

		xij = np.ndarray((3))
		for i in range(nesp):
			line = espot_file.readline().split()
			if not line:
				output_file.write("     premature end of potential file")
				sys.end()
			espi, xi, yi, zi = float(line[0]),float(line[1]),float(line[2]),float(line[3])
			ssvpot += wt2*espi*espi

			for j in range(ibeg[imol], iend[imol]+1):
				xij[0] = xi-crd[0][j]
				xij[1] = yi-crd[1][j]
				xij[2] = zi-crd[2][j]
				rij2 = xij[0]**2 + xij[1]**2 + xij[2]**2
				rij = math.sqrt(rij2)
				rij3 = rij*rij2

				j_idx = j-ibeg[imol]

				if ipol == 5 and igdm > 0:
					polj, radj = pol_dict[atype[j_idx]][0], pol_dict[atype[j_idx]][1]
					fe, ft, f0 = damp_facts(rij, 1, polj, 0, radj)
					print("i, j, fe, ft, f0",i, j, fe, ft, f0)
					rij /= f0
					rij3 /= fe

				dismatq[j_idx] = 1/rij
				if ipol > 0:
					for k in range(3):
						dismatdip[0][3*j_idx+k] = xij[k]/rij3

			if ipol > 0:
				dismatdip_AinvQ = np.matmul(dismatdip, AinvQ)

			# Diagonal elements of apot (charge)
			for j in range(ibeg[imol], iend[imol]+1):
				j_idx = j-ibeg[imol]

				if ipol > 0:
					matq[j_idx] = dismatq[j_idx] + dismatdip_AinvQ[0][j_idx]
				else:
					matq[j_idx] = dismatq[j_idx]

				mat_out[i][j] = matq[j_idx]
				bpot[j] += espi*matq[j_idx]
				apot[j][j] += matq[j_idx]**2

			# Off-diagonal elements of apot (charge)
			for j in range(ibeg[imol], iend[imol]):
				j_idx = j-ibeg[imol]
				for k in range(j+1, iend[imol]+1):
					k_idx = k-ibeg[imol]

					apot[j][k] += matq[j_idx]*matq[k_idx]

			if ipermdip > 0:
				dismatdip_AinvP_L2G = np.matmul(np.matmul(dismatdip, AinvP_I),L2G)

				# Diagonal elements of apot (dipole)
				for j in range(ibeg_p[imol], iend_p[imol]+1):
					j_idx1 = j-ibeg_p[imol]
					j_idx2 = j+iuniq
					matp[j_idx1] = dismatdip_AinvP_L2G[0][j_idx1]

					mat_out[i][j_idx2] = matp[j_idx1]
					bpot[j_idx2] += espi*matp[j_idx1]
					apot[j_idx2][j_idx2] += matp[j_idx1]**2

				# Off-diagonal elements of apot (dipole)
				for j in range(ibeg_p[imol], iend_p[imol]):
					j_idx1 = j-ibeg_p[imol]
					j_idx2 = j+iuniq
					for k in range(j+1, iend_p[imol]+1):
						k_idx1 = k-ibeg_p[imol]
						k_idx2 = k+iuniq
	
						apot[j_idx2][k_idx2] += matp[j_idx1]*matp[k_idx1]

				# Off-diagonal elements of apot (charge * dipole)
				for j in range(ibeg[imol], iend[imol]+1):
					j_idx = j-ibeg[imol]
					for k in range(ibeg_p[imol], iend_p[imol]+1):
						k_idx1 = k-ibeg_p[imol]
						k_idx2 = k+iuniq

						apot[j][k_idx2] += matq[j_idx]*matp[k_idx1]

		for j in range(ibeg[imol], iend[imol]+1):
			bwt[j] = wt2*bpot[j]
			awt[j][j] = wt2*apot[j][j]
			for k in range(j+1, iend[imol]+1):
				awt[j][k] = wt2*apot[j][k]

		if ipermdip > 0:
			for j in range(ibeg_p[imol], iend_p[imol]+1):
				j_idx = j+iuniq
				bwt[j_idx] = wt2*bpot[j_idx]
				awt[j_idx][j_idx] = wt2*apot[j_idx][j_idx]
				for k in range(j+1, iend_p[imol]+1):
					k_idx = k+iuniq
					awt[j_idx][k_idx] = wt2*apot[j_idx][k_idx]

			for j in range(ibeg[imol], iend[imol]+1):
				for k in range(ibeg_p[imol], iend_p[imol]+1):
					k_idx = k+iuniq
					awt[j][k_idx] = wt2*apot[j][k_idx]

	# symmetrize the potenitial and weighted potential matrices
	for j in range(iuniq+iuniq_p-1):
		for k in range(j+1, iuniq+iuniq_p):
			awt[k][j] = awt[j][k]
			apot[k][j] = apot[j][k]
	print()
	print("mat_out")
	print(mat_out)
	print()
	print("awt")
	print(awt)
	print("bwt")
	print(bwt)

	espot_file.close()

def bld_AinvQP(imol, atype):
	##################################################################################
	# A: Applequist matrix (Relay matrix)
	# Ainv: Inverse of Applequist matrix. This is the molecular polarizability matrix.
	# QFld: Charge field matrix
	# PFld: Dipole field matrix
	# AinvQ: Ainv * QFld
	# AinvP: Ainv * PFld
	# L2G: Matrix that converts permanent dipoles from local frames to global frames.
	##################################################################################

	# Initialize A, QFld (and PFld)
	natm = iend[imol]-ibeg[imol]+1
	A = np.zeros((3*natm, 3*natm))
	QFld = np.zeros((3*natm, natm))
	xij = np.ndarray((3)) # local variable, stores the distance elements between atoms i and j

	# Fill the diagonal matrices of A
	for i in range(natm):
		for j in range(3):
			idx = 3*i+j
			poli = pol_dict[atype[i]][0]
			A[idx][idx] = 1/poli

	if ipermdip > 0:
		npermdip = iend_p[imol]-ibeg_p[imol]+1
		PFld = np.zeros((3*natm, 3*natm))
		L2G = np.zeros((3*natm, npermdip))

	for i in range(natm-1):
		i_idx = ibeg[imol]+i
		for j in range(i+1, natm):
			j_idx = ibeg[imol]+j

			for k in range(3):
				xij[k] = crd[k][i_idx]-crd[k][j_idx]   # !!! from j to i

			rij2 = xij[0]**2 + xij[1]**2 + xij[2]**2
			rij = math.sqrt(rij2)
			rij3 = rij*rij2
			rij5 = rij2*rij3

			poli, polj = pol_dict[atype[i]][0], pol_dict[atype[j]][0]
			radi, radj = pol_dict[atype[i]][1], pol_dict[atype[j]][1]
			fe, ft, f0 = damp_facts(rij, poli, polj, radi, radj)
			print("i, j, fe, ft",i, j, fe, ft)

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
					PFld[idx1][idx2] = element2
					PFld[idx2][idx1] = element2

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
						PFld[idx1][idx2] = element
						PFld[idx3][idx4] = element
						PFld[idx2][idx1] = element
						PFld[idx4][idx3] = element

	print("A")
	print(A)
	print()
	if ipermdip > 0:
		# Fill L2G
		p_cnt = 0
		for i in range(natm):
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

		print("PFld before")
		print(PFld)
		print()
		print("L2G")
		print(L2G)
	print()
	print("QFld before")
	print(QFld)

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
	print()	
	print("QFld after12",QFld)
	print("PFld after12",PFld)

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
	print()
	print("QFld after13",QFld)
	print("PFld after13",PFld)

	Ainv = np.linalg.inv(A)
	print()
	print("Ainv",Ainv)

	AinvQ = np.matmul(Ainv, QFld)
	if ipermdip > 0:
		AinvP = np.matmul(Ainv, PFld)
		return AinvQ, AinvP, L2G
	else:
		return AinvQ

def damp_facts(r, poli, polj, radi, radj):
	fe, ft, f0 = 1, 1, 1  # Applequist

	if ipol == 2:  # Tinker-exponential
		coef1 = math.sqrt(radi/poli)
		coef2 = math.sqrt(radj/polj)
		v = r**3 * coef1 * coef2
		exp_v = math.exp(-v)
		fe = 1 - exp_v
		ft = fe - v * exp_v
		#f0 = exp_v

	elif ipol == 3:  # Exponential-Thole
		coef1 = math.sqrt(radi) / poli**(1/6)
		coef2 = math.sqrt(radj) / polj**(1/6)
		v = r * coef1 * coef2
		exp_v = math.exp(-v)
		v2 = v**2
		v3 = v * v2
		fe = 1 - (v2/2 + v + 1) * exp_v
		ft = fe - (v3/6) * exp_v
		#f0 = exp_v

	elif ipol == 4:  # linear Thole
		coef1 = math.sqrt(radi) * poli**(1/6)
		coef2 = math.sqrt(radj) * polj**(1/6)
		v = r / (coef1 * coef2)
		if v < 1:
			v3 = v**3
			v4 = v3*v
			fe = 4*v3 - 3*v4
			ft = v4
		#	f0 = v
		#else:
		#	f0 = 0

	elif ipol == 5:  # pGM
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
	global q0
	
	q0 = np.ndarray((iuniq))
	if ipermdip > 0:
		global p0
		p0 = np.ndarray((iuniq_p))

	if iqopt > 1:
		# replace initial charges q0 from qin if iqopt>1
		qin_file = open(qin, 'r')
		if ipermdip > 0:
			output_file.write("\n since iqopt>1, %4d new q0 values and %4d "%(iuniq, iuniq_p))
			output_file.write("\n new p0 values will be read from file %s"%qin)
		else:
			output_file.write("\n since iqopt>1, %4d new q0 values will "%iuniq)
			output_file.write("\n be read from file %s"%qin)

		q0_idx = 0
		if ipermdip > 0:
			p0_idx = 0

		# now read in replacement charges
		for imol in range(nmol):
			natm = iend[imol]-ibeg[imol]+1

			for i in range(natm+10):
				line = qin_file.readline().split()

			while line:
				print(line)
				q0[q0_idx] = float(line[3])
				q0_idx += 1
				line = qin_file.readline().split()

			if ipermdip > 0:
				for i in range(3):
					line = qin_file.readline().split()

				while line:
					print(line)
					p0[p0_idx] = float(line[4])
					p0_idx += 1
					line = qin_file.readline().split()

				for i in range(2*natm+5):
					line = qin_file.readline()
			else:
				for i in range(natm+2):
					line = qin_file.readline()

		qin_file.close()
		print("q0_idx = ", q0_idx)
		print("p0_idx = ", p0_idx)
		print("q0: ",q0)
		print("p0: ",p0)
	else:
		# setting initial charges/dipoles to 0; done if iqopt=1
		q0.fill(0.0)
		if ipermdip > 0:
			output_file.write("\n\n iqopt=1, all q0 and p0 values will be set to 0")
			p0.fill(0.0)
		else:
			output_file.write("\n\n iqopt=1, all q0 values will be set to 0")

def data_prep():
	# setup pointers for groups of charges based on "ivary" info
	#
	# called from Main

	#************************************************************************
	# begin section: set lists for combined and frozen charges
	#
	# ivary[i] = 0, it is a new charge center to be fitted
	# ivary[i] =+n, it is a charge center to be fitted with center n
	#                    (center n must be a previous center entered with
	#                    ivary[n] = 0
	# ivary[i] =-n, it is a frozen charge center to be kept at q0[i]
	#
	# Similarly:
	#
	# ivary_p[i] = 0, it is a new permanent dipole center to be fitted
	# ivary_p[i] =+n, it is a permanent dipole center to be fitted with 
	#                    center n (center n must be a previous center 
	#                    entered with ivary_p[n] = 0
	# ivary_p[i] =-n, it is a frozen permanent dipole center to be kept
	#                    at p0[i]
	#************************************************************************
	global iqpcntr, nqpl

	nqp = 0
	iqpcntr = np.zeros((iuniq+iuniq_p+nlgrng), dtype=int)
	for i in range(iuniq):
		if ivary[i] == 0:
			#nat += 1
			nqp += 1
			iqpcntr[i] = nqp
		elif ivary[i] > 0:
			iqpcntr[i] = iqpcntr[ivary[i]-1]

			if iqpcntr[i] > nqp:
				output_file.write("\n data_prep: charge equivalence input is screwy")
				sys.end()
		else:
			iqpcntr[i] = -1

	if nqp == 0:
		output_file.write("\n Warning: ALL charges are frozen!!!")
	else:
		output_file.write("\n\n\n Number of unique UNfrozen charge centers=%5d"%nqp)
		nq = nqp

	if ipermdip > 0:
		for i in range(iuniq_p):
			if ivary_p[i] == 0:
				nqp += 1
				iqpcntr[iuniq+i] = nqp
			elif ivary_p[i] > 0:
				iqpcntr[iuniq+i] = iqpcntr[iuniq+ivary_p[i]-1]

				if iqpcntr[iuniq+i] > nqp:
					output_file.write("\n data_prep: permanent dipole equivalence input is screwy")
					sys.end()
			else:
				iqpcntr[iuniq+i] = -1

		if (nqp-nq) == 0:
			output_file.write("\n Warning: ALL permanent dipoles are frozen!!!")
		else:
			output_file.write("\n Number of unique UNfrozen permanent dipoles=%5d"%(nqp-nq))

	# finish off list with Lagrange constraints
	for i in range(nlgrng):
		iqpcntr[iuniq+iuniq_p+i] = nqp+i+1

	# set nqpl to the total number of row elements (charges/permanent dipoles
	# to be independantly fit + constraints) in fitting matrix
	nqpl = nqp + nlgrng
	print("iqpcntr",iqpcntr)
	print("nqp",nqp)
	print("nlgrng",nlgrng)

	# done adding Lagrange constraints to elements list

	# read in charges must now be averaged
	# a posteriori averaging of replacement charges according
	# to current ivary charge-combining pointers
	#if iqopt == 3:
	#	for i in range(iuniq-1):
	#		qcntrs = q0[i]
	#		tmpctr = 1.0
	#		for j in range(i+1, iuniq):
	#			if ivary[j] == i:
	#				qcntrs += q0[j]
	#				tmpctr += 1.0
#
	#		if tmpctr > 0.99:
	#			qcntrs /= tmpctr
	#			q0[i] = qcntrs
	#			for j in range(i+1,iuniq):
	#				if ivary[j] == i:
	#					q0[j] = qcntrs

def charge_opt():
	# driver for the charge determinization/optimizaton
	#
	# called from Main
	global irstrnt, awork, bwork

	qold = np.zeros((iuniq))   # local variable

	irsave = 0 # local variable
	nitern = 0 # local variable

	# qtol & maxit are criteria for convergence & maximum iterations for
	# the non-linear optimizations.
	qtol = 0.000001
	maxit = 24

	# only on first pass through this function (indicated by nitern= 0),
	# if irstrnt > 0, transfer irstrnt to irsave and reset irstrnt to 0,
	# in order to get an initial guess using a harmonic constraint.  This is
	# done so restraint function rstran() will use a harmonic restraint.
	if irstrnt > 0:
		irsave = irstrnt
		irstrnt = 0
		output_file.write("\n\n Non-linear optimization requested.")

	# now go do a "harmonic restraint" run, restraint= qwt(qcal[i]-q0[i])**2
	# -- loop to convergence

	while nitern < maxit:
		matbld()
	
		# solve (Ax = b) where A and b are input, x is solution
		#              awork x = bwork
		# 
		# the solution "x" is returned in "b" (bwork)
		# 
		# -- condition the matrix diagonal to avoid DGETRF() detectingNN (wrapped in lu_factor)
		#    singularity
		for jn in range(nqpl):
			if abs(awork[jn][jn]) < 1.0E-10:
				awork[jn][jn] = 1.0E-10

		print()
		print("awork before", nitern)
		print(awork)
		print("bwork before", nitern)
		print(bwork)
	
		awork, piv = lu_factor(awork, overwrite_a=True)
		bwork = lu_solve((awork, piv), bwork, overwrite_b=True)

		print()
		print("awork after", nitern)
		print(awork)
		print("bwork after", nitern)
		print(bwork)
		print("irstrnt after",irstrnt)
	
		# -- copy solution vector "bwork" to 'calculated charges' vector
		#    qcal and pcal
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
	
		# -- a quick check from rstrn: if irstrnt is now negative,
		#    there are no restraints because no qwtval(i) > 0.1e-10,
		#    so reset irsave= 0
		if irstrnt < 0:
			irsave = 0
			output_file.write("\n\n WARNING: Restraints were requested, but the restraint weights were all zero\n\n")
	
		# -- we're finished if it's only a "harmonic restraint" run,
		#    but if it's a non-linear optimization (irsave>0)...
		#    we've only just begun (i.e. we have our initial guess)
		if irsave <= 0:
			return
		else:
			# -- it's a non-linear optimization: reset irstrnt (to now
			#    calculate the proper non-linear restraint derivatives
			#    in routine rstran)
			irstrnt = irsave
	
		# -- begin iterative optimization loop with comparison of
		#    old & new charges; calculate the convergence and replace
		#    the old charges with the new
		qchnge = 0.0
		for i in range(iuniq):
			qdiff = qcal[i] - qold[i]
			qchnge += qdiff*qdiff
			qold[i] = qcal[i]
	
		qchnge = math.sqrt(qchnge)/iuniq
		output_file.write("\n qchnge ={:20.10E}".format(qchnge))
	
		# -- if this is less than qtol then we're done
		if qchnge < qtol and nitern > 1:
			output_file.write("\n\n Convergence in%5d iterations\n\n"%nitern)
			return

		# loop again
		nitern += 1

	output_file.write("\n after %5d iterations, no convergence!"%maxit)

def matbld():
	# called from "chgopt"
	#
	# build up matrices for LU decomposition:
	#
	#   stage 1: copy weighted matrices awt and bwt to work arrays awork and bwork
	#            (which are destroyed in the LU decomp & back subst)
	#
	#   stage 2: if charge restraints are to be included,
	#            then modify awork and bwork appropriately
	global a, b, awork, bwork

	a = np.zeros((iuniq+iuniq_p+nlgrng,iuniq+iuniq_p+nlgrng))
	b = np.zeros((iuniq+iuniq_p+nlgrng))

	for k in range(iuniq+iuniq_p):
		b[k] = bwt[k]
		for j in range(iuniq+iuniq_p):
			a[j][k] = awt[j][k]

	# fill in the final columns & rows of A with the Lagrange
	# constraints which keep the charge on groups of atoms to a
	# constant
	#
	# note index counters!
	for i in range(nlgrng):
		b[iuniq+iuniq_p+i] = grpchg[i]
		for j in range(iuniq+iuniq_p+nlgrng):
			a[iuniq+iuniq_p+i][j] = lgrcnt[i][j]
			a[j][iuniq+iuniq_p+i] = lgrcnt[i][j]

	print("a")
	print(a)
	print()
	print("b")
	print(b)

	# add restraint to initial charge q0[i]:
	rstran()

	# build awork and bwork based on "combined and frozen centers" info:
	#
	# 1) frozen centers do not appear in the matrix of fitted charges
	# 2) combined centers appear as one single charge center for fitting
	#
	# first, since we accumulate values, zero out awork & bwork up to nqpl
	# (the independant + contraint number):
	awork = np.zeros((nqpl,nqpl))
	bwork = np.zeros((nqpl))

	# loop over all centers, building awork & bwork from A and
	# B based on iqpcntr: for each center, iqpcntr[i] dictates which of
	# the fitted charges it is and therefore where it goes in the matrices.
	# If iqpcntr[j] < 1, this center is a frozen charge and it is skipped as
	# far as forming a row in awork, and its esp contribution is subtracted
	# from bwork to take care of it's awork jth column-element for each i.
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
	# routine to assign the retraint weights
	# to the diagonal of A and to B
	#
	# called from "matbld"

	#----------------------------------------------------------------------
	# two kinds of restraint are available:
	#
	# a) a harmonic restraint to the initial charge.  Fine as long as there
	#  aren't any large charges that SHOULD be large... these really feel a
	#  strong force if they are restrained to a low value.
	#
	# b) a hyperbolic restraint to a charge of 0.  This gets asymptotic at
	#  "large" values, so "large" charges aren't pulled down any stronger
	#  than some (reasonable) limiting force.  This is a non-linear
	#  weighting function, so the fit procedure is iterative.
	#
	# other options for restraints to initial charge q0[i]:
	# if requested, restrain the charges by modifying the sum-of-squares
	# cost function derivative.  The scheme for doing this is as follows:
	#
	# if control variable ihfree > 0, let hydrogen charges float free
	#                                   (i.e. reset their qwtval to 0.0).
	#
	#-----------------------------------------------------------------------
	global irstrnt, qwtval
	qwtval = np.ndarray((iuniq))
	qwtval.fill(qwt)
	print("qwtval before ",qwtval)

	if ipermdip > 0:
		global pwtval
		pwtval = np.ndarray((iuniq_p))
		pwtval.fill(pwt)
		print("pwtval before ",pwtval)

	for i in range(iuniq):
		if ihfree > 0 and izan[i] == 1:
			qwtval[i] = 0.0

		if irstrnt == 0:
			a[i][i] += qwtval[i]

			# q0 has the initial and/or frozen charge
			b[i] += qwtval[i]*q0[i]
		elif irstrnt > 0 and qwtval[i] > 0.1E-10:
			# use analytic gradient of the hyperbola

			# qcal has the current (calculated) charge
			qwtval[i] = qwt/math.sqrt(qcal[i]*qcal[i] + 0.01)
			a[i][i] += qwtval[i]

	for i in range(iuniq_p):
		i_idx = iuniq+i
		if ihfree > 0 and izan_p[i] == 1:
			pwtval[i] = 0.0

		if irstrnt == 0:
			a[i_idx][i_idx] += pwtval[i]

			# p0 has the initial and/or frozen dipole
			b[i_idx] += pwtval[i]*p0[i]
		elif irstrnt > 0 and pwtval[i] > 0.1E-10:
			# use analytic gradient of the hyperbola

			# pcal has the current (calculated) dipol
			pwtval[i] = pwt/math.sqrt(pcal[i]*pcal[i] + 0.01)
			a[i_idx][i_idx] += pwtval[i]

	print("qwtval after ",qwtval)
	print("pwtval after ",pwtval)

	# if all qwtval[i] and pwtval[i] are 0.0, no restraints so reset irstrnt= -1
	for i in range(iuniq):
		if qwtval[i] > 0.1E-10:
			return
	for i in range(iuniq_p):
		if pwtval[i] > 0.1E-10:
			return

	irstrnt = -1

def calc_dip():
	# function to calculate the dipole moments
	global dipol_mol, dipmom_mol, dipind, dipperm, dipindperm

	dipol_mol = np.zeros((3,nmol))
	dipmom_mol = np.zeros((nmol))
	dipind = np.zeros((3,iuniq))
	dipperm = np.zeros((3,iuniq))
	dipindperm = np.zeros((3,iuniq))

	# calculate the induced dipole and permanent dipole on each atom (global frame)
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

	# calculate molecular dipole moment and quadrapole moment
	for imol in range(nmol):
		for i in range(ibeg[imol], iend[imol]+1):
			for j in range(3):
				dipindperm[j][i] = dipind[j][i] + dipperm[j][i]
				dipol_mol[j][imol] += qcal[i]*crd[j][i] + dipindperm[j][i]

		dipmom_mol[imol] = math.sqrt(dipol_mol[0][imol]**2 + dipol_mol[1][imol]**2 + dipol_mol[2][imol]**2)

	# convert dipoles from a.u. to debyes
	dipol_mol *= au2D
	dipmom_mol *= au2D

def reornt():
	##############################################
	# translates molecule to center of mass and  #
	# reorients it along principal axes          #
	# in preparation for dipole and quadrupole   #
	# moment calculation.                        #
	##############################################

	# ---- ATOMIC WEIGHT ARRAY FOR CENTER OF MASS ----
	#
	#  20 elements were originally handled: H(1) - Ca(20)
	# 103 elements are now considered:      H(1) - Lr(103)
	#
	# Atomic weight from The Merck Index - Thirteeth edition
	# Merck & Co., INC., Whitehouse Station, NJ, 2001
	#
	# F.-Y. Dupradeau & P. Cieplak
	# http://q4md-forcefieldtools.org/
	global co, cmas_mol, dipol_mol_com, dipindperm_com

	cmas_mol = np.ndarray((3,nmol))
	dipol_mol_com = np.zeros((3,nmol))
	co = np.ndarray((3,iuniq))
	dipindperm_com = np.zeros((3,iuniq))

	wt = [1.0079,4.0026,
	6.9410,9.0122,10.8110,12.0107,
	14.0067,15.9994,18.9984,20.1797,
	22.9898,24.3050,26.9815,28.0855,
	30.9738,32.0650,35.4530,39.9480,
	39.0983,40.0780,44.9559,47.8670,
	50.9415,51.9961,54.9380,55.8450,
	58.9332,58.6934,63.5460,65.3900,
	69.7230,72.6400,74.9216,78.9600,
	79.9040,83.8000,85.4678,87.6200,
	88.9058,91.2240,92.9064,95.9400,
	97.9072,101.0700,102.9055,106.4200,
	107.8682,112.4110,114.8180,118.7100,
	121.7600,127.6000,126.9045,131.2930,
	132.9054,137.3270,138.9055,140.1160,
	140.9076,144.2400,144.9127,150.3600,
	151.9640,157.2500,158.9253,162.5000,
	164.9303,167.2590,168.9342,173.0400,
	174.9670,178.4900,180.9479,183.8400,
	186.2070,190.2300,192.2170,195.0780,
	196.9665,200.5900,204.3833,207.2000,
	208.9804,208.9824,209.9871,222.0176,
	223.0197,226.0254,227.0277,232.0381,
	231.0359,238.0289,237.0482,244.0642,
	243.0614,247.0704,247.0703,251.0796,
	252.0830,257.0951,258.0984,259.1010,
	262.1097]

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
	
	# convert center of mass from bohr to angstroms
	cmas_mol *= au2A

def cmass(wt, imol):
	# called from reornt()

	# THIS FUNCTION CALCULATES THE CENTER OF MASS OF THE MOLECULE.
	sumx, sumy, sumz, Sum = 0.0, 0.0, 0.0, 0.0
	for i in range(ibeg[imol], iend[imol]+1):
		idx = izan[i]-1
		sumx += crd[0][i] * wt[idx]
		sumy += crd[1][i] * wt[idx]
		sumz += crd[2][i] * wt[idx]
		Sum += wt[idx]
	return sumx/Sum, sumy/Sum, sumz/Sum

def cmove(xc, yc, zc, imol):
	# called from reornt()

	# THIS FUNCTION MOVES THE ORIGIN TO THE CENTER OF MASS.
	for i in range(ibeg[imol], iend[imol]+1):
		co[0][i] = crd[0][i] - xc
		co[1][i] = crd[1][i] - yc
		co[2][i] = crd[2][i] - zc

def momin(wt, imol):
	# called from reornt()

	# THIS FUNCTION CALCULATES THE MOMENTS OF INERTIA TENSOR AND
	# THE PRINCIPAL AXES OF ROTATION OF THE MOLECULE.  IT THEN
	# REORIENTS THE MOLECULE ALONG THE PRINCIPAL AXES OF
	# ROTATION (WITH THE ORIGIN AT THE CENTER OF MASS) IN
	# PREPARATION FOR CALCULATION OF QUADRAPOLE MOMENT
	# COMPONENTS.
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
		syz += wt[idx]*co[1][i]*co[2][i]  # need confirm
	ain[0][0] = sxx
	ain[1][1] = syy
	ain[2][2] = szz
	ain[0][1] = -sxy
	ain[0][2] = -sxz
	ain[1][2] = -syz
	ain[1][0] = ain[0][1]
	ain[2][0] = ain[0][2]
	ain[2][1] = ain[1][2]
	print("ain before diagm",ain)

	#  ----- CALCULATE PRINCIPAL AXES OF INERTIA -----
	eigvals, s = np.linalg.eig(ain)
	print("eigvals after diagm",eigvals)
	print("s before sort",s)

	for i in range(2):
		for j in range(i+1,3):
			if eigvals[i] > eigvals[j]:
				eigvals[i], eigvals[j] = eigvals[j], eigvals[i]
				for k in range(3):
					s[k][i], s[k][j] = s[k][j], s[k][i]

	print("eigvals after sort",eigvals)
	print("s after sort",s)
	return s

def mominrot(s, vecin, vecout, imol):
	# called from reornt()

	# THIS FUNCTION ROTATE MOLECULES/DIPOLES ALONG THE PRINCIPLE AXES OF MOMENT OF INERTIA.
	if len(vecin[0]) == nmol:
		xs = vecin[0][imol]*s[0][0] + vecin[1][imol]*s[1][0] + vecin[2][imol]*s[2][0]
		ys = vecin[0][imol]*s[0][1] + vecin[1][imol]*s[1][1] + vecin[2][imol]*s[2][1]
		zs = vecin[0][imol]*s[0][2] + vecin[1][imol]*s[1][2] + vecin[2][imol]*s[2][2]
		vecout[0][imol] = xs
		vecout[1][imol] = ys
		vecout[2][imol] = zs
	else:
		for i in range(ibeg[imol], iend[imol]+1):
			xs = vecin[0][i]*s[0][0] + vecin[1][i]*s[1][0] + vecin[2][i]*s[2][0]
			ys = vecin[0][i]*s[0][1] + vecin[1][i]*s[1][1] + vecin[2][i]*s[2][1]
			zs = vecin[0][i]*s[0][2] + vecin[1][i]*s[1][2] + vecin[2][i]*s[2][2]
			vecout[0][i] = xs
			vecout[1][i] = ys
			vecout[2][i] = zs

def calc_quad(quad, coord, dip):
	# function to calculate the molecular quadrupole moments
	for imol in range(nmol):
		for i in range(ibeg[imol], iend[imol]+1):

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

	# convert quadrupoles to debye*angstroms
	quad *= au2D*au2A

def wrt_qout():
	# called from Main

	qout_file = open(qout, 'w')
	qout_file.write("All values are reported in atomic units\n")

	ioff = 0
	for imol in range(nmol):
		qout_file.write("%FLAG TITLE\n")
		qout_file.write(" molecule %d: %s"%(imol+1, title[imol]))
		qout_file.write("\n%FLAG ATOM CRD: I4,3E16.7\n")
		qout_file.write(" atm.no      X               Y               Z\n")
		for i in range(ibeg[imol], iend[imol]+1):
			qout_file.write("{:4d}{:16.7E}{:16.7E}{:16.7E}\n".format(i+1, crd[0][i], crd[1][i], crd[2][i]))
		qout_file.write("\n%FLAG ATOM CHRG: 2(I4,X7),I4,X2,E16.7\n")
		qout_file.write(" atm.no   element.no   ivary      q(opt)\n")
		for i in range(ibeg[imol], iend[imol]+1):
			qout_file.write("{:4d}       {:4d}       {:4d}  {:16.7E}\n".format(i+1, izan[i], ivary[i], qcal[i]))

		if ipol > 0:
			if ipermdip > 0:
				qout_file.write("\n%FLAG PERM DIP LOCAL: 3(I4,X5),I4,X2,E16.7\n")
				qout_file.write(" dip.no   atm.no   ref.no   ivary      p(opt)\n")
				#for atm1 in sorted(dict12_list[imol]):
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

				qout_file.write("\n%FLAG PERM DIP GLOBAL: I4,3E16.7\n")
				qout_file.write(" atm.no      X               Y               Z\n")
				for i in range(ibeg[imol], iend[imol]+1):
					qout_file.write("{:4d}{:16.7E}{:16.7E}{:16.7E}\n".format(i+1, dipperm[0][i], dipperm[1][i], dipperm[2][i]))

			qout_file.write("\n%FLAG IND DIP GLOBAL: I4,3E16.7\n")
			qout_file.write(" atm.no      X               Y               Z\n")
			for i in range(ibeg[imol], iend[imol]+1):
				qout_file.write("{:4d}{:16.7E}{:16.7E}{:16.7E}\n".format(i+1, dipind[0][i], dipind[1][i], dipind[2][i]))
		qout_file.write("\n")

	qout_file.close()

def evlchi():
	# called from Main
	#
	# Evaluate chi-square for linear function espclci = sum_j(qj*termij),
	# where j = number of terms, qj is the coefficient to termij, and
	# chi-square is the merit function: chi-square = sum_i((espi-espclci)**2),
	# where i is the number of data points for which esp is known.

	#cross = 0.0   # local variable
	#ssyclc = 0.0  # local variable

	#for j in range(iuniq):
	#	cross += qcal[j]*bpot[j]
	#	for k in range(iuniq):
	#		ssyclc += qcal[j]*qcal[k]*apot[j][k]
	#if ipermdip > 0:
	#	for j in range(iuniq_p):
	#		cross += pcal[j]*bpot[iuniq+j]
	#		for k in range(iuniq_p):
	#			ssyclc += pcal[j]*pcal[k]*apot[iuniq+j][iuniq+k]

	#chipot = ssvpot - 2.0*cross + ssyclc
	chipot = 0.0

	# read in the electrostatic potential points used in the fitting,
	# calculate esp using existing charges, and write out both esp's & residual

	# open the file containing the qm esp points & read in the no. of points
	espot_file = open(espot, 'r')

	if ioutopt == 1:
		esout_file = open(esout, 'w')

	for imol in range(nmol):
		chipot_mol = 0.0
		line = espot_file.readline().split()
		natm, nesp = int(line[0]), int(line[1])
		if ioutopt == 1:
			esout_file.write("molecule %d: %s"%(imol+1, title[imol]))
			esout_file.write("nesp: %d   natm: %d\n"%(nesp, natm))
			esout_file.write("      X         Y         Z        esp_qm      esp_clc       diff\n")
	
		for i in range(natm):
			line = espot_file.readline()
	
		for i in range(nesp):
			line = espot_file.readline().split()
			if not line:
				output_file.write("\n unexpected eof in %s"%espot)
				sys.exit()
			espqmi, xi, yi, zi = float(line[0]),float(line[1]),float(line[2]),float(line[3])
			espclc = 0.0
	
			for j in range(ibeg[imol], iend[imol]+1):
				espclc += mat_out[i][j]*qcal[j]

			if ipermdip > 0:
				for j in range(ibeg_p[imol], iend_p[imol]+1):
					espclc += mat_out[i][iuniq+j]*pcal[j]
	
			vresid = espqmi - espclc
			chipot_mol += vresid**2

			if ioutopt == 1:
				# write the coords, qm esps, calculated esps and residuals
				esout_file.write("%10.5f%10.5f%10.5f%12.5f%12.5f%12.5f\n"%(xi,yi,zi,espqmi,espclc,vresid))

		chipot += chipot_mol

		if ioutopt == 1:
			rmse_mol = math.sqrt(chipot_mol/nesp)
			esout_file.write("molecule {:d}  rss: {:13.7E}  rmse: {:13.7E}\n\n".format(imol+1,chipot_mol,rmse_mol))

	espot_file.close()
	if ioutopt == 1:
		esout_file.close()

	rmse = math.sqrt(chipot/tot_nesp)
	rrmse = rmse/math.sqrt(ssvpot)

	print("chipot: ",chipot)
	return chipot, rmse, rrmse

def wrt_out():
	# called from Main
	global rmse, rrmse

	# ---- print the optimized charges and coordinates ----
	for imol in range(nmol):
		output_file.write(" molecule %d: %s"%(imol+1, title[imol]))

	# print the charges
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

	# calculate residuals sum-of-squares (chi-square) for the esp's
	chipot,rmse,rrmse = evlchi()

	# now write all these stuff out
	output_file.write("\n\n        Statistics of the fitting:")
	output_file.write("\n  The initial sum of squares (ssvpot)  {:13.7E}".format(ssvpot))
	output_file.write("\n  The residual sum of squares (RSS)    {:13.7E}".format(chipot))
	output_file.write("\n  The root-mean-squared error (RMSE)   {:13.7E}".format(rmse))
	output_file.write("\n  The relative RMSE (RRMSE)            {:13.7E}".format(rrmse))

	# ----- print the dipole, quadrupole and center of mass ----
	output_file.write("\n\n Center of Mass (Angst.):")
	for imol in range(nmol):
		output_file.write("\n #MOL          X          Y          Z")
		output_file.write("\n %3d     %10.5f %10.5f %10.5f"%(imol+1,cmas_mol[0][imol], cmas_mol[1][imol], cmas_mol[2][imol]))
	output_file.write("\n\n Dipole (Debye):")
	for imol in range(nmol):
		output_file.write("\n #MOL         D          Dx         Dy         Dz")
		output_file.write("\n %3d     %10.5f %10.5f %10.5f %10.5f"%(imol+1,dipmom_mol[imol],dipol_mol[0][imol],dipol_mol[1][imol],dipol_mol[2][imol]))
	output_file.write("\n\n Quadrupole (Debye*Angst.):")
	for imol in range(nmol):
		output_file.write("\n #MOL          X          Y          Z")
		output_file.write("\n %3d    X %10.5f"%(imol+1,quad_mol[0][imol]))
		output_file.write("\n        Y %10.5f %10.5f"%(quad_mol[3][imol],quad_mol[1][imol]))
		output_file.write("\n        Z %10.5f %10.5f %10.5f"%(quad_mol[4][imol], quad_mol[5][imol], quad_mol[2][imol]))
	output_file.write("\n\n Dipole Reoriented (Debye):")
	for imol in range(nmol):
		output_file.write("\n #MOL         D          Dx         Dy         Dz")
		output_file.write("\n %3d     %10.5f %10.5f %10.5f %10.5f"%(imol+1,dipmom_mol[imol],dipol_mol_com[0][imol],dipol_mol_com[1][imol],dipol_mol_com[2][imol]))
	output_file.write("\n\n Quadrupole Reoriented (Debye*Angst.):")
	for imol in range(nmol):
		output_file.write("\n #MOL          X          Y          Z")
		output_file.write("\n %3d    X %10.5f"%(imol+1,quad_mol_com[0][imol]))
		output_file.write("\n        Y %10.5f %10.5f"%(quad_mol_com[3][imol],quad_mol_com[1][imol]))
		output_file.write("\n        Z %10.5f %10.5f %10.5f"%(quad_mol_com[4][imol], quad_mol_com[5][imol], quad_mol_com[2][imol]))

###### get the file names ######
Input,output,qin,polariz,qout,espot,esout = file_in()

input_file = open(Input, 'r')
output_file = open(output, 'w')

###### read the atomic centers and q0's, then read the potential inf ######
read_in()

# if nmol > 1, this is a multiple molecule run
# otherwise it is a single molecule run
if nmol > 1:
	mult_mol()
else:
	sing_mol()

input_file.close()

if ipol > 0:
	read_pol_dict()

###### read in the qm esp, forming the matrices apot(awt) and bpot(bwt)
matpot()

###### initialize q0 according to iqopt
init_q0_p0()

###### process the input (freezing, equivalencing charges)
data_prep()

qcal = np.zeros((iuniq))
if ipermdip > 0:
	pcal = np.zeros((iuniq_p))

if irstrnt == 2:
	# if irstrnt= 2 then we just want to compare esp's to q0's
	for k in range(iuniq):
		qcal[k] = q0[k]
	qwt = 0.0

	if ipermdip > 0:
		for k in range(iuniq_p):
			pcal[k] = p0[k]
		pwt = 0.0
else:
	# do the charge fitting
	charge_opt()

# now calculate dipole moments
calc_dip()

###### center & reorient molecule in preparation for dipole & quadrupole ######
reornt()

# now calculate quadrupole moments
quad_mol = np.zeros((6,nmol))
quad_mol_com = np.zeros((6,nmol))
calc_quad(quad_mol, crd, dipindperm)
calc_quad(quad_mol_com, co, dipindperm_com)

# now write charges & permanent dipoles
wrt_qout()

# now calculate and print sum-of-squares, rmse, and rrmse
wrt_out()

output_file.close()


