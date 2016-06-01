#def energySet():
from rosetta import *
init()
from toolbox import *
import io
import random
import math
aminoacids = ["PHE","ALA","ASN","ASP","CYS","GLN","GLU","ILE","LEU","PRO","SER","THR","TYR","VAL"]
oneletter = ["F","A","N","D","C","Q","E","I","L","P","S","T","Y","V"]
#oneletter = ["F","A","N","D","C","Q"]
#aminoacids = ["PHE","ALA","ASN","ASP","CYS","GLN"]
scorefxn = get_fa_scorefxn()
domain = pose_from_pdb("chainAA.pdb") #PHE domain
seqs = ["DLLFGIAVLK","DLTKLGEVGK","DLTKVGHIGK","DHESDVGITK","DAQDLGVVDK","DAWHFGGVDK","DGFFLGVVYK","DAWFLGNVVK","DMENLGLINK","DAWTIAAVCK","DIFLLGLLCK","DFQLLGVAVK","DVWHLSLIDK","DALVTGAVVK"]
positions = [235,236,239,278,299,301,322,330,331,517]
boundScore = [] #list of energy scores for binding between domain + its natural aminoacid
unboundScore = [] #list of energy scores for binding between domain + unnatural aminoacid (this will be considered poor/no binding)

from toolbox import generate_resfile_from_pdb
generate_resfile_from_pdb("chainAA.pdb","chainAA.resfile")

amIndex = 0 #aminoacid + seqs index
for currentAm in aminoacids:

	if currentAm == "PHE":
		print "**SCORING NATIVE PHE ENERGY**"
		boundScore.append(scorefxn(domain))
	else:
		print "**WORKING ON AMINOACID: %s ** " % currentAm
		
		#mutate to specificity code
		pdbPosn = domain.pdb_info().pdb2pose('A',positions[0])
		mutDom = mutate_residue(domain,pdbPosn,seqs[0][0],8.0)
		
		for posn in range(1,len(positions)):
		
			pdbPosn = domain.pdb_info().pdb2pose('A',positions[posn])
			mutDom = mutate_residue(mutDom,pdbPosn,seqs[amIndex][posn],8.0) 
		print "**DOMAIN MUTATED FOR %s SPECIFICITY**" % currentAm
		
		#variables for pdb filenames
		pdbname = currentAm + "specific.pdb"
		mutDom.dump_pdb(pdbname) #creates file for new specific domain (post mutation)
		
		'''PROPER BINDING VALUES'''
		print "**READING FROM FILE: chainAA.resfile**"
		with open("chainAA.resfile","r+") as resFile:
			data = resFile.readlines()
			end = len(data) - 1
			data[end] = " 566 A PIKAA %s \n" % oneletter[amIndex]
		
		resFilename = currentAm + ".resfile"
		print "**WRITING TO NEW FILE: %s **" % resFilename
		with open(resFilename,"w") as newResFile:
			newResFile.writelines(data)
		
		print "**SCORING BINDING ENERGY**"
		task_design = TaskFactory.create_packer_task(mutDom)
		parse_resfile(mutDom,task_design,resFilename)
		packMover = PackRotamersMover(scorefxn,task_design)
		packMover.apply(mutDom)
		boundScore.append(scorefxn(mutDom))
		
		#generate .params and .pdb files for current aminoacid
		##structName = currentAm + ".sdf"
		##run molfile_to_params.py structName -n currentAm
		#filename = currentAm + ".params"
		#print "**GENERATING RESISET FOR %s **" % currentAm
		#resiSet = generate_nonstandard_residue_set([filename])
		
		
		# #adds ligand (current aminoacid) lines to end of .pdb file, as per the tutorial
		# print "**ADDING %s TO DOMAIN'S .PDB FILE**" % currentAm
		# aminame = currentAm + "_0001.pdb"
		# mutPdb = open(pdbname,"a+")
		# aminPdb = open(aminame,"r+")
		# amin = aminPdb.read() #check if this works
		# mutPdb.write(amin)
		# mutPdb.close()
		# aminPdb.close()
		
		# print "**SCORING BINDING**"
		# #create domain+ligand pose, score it, save score in boundscore list
		# chm = rosetta.core.chemical.ChemicalManager.get_instance()
		# rts = chm.residue_type_set('fa_standard')
		# amSet = rosetta.core.conformation.ResidueFactory.create_residue(rts.name_map(currentAm))
		# boundDom = Pose()
		# pose_from_pdb(boundDom,rts,pdbname)
		# boundScore.append(scorefxn(boundDom))
	
	'''POOR BINDING VALUES'''
	#generate score for non-binding structures
	nonBinding = list(aminoacids)
	nonBinding.remove(currentAm)
	onebadletter = list(oneletter)
	onebadletter.remove(oneletter[amIndex])
	if nonBinding != []:
		if currentAm == "PHE":
			mutDom = Pose()
			mutDom = domain
		badIndex = 0
		for badAm in nonBinding:
		
			print "**READING FROM FILE: chainAA.resfile**"
			with open("chainAA.resfile","r+") as resFile:
				data = resFile.readlines()
				end = len(data) - 1
				data[end] = " 566 A PIKAA %s \n" % onebadletter[badIndex]
		
			badResFilename = badAm + ".resfile"
			print "**WRITING TO NEW FILE: %s **" % badResFilename
			with open(badResFilename,"w") as newResFile:
				newResFile.writelines(data)
			
			print "**SCORING NONBINDING ENERGY OF %s DOMAIN WITH %s" % (currentAm, badAm)
			task_design = TaskFactory.create_packer_task(mutDom)
			parse_resfile(mutDom,task_design,badResFilename)
			packMover = PackRotamersMover(scorefxn,task_design)
			packMover.apply(mutDom)
			unboundScore.append(scorefxn(mutDom))
			#print "**GENERATE RESISET FOR NONBINDING %s **" % badAm
			#generate .params and .pdb files for current aminoacids
			##badStruct = badAm + ".sdf"
			##run molfile_to_params.py structName -n badAm
			#badFilename = badAm + ".params"
			#badResiSet = generate_nonstandard_residue_set([badFilename])
			
			#print "**ADDING %s TO DOMAIN'S .PDB**" % badAm
			#adds ligand (current aminoacid) lines to end of .pdb file, as per the tutorial
			# badName = badAm + "_0001.pdb"
			# badMutPdb = open(pdbname,"a+")
			# badAminPdb = open(badName,"r+")
			# badAmin = badAminPdb.read() #check if this works
			# badMutPdb.write(amin)
			# badMutPdb.close()
			# badAminPdb.close()
			
			# print "**SCORING NONBINDING**"
			# #create domain+ligand pose, score it, save score in boundscore list
			# badchm = rosetta.core.chemical.ChemicalManager.get_instance()
			# badrts = badchm.residue_type_set('fa_standard')
			# badamSet = rosetta.core.conformation.ResidueFactory.create_residue(badrts.name_map(badAm))
			# unboundDom = Pose()
			# pose_from_pdb(unboundDom,badrts,pdbname)
			# unboundScore.append(scorefxn(unboundDom))
			badIndex = badIndex + 1
			
	amIndex = amIndex + 1
	
'''GENERATE TST AND TRN SETS'''
print "**SHUFFLING AND GENERATING TRN AND TST SETS**"
allSet = boundScore + unboundScore
random.shuffle(allSet)
trnNum = int(math.floor((2* len(allSet))/3))
trnSet = zip(*[iter(allSet)]*trnNum)
allSet.reverse()
temp = zip(*[iter(allSet)] * (len(allSet)-trnNum))
tstSet = list(temp[0])
tstSet.reverse()
sets = (trnSet,tstSet)
print "***DONE*** \n"
print "***boundScore/unboundScore: original sorted values*** \n ***tstSet/trnSet: shuffled + divided at trnNum*** \n ***sets: combined***"
#return sets
'''Write to csv files for analysis'''
with open("binding.csv","w") as out:
	csv.writer(out).writerow(boundScore)
	
with open("nonbinding.csv","w") as out2:
	csv.writer(out2).writerow(unboundScore)
			
