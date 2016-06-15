#def energySet(): 
'''For this to work you need to install this script and your A domain (chainAA.pdb) in the same folder.
It will create one .pdb file and one modified .resfile for each A domain mutated as per Stachelhaus code,
and output to binding.csv and nonbinding.csv. '''
from rosetta import *
init()
from toolbox import *
import csv
import io
import random
import math
import D070_Refinement
from D070_Refinement import sample_refinement
aminoacids = ["PHE","ALA","ASN","ASP","CYS","GLN","GLU","ILE","LEU","PRO","SER","THR","TYR","VAL"]
oneletter = {'ALA': 'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS_D':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
#Used for testing:
#aminoacids = ["PHE","ALA","ASN","ASP","CYS","PRO"]
scorefxn = get_fa_scorefxn()
domain = pose_from_pdb("chainAA.pdb") #PHE domain
mutDom = Pose()
mutDom.assign(domain)
positions = [235,236,239,278,299,301,322,330,331,517]
seqs = {"PHE":"DAWTIAAVCK","ALA":"DLLFGIAVLK","ASN":"DLTKLGEVGK","ASP":"DLTKVGHIGK","CYS":"DHESDVGITK","GLN":"DAQDLGVVDK","GLU":"DAWHFGGVDK","ILE":"DGFFLGVVYK","LEU":"DAWFLGNVVK","PRO":"DVQLIAHVVK","SER":"DVWHLSLIDK","THR":"DFWNIGMVHK","TYR":"DGTITAEVAK","VAL":"DAFWIGGTFK"}
boundScore = [] #list of energy scores for binding between domain + its natural aminoacid
unboundScore = [] #list of energy scores for binding between domain + unnatural aminoacid (this will be considered poor/no binding)

from toolbox import generate_resfile_from_pdb
generate_resfile_from_pdb("chainAA.pdb","chainAA.resfile")

amIndex = 0 #aminoacid + seqs index

for currentAm in aminoacids:
	data = []
	resFilename =""
	if currentAm == "PHE":
		print "**SCORING NATIVE PHE ENERGY**"
		values = sample_refinement("chainAA.pdb",jobs = 5, job_output = "PHE_binding_job")
		boundScore.append((sum(values))/5)
	else:
		print "**WORKING ON AMINOACID: %s ** " % currentAm
		
		'''BINDING'''
		#changes ligand
		
		print "**READING FROM FILE: chainAA.resfile**"
		with open("chainAA.resfile","r+") as resFile:
			data = resFile.readlines()
			end = len(data) - 1
			data[end] = " 566 A PIKAA %s \n" % oneletter[currentAm]
		
		resFilename = currentAm + ".resfile"
		print "**WRITING TO NEW FILE: %s **" % resFilename
		with open(resFilename,"w") as newResFile:
			newResFile.writelines(data)
		
		print "**SCORING BINDING ENERGY**"
		#repacks
		task_design = TaskFactory.create_packer_task(mutDom)
		parse_resfile(mutDom,task_design,resFilename)
		packMover = PackRotamersMover(scorefxn,task_design)
		packMover.apply(mutDom)
		# minMover = MinMover()
		# mmap = MoveMap()
		# ligand = mutDom.pdb_info().pdb2pose('A',566)
		# mmap.set_bb_true_range(ligand,ligand)
		# minMover.movemap(mmap)
		# minMover.score_function(scorefxn)
		# minMover.apply(mutDom)
		
		#mutate to specificity code
		# pdbPosn = domain.pdb_info().pdb2pose('A',positions[0])
		# mutDom = mutate_residue(domain,pdbPosn,seqs[currentAm][0],8.0)
		# task_pack = standard_packer_task(mutDom)
		# task_pack.restrict_to_repacking()
		# task_pack.temporarily_fix_everything()
		# task_pack.temporarily_set_pack_residue(pdbPosn,True)
		# pack_mover = PackRotamersMover(scorefxn,task_pack)
		# pack_mover.apply(mutDom)
		
		#Check that the mutation works! Go through one of the structures,
		#check that residues noted in [positions] actually are what they're
		#supposed to be.
		#pdbPosn = pose.pdb_info().pdb2pose('A',position) -> gives you the pdb position number of the residue, may differ!
		#print pose.residue(pdbPosn).name() -> gives name of residue at position.
		for posn in range(len(positions)):
			#repacks for each mutation
			pdbPosn = domain.pdb_info().pdb2pose('A',positions[posn])
			mutDom = mutate_residue(mutDom,pdbPosn,seqs[currentAm][posn],8.0) 
			# task_pack = standard_packer_task(mutDom)
			# task_pack.restrict_to_repacking()
			# task_pack.temporarily_fix_everything()
			# task_pack.temporarily_set_pack_residue(pdbPosn,True)
			# pack_mover = PackRotamersMover(scorefxn,task_pack)
			# pack_mover.apply(mutDom)
			# mmap.set_bb_true_range(pdbPosn,pdbPosn)
			# minMover.movemap(mmap)
			# minMover.score_function(scorefxn)
			# minMover.apply(mutDom)
		print "**DOMAIN MUTATED FOR %s SPECIFICITY**" % currentAm
		# mmap.set_bb(True)
		# minMover.movemap(mmap)
		# minMover.score_function(scorefxn)
		# minMover.apply(mutDom)
		#variables for pdb filenames
		pdbname = currentAm + "specific.pdb"
		jobname = currentAm + "_binding_job"
		#Check mutations from these files produced here!
		mutDom.dump_pdb(pdbname) #creates file for new specific domain (post mutation)
		values = sample_refinement(pdbname,jobs = 5, job_output = jobname)
		
		boundScore.append((sum(values))/5)
		
		
	'''NONBINDING'''
	# #generate score for non-binding structures
	nonBinding = list(aminoacids)
	nonBinding.remove(currentAm)
	if nonBinding != []:
		if currentAm == "PHE":
			mutDom.assign(domain)
		badIndex = 0
		for badAm in nonBinding:
			#changes ligand
			print "**READING FROM FILE: chainAA.resfile**"
			with open("chainAA.resfile","r+") as resFile:
				data = resFile.readlines()
				end = len(data) - 1
				data[end] = " 566 A PIKAA %s \n" % oneletter[badAm]
		
			badResFilename = badAm + ".resfile"
			print "**WRITING TO NEW FILE: %s **" % badResFilename
			with open(badResFilename,"w") as newResFile:
				newResFile.writelines(data)
			
			print "**SCORING NONBINDING ENERGY OF %s DOMAIN WITH %s" % (currentAm, badAm)
			task_design = TaskFactory.create_packer_task(mutDom)
			parse_resfile(mutDom,task_design,badResFilename)
			packMover = PackRotamersMover(scorefxn,task_design)
			packMover.apply(mutDom)
			badPdbname = currentAm + "with" + badAm + ".pdb"
			badJobname = currentAm + "_nonbinding_" + badAm + "_job"
			mutDom.dump_pdb(badPdbname)
			badValues = sample_refinement(badPdbname,jobs = 5, job_output = badJobname)
		
			unboundScore.append((sum(badValues))/5)
			
			badIndex = badIndex + 1
			
	amIndex = amIndex + 1
# domain.dump_pdb("PHEspecific.pdb")
# mutatedList = []
# index = 0

# for k in aminoacids:
	# filename = k +"specific.pdb"
	# pose = pose_from_pdb(filename)
	# newSeqs = []
	# for m in range(len(positions)):
		# pdbPosn = domain.pdb_info().pdb2pose('A',positions[m])
		# newSeqs.append(dict[pose.residue(pdbPosn).name()])
	# mutatedList.append(newSeqs)
	# index = index + 1
		
		
	
'''GENERATE TST AND TRN SETS'''
#print "**SHUFFLING AND GENERATING TRN AND TST SETS**"
'''code below used to divide into trn/tst'''
# allSet = boundScore + unboundScore
# random.shuffle(allSet)
# trnNum = int(math.floor((2* len(allSet))/3))
# trnSet = zip(*[iter(allSet)]*trnNum)
# allSet.reverse()
# temp = zip(*[iter(allSet)] * (len(allSet)-trnNum))
# tstSet = list(temp[0])
# tstSet.reverse()
# sets = (trnSet,tstSet)
print "***DONE*** \n"
#print "***boundScore/unboundScore: original sorted values*** \n ***tstSet/trnSet: shuffled + divided at trnNum*** \n ***sets: combined***"
#return sets
with open("binding.csv","w") as out:
	csv.writer(out).writerow(boundScore)
	
with open("nonbinding.csv","w") as out2:
	csv.writer(out2).writerow(unboundScore)
			
