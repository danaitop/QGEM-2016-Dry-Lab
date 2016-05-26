def trnSet():
	import io
	aminoacids = ['A','N','D','C','Q','E','I','L','O','F','P','S','T','Y','V']
	scorefxn = get_fa_scorefxn()
	domain = pose_from_pdb("chainA.clean.pdb") #empty PHE domain
	seqs = ["DLLFGIAVLK","DLTKLGEVGK","DLTKVGHIGK","DHESDVGITK","DAQDLGVVDK","DAWHFGGVDK","DGFFLGVVYK","DAWFLGNVVK","DMENLGLINK","DAWTIAAVCK","DIFLLGLLCK","DFQLLGVAVK","DVWHLSLIDK","DALVTGAVVK"]
	positions = [235,236,239,278,299,301,322,330,331,517]
	boundScore = [] #list of energy scores for binding between domain + its natural aminoacid
	unboundScore = [] #list of energy scores for binding between domain + unnatural aminoacid (this will be considered poor/no binding)
	
	amIndex = 0 #aminoacid + seqs index
	for currentAm in aminoacids:
		#generate .params and .pdb files for current aminoacids
		structName = currentAm + ".sdf"
		run molfile_to_params.py structName -n currentAm
		filename = currentAm + ".params"
		resiSet = generate_nonstandard_residue_set(filename)
		
		#mutate to specificity code
		pdbPosn = domain.pdb_info().pdb2pose('A',positions(0))
		mutDom = mutate_residue(domain,pdbPosn,seqs(0)(0),8.0)
		
		for posn in range(1,len(positions)):
		
			pdbPosn = domain.pdb_info().pdb2pose('A',positions(posn))
			
			mutDom = mutate_residue(mutDom,pdbPosn,seqs(amIndex)(posn),8.0) 
			
		#variables for pdb filenames
		pdbname = currentAm + "bound.pdb"
		aminame = currentAm + "_0001.pdb"
		mutDom.dump_pdb(pdbname) #creates file for new specific domain (post mutation)
		
		#adds ligand (current aminoacid) lines to end of .pdb file, as per the tutorial
		file mutPdb = open(pdbname,a+,1)
		file aminPdb = open(aminame,r+,1)
		amin = aminPdb.read() #check if this works
		mutPdb.write(amin)
		mutPdb.close()
		aminPdb.close()
		
		#create domain+ligand pose, score it, save score in boundscore list
		boundDom = Pose()
		pose_from_pdb(boundDom,resiSet,pdbname)
		
		boundscore(amIndex) = scorefxn(boundDom)
		
		#generate score for non-binding structures
		nonBinding = aminoacids.remove(currentAm)
		for n in range (len(nonBinding))
			
