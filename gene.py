import math
import re

global fpname
global rpname

pd_full=0
pd_extensible=1
oligo_conc = 200.00
global pd_temperature
mg_conc=0.7
dntp_conc=0.2
monovalent_cation_conc=50
pd_extensible = 1

global numlnc
numlnc=0
global repeat
repeat = 3

global prl
prl = 0
global pfl
pfl = 0
global pl
pl = 0
global dir_name
global score,rprimer_r,pos
global primer_hash,bind_string
global score_sort,d
global bind_string
global rating_hash

global bind_string 

global result


global oligo_dG
oligo_dG = {'initC': 0.98, 'initG':0.98, 'initA':1.03,	'initT':1.03}

	## deltaH (cal/K.mol)
oligo_dH={'AA': -7.9, 'TT': -7.9, 'AT': -7.2, 'TA': -7.2, 'CA': -8.5, 'TG': -8.5, 'GT': -8.4,
              'AC': -8.4, 'CT': -7.8, 'AG': -7.8, 'GA': -8.2, 'TC': -8.2, 'CG': -10.6, 'GC': -9.8,
              'GG': -8.0, 'CC': -8.0, 'initC': 0.1, 'initG': 0.1, 'initA': 2.3, 'initT': 2.3
       }
	
oligo_dH_full={'AATT': -7.9, 'TTAA': -7.9, 'ATTA': -7.2, 'TAAT': -7.2, 'CAGT': -8.5, 'TGAC': -8.5, 'GTCA': -8.4,
                  'ACTG': -8.4, 'CTGA': -7.8, 'AGTC': -7.8, 'GACT': -8.2, 'TCAG': -8.2, 'CGGC': -10.6,	'GCCG': -9.8,
                  'GGCC': -8.0, 'CCGG': -8.0, 'initC': 0.1, 'initG': 0.1, 'initA': 2.3, 'initT': 2.3,
		# Like pair mismatches
                   'AATA': 1.2, 'ATAA': 1.2, 'CAGA': -0.9,'AGAC': -0.9, 'GACA': -2.9, 'ACAG': -2.9, 'TAAA': 4.7,
                   'AAAT': 4.7, 'ACTC': 0.0, 'CTCA': 0.0, 'CCGC': -1.5, 'CGCC': -1.5, 'GCCC': 3.6, 'CCCG': 3.6,
                   'TCAC': 6.1, 'CACT': 6.1, 'AGTG': -3.1,'GTGA': -3.1, 'CGGG': -4.9, 'GGGC': -4.9, 'GGCG': -6.0,
                   'GCGG': -6.0,'TGAG': 1.6, 'GAGT': 1.6, 'ATTT': -2.7, 'TTTA': -2.7, 'CTGT': -5.0, 'TGTC': -5.0,
                   'GTCT': -2.2,'TCTG': -2.2,'TTAT': 0.2, 'TATT': 0.2,
		# G.T mismatches 
                 'AGTT': 1.0, 'TTGA': 1.0, 'ATTG': -2.5, 'GTTA': -2.5,'CGGT': -4.1, 'TGGC': -4.1, 'CTGG': -2.8, 'GGTC': -2.8,
                 'GGCT': 3.3, 'TCGG': 3.3, 'GGTT': 5.8,  'TTGG': 5.8, 'GTCG': -4.4, 'GCTG': -4.4, 'GTTG': 4.1,  'GTTG': 4.1,
                 'TGAT': -0.1,'TAGT': -0.1,'TGGT': -1.4, 'TGGT': -1.4,'TTAG': -1.3, 'GATT': -1.3, 
		# G.A mismatches 
		 'AATG': -0.6, 'GTAA': -0.6, 'AGTA': -0.7, 'ATGA': -0.7, 'CAGG': -0.7, 'GGAC': -0.7, 'CGGA': -4.0, 'AGGC': -4.0,
		 'GACG': -0.6, 'GCAG': -0.6, 'GGCA': 0.5,  'ACGG': 0.5,  'TAAG': 0.7,  'GAAT': 0.7,  'TGAA': 3.0,  'AAGT': 3.0,
		# C.T mismatches 
		 'ACTT': 0.7, 'TTCA': 0.7, 'ATTC': -1.2, 'CTTA': -1.2, 'CCGT': -0.8, 'TGCC': -0.8, 'CTGC': -1.5, 'CGTC': -1.5,
                 'GCCT': 2.3,'TCCG': 2.3, 'GTCC': 5.2, 'CCTG': 5.2,  'TCAT': 1.2,  'TACT': 1.2, 'TTAC': 1.0, 'CATT': 1.0, 
		# A.C mismatches 
		 'AATC': 2.3, 'CTAA': 2.3, 'ACTA': 5.3, 'ATCA': 5.3, 'CAGC': 1.9, 'CGAC': 1.9, 'CCGA': 0.6, 'AGCC': 0.6,
                 'GACC': 5.2, 'CCAG': 5.2, 'GCCA': -0.7,'ACCG': -0.7,'TAAC': 3.4, 'CAAT': 3.4, 'TCAA': 7.6, 'AACT': 7.6
                 }
	
	## deltaS (cal/K.mol) #	
oligo_dS={'AA': -22.2, 'TT': -22.2, 'AT': -20.4, 'TA': -21.3, 'CA': -22.7, 'TG': -22.7, 'GT': -22.4, 'AC': -22.4,
              'CT': -21.0, 'AG': -21.0, 'GA': -22.2, 'TC': -22.2, 'CG': -27.2, 'GC': -24.4, 'GG': -19.9, 'CC': -19.9,
              'initC': -2.8,'initG': -2.8, 'initA': 4.1, 'initT': 4.1, 'sym': -1.4
	}
	
oligo_dS_full={'AATT': -22.2, 'TTAA': -22.2, 'ATTA': -20.4, 'TAAT': -21.3, 'CAGT': -22.7, 'TGAC': -22.7, 'GTCA': -22.4,
                    'ACTG': -22.4, 'CTGA': -21.0, 'AGTC': -21.0, 'GACT': -22.2,	'TCAG': -22.2, 'CGGC': -27.2, 'GCCG': -24.4,
                    'GGCC': -19.9, 'CCGG': -19.9, 'initC': -2.8, 'initG': -2.8, 'initA': 4.1,  'initT': 4.1,  'sym': -1.4,		
		# Like pair mismatches	
		 'AATA': 1.7, 'ATAA': 1.7, 'CAGA': -4.2, 'AGAC': -4.2, 'GACA': -9.8, 'ACAG': -9.8, 'TAAA': 12.9, 'AAAT': 12.9, 
		 'ACTC': -4.4,'CTCA': -4.4,'CCGC': -7.2, 'CGCC': -7.2, 'GCCC': 8.9,  'CCCG': 8.9,  'TCAC': 16.4, 'CACT': 16.4, 
		 'AGTG': -9.5,'GTGA': -9.5, 'CGGG': -15.3,'GGGC': -15.3, 'GGCG': -15.8,'GCGG': -15.8, 'TGAG': 3.6, 'GAGT': 3.6, 
		 'ATTT': -10.8, 'TTTA': -10.8,'CTGT': -15.8,'TGTC': -15.8,'GTCT': -8.4 ,'TCTG': -8.4, 'TTAT': -1.5,'TATT': -1.5,
		# G.T mismatches
		 'AGTT': 0.9, 'TTGA': 0.9, 'ATTG': -8.3, 'GTTA': -8.3, 'CGGT': -11.7, 'TGGC': -11.7, 'CTGG': -8.0, 'GGTC': -8.0,
		 'GGCT': 10.4,'TCGG': 10.4,'GGTT': 16.3, 'TTGG': 16.3, 'GTCG': -12.3, 'GCTG': -12.3, 'GTTG': 9.5,  'GTTG': 9.5,
		 'TGAT': -1.7,'TAGT': -1.7,'TGGT': -6.2, 'TGGT': -6.2, 'TTAG': -5.3,   'GATT': -5.3, 
		# G.A mismatches
		 'AATG': -2.3, 'GTAA': -2.3, 'AGTA': -2.3, 'ATGA': -2.3, 'CAGG': -2.3, 'GGAC': -2.3, 'CGGA': -13.2, 'AGGC': -13.2,
                 'GACG': -1.0, 'GCAG': -1.0, 'GGCA': 3.2, 'ACGG': 3.2,	'TAAG': 0.7, 'GAAT': 0.7, 'TGAA': 7.4, 'AAGT': 7.4, 
		# C.T mismatches
		 'ACTT': 0.2, 'TTCA': 0.2,'ATTC': -6.2,	'CTTA': -6.2, 'CCGT': -4.5, 'TGCC': -4.5, 'CTGC': -6.1,	'CGTC': -6.1,
		 'GCCT': 5.4, 'TCCG': 5.4,'GTCC': 13.5,	'CCTG': 13.5, 'TCAT': 0.7,  'TACT': 0.7, 'TTAC': 0.7, 	'CATT': 0.7, 

		# A.C mismatches		
		 'AATC': 4.6, 'CTAA': 4.6, 'ACTA': 14.6, 'ATCA': 14.6, 'CAGC': 3.7, 'CGAC': 3.7, 'CCGA': -0.6,	'AGCC': -0.6,
		 'GACC': 14.2,'CCAG': 14.2,'GCCA': -3.8, 'ACCG': -3.8, 'TAAC': 8.0,'CAAT': 8.0,  'TCAA': 20.2, 	'AACT': 20.2 
	}
	
	
	# Genetic code hash
genetic_code={'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
		  'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
		  'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
		  'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
		  'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
		  'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
	}

def recalculate_dG():
    global score
    global primer_hash
    global score_sort
    global bind_string
    global rating_hash
    global bind_string 
    global result
    global oligo_dG
    # because dG = dH - TdS, and dS is dependent on the salt concentration ...
    temperature = pd_temperature

    if mg_conc > dntp_conc:
        salt_correction = math.sqrt(mg_conc - dntp_conc)
    else:
        salt_correction = 0

    na_eq = (monovalent_cation_conc + 120 * salt_correction)/1000.00

    # the length of each NN dimer is 2, therefore the modifier is 1
    entropy_adjust = (0.368 * math.log(na_eq))
    for key in oligo_dH_full:
        if 'init' in key:
            continue
        dS = oligo_dS_full[key] + entropy_adjust
        dG = oligo_dH_full[key]-((273.15 + temperature)*(dS/1000.00))
        oligo_dG[key] = dG

def check_degenerate(seq, l_or_m):
    flag = False
    if l_or_m == 'l':
        for x in seq.strip().upper():
            if x not in ('A','G','T','C'):
                print "One of your lncRNA sequences has a degenerate or non-DNA character,please check\n"
                flag = True
                break
    if l_or_m== 'm':
        for x in seq.strip().upper():
            if x not in ('A','G','T','C'):
                print "One of your mRNA sequences has a degenerate or non-DNA character,please check\n"
                flag = True
                break
    return flag
        
def tm(rna):
    primer = rna
    primer = primer.upper().strip()
    primer_len = len(primer)
    i = 0

    # calculate deltaH and deltaS 
    deltaH = 0
    deltaS = 0
    for i in range(primer_len-1):
        nn = primer[i:i+2]
        deltaH += oligo_dH[nn]
        deltaS += oligo_dS[nn]
    # initial term correction
    initterm = "init" + primer[0:1]
    deltaH += oligo_dH[initterm]
    deltaS += oligo_dS[initterm]
    
    endterm = "init" + primer[-1:]
    deltaH += oligo_dH[endterm]
    deltaS += oligo_dS[endterm]

    # Tm at 1M NaCl
    # tm= (deltaH * 1000) / (deltaS + (1.987 * log(oligo_conc / 4))) - 273.15

    # correct for salt concentration on deltaS

    # Big problems if [dNTPs] > [Mg++] !!  This is a quick fix ...
    if mg_conc > dntp_conc:
        salt_correction = math.sqrt(mg_conc - dntp_conc)
    else:
        salt_correction = 0
    
    na_eq = (monovalent_cation_conc + 120 * salt_correction)/1000.00
    deltaS += (0.368 * (primer_len - 1) * math.log(na_eq))
    oligo_conc_mols = oligo_conc / 1000000000.00
    
    # Salt corrected Tm
    # NB - for PCR I'm assuming for the moment that the [strand target] << [oligo]
    # and that therefore the C(t) correction term approx equals [oligo]
    corrected_tm=((deltaH * 1000.00) / (float(deltaS) + (1.987 * math.log(float(oligo_conc_mols)/4)))) - 273.15 

    return (corrected_tm, deltaH, deltaS)

def gc(rna):
    rna = rna.strip().upper()
    gc = 0
    countgc = 0
    for x in rna:
        if x == 'G' or x == 'C':
            countgc += 1
    counttotal = len(rna)
    gc = float(countgc)/float(counttotal)*100.00
    return gc

def complement(seq):
    new_seq = []
    for x in seq.upper().strip():
        if x == 'A':
            new_seq.append('T')
        if x == 'T':
            new_seq.append('A')
        if x == 'C':
            new_seq.append('G')
        if x == 'G':
            new_seq.append('C')
    return str(''.join(new_seq))

def draw_dimer(primer_f, primer_r, pos,pic, output):
    rprimer_r=primer_r[::-1]
    dimer_binding=""
    pr_space=""
    fspace=""
    rspace=""
    tang=0
    shan=0
    tianjin=0
    beijing=0
    hongqiao=0
    nankai=0
    heping=0
    langfang=0
    global prl,pl,pfl
    global score
    global primer_hash
    global score_sort
    global bind_string
    global rating_hash
    global bind_string 
    global result
    global oligo_dG
    
    rspace_len = pfl - pos - 1
    fspace_len = 0
    
    if pos >= pfl-1:
        fspace_len = pos - pfl + 1
        rspace_len = 0

    bind_len = rspace_len
    pos_first = -1
    pos_end = -1
    for j in range(pos+1):
        if bind_string[pos][j:j+1] == '1':
            dimer_binding=dimer_binding+"|"
            if pos_first == -1:
                pos_first = j
            pos_end = j
        elif bind_string[pos][j:j+1] == '0':
            dimer_binding=dimer_binding+"x"
        elif bind_string[pos][j:j+1] == '3':
            dimer_binding=dimer_binding+"."
            pos_end = j
        else:
            dimer_binding=dimer_binding+" "

    fspace = " " * fspace_len
    rspace = " " * rspace_len
    bind_space = " " * bind_len

    
    fstart = rspace_len + pos_first + 1 - fspace_len
    fend= rspace_len + pos_end + 1 - fspace_len

    rstart = pos_first + 1
    rend = pos_end + 1
    
    result_str = fspace+"5' "+primer_f+" 3'\n"+ \
             bind_space+"   "+dimer_binding+"\n"+ \
                 rspace+"3' "+rprimer_r+" 5'\n\n"
    
    global dir_name
    logFile = open(dir_name+"/"+output,'a')
    logFile.write(str(fstart)+'\t'+str(fend)+'\t'+str(rstart)+'\t'+str(rend)+'\n')
    logFile.close()
 
    if pic == 'T' or pic == 't':
        logFile = open(dir_name+"/"+output,'a')
        logFile.write(result_str)
        logFile.close()

def primer_dimer(primer_f, primer_r, pd_full):
    global prl,pl,pfl
    global score,rprimer_r,bind_string
    global primer_hash
    primer_hash = {}
    primer_hash['a'] = []
    primer_hash['g'] = []
    primer_hash['c'] = []
    primer_hash['t'] = []
    global score_sort
    global bind_string,rating_hash
    bind_string = {}
    rating_hash = {} 
    global result
    global oligo_dG
    
    if primer_f =='' or primer_r == '':
        return
    	
    # pl = greatest length
    pfl=len(primer_f)
    prl=len(primer_r)
    pl =  (pfl if pfl > prl else prl)
    rcompr = (complement(primer_r))[::-1]
    rcomprlc = rcompr.lower()
    fprimer_r = primer_f[::-1].lower()
    rprimer_r = primer_r[::-1]
	
    # create a binding array for each of the four bases
    l = 0
    for l in range(pfl):
        mbase = fprimer_r[l:l+1]
        for k in ('a','g','c','t'):
            if k == mbase:
                primer_hash[k].append(1)
            else:
                primer_hash[k].append(0)
     
    # create the primer matrix
    primer_comp = {}
    for k in range(prl):
        primer_comp[k] = []
        primer_comp[k]=primer_hash[rcomprlc[k:k+1]]
    		
    # read each combination from the matrix, calculate dG for each dimer
    pd_len = (pfl+prl-1 if pd_full > 0 else pl-2)
    
    score = {}
    for k in range(pd_len+1):
        score[k]=0
        bind = ''
        score_p = 0
        # extensible primer short-circuit - ignore all primers that will
	# not create extensible (i.e. amplifiable) dimers
        start = pfl-1 if k > pfl-1 else k
        end = prl-1 if k > prl-1 else k
        if pd_extensible and not pd_full:
            if primer_comp[0][start] == 1:
                pass
            else:
                continue
            if primer_comp[end][start-k] == 1:
                pass
            else:
                continue
		
	    # read the binding data
        for l in range(prl):
            if (k-l) < pfl:
                if (k-l) >= 0:
                    bind = bind + str(primer_comp[l][k-l])
                else:
                    break
            else:
                # spacer
                bind = bind + "2"   		
	# Single matched bases surrounded by mismatches are unstable,
	# so we remove them with the regexp (look ahead is needed otherwise
	# strings of consecutive match/mismatches are not caught)
        p = '01(?=[^1])'
        bind = re.sub(p, "03", bind) 
	
	# Short circuit if there's nothing to bind
        if '1' in bind:
            pass
        else:
            continue
	
	# Find start and end of similarity
	# my (pb_init,pb_end)
        pb_init = -1
        
        for l in range(len(bind)):
            # at first I tried finding the initiating terminal bases with
            # regexps, but that was much slower ...
            if bind[l:l+1]=="1":
                if pb_init == -1:
                    pb_init = l
                pb_end = l
				
        if pb_init != -1:
            # deltaG calculation
            for l in range(pb_init, pb_end):
                if bind[l:l+2] == "00" or bind[l:l+2] == "30" or bind[l:l+2] == "03":
                    continue
                if bind[l:l+1] == "2":
                    continue
                score_p += oligo_dG[primer_f[pfl-k+l-1:pfl-k+l+1] + rprimer_r[l:l+2]]
 
	    # init term corrections
            initterm="init" + rprimer_r[pb_init:pb_init+1]
            score_p += oligo_dG[initterm]
	    
            endterm="init" + rprimer_r[pb_end:pb_end+1]
            score_p += oligo_dG[endterm]

            # add to the hash ...
            score[k]=round(float(score_p), 2)
            bind_string[k]=bind
            rating_hash[score[k]]=k
    
    # sort the dimers to give the most stable:
    score_list = []
    for score_key in score:
        score_list.append(score[score_key])
    score_sort = sorted(score_list)
    # Returns the most stable dimer
##    print score_sort[0]
    return score_sort[0]

def get_tm(fprimer, rprimer, pd_full, pic, output):
    global dir_name
    global prl,pl,pfl,numlnc,repeat
    global score,d,pos
    global primer_hash
    global score_sort
    global bind_string
    global rating_hash
    global bind_string
    global result
    global oligo_dG
    
    fprimer = fprimer.upper().strip()
    rprimer = rprimer.upper().strip()
    
    oligo_conc_mols = oligo_conc / 1000000000.00   
    if (fprimer != '') and (not check_degenerate(fprimer,'l')):
        (fprimer_tm, deltaH, deltaS) = tm(fprimer)
        fprimer_gc = int(gc(fprimer))
        fprimer_tm = round(float(fprimer_tm), 2)
        deltaG = deltaH-((273.15 + pd_temperature)*(float(deltaS)/1000.00))
        fprimer_ds = round(float(deltaS), 2)
        fprimer_dh = round(float(deltaH), 2)
        fprimer_dg = round(float(deltaG), 2)
        fprimer_len = len(fprimer)
    else:
        exit
    
    if (rprimer != '') and (not check_degenerate(rprimer,'m')):
        (rprimer_tm, deltaH, deltaS) = tm(rprimer)
        rprimer_gc = int(gc(rprimer))
        rprimer_tm = round(float(rprimer_tm), 2)
        deltaG = deltaH-((273.15 + pd_temperature)*(float(deltaS)/1000.00))
        rprimer_ds = round(float(deltaS), 2)
        rprimer_dh = round(float(deltaH), 2)
        rprimer_dg = round(float(deltaG), 2)
        rprimer_len = len(rprimer)
    else:
        exit
    
####    print "lncRNA is:" + fpname + "\tlength:" + str(fprimer_len) + "\tmRNA is:" + rpname + "\tlength:" + str(rprimer_len) + "\n"
##
##    if rprimer_len >= 9000:
##        step= int(rprimer_len/5)
####        print "%%%%%%%%%%%%%Step:" + step + "\n"
##        zch_1 = rprimer[0:step]
##        primer_dimer(fprimer,zch_1,1)
##        pos_1=rating_hash[score_sort[0]]
##        score_1=score_sort[0]
##
##        start = step - fprimer_len
##        end = 2*step
##        zch_2 = rprimer[start:start+end]
##        primer_dimer(fprimer,zch_2,1)
##        pos_2=rating_hash[score_sort[0]]
##        score_2=score_sort[0]
##
##        start=2*step-fprimer_len
##        end=3*step
##        zch_3=rprimer[start:start+end]
##        primer_dimer(fprimer,zch_3,1)
##        pos_3=rating_hash[score_sort[0]]
##        score_3=score_sort[0]
##
##        start=3*step-fprimer_len
##        end=4*step
##        zch_4=rprimer[start:start+end]
##        primer_dimer(fprimer,zch_4,1)
##        pos_4=rating_hash[score_sort[0]]
##        score_4=score_sort[0]
##
##        start=4*step-fprimer_len
##        zch_5=rprimer[start,start+rprimer_len]
##        primer_dimer(fprimer,zch_5,1)
##        pos_5=rating_hash[score_sort[0]]
##        score_5=score_sort[0]
##
##        if score_1<=score_2 and score_1<=score_3 and score_1<=score_4 and score_1<=score_5:
##            score_sort[0]=score_1  
##            pos= pos_1  
##        elif score_2<=score_1 and score_2<=score_3 and score_2<=score_4 and score_2<=score_5:
##            score_sort[0]=score_2 
##            pos= pos_2   
##        elif score_3<=score_1 and score_3<=score_2 and score_3<=score_4 and score_3<=score_5:
##            score_sort[0]=score_3
##            pos= pos_3   
##        elif score_4<=score_1 and score_4<=score_3 and score_4<=score_2 and score_4<=score_5:
##            score_sort[0]=score_4
##            pos= pos_4   
##        else:
##            score_sort[0]=score_5
##            pos= pos_5
##
##        if fprimer_len >= rprimer_len:
##            min_leng=rprimer_len
##        else:
##            min_leng=fprimer_len
##
##        sdelta = score_sort[0]/float(min_leng)
##
##        if score_sort[0] < 0:
##            if sdelta<=d or sdelta > d:
##                logFile = open("/var/www/html/lrt/jobs/"+dir_name+"/"+output,'a')
##                logFile.write("Query:"+fpname+"\tLength:"+str(fprimer_len)+"\t")
##                logFile.write("Target:"+rpname+"\tLength:"+str(rprimer_len)+"\t")
##                logFile.write("dG:"+str(score_sort[0])+"\tStandard deltaG:"+str(sdelta)+"\n")
##                logFile.close()
##                if pic == 'T' or pic == 't':
##                    draw_dimer(fprimer, rprimer, pos,pic, output)
##            numlnc += 1	
##        return
    
    if rprimer and not check_degenerate(rprimer,'m'):
        primer_dimer(fprimer,rprimer,1)
        score_1=score_sort[0]
        pos_a=rating_hash[score_sort[0]]

        score_sort[0] = score_1
        pos=pos_a
        if fprimer_len>=rprimer_len:
            min_leng=rprimer_len
        else:
            min_leng=fprimer_len
        sdelta=score_sort[0]/min_leng

        if score_sort[0]<0:
            if sdelta<=d or sdelta > d:
                logFile = open(dir_name+"/"+output,'a')
                logFile.write(fpname+ "\t"+ str(fprimer_len)+"\t")
                logFile.write(rpname+ "\t"+str(rprimer_len)+"\t")
                logFile.write(str(score_sort[0])+"\t"+str(sdelta)+"\t")
                logFile.close()
                draw_dimer(fprimer,rprimer,pos,pic,output)
            numlnc += 1

	
def Calculate_G(Dir, lnc_name, fprimer, m_name, rprimer, dG, pic, output, Tempa):
    global fpname,rpname,d, dir_name,pd_temperature
    fpname=lnc_name
    rpname=m_name
    dir_name = Dir
    d = dG
    pd_temperature=float(Tempa)
    recalculate_dG()
    x = get_tm(fprimer, rprimer, pd_full, pic, output)
    return x




