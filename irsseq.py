# coding: UTF-8
import sys
import numpy as np
from Bio import Seq,SeqIO

if len(sys.argv) != 5:
    print("Usage: python3 gtfbamcount.py [bamcount file list: control] [bamcount file list: treated] [genome FASTA file name] [output file name]")
    print("Bamcount filenames are delimited by commas without spaces.")
    print("Example: python3 gtfbamcount.py wt1.bam.txt,wt2.bam.txt,wt3.bam.txt ko1.bam.txt,ko2.bam.txt,ko3.bam.txt genome.fa outfile.txt")
    exit(1)
    
ctrlist = sys.argv[1].split(",")
trelist = sys.argv[2].split(",")
fastafi = sys.argv[3]
outf = sys.argv[4]

def readfiles(ctrlist:str,trelist:str) -> tuple: #input filenames delimited by commas, return tuple of file content list [replicates: introns: contents]
    clist = []
    for fn in ctrlist:
        ctfi = []
        with open(fn) as ctlifh:
            for ctli in ctlifh:
                ct = ctli[:-1].split("\t")
                if len(ct) == 15 and ct[0] != "CHR":
                    ctfi.append(ct)
        clist.append(ctfi)
        
    tlist = []
    for fn in trelist:
        ctfi = []
        with open(fn) as ctlifh:
            for ctli in ctlifh:
                ct = ctli[:-1].split("\t")
                if len(ct) == 15 and ct[0] != "CHR":
                    ctfi.append(ct)
        tlist.append(ctfi)
    print("Count files read!")
    return (clist,tlist)

def getfastadic(fname:str) -> dict: #return dict(sequences))
    chrdict = {}
    for farec in SeqIO.parse(fname,"fasta"):
        chrdict[farec.id] = str(farec.seq)
    
    print("Fasta file read!")
    return chrdict


def dupcor(corlist): #count the frequency of alternative sites and return the total number and the top three
    altc = 0 #number of alternative reads
    sdict = {}
    for cl in corlist:
        if cl != "":
            altc += 1
            try:
                sdict[str(cl)] += 1
            except KeyError:
                sdict[str(cl)] = 1
    if altc == 0:
        return ""
    rlis = [] #[0]=altc, from [1:] is list of tuples (coordinate, number)
    for lname, lnumb in list(sdict.items()):
        rlis.append([lname,float(lnumb)/altc])
    rlis.sort(key=lambda x:x[1],reverse=True)
    return rlis



def getcalccounts(clist, tlist, fasdic):
    if min([len(i) for i in clist]+[len(i) for i in tlist]) != max([len(j) for j in clist]+[len(j) for j in tlist]):
        print("Number of introns not equal")
        input()
        
    retlist = []
    finlist = []

    outcols = ["ID","CHR","INTRST","INTREN","STRAND","GENEID","GENENAME","TRANID","INTRNUM", \
    "s5p","s3p","minreads", \
    "conirsme","treirsme","conirssd","treirssd","irsdiff","irszscore", \
    "n5pe3","n5pe2","n5pe1","n5pi4","seq20i","seq5p15","seq3p23","conpcsme","trepcsme"]

#    outcols = ["ID","CHR","INTRST","INTREN","STRAND","GENEID","GENENAME","TRANID","EXNUM", \
#    "s5p","s3p","minreads", \
#    "conirsme","treirsme","conirssd","treirssd","irsdiff","irszscore", \
#    "cona5sme","trea5sme","cona5ssd","trea5ssd","a5sdiff","a5szscore", \
#    "cona3sme","trea3sme","cona3ssd","trea3ssd","a3sdiff","a3szscore", \
#    "cona5pos","trea5pos","cona3pos","trea3pos","n5pe3","n5pe2","n5pe1","n5pi4","seq20i","seq5p15","seq3p23","seqalt5","seqalt3","conpseme","trepseme"]

    for i in ["eijr","iejr","csr","a5r","a3r"]:
        for j in range(len(clist)):
            outcols.append(i+"con"+str(j+1))
        for j in range(len(tlist)):
            outcols.append(i+"tre"+str(j+1))
    
    #for i in ["irs","a5s","a3s"]:
    for i in ["irs"]:
        for j in range(len(clist)):
            outcols.append(i+"_con"+str(j+1))
        for j in range(len(tlist)):
            outcols.append(i+"_tre"+str(j+1))

    listid = 1
    for i in range(len(clist[0])):
        if "_".join(clist[0][i][:4]) not in finlist:
            finlist.append("_".join(clist[0][i][:4]))
            
            sys.stdout.write("\r" + "Calculating " + "_".join(clist[0][i][:4]) + "            ")
            sys.stdout.flush()
            
            blist = {}
            
            
            blist["CHR"] = clist[0][i][0]
            blist["INTRST"] = int(clist[0][i][1])
            blist["INTREN"] = int(clist[0][i][2])
            blist["STRAND"] = clist[0][i][3]
            blist["GENEID"] = clist[0][i][4]
            blist["GENENAME"] = clist[0][i][5]
            blist["TRANID"] = clist[0][i][6]
            blist["INTRNUM"] = clist[0][i][7]
            
            
            
            istart = int(clist[0][i][1])
            iend = int(clist[0][i][2])
            istrand = clist[0][i][3]
            i20seq = fasdic[clist[0][i][0]][blist["INTRST"]-21:blist["INTREN"]+20]
            
            if istrand == "-":
                i20seq = Seq.reverse_complement(i20seq)
            
            blist["seq20i"] = i20seq
            blist["seq5p15"] = i20seq[14:29]
            blist["seq3p23"] = i20seq[-40:-17]
            blist["n5pe3"] = i20seq[17:18]
            blist["n5pe2"] = i20seq[18:19]
            blist["n5pe1"] = i20seq[19:20]
            blist["n5pi4"] = i20seq[23:24]
            
            for p,q in enumerate(["eijr","iejr","csr","a5r","a3r"]):
                for j in range(len(clist)):
                    blist[q+"con"+str(j+1)] = clist[j][i][p+8]
                for j in range(len(tlist)):
                    blist[q+"tre"+str(j+1)] = tlist[j][i][p+8]
            
            donseq = i20seq[20:22]
            accseq = i20seq[-22:-20]
            blist["s5p"] = donseq
            blist["s3p"] = accseq
            
            contot = [(float(j[i][8])+float(j[i][9]))/2+float(j[i][11])+float(j[i][12])+float(j[i][10]) for j in clist]
            tretot = [(float(j[i][8])+float(j[i][9]))/2+float(j[i][11])+float(j[i][12])+float(j[i][10]) for j in tlist]
            csrconavg = np.mean([float(j[i][10]) for j in clist])
            csrtreavg = np.mean([float(j[i][10]) for j in tlist])
            
            blist["minreads"] = min(contot+tretot)
            if min(contot+tretot) < 10 or (donseq != "GT" and donseq != "GC") or accseq != "AG" or max(csrconavg,csrtreavg) < 10:
                #retlist.append(blist)
                continue
                #if the intron is not GY-AG or has smaller read number than 10, skip
            
            #[0"CHR",1"INTRST",2"INTREN",3"STRAND",4"GENEID",5"GENENAME",6"TRANID",7"EXNUM",8"EIJ",9"IEJ",10"E1E2",11"ALTSP5",12"ALTSP3",13"ALT5LIS",14"ALT3LIS"] len=15
            
            conirs = [np.log2((float(j[i][8])+float(j[i][9])+0.1)/(2*float(j[i][10])+0.1)) for j in clist]
            treirs = [np.log2((float(j[i][8])+float(j[i][9])+0.1)/(2*float(j[i][10])+0.1)) for j in tlist]
            #cona5s = [np.log2((float(j[i][11])+0.1)/(float(j[i][10])+0.1)) for j in clist]
            #trea5s = [np.log2((float(j[i][11])+0.1)/(float(j[i][10])+0.1)) for j in tlist]
            #cona3s = [np.log2((float(j[i][12])+0.1)/(float(j[i][10])+0.1)) for j in clist]
            #trea3s = [np.log2((float(j[i][12])+0.1)/(float(j[i][10])+0.1)) for j in tlist]
            conpse = [100.0*((float(j[i][8])+float(j[i][9]))/2+float(j[i][11])+float(j[i][12]))/((float(j[i][8])+float(j[i][9]))/2+float(j[i][11])+float(j[i][12])+float(j[i][10])) for j in clist]
            trepse = [100.0*((float(j[i][8])+float(j[i][9]))/2+float(j[i][11])+float(j[i][12]))/((float(j[i][8])+float(j[i][9]))/2+float(j[i][11])+float(j[i][12])+float(j[i][10])) for j in tlist]
            
            #conlis = [conirs,cona5s,cona3s]
            #trelis = [treirs,trea5s,trea3s]
            conlis = [conirs]
            trelis = [treirs]
            
            #for pind,pnam in enumerate(["irs","a5s","a3s"]):
            for pind,pnam in enumerate(["irs"]):
                for j in range(len(clist)):
                    blist[pnam+"_con"+str(j+1)] = conlis[pind][j]
                for j in range(len(tlist)):
                    blist[pnam+"_tre"+str(j+1)] = trelis[pind][j]
            
            #cona5poss = dupcor((",".join([j[i][13] for j in clist])).split(","))
            #trea5poss = dupcor((",".join([j[i][13] for j in tlist])).split(","))
            #cona3poss = dupcor((",".join([j[i][14] for j in clist])).split(","))
            #trea3poss = dupcor((",".join([j[i][14] for j in tlist])).split(","))
            
            #blist["cona5pos"] = "_".join([":".join([str(s) for s in k]) for k in cona5poss])
            #blist["trea5pos"] = "_".join([":".join([str(s) for s in k]) for k in trea5poss])
            #blist["cona3pos"] = "_".join([":".join([str(s) for s in k]) for k in cona3poss])
            #blist["trea3pos"] = "_".join([":".join([str(s) for s in k]) for k in trea3poss])
            
            blist["conirsme"] = np.mean(conirs)
            blist["treirsme"] = np.mean(treirs)
            blist["conirssd"] = np.std(conirs, ddof=1)
            blist["treirssd"] = np.std(treirs, ddof=1)
            blist["irsdiff"]  = blist["treirsme"] - blist["conirsme"]
            blist["irszscore"] = (blist["treirsme"] - blist["conirsme"]) / np.sqrt((blist["conirssd"]**2)/len(clist) + (blist["treirssd"]**2)/len(tlist))
            
            #blist["cona5sme"] = np.mean(cona5s)
            #blist["trea5sme"] = np.mean(trea5s)
            #blist["cona5ssd"] = np.std(cona5s, ddof=1)
            #blist["trea5ssd"] = np.std(trea5s, ddof=1)
            #blist["a5sdiff"]  = blist["trea5sme"] - blist["cona5sme"]
            #blist["a5szscore"] = (blist["trea5sme"] - blist["cona5sme"]) / np.sqrt((blist["cona5ssd"]**2)/len(clist) + (blist["trea5ssd"]**2)/len(tlist))
            
            #blist["cona3sme"] = np.mean(cona3s)
            #blist["trea3sme"] = np.mean(trea3s)
            #blist["cona3ssd"] = np.std(cona3s, ddof=1)
            #blist["trea3ssd"] = np.std(trea3s, ddof=1)
            #blist["a3sdiff"]  = blist["trea3sme"] - blist["cona3sme"]
            #blist["a3szscore"] = (blist["trea3sme"] - blist["cona3sme"]) / np.sqrt((blist["cona3ssd"]**2)/len(clist) + (blist["trea3ssd"]**2)/len(tlist))
            
            blist["conpcsme"] = np.mean(conpse)
            blist["trepcsme"] = np.mean(trepse)
            
            #if blist["cona5sme"] < blist["trea5sme"] and 10 < min([int(j[i][11]) for j in tlist]):
            #    chancoord = int(trea5poss[0][0])
            #    if istrand == "+":
            #        pseq=fasdic[clist[0][i][0]][chancoord-6:chancoord+9]
            #    else:
            #        pseq=str(Seq.reverse_complement(fasdic[clist[0][i][0]][chancoord-9:chancoord+6]))
            #    blist["seqalt5"] = pseq
            #elif blist["cona5sme"] > blist["trea5sme"] and 10 < min([int(j[i][11]) for j in clist]):
            #    chancoord = int(cona5poss[0][0])
            #    if istrand == "+":
            #        pseq=fasdic[clist[0][i][0]][chancoord-6:chancoord+9]
            #    else:
            #        pseq=str(Seq.reverse_complement(fasdic[clist[0][i][0]][chancoord-9:chancoord+6]))
            #    blist["seqalt5"] = pseq
            #else:
            #    blist["seqalt5"] = ""
            
            #if blist["cona3sme"] < blist["trea3sme"] and 10 < min([int(j[i][12]) for j in tlist]):
            #    chancoord = int(trea3poss[0][0])
            #    if istrand == "+":
            #        pseq=fasdic[clist[0][i][0]][chancoord-6:chancoord+9]
            #    else:
            #        pseq=str(Seq.reverse_complement(fasdic[clist[0][i][0]][chancoord-9:chancoord+6]))
            #    blist["seqalt3"] = pseq
            #    
            #elif blist["cona3sme"] > blist["trea3sme"] and 10 < min([int(j[i][12]) for j in clist]):
            #    chancoord = int(cona3poss[0][0])
            #    if istrand == "+":
            #        pseq=fasdic[clist[0][i][0]][chancoord-6:chancoord+9]
            #    else:
            #        pseq=str(Seq.reverse_complement(fasdic[clist[0][i][0]][chancoord-9:chancoord+6]))
            #    blist["seqalt3"] = pseq
                
            #else:
            #    blist["seqalt3"] = ""
            
            blist["ID"] = listid
            listid += 1
            
            retlist.append(blist)
    
    print("")
    
    outlist = [outcols]
    for n in retlist:
        outlist.append([n[m] if m in n else "0" for m in outcols])
    
    return (outlist,retlist)

def savetxt(nda:list, fname:str):
    with open(fname, "w") as wf:
        #for i, n in enumerate(nda):
        for n in nda:
            wf.write("\t".join([str(i) for i in n])+"\n")


(clist, tlist) = readfiles(ctrlist, trelist)
fasdic = getfastadic(fastafi)
(calccounts, caldic) = getcalccounts(clist, tlist, fasdic)
savetxt(calccounts, outf)







