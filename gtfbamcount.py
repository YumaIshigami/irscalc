import sys, pysam, datetime

if len(sys.argv) != 4:
    print("Usage: python3 gtfbamcount.py [BAM file name] [GTF file name] [output bamcount file name]")
    exit(1)
    
inbam = pysam.AlignmentFile(sys.argv[1],"rb")
ingtf = sys.argv[2]
outf = sys.argv[3]

def countintronic(chromosome, strand, exon1end, exon2start):
    eij = 0
    iej = 0
    e1e2 = 0
    altsp5 = 0
    altsp3 = 0
    alt5lis = []
    alt3lis = []
    rnames = []
    for infe in inbam.fetch(chromosome, exon1end-1, exon2start):
        if infe.is_proper_pair and not infe.is_unmapped and infe.query_name not in rnames:
            rnames.append(infe.query_name)
            if (strand == "+" and ((infe.is_read1 and infe.is_reverse) or (infe.is_read2 and not(infe.is_reverse)))) or \
                    (strand == "-" and ((infe.is_read1 and not(infe.is_reverse)) or (infe.is_read2 and infe.is_reverse))):
                if infe.get_overlap(exon1end-1,exon1end+1) == 2:
                    eij += 1
                if infe.get_overlap(exon2start-2,exon2start) == 2:
                    iej += 1
                curloc = infe.reference_start #current location
                for incig in infe.cigar:
                    if incig[0] == 3: # if there is an intron gap
                        if curloc == exon1end and curloc+incig[1] == exon2start-1:
                            e1e2 += 1
                        elif curloc == exon1end and curloc+incig[1] != exon2start-1:
                            if strand == "+":
                                altsp3 += 1
                                alt3lis.append(str(curloc+incig[1]))
                            elif strand == "-":
                                altsp5 += 1
                                alt5lis.append(str(curloc+incig[1]))
                        elif curloc != exon1end and curloc+incig[1] == exon2start-1:
                            if strand == "+":
                                altsp5 += 1
                                alt5lis.append(str(curloc))
                            elif strand == "-":
                                altsp3 += 1
                                alt3lis.append(str(curloc))
                    curloc += incig[1]
                    
    if (e1e2+(eij+iej)/2.0+altsp5+altsp3) < 10:
        return ["-1","-1","-1","-1","-1","",""]
    return [str(eij),str(iej),str(e1e2),str(altsp5),str(altsp3),",".join(alt5lis),",".join(alt3lis)]
    

if __name__ == "__main__":
    d1 = datetime.datetime.now() #start counting time
    genename = ""
    exonnum = 0 #exon number
    intlist = [["CHR","INTRST","INTREN","STRAND","GENEID","GENENAME","TRANID","EXNUM","EIJ","IEJ","E1E2","ALTSP5","ALTSP3","ALT5LIS","ALT3LIS"]] #column names
    coordnow = -1 #current 3' end of the exon
    gtflen = sum([1 for _ in open(ingtf)])
    print("Lines read!")
    
    with open(ingtf) as inputgtffile:
        for i,gline in enumerate(inputgtffile):
            gl = gline[:-1].split("\t")
            if len(gl) == 9:
                attr = gl[8].split("\"")
                if gl[2] == "transcript":
                    geneid = attr[1]
                    if "; gene_name " in attr:
                        genename = attr[attr.index("; gene_name ")+1]
                    else:
                        genename = geneid
                    tranid = attr[attr.index("; transcript_id ")+1]
                    coordnow = -1
                    sys.stdout.write("\r" + sys.argv[1] + ": " + str(round(i/gtflen*100,2)) + "%" + "; checking:chr" + gl[0] + ": " + gl[3] + ": " + genename + "               ")
                    sys.stdout.flush()
                if gl[2] == "exon":
                    if gl[6] == "+":
                        if coordnow == -1:
                            coordnow = int(gl[4])
                        else:
                            tmpl = [gl[0],str(coordnow+1),str(int(gl[3])-1),gl[6],geneid,genename,tranid,str(int(attr[5])-1)]
                            tmpl.extend(countintronic(gl[0],gl[6],coordnow,int(gl[3])))
                            intlist.append(tmpl)
                            coordnow = int(gl[4])
                    elif gl[6] == "-":
                        if coordnow == -1:
                            coordnow = int(gl[3])
                        else:
                            tmpl = [gl[0],str(int(gl[4])+1),str(coordnow-1),gl[6],geneid,genename,tranid,str(int(attr[5])-1)]
                            #print tmpl
                            tmpl.extend(countintronic(gl[0],gl[6],int(gl[4]),coordnow))
                            intlist.append(tmpl)
                            coordnow = int(gl[3])
                            
                    else:
                        print("error in exon of gtf file")
                        input()
    
    with open(outf,"w") as wf:
        for tl in intlist:
            try:
                wf.write("\t".join(tl)+"\n")
            except TypeError:
                print(tl)
                input()
    
    print("")
    d2=datetime.datetime.now()
    dt=d2-d1
    print(sys.argv[1] + " count complete with " + str(dt.seconds/60.0) + " minutes!")
