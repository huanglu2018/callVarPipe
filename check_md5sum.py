def check_md5sum(fq1, fq2, md5_re, md5_log, md5_golden):

    def readTwoColTxtAsDic(TwoColTxt):
        with open(TwoColTxt,"r",encoding="utf-8") as f:
            my_dict={}
            for line in f.readlines():
                line=line.strip('\n')
                b=line.split('  ',1)
                fq=b[1].rsplit('/',1)[-1]
                my_dict[fq]=b[0]
        return my_dict

    import os,subprocess
    fqId=os.path.basename(fq1).split("_1.fq.gz")[0]
    fqFile1=os.path.basename(fq1)
    fqFile2=os.path.basename(fq2)
    if not os.path.exists(md5_golden):
        with open(md5_log, 'w') as f1:
            f1.write("ERROR1: "+ fqId + " got no golden md5 file, check with manufacturer !!!")
        if not os.path.exists(md5_re):
            cmd="cd $(dirname " + fq1 + ") ; md5sum " + fqFile1 + " " + fqFile2 + " > " + md5_re
            p=subprocess.Popen(cmd,shell=True)
            return_code=p.wait()
    else:
        if not os.path.exists(md5_re):
            cmd="cd $(dirname " + fq1 + ") ; md5sum " + fqFile1 + " " + fqFile2 + " > " + md5_re
            p=subprocess.Popen(cmd,shell=True)
            return_code=p.wait()
        else:
            nrow=len(open(md5_re,"r").readlines())
            if nrow != 2:
                print("incomplete md5_re1, re calculating...  "+ md5_re)
                cmd="cd $(dirname " + fq1 + ") ; md5sum " + fqFile1 + " " + fqFile2 + " > " + md5_re
                p=subprocess.Popen(cmd,shell=True)
                return_code=p.wait()
        dicGold=readTwoColTxtAsDic(md5_golden)
        dicRe1=readTwoColTxtAsDic(md5_re)
        with open(md5_log, 'w') as f2:
            if (dicGold[fqFile1] == dicRe1[fqFile1]) and (dicGold[fqFile2] == dicRe1[fqFile2]):
                f2.write("PASS: "+ fqId + " md5 validation succeed.")
            else:
                f2.write("ERROR2: "+ fqId + " md5 validation failed.")
check_md5sum(snakemake.input[0], snakemake.input[1], snakemake.output[0], snakemake.output[1], snakemake.params[0])
