# -*- coding: utf8 -*-
######################################################
## Author:  
## Date:    2018-06-27
## Place:   PUMC&CAMS-IMPLAD
## Version: 2.02
######################################################

from gene import *
import os, sys

def check_seq(seq, l_or_m):
    flag = False
    if l_or_m == 'l':
        for x in seq.strip().upper():
            if x not in ('A','G','T','C'):
                print("ERROR: One of your Query (lncRNA) sequences has a degenerate or non-DNA character, please check\n")
                return False
    if l_or_m== 'm':
        for x in seq.strip().upper():
            if x not in ('A','G','T','C'):
                print("ERROR: One of your Subject (mRNA) sequences has a degenerate or non-DNA character, please check\n")
                return False
    return True

def RNA_Target(Dir, RNAfile, mRNAfile, deltaG, align, outFile, Tempera):
    
    seq_q = open(Dir + "/seq_temp_qq.txt", 'w')
    fp = open(RNAfile, 'r')
    first_line = ""
    seq_line = ""
    q_num = 0
    i = 0
    fp_lines = fp.readlines()
    for line in fp_lines:
        if line[0:1]=='>':
            if i == 0:
                i = i+1
            else:
                seq_q.write(first_line+"\t"+seq_line+'\n')
            first_line = line.strip()
            seq_line = ""
            q_num = q_num + 1
        else:
            seq_line = seq_line + line.strip()
    seq_q.write(first_line+"\t"+seq_line)
    fp.close()
    seq_q.close()

    seq_s = open(Dir + "/seq_temp_ss.txt", 'w')
    fp = open(mRNAfile, 'r')
    first_line = ""
    seq_line = ""
    s_num = 0
    i = 0
    fp_lines = fp.readlines()
    for line in fp_lines:
        if line[0:1]=='>':
            if i == 0:
                i = i+1
            else:
                seq_s.write(first_line+"\t"+seq_line+'\n')
            first_line = line.strip()
            seq_line = ""
            s_num = s_num + 1
        else:
            seq_line = seq_line + line.strip()
    seq_s.write(first_line+"\t"+seq_line)
    fp.close()
    seq_s.close()
    
    lncRNAFile="seq_temp_qq.txt"
    mRNAFile="seq_temp_ss.txt"
    deltaG=float(deltaG)
    Tempa=float(Tempera)
    isPic=align

    lncRNAs = {}
    mRNAs = {}

    outFile_fp = open(Dir + "/" + outFile, 'w')

    sTitle='lncRNA\tlncRNA_Length\tmRNA\tmRNA_Length\tdG\tndG\tStart_Position_lncRNA\tEnd_Position_lncRNA\tStart_Position_mRNA\tEnd_Position_mRNA\n'
    outFile_fp.write(sTitle)
    outFile_fp.close()
    
    lncRNAFile = open(Dir+"/"+lncRNAFile, 'r')
    for lnc_line in lncRNAFile.readlines():
        
        lncTemp = lnc_line.split('\t')
        lncRNA_Name = lncTemp[0]
        lncRNA_Seq = lncTemp[1]

        if lncRNA_Name == '' or lncRNA_Seq == '':
            return 'The format of lncRNA file is error, Please check.\n '
              
        lncRNA_Seq = lncRNA_Seq.upper()
        lncRNA_Seq.replace(' ', '')
        lncRNA_Seq.replace('\n', '')
        lncRNA_Seq.replace('U', 'T')

        lncRNAs[lncRNA_Name] = lncRNA_Seq
    lncRNAFile.close()
        
    mRNAFile = open(Dir+"/"+mRNAFile, 'r')
    for m_line in mRNAFile.readlines():
        mTemp = m_line.split('\t')
        mRNA_Name = mTemp[0]
        mRNA_Seq = mTemp[1]

        if mRNA_Name == '' or mRNA_Seq == '':
            return 'The format of mRNA file is error, Please check.\n '
        
        mRNA = mRNA_Seq.upper()
        mRNA_Seq.replace(' ', '')
        mRNA_Seq.replace('\n', '')
        mRNA_Seq.replace('U', 'T')

        mRNAs[mRNA_Name] = mRNA_Seq
    mRNAFile.close()

    sum_num = q_num * s_num
    now_num = 0
    percent = 1
    for lnc_name in lncRNAs:
        for m_name in mRNAs:
            if(check_seq(lncRNAs[lnc_name],'l')) and (check_seq(mRNAs[m_name],'m')):
                y = Calculate_G(Dir,lnc_name, lncRNAs[lnc_name], m_name, mRNAs[m_name], deltaG, isPic, outFile, Tempa)
                now_num = now_num + 1
                percent = round(float(float(now_num)/float(sum_num)), 2)*100.00
                print("%s" % (percent))
            else:
                break
    line_num = 0
    temp = []
    result={}
    res_sort = {}

    resilt_table = open(Dir+"/result_table.txt", 'w') 
    result_fp = open(Dir+"/"+outFile, 'r')
    result_line = result_fp.readline()
    result_line = result_fp.readline()
    while result_line:
        line_num = line_num + 1
        if line_num%5 == 1:
            vs_title = result_line.strip('\n')
            res_sort[vs_title.split('\t')[0]+vs_title.split('\t')[2]] = vs_title.split('\t')[5]
##            print(str(vs_title.split('\t')[5])+'\n')
##            print(str(line_num))
            temp.append(vs_title)
        if line_num%5 == 2:
            vs_l = result_line.rstrip()
            temp.append(vs_l)
        if line_num%5 == 3:
            vs_align = result_line.rstrip()
            temp.append(vs_align)
        if line_num%5 == 4:
            vs_m = result_line.rstrip()
            temp.append(vs_m)
        if line_num%5 == 0:
            result[vs_title.split('\t')[0]+vs_title.split('\t')[2]]=temp
##            print(str(line_num))
##            print(str(result[line_num]))
            temp = []   
        result_line = result_fp.readline()
    if line_num == 0:
        return "There is no result record"
    else:
        resilt_table.write('<table>')
        resilt_table.write('<tr><th>Query</th><th>Q_Len</th><th>Subjct</th><th>S_Len</th><th>&Delta;G</th><th>&Delta;G/n</th></tr>')
        temp_sort = sorted(res_sort.items(), key=lambda d: d[1],reverse = True)
        for val in range(len(temp_sort)):
            if float(result[temp_sort[val][0]][0].split('\t')[5]) > deltaG:
                resilt_table.write("<tr class='red' id='"+str(val)+"' onclick='ex_cla(this.id)'>")
            else:
                resilt_table.write("<tr class='green' id='"+str(val)+"' onclick='ex_cla(this.id)'>")
            resilt_table.write("<td>"+result[temp_sort[val][0]][0].split('\t')[0].replace('>','')+"</td>")
            resilt_table.write("<td>"+result[temp_sort[val][0]][0].split('\t')[1]+"</td>")
            resilt_table.write("<td>"+result[temp_sort[val][0]][0].split('\t')[2].replace('>','')+"</td>")
            resilt_table.write("<td>"+result[temp_sort[val][0]][0].split('\t')[3]+"</td>")
            resilt_table.write("<td>"+result[temp_sort[val][0]][0].split('\t')[4]+"</td>")
            resilt_table.write("<td>"+result[temp_sort[val][0]][0].split('\t')[5]+"</td>")
            resilt_table.write("</tr>")
            index_c = result[temp_sort[val][0]][2].count(' ')
            sum_len = len(result[temp_sort[val][0]][2].strip())
            start = index_c
            l = 0
            if result[temp_sort[val][0]][1][0:1]==' ' or result[temp_sort[val][0]][3][0:1]=='3':
                index_p = 1
                index_r = int(index_c)-2
                while l < sum_len:
                    resilt_table.write("<tr class='align' name='"+str(val)+"' style='display:none;'>")
                    resilt_table.write("<td>Query&nbsp;&nbsp;&nbsp;&nbsp;5'  <br><br>Subjct&nbsp;&nbsp;&nbsp;&nbsp;3'  </td>")
                    resilt_table.write("<td><p class='align_p'>"+str(index_p)+"<br><br>"+str(index_r)+"</p></td>")
                    resilt_table.write("<td colspan = 3><p class='align_p'>")
                    resilt_table.write(result[temp_sort[val][0]][1][start+l:start+l+100]+"<br>")
                    resilt_table.write(result[temp_sort[val][0]][2][start+l:start+l+100]+"<br>")
                    resilt_table.write(result[temp_sort[val][0]][3][start+l:start+l+100])
                    resilt_table.write("</p></td>")
                    resilt_table.write("<td>3'<br><br>5'</td>")
                    resilt_table.write("</tr>")
                    l = l + 100
                    index_p = index_p + 100
                    index_r = index_r + 100
            if result[temp_sort[val][0]][3][0:1]==' ' or result[temp_sort[val][0]][1][0:1]=='5':
                index_p = int(index_c)-2
                index_r = 1
                while l < sum_len:
                    resilt_table.write("<tr class='align' name='"+str(val)+"' style='display:none;'>")
                    resilt_table.write("<td>Query&nbsp;&nbsp;&nbsp;&nbsp;5'  <br><br>Subjct&nbsp;&nbsp;&nbsp;&nbsp;3'  </td>")
                    resilt_table.write("<td><p class='align_p'>"+str(index_p)+"<br><br>"+str(index_r)+"</p></td>")
                    resilt_table.write("<td colspan = 3><p class='align_p'>")
                    resilt_table.write(result[temp_sort[val][0]][1][start+l:start+l+100]+"<br>")
                    resilt_table.write(result[temp_sort[val][0]][2][start+l:start+l+100]+"<br>")
                    resilt_table.write(result[temp_sort[val][0]][3][start+l:start+l+100])
                    resilt_table.write("</p></td>")
                    resilt_table.write("<td>3'<br><br>5'</td>")
                    resilt_table.write("</tr>")
                    l = l + 100
                    index_p = index_p + 100
                    index_r = index_r + 100
        resilt_table.write("</table>")
    resilt_table.close()

    table_str = "\n".join(open(Dir+"/result_table.txt", 'r').readlines()); 

    resultHTML_fp = open("result.html","w")
    template_fp = open("js_css/template.html", 'r')
    for line_i in template_fp.readlines():
        if line_i.strip() == "<!--RESULT TABLE-->":
            resultHTML_fp.write(table_str)
        else:
            resultHTML_fp.write(line_i)
    template_fp.close()
    resultHTML_fp.close()




def main():
    ## Parametars
    Dir = "output"
    align = "T"
    outFile = "result.txt"
    Tempera = 37.0
    deltaG = -0.05
    argv_ls = ['-s', '-t', '-q', '-g']

    if not os.path.exists("output"):
        os.makedirs("output")
    
    print sys.argv
    i = 1
    for argv_i in sys.argv[1::2]:
        if argv_i not in argv_ls:
            print "ERROR: '" + argv_i + "' is invalid !"
        else:
            if argv_i == '-q':
                RNAfile = sys.argv[i+1]
            if argv_i == '-s':
                mRNAfile = sys.argv[i+1]
            if argv_i == '-t':
                Tempera = sys.argv[i+1]
            if argv_i == '-g':
                deltaG = sys.argv[i+1]
        i += 2

    RNA_Target(Dir, RNAfile, mRNAfile, deltaG, align, outFile, Tempera)

if __name__ == '__main__':
    main()
