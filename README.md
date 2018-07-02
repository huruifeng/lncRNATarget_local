# lncRNATarget_local
## lncRNATargets: A platform for lncRNA target prediction based on nucleic acid thermodynamics.

This tool is for the publication of "lncRNATargets: A platform for lncRNA target prediction based on nucleic acid thermodynamics."
## Due to the instable of the server, I bulid this local-run version.

### Usage:

```bash
python rnatarget.py -t  37  -s [subject_file.fa/.txt/.fasta] -q [query_file.fa/.txt/.fasta] -g -0.05
```

-s : the Subject sequences file (Required);<br>
-q : the Query sequences file (Required);<br>
-t : the Temperature value in degree Celsius (Default: 37);<br>
-g : the threshold value for deltaG/n (Default: -0.05); <br>


If this tool can have any help for you, please cite

`Hu, R. and Sun, X., 2016. lncRNATargets: a platform for lncRNA target prediction based on nucleic acid thermodynamics. Journal of bioinformatics and computational biology, 14(04), p.1650016.`
