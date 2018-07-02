# lncRNATarget_local
##lncRNATargets: A platform for lncRNA target prediction based on nucleic acid thermodynamics.

This tool is for the publication of "lncRNATargets: A platform for lncRNA target prediction based on nucleic acid thermodynamics."

## Due to the instable of the server, I bulid this local-run version.

### Usage:

```bash
python rnatarget.py -t  37  -s [subject_file.fa/.txt/.fasta] -q [query_file.fa/.txt/.fasta] -g -0.05
```

-t : the Temperature value in degree Celsius (Default: 37);<br>
-s : the Subject sequences file;<br>
-q : the Query sequences file;<br>
-g : the threshold value for deltaG/n ; <br>
