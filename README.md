## ProxLogPRF: A Proximity-based Log-logistic Feedback Model for Pseudo-relevance Feedback

This is the official repository of the manuscript "ProxLogPRF: A Proximity-based Log-logistic Feedback Model for Pseudo-relevance Feedback" 
submitted to Information Processing & Management (IP&M).

### Updates
- Jan 23, 2022: the project has been released to GitHub.

### Requirements
- java 1.7 -- development environment
- some necessary *.jar packages -- we include them in the 'lib' folder.
- pandas -- used by trecEval.py

### Experimentation instructions
- Step 1: create the index for each data collection (e.g. WT2G).
    ```
    $ javac -cp "lib/*:" ./*.java                                             # compile
    $ java -cp "lib/*:" index.IndexTREC -docs datasets/WT2G/ -data WT2G       # create index for WT2G
    ```

- Step 2: retrieve the documents using each model (e.g. BM25). ps: Step 2 generates the ranking results and outputs 
them to a file named *-report.txt under 'ProxLogPRF/result'. In addition, the *.txt files containing the metric results 
(i.e. MAP and P@k) will also be generated under 'ProxLogPRF/result'.
    ```
    $ javac -cp "lib/*:" ./*.java                                             # compile
    $ java -cp "lib/*:" ProxPRF.BM25 -k1 1.2 -b 0.35                          # use BM25 model to retrieve documents
    $ java -cp "lib/*:" ProxPRF.BM25 -h                                       # use this command to check arguments usages
    ```

- Step 3: evaluate the retrieval model
We use trecEval.py to evaluate the model performance via MAP, P@k, nDCG and nDCG@k.
    ```
    $ python trecEval.py result/BM25/WT2G-BM25-1.2-0.35-report.txt query-judge/qrels.WT2G result/WT2G-BM25-1.2-0.35.xls
    ```

### Package structure
- analyzer
    - MyStopAndStemmingAnalyzer.java: stopwords removal and stemming
- common
    - ByWeightComparator.java -- numerical comparator.
    - MyQQParser.java -- simplistic quality query parser.
    - MyTrecParser.java -- TREC document analyzer
    - QualityStats.java -- compute the results (MAP, P@k and MRR) of quality benchmark run for a single query or for a set of queries.
    - StaTools.java -- implementation on some basic statistical functions
- datasets -- directory to the data collections
- index
    - IndexTREC.java -- create index for data collections
- indices -- directory to the files containing index of each data collection
- lib -- directory to all the *.jar packages used for the project
- models -- directory to all the retrieval models (i.e. BM25, DLM, LL, LLPRF, LLEXPStar (LL+EXP*), PRoc2, PRoc3 and ProxLogPRF)
- query-judge -- directory to all the query topics
- result -- directory to the experimental results
- stopwords.txt -- stopwords used in our experiments
- trecEval.py -- evaluation metrics


### Data collections
We tested baselines, SOTA proximity-based PRF models and our model variants on eight standard TREC collections, 
namely AP (Associated Press 1988-90), DISK1&2, DISK4&5, ROBUST04 (TREC Robust Track 2004), WSJ (Wall Street Journal), 
WT2G (TREC Web Track 2000), WT10G (TREC Web Track 2001- 2002) and GOV2. Note that AP, DISK1&2, DISK4&5, ROBUST04 and 
WSJ are popular newswire collections where noise is rare, while WT2G, WT10G and GOV2 are collections consisting of web 
documents with inherent noises.
- [AP](http://opus.nlpl.eu/OpenSubtitles-v2018.php).
- [DISK1\&2](http://opus.nlpl.eu/OpenSubtitles-v2018.php).
- [DISK4\&5](http://yanran.li/dailydialog.html).
- [ROBUST04](http://opus.nlpl.eu/OpenSubtitles-v2018.php).
- [WSJ](http://opus.nlpl.eu/OpenSubtitles-v2018.php).
- [WT2G](http://opus.nlpl.eu/OpenSubtitles-v2018.php).
- [WT10G](http://opus.nlpl.eu/OpenSubtitles-v2018.php).
- [GOV](http://opus.nlpl.eu/OpenSubtitles-v2018.php).


### Acknowledgments
This research is supported by the [Natural Sciences and Engineering Research Council (NSERC) of Canada](https://www.nserc-crsng.gc.ca/index_eng.asp), 
the [York Research Chairs (YRC) program](https://www.yorku.ca/research/york-research-chairs/),
[NSERC CREATE award](https://www.nserc-crsng.gc.ca/Professors-Professeurs/Grants-Subs/CREATEResults-ResultatsFONCER_eng.asp?Year=2015) 
and an [ORF-RE (Ontario Research Fund Research Excellence award](https://www.ontario.ca/page/ontario-research-fund-research-excellence) 
in [BRAIN Alliance](https://brainalliance.ca/en).
