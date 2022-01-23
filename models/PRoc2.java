package ProxPRF;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Map.Entry;
import org.apache.lucene.benchmark.quality.Judge;
import org.apache.lucene.benchmark.quality.QualityQuery;
import org.apache.lucene.benchmark.quality.QualityQueryParser;
import org.apache.lucene.benchmark.quality.trec.TrecJudge;
import org.apache.lucene.benchmark.quality.trec.TrecTopicsReader;
import org.apache.lucene.benchmark.quality.utils.DocNameExtractor;
import org.apache.lucene.benchmark.quality.utils.SubmissionReport;
import org.apache.lucene.index.AtomicReader;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.DocsAndPositionsEnum;
import org.apache.lucene.index.DocsEnum;
import org.apache.lucene.index.Fields;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.MultiFields;
import org.apache.lucene.index.SlowCompositeReaderWrapper;
import org.apache.lucene.index.Term;
import org.apache.lucene.index.Terms;
import org.apache.lucene.index.TermsEnum;
import org.apache.lucene.search.DocIdSetIterator;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.TopDocs;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.FSDirectory;
import org.apache.lucene.util.BytesRef;

import common.ByWeightComparator;
import common.MyQQParser;
import common.QualityStats;
import common.StaTools;

public class PRoc2 {

	QualityQuery[] qualityQueries;
	AtomicReader atomicReader;

	QualityQueryParser qqParser;	

	IndexReader reader;

	IndexSearcher searcher;

	String docNameField = "docno";

	String docContentField = "contents";

	PrintWriter qualityLog;

	SubmissionReport submitRep;

	Judge judge;

	int maxResults = 1000;

	private static int totalDocNum;
	private static double avg_doc_length = 0.0;

	private static HashMap<String, Long> term_total_freq_map = new HashMap<String, Long>();
	private static HashMap<String, Integer> term_doc_freq_map = new HashMap<String, Integer>();
	private static HashMap<Integer, Integer> doc_length_map = new HashMap<Integer, Integer>();
	private static HashMap<Integer, Double> doc_avg_tf_map = new HashMap<Integer, Double>();

	private int PRFDocNum = 10;		
	private int PRFTermNum = 50;	
	private double alpha = 1.0;		
	private double beta = 1.0;

	private double k1 = 1.2;
	private double b = 0.35;
	private double k3 = 8.0;
	private double sigma = 50;

	public PRoc2() {
	}

	public PRoc2(double k1, double b) {
		this.k1 = k1;
		this.b = b;
	}
	
	public PRoc2(double k1, double b, double sigma, int PRFDocNum, int PRFTermNum, double alpha, double beta) {
		this.k1 = k1;
		this.b = b;
		this.sigma = sigma;
		this.PRFDocNum = PRFDocNum;
		this.PRFTermNum = PRFTermNum;
		this.alpha = alpha;
		this.beta = beta;
	}

	/**
	 * PRoc2 model in Proximity-based Rocchio’s Model for Pseudo Relevance Feedback 
	 * @param args
	 * @throws Exception
	 * @author JFH
	 */
	public static void main(String[] args) throws Exception {

		String usage =
				"Usage: [-index indexFile] [-topics topicsFile] [-qrels qrelsFile] [-data dataSetName] [-k1 k1] [-b b] "
				+ "[-s simgma]  [-fd PRFDocNum] [-ft PRFTermNum] [-alpha alpha] [-beta beta]";
		if (args.length > 0 && ("-h".equals(args[0]) || "-help".equals(args[0]))) {
			System.out.println(usage);
			System.exit(0);
		}
		
		String indexFile = "indices/index_WT2G";
		String topicsFile = "query-judge/topics.WT2G";
		String qrelsFile = "query-judge/qrels.WT2G";
		String dataSetName = "WT2G";
		
		double k1 = 1.2;
		double b = 0.35;
		double sigma = 50;
		int PRFDocNum = 10;		
		int PRFTermNum = 50;	
		double alpha = 1.0;		
		double beta = 1.0;		
		
		for(int i = 0;i < args.length;i++) {
			if ("-index".equals(args[i])) {
				indexFile = args[i+1];
				i++;
			} else if ("-topics".equals(args[i])) {
				topicsFile = args[i+1];
				i++;
			} else if ("-qrels".equals(args[i])) {
				qrelsFile = args[i+1];
				i++;
			} else if ("-data".equals(args[i])) {
				dataSetName = args[i+1];
				i++;
			} else if ("-k1".equals(args[i])) {
				k1 =  Float.parseFloat(args[i+1]);
				i++;
			} else if ("-b".equals(args[i])) {
				b =  Float.parseFloat(args[i+1]);
				i++;
			} else if ("-s".equals(args[i])) {
				sigma =  Float.parseFloat(args[i+1]);
				i++;
			} else if ("-fd".equals(args[i])) {
				PRFDocNum =  Integer.parseInt(args[i+1]);
				i++;
			}  else if ("-ft".equals(args[i])) {
				PRFTermNum =  Integer.parseInt(args[i+1]);
				i++;
			} else if ("-alpha".equals(args[i])) {
				alpha =  Float.parseFloat(args[i+1]);
				i++;
			} else if ("-beta".equals(args[i])) {
				beta =  Float.parseFloat(args[i+1]);
				i++;
			}   
		}

		final PRoc2 myproc2 = new PRoc2(k1, b, sigma, PRFDocNum, PRFTermNum, alpha, beta);
		final String resultFile =  String.format("result/PRoc2/%s-PRoc2-%.1f-%.2f-%.0f-%d-%d-%.1f-%.1f.txt", dataSetName, myproc2.k1, myproc2.b, myproc2.sigma, myproc2.PRFDocNum, myproc2.PRFTermNum, myproc2.alpha, myproc2.beta);
		final String reportFile = String.format("result/PRoc2/%s-PRoc2-%.1f-%.2f-%.0f-%d-%d-%.1f-%.1f-report.txt", dataSetName, myproc2.k1, myproc2.b, myproc2.sigma, myproc2.PRFDocNum, myproc2.PRFTermNum, myproc2.alpha, myproc2.beta); 
		myproc2.Run(topicsFile, qrelsFile, indexFile, resultFile, reportFile); 

	}

	public void Run(String topicsFile, String qrelsFile, String indexFile, String resultFile, String reportFile) throws Exception {
		File result = new File(resultFile);

		if (result.exists()) {
			System.err.println(String.format("[%s] exists, pass!!!", result.getName()));
			return;
		} else {
			result.getParentFile().mkdirs();
			System.out.println(String.format("--------------------%s--------------------", result.getName().replace(".txt", "")));
		}

		Directory directory = FSDirectory.open(new File(indexFile));
		reader = DirectoryReader.open(directory);
		searcher = new IndexSearcher(reader);
		atomicReader = SlowCompositeReaderWrapper.wrap(reader);

		totalDocNum = reader.numDocs();

		TrecTopicsReader qReader = new TrecTopicsReader();
		qualityQueries = qReader.readQueries(new BufferedReader(new FileReader(new File(topicsFile))));

		qqParser = new MyQQParser("title", "contents");

		judge = new TrecJudge(new BufferedReader(new FileReader(new File(qrelsFile))));
		judge.validateData(qualityQueries, qualityLog); 

		qualityLog = new PrintWriter(new File(resultFile), "UTF-8");

		submitRep = new SubmissionReport(new PrintWriter(new File(reportFile), "UTF-8"), "TEST");

		execute( );

		directory.close();
		qualityLog.close();
	}


	public void termStats() throws Exception {

		Fields fields = MultiFields.getFields(reader);
		Terms terms = fields.terms(docContentField);

		TermsEnum iterator = terms.iterator(null);
		BytesRef byteRef = null;

		while ((byteRef = iterator.next()) != null) {
			String term = new String(byteRef.bytes, byteRef.offset, byteRef.length);
			term_total_freq_map.put(term, iterator.totalTermFreq());	
			term_doc_freq_map.put(term, iterator.docFreq());	
		}

	}

	public void docStats() throws Exception {	

		long totalDocLength = 0;

		for (int j = 0; j < reader.numDocs(); j++) {

			int docLen = 0;		
			int term_num = 0;	

			Terms terms = reader.getTermVector(j, docContentField);		

			if (terms != null && terms.size() > 0) {

				TermsEnum termsEnum = terms.iterator(null);		


				while ((termsEnum.next()) != null) {
					int freq = (int) termsEnum.totalTermFreq();
					docLen += freq;
					term_num++;
				}
			}

			totalDocLength += docLen;
			doc_length_map.put(j, docLen);

			double avg_tf = (term_num == 0) ? 0 : ((double) docLen) / term_num;
			doc_avg_tf_map.put(j, avg_tf);
		}
		avg_doc_length = ((double) totalDocLength) / totalDocNum;
	}

	public static double log2(double n) {
		return (Math.log(n) / Math.log(2));
	}

	public void execute( ) throws Exception {
		termStats();
		docStats();
		search();
	}

	public void search() throws Exception {

		QualityStats stats[] = new QualityStats[qualityQueries.length];

		for (int i = 0; i < qualityQueries.length; i++) {

			QualityQuery qq = qualityQueries[i];
			Query query = qqParser.parse(qq);
			HashSet<Term> term_set = new HashSet<Term>();
			query.extractTerms(term_set);	

			ArrayList<Term> term_list = new ArrayList<Term>();

			for (Term term : term_set) {
				DocsEnum docsEnum = atomicReader.termDocsEnum(term);	
				if(docsEnum!=null){
					term_list.add(term);
				}
			}
			HashMap<String, Double> queryVector = generateQueryVectorFromQuery(query);

			ScoreDoc[] firstRetrievalScoreDocs = retrievalByBM25(queryVector);

			ScoreDoc[] PRFScoreDocs = new ScoreDoc[Math.min(PRFDocNum, firstRetrievalScoreDocs.length)];

			for (int j = 0; j < PRFScoreDocs.length; j++) {
				PRFScoreDocs[j] = firstRetrievalScoreDocs[j];
			}

			Set<Integer> doc_set = get_doc_set(PRFScoreDocs);
			ArrayList<Term> term_list1 = get_term_list(PRFScoreDocs);

			HashMap<String, HashMap<Integer, ArrayList<Integer>>> within_query_freq_map = get_term_position(term_list1, doc_set); 

			int[][]freq_array = get_freq_array(term_list);
			
			HashMap<String, Double> topTermWeightVector = getTermProximityVectorFromPRFScoreDocs(term_list, freq_array, within_query_freq_map, PRFScoreDocs);

			HashMap<String, Double> newQueryVector = generateNewQueryVector(queryVector, topTermWeightVector);

			ScoreDoc[] secondRetrievalScoreDocs = retrievalByBM25(newQueryVector);

			CreateResult(secondRetrievalScoreDocs, qualityQueries, i, query, stats);
		}

		WriteResult(stats);
	}

	/**
	 * 
	 * @param query 
	 * @return queryVector
	 * @throws IOException 
	 */
	public HashMap<String, Double> generateQueryVectorFromQuery(Query query) throws IOException {
		HashMap<String, Double> queryVector = new HashMap<String, Double>();
		HashSet<Term> termSet = new HashSet<Term>();
		query.extractTerms(termSet);
		for (Term term : termSet) {
			DocsEnum docsEnum = atomicReader.termDocsEnum(term); 
			if (docsEnum != null) {
				queryVector.put(term.text(), 1.0);
			}
		}
		return queryVector;
	}

	/**
	 * @param queryVector 
	 * @return scoreDocs
	 * @throws IOException
	 */

	public ScoreDoc[] retrievalByBM25(HashMap<String, Double> queryVector) throws IOException {

		Set<Integer> queryRelDocIDSet = new HashSet<Integer>();

		int[][] termDocFreq = new int[queryVector.size()][totalDocNum];

		ArrayList<String> termListOfQuery = new ArrayList<String>();

		for (String key : queryVector.keySet()) {

			Term term = new Term(docContentField, key);
			DocsEnum docsEnum = atomicReader.termDocsEnum(term); 

			if (docsEnum != null) {
				termListOfQuery.add(key);
				while (docsEnum.nextDoc() != DocIdSetIterator.NO_MORE_DOCS) {
					int docID = docsEnum.docID();
					queryRelDocIDSet.add(docID);
					termDocFreq[termListOfQuery.size() - 1][docID] = docsEnum.freq();
				}
			}
		}

		ScoreDoc[] scoreDocs = new ScoreDoc[queryRelDocIDSet.size()];

		int scoreDocIndex = 0;

		for (int i:queryRelDocIDSet) {

			double docScore = 0.0;

			int docLength = doc_length_map.get(i);

			double K = k1 * ((1 - b) + b * docLength / avg_doc_length);

			for (int j = 0; j < termListOfQuery.size(); j++) {

				double qtf = queryVector.get(termListOfQuery.get(j)); 
				double QTF = (k3 + 1.0) * qtf / (k3 + qtf);

				int df = term_doc_freq_map.get(termListOfQuery.get(j));
				int tf = termDocFreq[j][i];
				double TF = (k1 + 1) * tf / (K + tf);
				double IDF = log2((totalDocNum - df + 0.5) / (df + 0.5));

				docScore += TF * IDF * QTF;
			}

			scoreDocs[scoreDocIndex++] = new ScoreDoc(i, (float) docScore);
		}

		Arrays.sort(scoreDocs, new ByWeightComparator());

		return scoreDocs;
	}

	public void CreateResult(ScoreDoc[] scoreDocs, QualityQuery[] qualityQueries, int i, Query query, QualityStats[] stats) throws IOException {

		int maxResultNum = Math.min(maxResults, scoreDocs.length);

		ScoreDoc[] resultDocs = new ScoreDoc[maxResultNum];

		System.arraycopy(scoreDocs, 0, resultDocs, 0, maxResultNum); 

		TopDocs td = new TopDocs(maxResultNum, resultDocs, (float) resultDocs[0].score);

		stats[i] = analyzeQueryResults(qualityQueries[i], query, td, judge, qualityLog, 1);

		submitRep.report(qualityQueries[i], td, docNameField, searcher);
		submitRep.flush();
	}

	private QualityStats analyzeQueryResults(QualityQuery qq, Query q, TopDocs td, Judge judge, PrintWriter logger, long searchTime) throws IOException {

		QualityStats qualityStats = new QualityStats(judge.maxRecall(qq), searchTime);

		long t1 = System.currentTimeMillis();

		ScoreDoc[] scoreDocs = td.scoreDocs;

		DocNameExtractor xt = new DocNameExtractor(docNameField);

		for (int i = 0; i < scoreDocs.length; i++) {

			String docName = xt.docName(searcher, scoreDocs[i].doc);

			long docNameExtractTime = System.currentTimeMillis() - t1;

			t1 = System.currentTimeMillis();

			boolean isRelevant = judge.isRelevant(docName, qq);

			qualityStats.addResult(i + 1, isRelevant, docNameExtractTime);
		}

		if (logger != null) {
			logger.println(qq.getQueryID() + "  -  " + q);
			qualityStats.log(qq.getQueryID() + " Stats:", 1, logger, "  ");
		}

		return qualityStats;
	}


	public void WriteResult(QualityStats[] stats) {
		QualityStats avg = QualityStats.average(stats);
		avg.log("SUMMARY", 2, qualityLog, "  ");
	}
	
	
	/**
	 * 
	 * @param A 
	 * @param B 
	 * @return A+B
	 */
	public static HashMap<String, Double> addVector(HashMap<String, Double> A, HashMap<String, Double> B) {

		HashMap<String, Double> C = new HashMap<String, Double>();

		for (Entry<String, Double> a : A.entrySet()) {
			if (B.containsKey(a.getKey())) {
				C.put(a.getKey(), (a.getValue() + B.get(a.getKey())));
			} else {
				C.put(a.getKey(), a.getValue());
			}
		}

		for (Entry<String, Double> b : B.entrySet()) {
			if (!C.containsKey(b.getKey())) {
				C.put(b.getKey(), b.getValue());
			}
		}

		return C;
	}
	
	/**
	 * @param queryVector 
	 * @param topTermWeightVector
	 * @return newQueryVector
	 */
	public HashMap<String, Double> generateNewQueryVector(HashMap<String, Double> queryVector, HashMap<String, Double> topTermWeightVector) {

		HashMap<String, Double> newQueryVector = new HashMap<String, Double>();

		HashMap<String, Double> alphaQueryVector = new HashMap<String, Double>();

		for (String term : queryVector.keySet()) {
			double weight = queryVector.get(term); 
			alphaQueryVector.put(term, alpha * weight);
		}

		HashMap<String, Double> expansionVector = new HashMap<String, Double>();

		for (String term : topTermWeightVector.keySet()) {
			double weight = topTermWeightVector.get(term);
			expansionVector.put(term, beta * weight);
		}

		newQueryVector = addVector(alphaQueryVector, expansionVector);

		return newQueryVector;
	}
	
	/**
	 * @param ScoreDocs
	 * @return doc_set
	 * @throws IOException
	 */
	public 	Set<Integer> get_doc_set(ScoreDoc[] ScoreDocs){

		Set<Integer> doc_set = new HashSet<Integer>();              
		for(int doc_i=0;doc_i < ScoreDocs.length ;doc_i++)
		{
			doc_set.add(ScoreDocs[doc_i].doc);
		}
		return doc_set;
	}

	/**
	 * @param ScoreDocs 
	 * @return term_list
	 * @throws IOException
	 */
	public ArrayList<Term> get_term_list(ScoreDoc[] ScoreDocs) throws IOException{
		ArrayList<Term> term_list = new ArrayList<Term>();
		for(int n=0;n<ScoreDocs.length;n++)
		{
			int docNum=ScoreDocs[n].doc;                             
			Terms terms = reader.getTermVector(docNum, docContentField);	          
			TermsEnum iterator = terms.iterator(null);                                  
			BytesRef byteRef = null;
			while ((byteRef = iterator.next()) != null) {  
				String term_text = new String(byteRef.bytes, byteRef.offset, byteRef.length);
				Term term = new Term(docContentField,term_text );
				if(!term_list.contains(term))
					term_list.add(term);                                                
			}
		}
		return term_list;
	}
	
	
	/**
	 * @param term_list 
	 * @param doc_set 
	 * @return within_query_freq_map
	 * @throws IOException
	 */
	public HashMap<String, HashMap<Integer, ArrayList<Integer>>> get_term_position(ArrayList<Term> term_list, Set<Integer> doc_set) throws Exception {
		HashMap<String, HashMap<Integer, ArrayList<Integer>>> within_query_freq_map = new HashMap<String, HashMap<Integer, ArrayList<Integer>>>(); 
		AtomicReader atomicReader = SlowCompositeReaderWrapper.wrap(reader);
		for (int k = 0; k < term_list.size(); k++) {
			Term term = term_list.get(k);                                 
			String termText = term.text();
			DocsAndPositionsEnum docsAndPositionsEnum = atomicReader.termPositionsEnum(term);   
			if (docsAndPositionsEnum == null) {
				continue;
			}             
			int doc_id;
			HashMap<Integer, ArrayList<Integer>> Docid_position_map = new HashMap<Integer, ArrayList<Integer>>();     
			while ((doc_id = docsAndPositionsEnum.nextDoc()) != DocIdSetIterator.NO_MORE_DOCS) {
				if(!doc_set.contains(doc_id)) continue; 
				int freq = docsAndPositionsEnum.freq();       
				int position;
				ArrayList<Integer> query_terms_position = new ArrayList<Integer>();  
				for (int j = 0; j < freq; j++) {
					position = docsAndPositionsEnum.nextPosition();
					query_terms_position.add(position);                         
				}
				Docid_position_map.put(doc_id, query_terms_position);       
			}
			within_query_freq_map.put(termText, Docid_position_map);     
		}
		return within_query_freq_map;
	}	
	
	
	/**
	 * @param term_list 
	 * @return freq_array
	 * @throws IOException
	 */
	public int[][] get_freq_array (ArrayList<Term>term_list) throws IOException{
		int[][] freq_array = new int[term_list.size()][totalDocNum];
		for (int k = 0; k < term_list.size(); k++) {
			Term term = term_list.get(k);
			DocsEnum docsEnum = atomicReader.termDocsEnum(term);
			while (docsEnum.nextDoc() != DocIdSetIterator.NO_MORE_DOCS) {
				int doc_id = docsEnum.docID();

				freq_array[k][doc_id] = docsEnum.freq();
			}
		}
		return freq_array;
	}
	
	/**
	 * @param term_list, doc_set, freq_array, within_query_freq_map, PRFScoreDocs 
	 * @return topTermProximityVector
	 * @throws Exception 
	 */
	public HashMap<String, Double> getTermProximityVectorFromPRFScoreDocs(ArrayList<Term> term_list, int[][]freq_array, HashMap<String, HashMap<Integer, ArrayList<Integer>>> within_query_freq_map, ScoreDoc[] PRFScoreDocs) throws Exception {

		HashMap<String, Double> topTermProximityVector = new HashMap<String, Double>();

		ArrayList<HashMap<String, Double>> docVectors = new ArrayList<HashMap<String, Double>>();
		double cl = avg_doc_length * totalDocNum;
		for (int n = 0; n < PRFScoreDocs.length; n++) {                                               
			HashMap<String,Double> docVector = new HashMap<String,Double>();
			int j = PRFScoreDocs[n].doc;  
			double doc_length = doc_length_map.get(j);
			ArrayList<Term> term_list_tmp=new ArrayList<Term>();
			Terms terms = reader.getTermVector(j, docContentField);	    			
			TermsEnum iterator = terms.iterator(null);                                  
			BytesRef byteRef = null;
			while ((byteRef = iterator.next()) != null) {  
				String term_text = new String(byteRef.bytes, byteRef.offset, byteRef.length);
				Term term = new Term(docContentField,term_text );
				if(!term_list_tmp.contains(term))
					term_list_tmp.add(term);                                        
			}

			for (int ti = 0; ti < term_list_tmp.size(); ti++) {                     
				int tftiD = within_query_freq_map.get(term_list_tmp.get(ti).text()).get(j).size();    
				double ctf = term_total_freq_map.get(term_list_tmp.get(ti).text());
				double  ptf = 0;
				double score = 0;
				for (int qj = 0; qj < term_list.size(); qj++) {                     
					double tftiqjD = 0.0;
					int tfqjD = freq_array[qj][j];
					if (tftiD == 0 || tfqjD == 0) {
						continue;
					}

					int qjdf=term_doc_freq_map.get(term_list.get(qj).text());   
					double IDFqj = log2((totalDocNum - qjdf + 0.5) / (qjdf + 0.5)); 

					for (int tk = 0; tk < tftiD; tk++) {                  
						int positionkti = within_query_freq_map.get(term_list_tmp.get(ti).text()).get(j).get(tk);
						for (int qk = 0; qk < tfqjD; qk++) {
							int positionkqj = within_query_freq_map.get(term_list.get(qj).text()).get(j).get(qk);
							double dist = Math.abs(positionkti - positionkqj)/2;
							double kerneltemp = Math.exp(-dist*dist/(2.0*sigma*sigma)); 
							tftiqjD += kerneltemp;
						}     
					}
					double TFtiqjD = tftiqjD;
					ptf += TFtiqjD * IDFqj;
				}

				double prob = (ptf/doc_length);

				score = (prob == 0)? 0:prob*log2(prob/(ctf/cl));

				docVector.put(term_list_tmp.get(ti).text(), score);    
			}

			docVectors.add(docVector);                                              
		}

		HashMap<String, Double> sumOfDocVectors = new HashMap<String, Double>();

		for (int i = 0; i < docVectors.size(); i++) {
			sumOfDocVectors = addVector(sumOfDocVectors, docVectors.get(i));
		}			

		List<Entry<String, Double>> tempList = new ArrayList<Entry<String, Double>>(sumOfDocVectors.entrySet());

		Collections.sort(tempList, new Comparator<Entry<String, Double>>() {
			@Override
			public int compare(Entry<String, Double> o1, Entry<String, Double> o2) {
				return -o1.getValue().compareTo(o2.getValue());
			}
		});

		for (int i = 0; i < Math.min(PRFTermNum, tempList.size()); i++) {
			topTermProximityVector.put(tempList.get(i).getKey(), tempList.get(i).getValue()/PRFScoreDocs.length);
		}

		return StaTools.MAXNormalizedarray(topTermProximityVector);
	}
	
	
	/**
	 * @param A 
	 * @param B 
	 * @return A * B
	 */
	public static HashMap<String, Double> multipleVector(HashMap<String, Double> A, HashMap<String, Double> B) {

		HashMap<String, Double> C = new HashMap<String, Double>();

		for (Entry<String, Double> a : A.entrySet()) {
			if (B.containsKey(a.getKey())&B.get(a.getKey())!=0) {
				C.put(a.getKey(), (a.getValue() * B.get(a.getKey())));
			} 
		}
		return C;
	}

}

