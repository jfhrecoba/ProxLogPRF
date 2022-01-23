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

public class LLPRF {

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

	private double c = 1.0;

	public LLPRF() {
	}

	public LLPRF(double c) {
		this.c = c;
	}
	
	public LLPRF(double c, int PRFDocNum, int PRFTermNum, double alpha, double beta) {
		this.c = c;
		this.PRFDocNum = PRFDocNum;
		this.PRFTermNum = PRFTermNum;
		this.alpha = alpha;
		this.beta = beta;
	}

	/**
	 * Log-Logistic based Pseudo Relevance Feedback 
	 * @param args
	 * @throws Exception
	 * @author JFH
	 */
	public static void main(String[] args) throws Exception {

		String usage =
				"Usage: [-index indexFile] [-topics topicsFile] [-qrels qrelsFile] [-data dataSetName] [-c c] "
				+ "[-fd PRFDocNum] [-ft PRFTermNum] [-alpha alpha] [-beta beta]";
		if (args.length > 0 && ("-h".equals(args[0]) || "-help".equals(args[0]))) {
			System.out.println(usage);
			System.exit(0);
		}
		
		String indexFile = "indices/index_WT2G";
		String topicsFile = "query-judge/topics.WT2G";
		String qrelsFile = "query-judge/qrels.WT2G";
		String dataSetName = "WT2G";
		
		double c = 1.0;
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
			} else if ("-c".equals(args[i])) {
				c =  Float.parseFloat(args[i+1]);
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

		final LLPRF myllprf = new LLPRF(c, PRFDocNum, PRFTermNum, alpha, beta);
		final String resultFile =  String.format("result/LLPRF/%s-LLPRF-%.2f-%d-%d-%.1f-%.1f.txt", dataSetName, myllprf.c, myllprf.PRFDocNum, myllprf.PRFTermNum, myllprf.alpha, myllprf.beta);
		final String reportFile = String.format("result/LLPRF/%s-LLPRF-%.2f-%d-%d-%.1f-%.1f-report.txt", dataSetName, myllprf.c, myllprf.PRFDocNum, myllprf.PRFTermNum, myllprf.alpha, myllprf.beta); 
		myllprf.Run(topicsFile, qrelsFile, indexFile, resultFile, reportFile); 

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

			HashMap<String, Double> queryVector = generateQueryVectorFromQuery(query);

			ScoreDoc[] firstRetrievalScoreDocs = retrievalByLLogistic(queryVector);

			ScoreDoc[] PRFScoreDocs = new ScoreDoc[Math.min(PRFDocNum, firstRetrievalScoreDocs.length)];

			for (int j = 0; j < PRFScoreDocs.length; j++) {
				PRFScoreDocs[j] = firstRetrievalScoreDocs[j];
			}
			
			HashMap<String, Double> topTermWeightVector = getTopTermWeightVectorFromPRFScoreDocs(PRFScoreDocs);
			
			HashMap<String, Double> newQueryVector = generateNewQueryVector(queryVector, topTermWeightVector);

			ScoreDoc[] secondRetrievalScoreDocs = retrievalByLLogistic(newQueryVector);
			
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

	public ScoreDoc[] retrievalByLLogistic(HashMap<String, Double> queryVector) throws IOException {

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

		for (int i = 0; i < totalDocNum; i++) {

			if (queryRelDocIDSet.contains(i)) {

				double docScore = 0.0;

				int docLength = doc_length_map.get(i);

				for (int j = 0; j < termListOfQuery.size(); j++) {

					double qtf = queryVector.get(termListOfQuery.get(j)); 

					double df = term_doc_freq_map.get(termListOfQuery.get(j));

					double tf = termDocFreq[j][i];
					double TF = tf * log2(1 + c * avg_doc_length / docLength);
					double lambda = 1.0 * df / totalDocNum; 
					double weight = -log2(lambda/(lambda+TF));

					docScore += weight * qtf;
				}

				scoreDocs[scoreDocIndex++] = new ScoreDoc(i, (float) docScore);
			}
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
	 * @param PRFScoreDocs 
	 * @return topTermWeightVector
	 * @throws IOException
	 */
	public HashMap<String, Double> getTopTermWeightVectorFromPRFScoreDocs(ScoreDoc[] PRFScoreDocs) throws IOException {

		HashMap<String, Double> topTermWeightVector = new HashMap<String, Double>();

		ArrayList<HashMap<String, Double>> docVectors = new ArrayList<HashMap<String, Double>>();

		for (int i = 0; i < PRFScoreDocs.length; i++) {

			HashMap<String, Double> docVector = new HashMap<String, Double>();

			int docID = PRFScoreDocs[i].doc;

			double dl = doc_length_map.get(docID);

			Terms terms = reader.getTermVector(docID, docContentField);
			TermsEnum iterator = terms.iterator(null);

			BytesRef byteRef = null;
			while ((byteRef = iterator.next()) != null) {

				String term = new String(byteRef.bytes, byteRef.offset, byteRef.length);
				int df = term_doc_freq_map.get(term);
				double tf = iterator.totalTermFreq();
				double TF = tf * log2(1.0+c*avg_doc_length/dl);
				double lambda = 1.0*df/totalDocNum; 
				double weight = -log2(lambda/(lambda+TF));

				docVector.put(term, weight);
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
			topTermWeightVector.put(tempList.get(i).getKey(), tempList.get(i).getValue());
		}

		return StaTools.MAXNormalizedarray(topTermWeightVector);
	}
	
	/**
	 * 
	 * @param A 
	 * @param B 
	 * @return A + B
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
			expansionVector.put(term, beta * (1.0 / PRFDocNum) * weight);
		}

		newQueryVector = addVector(alphaQueryVector, expansionVector);

		return newQueryVector;
	}

}

