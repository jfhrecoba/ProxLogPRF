package models;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
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

public class BM25 {

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


	private double k1 = 1.2;
	private double b = 0.35;
	final double k3 = 8.0;

	public BM25() {
	}

	public BM25(double k1, double b) {
		this.k1 = k1;
		this.b = b;
	}

	/**
	 * BM25
	 * @param args
	 * @throws Exception
	 * @author JFH
	 */
	public static void main(String[] args) throws Exception {

		String usage =
				"Usage: [-index indexFile] [-topics topicsFile] [-qrels qrelsFile] [-data dataSetName] [-k1 k1] [-b b]";
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
			}
		}

		final BM25 mybm25 = new BM25(k1, b);
		final String resultFile = String.format("result/BM25/%s-BM25-%.1f-%.2f.txt", dataSetName, mybm25.k1, mybm25.b);
		final String reportFile = String.format("result/BM25/%s-BM25-%.1f-%.2f-report.txt", dataSetName, mybm25.k1, mybm25.b);
		mybm25.Run(topicsFile, qrelsFile, indexFile, resultFile, reportFile); 

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
		search( );
	}

	public void search() throws Exception {

		QualityStats stats[] = new QualityStats[qualityQueries.length];

		for (int i = 0; i < qualityQueries.length; i++) {

			QualityQuery qq = qualityQueries[i];
			Query query = qqParser.parse(qq);	

			HashMap<String, Double> queryVector = generateQueryVectorFromQuery(query);

			ScoreDoc[] firstRetrievalScoreDocs = retrievalByBM25(queryVector);

			CreateResult(firstRetrievalScoreDocs, qualityQueries, i, query, stats);
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
}

