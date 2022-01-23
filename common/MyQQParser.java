package common;
import org.apache.lucene.benchmark.quality.QualityQuery;
import org.apache.lucene.benchmark.quality.QualityQueryParser;
import org.apache.lucene.queryparser.classic.ParseException;
import org.apache.lucene.queryparser.classic.QueryParser;
import org.apache.lucene.queryparser.classic.QueryParserBase;
import org.apache.lucene.search.BooleanClause;
import org.apache.lucene.search.BooleanQuery;
import org.apache.lucene.search.Query;
import org.apache.lucene.util.Version;
import org.apache.lucene.analysis.Analyzer;
import analyzer.MyStopAndStemmingAnalyzer;

/**
 * Simplistic quality query parser. A Lucene query is created by passing 
 * the value of the specified QualityQuery name-value pair(s) into 
 * a Lucene's QueryParser using StandardAnalyzer. */
public class MyQQParser implements QualityQueryParser {

	private String qqNames[];
	private String indexField;
	ThreadLocal<QueryParser> queryParser = new ThreadLocal<QueryParser>();

	/**
	 * Constructor of a simple qq parser.
	 * @param qqNames name-value pairs of quality query to use for creating the query
	 * @param indexField corresponding index field  
	 */
	public MyQQParser(String qqNames[], String indexField) {
		this.qqNames = qqNames;
		this.indexField = indexField;
	}

	/**
	 * Constructor of a simple qq parser.
	 * @param qqName name-value pair of quality query to use for creating the query
	 * @param indexField corresponding index field  
	 */
	public MyQQParser(String qqName, String indexField) {
		this(new String[] { qqName }, indexField);
	}

	/* (non-Javadoc)
	 * @see org.apache.lucene.benchmark.quality.QualityQueryParser#parse(org.apache.lucene.benchmark.quality.QualityQuery)
	 */
	@Override
	public Query parse(QualityQuery qq) throws ParseException {

		QueryParser qp = queryParser.get();

		if (qp==null) {
			Analyzer analyzer = new MyStopAndStemmingAnalyzer();

			qp = new QueryParser(Version.LUCENE_41, indexField, analyzer);
			queryParser.set(qp);
		}

		BooleanQuery bq = new BooleanQuery();

		for (int i = 0; i < qqNames.length; i++)
			bq.add(qp.parse(QueryParserBase.escape(qq.getValue(qqNames[i]))), BooleanClause.Occur.SHOULD);

		return bq;
	}
}
