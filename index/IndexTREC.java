package index;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;
import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.index.IndexWriterConfig.OpenMode;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.FSDirectory;
import org.apache.lucene.util.Version;
import org.apache.lucene.benchmark.byTask.feeds.DocData;
import org.apache.lucene.benchmark.byTask.feeds.NoMoreDataException;
import org.apache.lucene.benchmark.byTask.feeds.TrecContentSource;
import org.apache.lucene.benchmark.byTask.utils.Config;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.Field;
import org.apache.lucene.document.FieldType;
import org.apache.lucene.document.StringField;
import org.apache.lucene.document.TextField;
import analyzer.MyStopAndStemmingAnalyzer;

public class IndexTREC {

	private IndexTREC() {
	}

	public static void main(String[] args) throws Exception {
		String usage = "Usage: [-doc docsDir] [-data dataSetName]";
		if (args.length > 0 && ("-h".equals(args[0]) || "-help".equals(args[0]))) {
			System.out.println(usage);
			System.exit(0);
		}

		String docsDir = "./datasets/WT2G/";
		String dataSetName = "WT2G";

		for(int i = 0;i < args.length;i++) {
			if ("-doc".equals(args[i])) {
				docsDir = args[i+1];
				i++;
			} else if ("-data".equals(args[i])) {
				dataSetName = args[i+1];
				i++;
			}
		}
		
		String indexDir = "indexs/index_" + dataSetName;
		createIndex(dataSetName, docsDir, indexDir) ;

	}

	/**
	 * 
	 * 
	 * @param dataSetName
	 * @param docsDir
	 * @param indexDir
	 * @return 
	 * @return
	 * @throws IOException
	 */
	public static void createIndex(String dataSetName, String docsDir, String indexDir) throws IOException {
		Directory dir = FSDirectory.open(new File(indexDir));
		File index = new File(indexDir);

		if (index.exists()) {
			System.err.println(String.format("[%s] exists, pass!!!", index.getName()));
			return;
		} else {
			index.getParentFile().mkdirs();
			System.out.println(String.format("--------------------%s--------------------", index.getName().replace(".txt", "")));
		}
		
		Analyzer analyzer = new MyStopAndStemmingAnalyzer();

		IndexWriterConfig iwc = new IndexWriterConfig(Version.LUCENE_41, analyzer);
		iwc.setOpenMode(OpenMode.CREATE);
		iwc.setRAMBufferSizeMB(256.0);

		IndexWriter writer = new IndexWriter(dir, iwc);

		Properties props = new Properties();
		props.setProperty("print.props", "false");
		props.setProperty("content.source.verbose", "false");
		props.setProperty("content.source.excludeIteration", "true");
		props.setProperty("doc.maker.forever", "false");

		props.setProperty("work.dir", "."); 
		props.setProperty("docs.dir", docsDir);

		props.setProperty("trec.doc.parser", common.MyTrecParser.class.getName());
		props.setProperty("content.source.forever", "false");

		TrecContentSource tcs = new TrecContentSource();
		tcs.setConfig(new Config(props));
		tcs.resetInputs();

		int n = 0;

		DocData dd = new DocData();

		while (true) {

			try {
				dd = tcs.getNextDocData(dd);
			} catch (NoMoreDataException e) {
				break;
			}

			Document doc = new Document();
			doc.add(new StringField("docno", dd.getName(), Field.Store.YES));

			FieldType textWithTermVectors = new FieldType(TextField.TYPE_STORED);
			textWithTermVectors.setIndexed(true);
			textWithTermVectors.setStoreTermVectors(true);

			doc.add(new Field("contents", dd.getBody(), textWithTermVectors));


			if (hasAnalyzedTerms(dd.getBody(), analyzer)) {
				writer.addDocument(doc);
				n++;
			}

			if (n % 10000 == 0) {
				System.out.println(String.format("========== %d documents indexed ==========", n));
			}
		}

		System.out.println(String.format("\n----------The total number of documents containing at least one term. --> %d ----------", n));

		tcs.close();
		writer.forceMerge(1);
		writer.close();
	}

	/**
	 * 
	 * 
	 * @param writer
	 * @param analyzer
	 * @param text
	 * @return
	 * @throws IOException
	 */
	public static boolean hasAnalyzedTerms(String text, Analyzer analyzer) throws IOException {
		ArrayList<String> terms = MyStopAndStemmingAnalyzer.getTermList(analyzer, text);
		return terms.size() != 0;
	}
}
