/* This file is used for preprocessing the documents by stopwords removal and stemming. */
package analyzer;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.Reader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.TokenFilter;
import org.apache.lucene.analysis.TokenStream;
import org.apache.lucene.analysis.Tokenizer;
import org.apache.lucene.analysis.core.LowerCaseFilter;
import org.apache.lucene.analysis.core.StopAnalyzer;
import org.apache.lucene.analysis.core.StopFilter;
import org.apache.lucene.analysis.en.PorterStemFilter;
import org.apache.lucene.analysis.standard.StandardTokenizer;
import org.apache.lucene.analysis.tokenattributes.CharTermAttribute;
import org.apache.lucene.analysis.util.CharArraySet;
import org.apache.lucene.util.Version;

public final class MyStopAndStemmingAnalyzer extends Analyzer {

	public static void main(String[] args) throws IOException {
	    System.out.println(getTermList(new MyStopAndStemmingAnalyzer(), "Iran's government is intensifying a birth control program _ despite opposition from radicals _ because the country's fast-growing population is imposing strains on a struggling economy."));
	}

	private static CharArraySet stopWords;

	static {
	    initStopWords();	// adding stopwords
	}

	/**
	 * Adding stopwords from stopwords.txt.
	 * If stopwords.txt doesn't exist, using the stopwords in Lucene.
	 */
	private static final void initStopWords() {

		LineNumberReader reader = null;

		Set<String> set = new HashSet<String>();

		try {
			reader = new LineNumberReader(new FileReader("./stopwords.txt"));

			String stopWord = null;

			while ((stopWord = reader.readLine()) != null) {
				set.add(stopWord.trim());
			}

			stopWords = new CharArraySet(Version.LUCENE_41, set, true);		

		} catch (FileNotFoundException e) {
			System.err.println("The file stopwords.txt doesn't exist, using the stopwords in Lucene！");
			stopWords = StopAnalyzer.ENGLISH_STOP_WORDS_SET;
		} catch (IOException e) {
			System.err.println("Can't read the stopwords.txt，using the stopwords in Lucene！");
			stopWords = StopAnalyzer.ENGLISH_STOP_WORDS_SET;
		} finally {
			try {
				if (reader != null) {
					reader.close();
					reader = null;
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	@Override
	protected TokenStreamComponents createComponents(String fieldName, Reader reader) {

		Tokenizer tokenizer = new StandardTokenizer(Version.LUCENE_41, reader);		

		TokenFilter lowerCaseFilter = new LowerCaseFilter(Version.LUCENE_41, tokenizer);
		TokenFilter stopFilter = new StopFilter(Version.LUCENE_41, lowerCaseFilter, stopWords);
		TokenFilter stemFilter = new PorterStemFilter(stopFilter);

		return new TokenStreamComponents(tokenizer, stemFilter);	
	}

	/**
	 * 
	 * @param analyzer 
	 * @param text 
	 * @return
	 * @throws IOException
	 */
	public static ArrayList<String> getTermList(Analyzer analyzer, String text) throws IOException {

		ArrayList<String> result = new ArrayList<String>();

		Reader reader = new StringReader(text);

		try {
			TokenStream ts = analyzer.tokenStream(null, reader);
			ts.reset(); 

			while(ts.incrementToken()) {
				CharTermAttribute ta = ts.getAttribute(CharTermAttribute.class);
				result.add(ta.toString());
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		return result;
	}
}
