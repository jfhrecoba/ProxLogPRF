package common;

import java.io.IOException;

import org.apache.lucene.benchmark.byTask.feeds.DocData;
import org.apache.lucene.benchmark.byTask.feeds.TrecContentSource;
import org.apache.lucene.benchmark.byTask.feeds.TrecDocParser;

public class MyTrecParser extends TrecDocParser {

	@Override
	public DocData parse(DocData docData, String name, TrecContentSource trecSrc, StringBuilder docBuf, ParsePathType pathType) throws IOException {

		int mark = 0; 

		docData.clear();
		docData.setName(name);

		String content1 = docBuf.toString().replaceAll("<DOCHDR>[\u0000-\uFFFF]+</DOCHDR>", "");
		String content2 = stripTags(content1, mark).toString();
		docData.setBody(content2);

		return docData;
	}
}
