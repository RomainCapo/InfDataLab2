package ch.heigvd.iict.dmg.labo2;

import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.StopwordAnalyzerBase;
import org.apache.lucene.analysis.core.WhitespaceAnalyzer;
import org.apache.lucene.analysis.en.EnglishAnalyzer;
import org.apache.lucene.analysis.standard.StandardAnalyzer;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

import javax.print.attribute.SetOfIntegerSyntax;

public class Evaluation {

	private static Analyzer analyzer = null;

	private static void readFile(String filename, Function<String, Void> parseLine) throws IOException {
		try (BufferedReader br = new BufferedReader(
				new InputStreamReader(new FileInputStream(filename), StandardCharsets.UTF_8))) {
			String line = br.readLine();
			while (line != null) {
				line = line.trim();
				if (!line.isEmpty()) {
					parseLine.apply(line);
				}
				line = br.readLine();
			}
		}
	}

	/*
	 * Reading CACM queries and creating a list of queries.
	 */
	private static List<String> readingQueries() throws IOException {
		final String QUERY_SEPARATOR = "\t";

		List<String> queries = new ArrayList<>();

		readFile("evaluation/query.txt", line -> {
			String[] query = line.split(QUERY_SEPARATOR);
			queries.add(query[1]);
			return null;
		});
		return queries;
	}

	/*
	 * Reading stopwords
	 */
	private static List<String> readingCommonWords() throws IOException {
		List<String> commonWords = new ArrayList<>();

		readFile("common_words.txt", line -> {
			commonWords.add(line);
			return null;
		});
		return commonWords;
	}

	/*
	 * Reading CACM qrels and creating a map that contains list of relevant
	 * documents per query.
	 */
	private static Map<Integer, List<Integer>> readingQrels() throws IOException {
		final String QREL_SEPARATOR = ";";
		final String DOC_SEPARATOR = ",";

		Map<Integer, List<Integer>> qrels = new HashMap<>();

		readFile("evaluation/qrels.txt", line -> {
			String[] qrel = line.split(QREL_SEPARATOR);
			int query = Integer.parseInt(qrel[0]);

			List<Integer> docs = qrels.get(query);
			if (docs == null) {
				docs = new ArrayList<>();
			}

			String[] docsArray = qrel[1].split(DOC_SEPARATOR);
			for (String doc : docsArray) {
				docs.add(Integer.parseInt(doc));
			}

			qrels.put(query, docs);
			return null;
		});
		return qrels;
	}

	public static List<Integer> intersectionBtwTwoList(List<Integer> a, List<Integer> b) {
		List<Integer> result = new ArrayList<Integer>(a);
		result.retainAll(new ArrayList<Integer>(b));
		return result;
	}

	public static List<Double> interpolatedPrecision(List<Double> unInterpolatedList)
	{
		List<Double> interpolatedList = new ArrayList<Double>();
		interpolatedList.add(unInterpolatedList.get(0));
		
		for(int i = 0; i < unInterpolatedList.size();i++)
		{
		interpolatedList.add(Collections.max(unInterpolatedList));
		unInterpolatedList.set(i, -1.0);
		}
		
		return interpolatedList;
	}
	
	public static void main(String[] args) throws IOException {
        ///
        /// Reading queries and queries relations files
        ///
        List<String> queries = readingQueries();
        System.out.println("Number of queries: " + queries.size());

        Map<Integer, List<Integer>> qrels = readingQrels();
        System.out.println("Number of qrels: " + qrels.size());

        double avgQrels = 0.0;
        for (int q : qrels.keySet()) {
            avgQrels += qrels.get(q).size();
        }
        avgQrels /= qrels.size();
        System.out.println("Average number of relevant docs per query: " + avgQrels);

        //TODO student: use this when doing the english analyzer + common words
        List<String> commonWords = readingCommonWords();

        ///
        ///  Part I - Select an analyzer
        ///
        // TODO student: compare Analyzers here i.e. change analyzer to
        // the asked analyzers once the metrics have been implemented
        analyzer = new WhitespaceAnalyzer();


        ///
        ///  Part I - Create the index
        ///
        Lab2Index lab2Index = new Lab2Index(analyzer);
        lab2Index.index("documents/cacm.txt");

        ///
        ///  Part II and III:
        ///  Execute the queries and assess the performance of the
        ///  selected analyzer using performance metrics like F-measure,
        ///  precision, recall,...
        ///
        //int totalDocs = 0; //1.a
        int totalRetrievedDocs = 0; //1.b
        int totalRelevantDocs = 0; //1.c
        int totalRetrievedRelevantDocs = 0; //1.d
        List<Double> AP = new ArrayList<Double>(); //2
        List<Double> RPrecisions = new ArrayList<Double>(); //4
        
        List<Double> precisons = new ArrayList<Double>(); //5
        
        //System.out.println(qrels);

        for(int queryNumber = 0; queryNumber < queries.size();queryNumber++)
        {
            List<Integer> queryResults = lab2Index.search(queries.get(queryNumber));
            List<Integer> qrelResults = qrels.get(queryNumber+1);
            
            if(qrelResults == null)
            {
            	qrelResults = new ArrayList<Integer>();
            }
            
            totalRetrievedDocs += queryResults.size();
            totalRelevantDocs += qrelResults.size();
            
            totalRetrievedRelevantDocs += intersectionBtwTwoList(queryResults, qrelResults).size();
            
            AP.add((double)totalRetrievedRelevantDocs/totalRelevantDocs);
            
            int R = qrelResults.size();
            List<Integer> rRelevantDocs = queryResults.stream().limit(R).collect(Collectors.toList());

            if(rRelevantDocs.size()!=0)
            {
                RPrecisions.add((double) (intersectionBtwTwoList(rRelevantDocs, qrelResults).size()/rRelevantDocs.size()));
            }
            else
            {
            	RPrecisions.add(0.0);
            }
        }
                
        double avgPrecision = totalRetrievedRelevantDocs/(double)totalRetrievedDocs; //1.e
        double avgRecall = totalRetrievedRelevantDocs/(double)totalRelevantDocs; //1.f
        double fMeasure = (2*avgRecall*avgPrecision)/(avgRecall + avgPrecision); //1.g
        
        double meanAveragePrecision = AP.stream().mapToDouble(Double::doubleValue).sum() / AP.size(); //3
        
        double avgRPrecision = RPrecisions.stream().mapToDouble(Double::doubleValue).sum() / RPrecisions.size(); //4
        
        
        // average precision at the 11 recall levels (0,0.1,0.2,...,1) over all queries
        
        double[] avgPrecisionAtRecallLevels = createZeroedRecalls(); //5

        // TODO student
        // compute the metrics asked in the instructions
        // you may want to call these methods to get:
        // -  The query results returned by Lucene i.e. computed/empirical
        //    documents retrieved
        //        List<Integer> queryResults = lab2Index.search(query);
        //
        // - The true query results from qrels file i.e. genuine documents
        //   returned matching a query
        //        List<Integer> qrelResults = qrels.get(queryNumber);




        ///
        ///  Part IV - Display the metrics
        ///

        //TODO student implement what is needed (i.e. the metrics) to be able
        // to display the results
        displayMetrics(totalRetrievedDocs, totalRelevantDocs,
                totalRetrievedRelevantDocs, avgPrecision, avgRecall, fMeasure,
                meanAveragePrecision, avgRPrecision,
                avgPrecisionAtRecallLevels);
    }

	private static void displayMetrics(int totalRetrievedDocs, int totalRelevantDocs, int totalRetrievedRelevantDocs,
			double avgPrecision, double avgRecall, double fMeasure, double meanAveragePrecision, double avgRPrecision,
			double[] avgPrecisionAtRecallLevels) {
		String analyzerName = analyzer.getClass().getSimpleName();
		if (analyzer instanceof StopwordAnalyzerBase) {
			analyzerName += " with set size " + ((StopwordAnalyzerBase) analyzer).getStopwordSet().size();
		}
		System.out.println(analyzerName);

		System.out.println("Number of retrieved documents: " + totalRetrievedDocs);
		System.out.println("Number of relevant documents: " + totalRelevantDocs);
		System.out.println("Number of relevant documents retrieved: " + totalRetrievedRelevantDocs);

		System.out.println("Average precision: " + avgPrecision);
		System.out.println("Average recall: " + avgRecall);

		System.out.println("F-measure: " + fMeasure);

		System.out.println("MAP: " + meanAveragePrecision);

		System.out.println("Average R-Precision: " + avgRPrecision);

		System.out.println("Average precision at recall levels: ");
		for (int i = 0; i < avgPrecisionAtRecallLevels.length; i++) {
			System.out.println(String.format("\t%s: %s", i, avgPrecisionAtRecallLevels[i]));
		}
	}

	private static double[] createZeroedRecalls() {
		double[] recalls = new double[11];
		Arrays.fill(recalls, 0.0);
		return recalls;
	}
}