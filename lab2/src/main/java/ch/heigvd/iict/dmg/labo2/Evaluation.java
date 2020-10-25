package ch.heigvd.iict.dmg.labo2;

import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.CharArraySet;
import org.apache.lucene.analysis.StopwordAnalyzerBase;
import org.apache.lucene.analysis.core.WhitespaceAnalyzer;
import org.apache.lucene.analysis.en.EnglishAnalyzer;
import org.apache.lucene.analysis.standard.StandardAnalyzer;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
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

	public static List<Double> interpolatedPrecision(List<Double> precisions, List<Double> recall) {
		List<Double> precisionLevel = new ArrayList<Double>();

		for (double i = 0.0; i < 1.0; i += 0.1) {
			double iRound = Math.round(i * 100.0) / 100.0;

			int j = 0;
			double recallValue = recall.get(0);
			int indexRecall = 0;

			while (recallValue < iRound && j < precisions.size()) {
				recallValue = recall.get(j);
				indexRecall = j;
				j++;
			}

			if (j >= precisions.size()) {
				precisionLevel.add(0.0);
			} else {
				List<Double> sublist = precisions.subList(indexRecall, precisions.size());
				precisionLevel.add(Collections.max(sublist));
			}
		}
		return precisionLevel;
	}

	private static List<Double>[] computePrecisionRecallAtk(List<Integer> queryResults, List<Integer> qrelResults) {
		List<Double> allPrecision = new ArrayList<Double>();
		List<Double> allRecall = new ArrayList<Double>();
		List<Double> precisionsForAP = new ArrayList<Double>();

		for (int i = 1; i <= queryResults.size(); i++) {
			List<Integer> subQuery = queryResults.stream().limit(i).collect(Collectors.toList());

			int retrievedDocs = subQuery.size();
			int relevantDocs = qrelResults.size();

			int retrievedRelevantDocs = intersectionBtwTwoList(subQuery, qrelResults).size();

			double precision = retrievedRelevantDocs / (double) retrievedDocs;
			allPrecision.add(precision);

			if (relevantDocs != 0)
				allRecall.add(retrievedRelevantDocs / (double) relevantDocs);
			else
				allRecall.add(0.0);

			if (qrelResults.contains(queryResults.get(i - 1))) {
				precisionsForAP.add(precision);
			}
		}

		List<Double>[] results = new ArrayList[3];
		results[0] = allPrecision;
		results[1] = allRecall;
		results[2] = precisionsForAP;

		return results;
	}

	public static void main(String[] args) throws IOException {
		/*
		 * List<Double> recall = new ArrayList<Double>(); recall.add(0.2);
		 * recall.add(0.2); recall.add(0.4); recall.add(0.4); recall.add(0.4);
		 * recall.add(0.6); recall.add(0.6); recall.add(0.6); recall.add(0.8);
		 * recall.add(0.8);
		 * 
		 * List<Double> precisions = new ArrayList<Double>(); precisions.add(1.0);
		 * precisions.add(0.5); precisions.add(0.67); precisions.add(0.5);
		 * precisions.add(0.4); precisions.add(0.5); precisions.add(0.43);
		 * precisions.add(0.38); precisions.add(0.44); precisions.add(0.5);
		 */

		// System.out.println(interpolatedPrecision(precisions, recall));

		/*
		 * List<Integer> relevant = new ArrayList<Integer>(); relevant.add(1);
		 * 
		 * List<Integer> retrived = new ArrayList<Integer>(); //retrived.add(2);
		 * //retrived.add(2); //retrived.add(2); //retrived.add(2);
		 * 
		 * System.out.println(Arrays.toString(computePrecisionRecallAtk(retrived,
		 * relevant)));
		 */

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

		// TODO student: use this when doing the english analyzer + common words
		List<String> commonWords = readingCommonWords();
		CharArraySet charArraySet = new CharArraySet(commonWords, true);

		///
		/// Part I - Select an analyzer
		///
		// TODO student: compare Analyzers here i.e. change analyzer to
		// the asked analyzers once the metrics have been implemented
		analyzer = new WhitespaceAnalyzer();
		//analyzer = new StandardAnalyzer();
		//analyzer = new EnglishAnalyzer();
		//analyzer = new EnglishAnalyzer(charArraySet);

		///
		/// Part I - Create the index
		///
		Lab2Index lab2Index = new Lab2Index(analyzer);
		lab2Index.index("documents/cacm.txt");

		///
		/// Part II and III:
		/// Execute the queries and assess the performance of the
		/// selected analyzer using performance metrics like F-measure,
		/// precision, recall,...
		///
		int totalRetrievedDocs = 0; // 1.b
		int totalRelevantDocs = 0; // 1.c

		List<Double> precisions = new ArrayList<Double>(); //1.e
		List<Double> recalls = new ArrayList<Double>(); //1.f
		List<Double> f1scores = new ArrayList<Double>();
		
		int totalRetrievedRelevantDocs = 0; // 1.d
		List<Double> AP = new ArrayList<Double>(); // 2
		List<Double> RPrecisions = new ArrayList<Double>(); // 4

		double[] avgPrecisionAtRecallLevels = createZeroedRecalls(); // Part 5

		int queryNumber = 0;
		for (; queryNumber < queries.size(); queryNumber++) {
			List<Integer> queryResults = lab2Index.search(queries.get(queryNumber));
			List<Integer> qrelResults = qrels.get(queryNumber + 1);

			if (qrelResults == null) {
				qrelResults = new ArrayList<Integer>();
			}

			totalRetrievedDocs += queryResults.size(); // 1.b
			totalRelevantDocs += qrelResults.size(); // 1.c
			
			totalRetrievedRelevantDocs += intersectionBtwTwoList(queryResults, qrelResults).size(); // 1.d
			
			// 1.e
			double precision = totalRetrievedRelevantDocs / (double) totalRetrievedDocs;
			precisions.add(precision);
			
			//1.f
			double recall = totalRetrievedRelevantDocs / (double) totalRelevantDocs;
			recalls.add(recall);
			
			//1.g
			f1scores.add((2*precision*recall)/(precision+recall));

			// Part 2
			List<Double>[] recallPrecisionAtk = computePrecisionRecallAtk(queryResults, qrelResults);
			List<Double> precisionForAp = recallPrecisionAtk[2];
			if(precisionForAp.size() != 0)
			{
				AP.add(precisionForAp.stream().mapToDouble(Double::doubleValue).sum() / precisionForAp.size());
			}
			else 
			{
				AP.add(0.0);
			}

			// Part 4
			int R = qrelResults.size();
			List<Integer> rRelevantDocs = queryResults.stream().limit(R).collect(Collectors.toList());

			if (rRelevantDocs.size() != 0) {
				RPrecisions.add(
						(double) (intersectionBtwTwoList(rRelevantDocs, qrelResults).size() / rRelevantDocs.size()));
			} else {
				RPrecisions.add(0.0);
			}

			// Part 5
			List<Double> allPrecision = recallPrecisionAtk[0];
			List<Double> allRecall = recallPrecisionAtk[1];

			List<Double> interpolated = interpolatedPrecision(allPrecision, allRecall);

			for (int i = 0; i < avgPrecisionAtRecallLevels.length; i++) {
				avgPrecisionAtRecallLevels[i] += interpolated.get(i);
			}
		}
		
		double avgPrecision = precisions.stream().mapToDouble(Double::doubleValue).sum() / precisions.size(); // 1.e
		double avgRecall = recalls.stream().mapToDouble(Double::doubleValue).sum() / recalls.size(); // 1.f
		double fMeasure = f1scores.stream().mapToDouble(Double::doubleValue).sum() / f1scores.size(); // 1.g

		double meanAveragePrecision = AP.stream().mapToDouble(Double::doubleValue).sum() / AP.size(); // Part 3

		double avgRPrecision = RPrecisions.stream().mapToDouble(Double::doubleValue).sum() / RPrecisions.size(); // Part 4
		
		// Part 5
		for (int i = 0; i < avgPrecisionAtRecallLevels.length; i++) {
			avgPrecisionAtRecallLevels[i] /= queryNumber;
		}
		// average precision at the 11 recall levels (0,0.1,0.2,...,1) over all queries

		// TODO student
		// compute the metrics asked in the instructions
		// you may want to call these methods to get:
		// - The query results returned by Lucene i.e. computed/empirical
		// documents retrieved
		// List<Integer> queryResults = lab2Index.search(query);
		//
		// - The true query results from qrels file i.e. genuine documents
		// returned matching a query
		// List<Integer> qrelResults = qrels.get(queryNumber);

		///
		/// Part IV - Display the metrics
		///

		// TODO student implement what is needed (i.e. the metrics) to be able
		// to display the results
		displayMetrics(totalRetrievedDocs, totalRelevantDocs, totalRetrievedRelevantDocs, avgPrecision, avgRecall,
				fMeasure, meanAveragePrecision, avgRPrecision, avgPrecisionAtRecallLevels);
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