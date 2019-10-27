package NER;

import gate.AnnotationSet;
import gate.Document;
import gate.Factory;
import gate.ProcessingResource;
import gate.Resource;
import gate.annotation.NodeImpl;
import gate.creole.AbstractLanguageAnalyser;
import gate.creole.ExecutionException;
import gate.creole.ResourceInstantiationException;
import gate.creole.metadata.CreoleParameter;
import gate.creole.metadata.CreoleResource;
import gate.creole.metadata.Optional;
import gate.creole.metadata.RunTime;
import gate.creole.ontology.InvalidURIException;
import gate.util.AnnotationDiffer;
import gate.util.InvalidOffsetException;
import gate.Annotation;
import gate.Corpus;
import gate.CorpusController;
import gate.AnnotationSet;
import gate.DocumentContent;
import gate.FeatureMap;
import gate.Gate;
import gate.Node;
import gate.SimpleAnnotation;
import gate.SimpleAnnotationSet;
import gate.Utils;
import gate.persist.PersistenceException;
import gate.persist.SerialDataStore;
import gate.util.*;
import gate.util.persistence.PersistenceManager;

import org.apache.commons.math3.distribution.BetaDistribution;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Properties;
import java.util.Random;
import java.util.Set;
import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeSet;
import java.util.UUID;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.net.URL;


import junit.framework.Test;

/**
 * This class illustrates how to do simple batch processing with GATE. It loads
 * an application from a .gapp file (created using "Save application state" in
 * the GATE GUI), and runs the contained application over one or more files. The
 * results are written out to txt files, either in GateXML format (all
 * annotation sets preserved, as in "save as XML" in the GUI), or with inline
 * format"). In this example, the output file names are simply the input file
 * names with ".out.xml" appended.
 *
 * To keep the example simple, we do not do any exception handling - any error
 * will cause the process to abort.
 */
public class Bandit {
	private static final String DS_DIR = "file:///Users/Zahraa/Documents/GATE/Datastore";
	private static final long TIMEOUT_BEFORE_GETTING_RESPONSE = 0;
	static List<String> annotationTypes = new ArrayList<String>();
	static String statsFile = null;
	private static int corpusIndex;
	String logFile = null;
	static int Numebr_Of_Systems = 8;
	static double[] perfMetric = new double[Numebr_Of_Systems];
	static double[] averageMetric = new double[Numebr_Of_Systems];
	static int[] numCalls = new int[Numebr_Of_Systems];

	/**
	 * The main entry point. First we parse the command line options (see
	 * usage() method for details), then we take all remaining command
	 * lineANNIE_with_defaults.gapp parameters to be file names to process. Each
	 * file is loaded, processed using the application and the results written
	 * to the output file (inputFile.out.xml).
	 */
	public static void main(String[] args) throws Exception {

		Gate.setGateHome(new File(
				"/Applications/GATE_Developer_8.1-SNAPSHOT-b5085"));
		Gate.setPluginsHome(new File(
				"/Applications/GATE_Developer_8.1-SNAPSHOT-b5085/plugins/"));
		Gate.setSiteConfigFile(new File(
				"/Applications/GATE_Developer_8.1-SNAPSHOT-b5085/gate.xml"));
		Properties props = System.getProperties();
		props.setProperty("gate.home", "http://gate.ac.uk/wiki/code-repository");
		parseCommandLine(args);
		// initialise GATE - this must be done before calling any GATE APIs
		Gate.init();

		//

		// Create a Corpus to use. We recycle the same Corpus object for each
		// iteration. The string parameter to newCorpus() is simply the
		// GATE-internal name to use for the corpus. It has no particular
		// significance.
		// System.out.println("1111");
		// Corpus corpus = Factory.newCorpus("BatchProcessApp Corpus");
		SerialDataStore sds = (SerialDataStore) Factory.openDataStore(
				"gate.persist.SerialDataStore", DS_DIR);

		if (corpusIndex != -1) {
			ExecuteOnCorpus(corpusIndex);
		} else {
			for (int k = 0; k < sds.getLrIds("gate.corpora.SerialCorpusImpl")
					.size(); k++) {
				ExecuteOnCorpus(k);
			}

		}

		System.out.println("All done");
	}

	private static void ExecuteOnCorpus(int corpusIndex2)
			throws PersistenceException, ResourceInstantiationException,
			IOException, ExecutionException, InterruptedException,
			InvalidOffsetException {
		SerialDataStore sds = (SerialDataStore) Factory.openDataStore(
				"gate.persist.SerialDataStore", DS_DIR);
		Out.prln("serial datastore ...");
		Object corpusID = sds.getLrIds("gate.corpora.SerialCorpusImpl").get(
				corpusIndex2);
		FeatureMap corpFeatures = Factory.newFeatureMap();
		corpFeatures.put(sds.LR_ID_FEATURE_NAME, corpusID);
		corpFeatures.put(sds.DATASTORE_FEATURE_NAME, sds);
		// tell the factory to load the Serial Corpus with the specified ID from
		// the specified datastore
		gate.Corpus corpus = (gate.Corpus) Factory.createResource(
				"gate.corpora.SerialCorpusImpl", corpFeatures);

		System.out.println("got the corpus: " + corpus.getName()
				+ " No. of documents: " + corpus.size());
		statsFile = "Corpus Name: " + corpus.getName()
				+ " \n No. of documents:" + corpus.size() + "\n";

		Long currentTime = System.nanoTime();
		String tableString = null;
		// process the files one by one
		// Writing each document xml
		final int TIMEOUT_BEFORE_GETTING_RESPONSE = 5000;
		Corpus tempCorpus = Factory.newCorpus("Temp Corpus");
		File keyGappFile = new File("KeyAnnotation.gapp");
		CorpusController keyApplication = (CorpusController) PersistenceManager
				.loadObjectFromFile(keyGappFile);
		System.out.println("Application state for Key annotation is loaded");

		// application.setCorpus(corpus);
		int[] index = new int[Numebr_Of_Systems]; 
		int nIteration=0; 
		int[] ar = new int[corpus.size()];
		int[] random = shuffleArray(ar);
		tableString = "Loop, index,  MF1_ANNIE, N-ANNIE, Theta-ANNIE,Lower bound,Upper bound";
		tableString += ",index,  MF1_Lingpipe, N-Lingpipe, Theta-Lingpipe,Lower bound,,Upper bound";
		tableString += ",index,  MF1_OpenNLP, N-OpenNLP, Theta-Lingpipe,Lower bound,Upper bound";
		tableString += ",index, MF1_Stanford, N-Stanford, Theta-Stanford,Lower bound,Upper bound";
		tableString += ", index, MF1_Alchemy, N-Alchemy, Theta-Alchemy,Lower bound,Upper bound";
		tableString += ", index, MF1_MC, N-MC, Theta-MC,Lower bound,Upper bound ";	
		tableString += ", index, MF1_Spotlight, N-Spotlight, Theta-Spotlight,Lower bound,Upper bound";
		tableString += ", index, MF1_Razor, N-Razor, Theta-Razor,Lower bound,Upper bound";
//		tableString += ", index, MF1_Semantria, N-Semantria, Theta-Semantria";
		
//		tableString += ", index, MF1_Calais, N-Calais, Theta-Calais";

		tableString += ", Selected ID, Selected Name, Posterior\n";
		Document doc = corpus.get(random[0]);
	

		doc.setPreserveOriginalContent(true);
		// docXMLString = doc.getName()
		// +
		// " Before any processing annotation set names size \n------------\n"
		// + doc.getAnnotationSetNames().size();
		// System.out.println(docXMLString);
		// doc.setMarkupAware(true);
		System.out.print("Processing document " + doc.getName() + "...");

		for (int i = 0; i < Numebr_Of_Systems; i++){
			numCalls[i] = 1;
			index[i]=1;
		}
		
		tempCorpus.clear();
		tempCorpus.add(doc);
		keyApplication.setCorpus(tempCorpus);
		if (!doc.getName().contains(".xml")) {
			System.out.println(doc.getName() + "\n");
			keyApplication.execute();
			System.out.println("done!");
		}
		keyApplication.cleanup();

//		 gappFile = new File("OpenCalais.gapp");
//		CorpusController OpenCalaisApp = (CorpusController) PersistenceManager
//				.loadObjectFromFile(gappFile);
//		System.out.println("Application state for OpenCalais is loaded");
//		OpenCalaisApp.setCorpus(tempCorpus);
//		OpenCalaisApp.execute();
//		System.out.println("Get the results");
//		perfMetric[9] = corpusQA(tempCorpus);
//		averageMetric[9] = perfMetric[9];
//		doc.removeAnnotationSet("Response");
//		
		// application.cleanup();
		
		
//		keyApplication.setCorpus(tempCorpus);
//		if (!doc.getName().contains(".xml")) {
//			System.out.println(doc.getName() + "\n");
//			keyApplication.execute();
//			System.out.println("done!");
//		}
//		keyApplication.cleanup();
		
//		gappFile = new File("Semantria_test.gapp");
//		CorpusController SemantriaApp = (CorpusController) PersistenceManager
//				.loadObjectFromFile(gappFile);
//		System.out.println("Application state for Semantria is loaded");
//		SemantriaApp.setCorpus(tempCorpus);
//		SemantriaApp.execute();
//		System.out.println("Get the results");
//		perfMetric[8] = corpusQA(tempCorpus);
//		averageMetric[8] = perfMetric[8];
//		doc.removeAnnotationSet("Response");
		
		

		gappFile = new File("Alchemy_2.gapp");
		CorpusController AlchemyApp = (CorpusController) PersistenceManager
				.loadObjectFromFile(gappFile);
		System.out.println("Application state for AlchemyAPI is loaded");
		AlchemyApp.setCorpus(tempCorpus);
		AlchemyApp.execute();
		System.out.println("Get the results");
		perfMetric[4] = corpusQA(tempCorpus);
		averageMetric[4] = perfMetric[4];
		doc.removeAnnotationSet("Response");

//		tempCorpus.clear();
//		tempCorpus.add(doc);
//		keyApplication.setCorpus(tempCorpus);
//		if (!doc.getName().contains(".xml")) {
//			System.out.println(doc.getName() + "\n");
//			keyApplication.execute();
//			System.out.println("done!");
//
//		}
//		keyApplication.cleanup();

		gappFile = new File("ANNIE.gapp");
		CorpusController ANNIEapp = (CorpusController) PersistenceManager
				.loadObjectFromFile(gappFile);
		System.out.println("Application state for ANNIE is loaded");
		ANNIEapp.setCorpus(tempCorpus);
		ANNIEapp.execute();
		System.out.println("Get the results");
		perfMetric[0] = corpusQA(tempCorpus);
		averageMetric[0] = perfMetric[0];
		doc.removeAnnotationSet("Response");

//		tempCorpus.add(doc);
//		keyApplication.setCorpus(tempCorpus);
//		if (!doc.getName().contains(".xml")) {
//			System.out.println(doc.getName() + "\n");
//			keyApplication.execute();
//			System.out.println("done!");
//		}
//		keyApplication.cleanup();



//		tempCorpus.add(doc);
//		keyApplication.setCorpus(tempCorpus);
//		if (!doc.getName().contains(".xml")) {
//			System.out.println(doc.getName() + "\n");
//			keyApplication.execute();
//			System.out.println("done!");
//		}
//		keyApplication.cleanup();

		gappFile = new File("Stanford.gapp");
		CorpusController StanfordApp = (CorpusController) PersistenceManager
				.loadObjectFromFile(gappFile);
		System.out.println("Application state for Stanford is loaded");
		StanfordApp.setCorpus(tempCorpus);
		StanfordApp.execute();
		System.out.println("Get the results");
		perfMetric[3] = corpusQA(tempCorpus);
		averageMetric[3] = perfMetric[3];
		doc.removeAnnotationSet("Response");
		
//		tempCorpus.add(doc);
//		keyApplication.setCorpus(tempCorpus);
//		if (!doc.getName().contains(".xml")) {
//			System.out.println(doc.getName() + "\n");
//			keyApplication.execute();
//			System.out.println("done!");
//		}
//		keyApplication.cleanup();
		CorpusController LingPipeApp = (CorpusController) PersistenceManager
				.loadObjectFromFile(gappFile);
		gappFile = new File("LingPipe.gapp");
		LingPipeApp = (CorpusController) PersistenceManager
				.loadObjectFromFile(gappFile);
		System.out.println("Application state for LingPipe is loaded");
		LingPipeApp.setCorpus(tempCorpus);
		LingPipeApp.execute();
		System.out.println("Get the results");
		perfMetric[1] = corpusQA(tempCorpus);
		averageMetric[1] = perfMetric[1];
		doc.removeAnnotationSet("Response");

//		tempCorpus.add(doc);
//		keyApplication.setCorpus(tempCorpus);
//		if (!doc.getName().contains(".xml")) {
//			System.out.println(doc.getName() + "\n");
//			keyApplication.execute();
//			System.out.println("done!");
//		}
//		keyApplication.cleanup();
		gappFile = new File("OpenNLP.gapp");
		CorpusController OpenNLPApp = (CorpusController) PersistenceManager
				.loadObjectFromFile(gappFile);
		System.out.println("Application state for OpenNLP is loaded");
		OpenNLPApp.setCorpus(tempCorpus);
		OpenNLPApp.execute();
		System.out.println("Get the results");
		perfMetric[2] = corpusQA(tempCorpus);
		averageMetric[2] = perfMetric[2];
		doc.removeAnnotationSet("Response");

		
		// application.cleanup();
		//
//		tempCorpus.add(doc);
//		keyApplication.setCorpus(tempCorpus);
//		if (!doc.getName().contains(".xml")) {
//			System.out.println(doc.getName() + "\n");
//			keyApplication.execute();
//			System.out.println("done!");
//		}
//		keyApplication.cleanup();

		gappFile = new File("MeaningCloud.gapp");
		CorpusController MCApp = (CorpusController) PersistenceManager
				.loadObjectFromFile(gappFile);
		System.out.println("Application state for MeaningCloud is loaded");
		MCApp.setCorpus(tempCorpus);
		MCApp.execute();
		System.out.println("Get the results");
		perfMetric[5] = corpusQA(tempCorpus);
		averageMetric[5] = perfMetric[5];
//		tempCorpus.clear();
		doc.removeAnnotationSet("Response");
		
		//
//		tempCorpus.add(doc);
//		keyApplication.setCorpus(tempCorpus);
//		if (!doc.getName().contains(".xml")) {
//			System.out.println(doc.getName() + "\n");
//			keyApplication.execute();
//			System.out.println("done!");
//		}
//		keyApplication.cleanup();

		gappFile = new File("Razor.gapp");
		CorpusController RazorAPP = (CorpusController) PersistenceManager
				.loadObjectFromFile(gappFile);
		System.out.println("Application state for Razor is loaded");
		RazorAPP.setCorpus(tempCorpus);
		RazorAPP.execute();
		System.out.println("Get the results");
		perfMetric[7] = corpusQA(tempCorpus);
		averageMetric[7] = perfMetric[7];
		doc.removeAnnotationSet("Response");
		
		//
//		tempCorpus.add(doc);
//		keyApplication.setCorpus(tempCorpus);
//		if (!doc.getName().contains(".xml")) {
//			System.out.println(doc.getName() + "\n");
//			keyApplication.execute();
//			System.out.println("done!");
//		}
//		keyApplication.cleanup();
		// application = null;
		gappFile = new File("Spotlight.gapp");
		CorpusController SpotlightAPP = (CorpusController) PersistenceManager
				.loadObjectFromFile(gappFile);
		System.out.println("Application state for Spotlight is loaded");
		SpotlightAPP.setCorpus(tempCorpus);
		SpotlightAPP.execute();
		System.out.println("Get the results");
		perfMetric[6] = corpusQA(tempCorpus);
		averageMetric[6] = perfMetric[6];
		doc.removeAnnotationSet("Response");
		
		//
// The main iteration is running until convergence  or reached the max number of documents. If converged, cont variable is false. 
		boolean cont = true;
		boolean firstRun = true;
		int iterations = 0;
		
		int maxindex= getMaxIndex(index); 
		while (cont && maxindex < corpus.size() - 1) {
			if (firstRun) {
				if (maxindex + 20 < corpus.size() - 1)
					iterations = 20;
				else
					iterations = corpus.size() - 1 - maxindex;
			} else {
				if (maxindex + 10 < corpus.size() - 1)
					iterations = 10;
				else
					iterations = corpus.size() - 1 - maxindex;
			}
			firstRun = false;
			for (int m = 0; m < iterations; m++) {
				System.out.println("Loop iteration " + m + " \n");
				tableString += nIteration  + "," ;
				// select theta for systems
				double[] theta = new double[Numebr_Of_Systems];
				for (int j = 0; j < Numebr_Of_Systems; j++) {
					// Beta distribution for each system
					double[] betaArr = beta(averageMetric[j], numCalls[j]);
					theta[j] = betaArr[0];
					System.out.println("\nAverageF1[" + j + "]= "
							+ averageMetric[j] + " N = " + numCalls[j] + "\n");
					System.out.println("Theta[" + j + "]= " + theta[j] + "Last document: "+ index[j]+"\n");
					
					tableString +=random[index[j]]+","+ averageMetric[j] + "," + numCalls[j] + ","
							+ theta[j] + "," + betaArr[1] + "," + betaArr[2]
							+ ",";
				}
// Selecting the winner
				double[] k = argmax(theta);
//				pointer= (int) k[0]; 
				String name = getSystemName((int) k[0]);
				double posterior = k[1];
				String shortName = name.substring(0, name.indexOf("."));
				System.out.println("Selected system  index= " + k[0]
						+ " name= " + shortName + "\n");

				tableString += k[0] + "," + shortName + "," + posterior + "\n";

				gappFile = new File(name);

				// Run selected tool on the next documet
				doc = corpus.get(random[index[(int)k[0]]]);
				index[(int)k[0]]++;
				doc.setPreserveOriginalContent(true);
				System.out
						.print("Processing document " + doc.getName() + "...");
				tempCorpus.clear();
				tempCorpus.add(doc);
				keyApplication.setCorpus(tempCorpus);
				if (!doc.getName().contains(".xml")) {
					System.out.println(doc.getName() + "\n");
					keyApplication.execute();
					System.out.println("done!");
				}
				keyApplication.cleanup();
				//Apply the winner on the next document
				int nn = (int) k[0];
			
				try {
					switch (nn) {
					case 0:
						System.out.println("K " + (int) k[0] + " ANNIE \n");
						ANNIEapp.setCorpus(tempCorpus);
						ANNIEapp.execute();
						// ANNIEapp.cleanup();
						break;

					case 1:
						System.out.println("K " + (int) k[0] + " LingPipe \n");
						LingPipeApp.setCorpus(tempCorpus);
						LingPipeApp.execute();
						// LingPipeApp.cleanup();
						break;

					case 2:
						System.out.println("K " + (int) k[0] + " OpenNLP \n");
						OpenNLPApp.setCorpus(tempCorpus);
						OpenNLPApp.execute();
						// OpenNLPApp.cleanup();
						break;

					case 3:
						System.out.println("K " + (int) k[0] + " Stanford \n");
						StanfordApp.setCorpus(tempCorpus);
						StanfordApp.execute();
						// StanfordApp.cleanup();
						break;

					case 4:
						System.out.println("K " + (int) k[0] + " Alchemy \n");
						AlchemyApp.setCorpus(tempCorpus);
						AlchemyApp.execute();
//						 StanfordApp.cleanup();
						break;

					case 5:
						System.out.println("K " + (int) k[0]
								+ " MeaningCloud \n");
						MCApp.setCorpus(tempCorpus);
						MCApp.execute();
						// StanfordApp.cleanup();
						break;

//					case 9:
//						System.out
//								.println("K " + (int) k[0] + " OpenCalais \n");
//						OpenCalaisApp.setCorpus(tempCorpus);
//						OpenCalaisApp.execute();
						// StanfordApp.cleanup();
//						break;

					case 7:
						System.out.println("K " + (int) k[0] + " Razor \n");
						RazorAPP.setCorpus(tempCorpus);
						RazorAPP.execute();
						// StanfordApp.cleanup();
						break;

//					case 8:
//						System.out.println("K " + (int) k[0] + " Semantria \n");
//							SemantriaApp.setCorpus(tempCorpus);
//							SemantriaApp.execute();
////						 StanfordApp.cleanup();
//						break;

					case 6:
						System.out.println("K " + (int) k[0] + " Spotlight \n");
						SpotlightAPP.setCorpus(tempCorpus);
						SpotlightAPP.execute();
						// StanfordApp.cleanup();
						break;

					}
				} catch (Exception e) {
					writeToFile(tableString, 0);
					e.printStackTrace();
				}

				System.out.println("Get the results");
				perfMetric[(int) k[0]] = corpusQA(tempCorpus);
				doc.removeAnnotationSet("Response");
				updateMetrics((int) k[0]);
				tempCorpus.clear();
				System.out.println("average performance for system "
						+ getSystemName((int) k[0]) +" = "+ averageMetric[(int) k[0]]
						+ " no of calls " + numCalls[(int) k[0]] + "\n");

			}
//Check Stopping criteria 
			double winner[] = convergenceCheck(averageMetric, numCalls);
			if (winner[0] == -1) {
				System.out
						.println("Confidence level below target for all systems");
			} else {

				System.out.println("Winner is: " + winner[0]
						+ " With probability: " + winner[1] + "\n");
				cont = false;
			}

		}

		System.out.println("Conclusion\n--------\n");
		System.out.println("Number of documents processed : " + index + "\n");
		for (int i = 0; i < Numebr_Of_Systems; i++) {
			System.out.println("average performance of " + getSystemName(i)
					+ " is " + averageMetric[i] + " no of calls " + numCalls[i]
					+ "\n");
			tableString += "\n" + i + "," + averageMetric[i] + ","
					+ numCalls[i];
			// writeToFile(tableString, 0);
		}
		writeToFile(tableString, 0);

	}

	private static int getMaxIndex(int[] index) {
		int max=0; 
		for (int i=0; i< index.length;i++){
			if(index[i]>max)
				max= index[i]; 
		}
		return max;
	}

	private static int[] shuffleArray(int[] ar) {

		for (int i = 0; i < ar.length; i++) {
			ar[i] = i;
		}
		Random rnd = new Random();
		for (int i = ar.length - 1; i > 0; i--) {
			int j = rnd.nextInt(i + 1);
			// Simple swap
			int a = ar[j];
			ar[j] = ar[i];
			ar[i] = a;
		}

		return ar;
	}

	private static double[] convergenceCheck(double[] avergae, int[] N)
			throws IOException {
		double theta[] = new double[Numebr_Of_Systems];
		double probability[] = new double[Numebr_Of_Systems];
		String stats = "";
		for (int i = 0; i < 1000; i++) {
			stats += i + ",";
			for (int j = 0; j < Numebr_Of_Systems; j++) {
				theta[j] = beta(averageMetric[j], numCalls[j])[0];
				stats += theta[j] + ",";
				// System.out.println("theta for system " + getSystemName(j)
				// + " in loop: " + j + " = " + theta[j] + "\n");
			}
			// System.out.println("Max theta for system: " + argmax(theta)[0]
			// + " With value: " + argmax(theta)[1] + "\n");
			stats += argmax(theta)[0] + "," + argmax(theta)[1] + "\n";
			probability[(int) argmax(theta)[0]]++;

		}
		for (int j = 0; j < Numebr_Of_Systems; j++) {
			probability[j] = probability[j] / 1000;
			System.out.println("Probability[" + j + "]=" + probability[j]
					+ "\n");
			stats += "Probability[" + j + "]=" + probability[j] + "\n";
		}

		System.out.println("selected system: " + argmax(probability)[0]
				+ " Of max value :" + argmax(probability)[1] + "\n");
		double winner[] = new double[2];
		//70% sure this is the winner
		if (argmax(probability)[1] > 0.7) {

			winner[0] = argmax(probability)[0];
			winner[1] = argmax(probability)[1];
			return winner;
		}
		winner[0] = -1;
		writeToFile(stats, 1);
		return winner;
	}

	private static void writeToFile(String statsFile, int i) throws IOException {
		long time = System.currentTimeMillis();
		String fileName = "";
		if (i == 0)
			fileName = "runStats_" + time + ".csv";
		else
			fileName = "probStats_" + time + ".csv";

		File outputFile = new File("Output/", fileName);
		FileOutputStream fos;
		try {
			fos = new FileOutputStream(outputFile);
			BufferedOutputStream bos = new BufferedOutputStream(fos);
			System.out.println("Writing to Stats file");
			OutputStreamWriter out = new OutputStreamWriter(bos);
			out.write(statsFile);
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		System.out.println("Done with writing to files");
		return;

	}

	private static double[] beta(double x, int n) {
		double term1 = x * n;
		double term2 = n * (1 - x);
		BetaDistribution beta = new BetaDistribution(term1 + 1, term2 + 1);
		double[] b = new double[3];
		b[0] = beta.sample();
		double delta = 0.25;
		b[1] = beta.inverseCumulativeProbability(delta);
		b[2] = beta.inverseCumulativeProbability(1 - delta);
		// System.out.println(" Beta sample : " + b);
		return b;
	}

	private static void updateMetrics(int k) {
		int oldN = numCalls[k];
		numCalls[k]++;
		double oldSum = averageMetric[k] * oldN;
		averageMetric[k] = (oldSum + perfMetric[k]) / numCalls[k];
	}

	private static String getSystemName(int k) {
		// TODO Auto-generated method stub

		switch (k) {

		case 0:
			return "ANNIE.gapp";
		case 1:
			return "LingPipe.gapp";
		case 2:
			return "OpenNLP.gapp";
		case 3:
			return "Stanford.gapp";
		case 4:
			return "Alchemy_2.gapp";
		case 5:
			return "MeaningCloud.gapp";
//		case 9:
//			return "OpenCalais.gapp";
		case 7:
			return "Razor.gapp";
//		case 8:
//			return "Semantria_test.gapp";
		case 6:
			return "Spotlight.gapp";
		default:
			return "";
		}

	}
// This class is to choose the winner among theta (the perfomance metrics for all systems)
	private static double[] argmax(double[] theta) {
		// TODO Auto-generated method stub
		int n = 0;
		double[] k = new double[2];
		double temp = theta[0];
		for (int i = 1; i < theta.length; i++) {
			if (theta[i] > temp) {
				temp = theta[i];
				n = i;
			}
		}
		k[0] = n;
		k[1] = temp;
		return k;
	}
// For evaluating teh performance o
	private static double corpusQA(Corpus corpus) throws IOException {

		final int FSCORE_MEASURES = 0;
		ArrayList<String> documentNames = new ArrayList<String>();
		TreeSet<String> types = new TreeSet<String>();
		Set<String> features = new HashSet<String>();
		List<Map<String, AnnotationDiffer>> differsByDocThenType = new ArrayList<Map<String, AnnotationDiffer>>();
		int measuresType = FSCORE_MEASURES;
		Object[] measures = new Object[] { "F1.0-score strict",
				"F1.0-score lenient", "F1.0-score averge" };
		String keySetName = "Key";
		String responseSetName = "Response";
		System.out.println("annotationScheme:" + annotationScheme + "\n");
		if (annotationScheme.toLowerCase().trim().contains("mentions")) {
			types.add("Entity");
		} else if (annotationScheme.toLowerCase().trim().contains("3types")) {
			types.add("person");
			types.add("location");
			types.add("org");
		} else if (annotationScheme.toLowerCase().trim().contains("4types")) {
			types.add("person");
			types.add("location");
			types.add("org");
			types.add("date");
		}

		annotationTypes = new ArrayList<String>();
		annotationTypes.addAll(types);
		statsFile += "Key set: " + keySetName + "\n";
		statsFile += "Response set: " + responseSetName + "\n";
		// for each document
		statsFile += "\n Statistics on corpus documents (per doc)\n";
		for (int row = 0; row < corpus.size(); row++) {
			boolean documentWasLoaded = corpus.isDocumentLoaded(row);
			Document document = (Document) corpus.get(row);
			documentNames.add(document.getName());
			Set<Annotation> keys = new HashSet<Annotation>();
			Set<Annotation> responses = new HashSet<Annotation>();
			// get annotations from selected annotation sets
			keys = document.getAnnotations(keySetName);
			responses = document.getAnnotations(responseSetName);
			if (!documentWasLoaded) { // in case of datastore
				corpus.unloadDocument(document);
				Factory.deleteResource(document);
			}

			// fscore document table
			if (measuresType == FSCORE_MEASURES) {
				HashMap<String, AnnotationDiffer> differsByType = new HashMap<String, AnnotationDiffer>();
				AnnotationDiffer differ;
				Set<Annotation> keysIter = new HashSet<Annotation>();
				Set<Annotation> responsesIter = new HashSet<Annotation>();
				for (String type : types) {
					if (!keys.isEmpty() && !types.isEmpty()) {
						keysIter = ((AnnotationSet) keys).get(type);
					}
					if (!responses.isEmpty() && !types.isEmpty()) {
						responsesIter = ((AnnotationSet) responses).get(type);
					}
					differ = new AnnotationDiffer();
					differ.setSignificantFeaturesSet(features);
					differ.calculateDiff(keysIter, responsesIter); // compare
					differsByType.put(type, differ);
				}

				differsByDocThenType.add(differsByType);
				differ = new AnnotationDiffer(differsByType.values());
				List<String> measuresRow = differ.getMeasuresRow(measures,
						documentNames.get(documentNames.size() - 1));
				// }

				statsFile += Arrays.deepToString(measuresRow.toArray()) + "\n";
				System.out.println(Arrays.deepToString(measuresRow.toArray()));

			}
		}
		statsFile += "\n Statistics on types (per type) \n";
		System.out.println("Corpus Stats");
		String[] corpusStats = calculateCorpusStats(differsByDocThenType);
		statsFile += corpusStats[0];
		String s = gappFile.getName().substring(0,
				gappFile.getName().indexOf("."));
		System.out.println(corpusStats[0]);
		return Double.valueOf(corpusStats[1]);
	}

	/**
	 * Parse command line options.
	 */
	private static void parseCommandLine(String[] args) throws Exception {
		int i;
		// iterate over all options (arguments starting with '-')
		for (i = 0; i < args.length && args[i].charAt(0) == '-'; i++) {
			switch (args[i].charAt(1)) {
			// -a type = write out annotations of type a.
			case 'a':
				if (annotTypesToWrite == null)
					annotTypesToWrite = new ArrayList();
				annotTypesToWrite.add(args[++i]);
				break;

			// -g gappFile = path to the saved application
			// case 'g':
			// gappFile = new File(args[++i]);
			// break;

			// -e encoding = character encoding for documents
			case 'e':
				encoding = args[++i];

			case 'c':
				corpusIndex = Integer.valueOf(args[++i]);
				break;
			case 's':
				annotationScheme = String.valueOf(args[++i]).toLowerCase()
						.trim();
				break;

			// default:
			// System.err.println("Unrecognised option " + args[i]);
			// usage();
			}
		}

	}

	/**
	 * Print a usage message and exit.
	 */
	private static final void usage() {
		System.err
				.println("Usage:\n"
						+ "   java sheffield.examples.BatchProcessApp -g <gappFile> [-e encoding]\n"
						+ "            [-a annotType] [-a annotType] file1 file2 ... fileN\n"
						+ "\n"
						+ "-g gappFile : (required) the path to the saved application state we are\n"
						+ "              to run over the given documents.  This application must be\n"
						+ "              a \"corpus pipeline\" or a \"conditional corpus pipeline\".\n"
						+ "\n"
						+ "-e encoding : (optional) the character encoding of the source documents.\n"
						+ "              If not specified, the platform default encoding (currently\n"
						+ "              \""
						+ "-s Annotation Scheme : (optional) mentions for only entity mentions"
						+ " \n 3types for person, location and org"
						+ " \n 4types for person, location, org and date"
						+ "              \""
						+ System.getProperty("file.encoding")
						+ "\") is assumed.\n"
						+ "\n"
						+ "-a type     : (optional) write out just the annotations of this type as\n"
						+ "              inline XML tags.  Multiple -a options are allowed, and\n"
						+ "              annotations of all the specified types will be output.\n"
						+ "              This is the equivalent of \"save preserving format\" in the\n"
						+ "              GATE GUI.  If no -a option is given the whole of each\n"
						+ "              processed document will be output as GateXML (the equivalent\n"
						+ "              of \"save as XML\").");

		System.exit(1);
	}

	/** Index of the first non-option argument on the command line. */
	private static int firstFile = 0;

	/** Path to the saved application file. */
	private static File gappFile = null;

	/**
	 * List of annotation types to write out. If null, write everything as
	 * GateXML.
	 */
	private static List annotTypesToWrite = null;

	/**
	 * The character encoding to use when loading the docments. If null, the
	 * platform default encoding is used.
	 */
	private static String encoding = null;
	private static String annotationScheme = null;

	private static String[] calculateCorpusStats(
			List<Map<String, AnnotationDiffer>> differsByDocThenType) {

		// annotation types found in the document
		String[] typesNames = new String[annotationTypes.size() + 2];

		// column names used in the html table
		String[] colnames = { "Annotation Type", "Match", "Only in Key",
				"Only in Response", "Overlap", "Rec.B/A", "Prec.B/A",
				"f1.0-strict", "Rec.B/A", "Prec.B/A", "f1.0-lenient",
				"Rec.B/A", "Prec.B/A", "f1.0-Average" };

		// last two rows used for macro and micro averages
		double[][] vals = new double[annotationTypes.size() + 2][13];

		// one annotation type at a time
		for (int rowIndex = 0; rowIndex < annotationTypes.size(); rowIndex++) {
			// get the counts and measures for the current document/row
			String type = annotationTypes.get(rowIndex);

			// by iterating over all documents, obtain differs created for the
			// type
			// under consideration
			ArrayList<AnnotationDiffer> differs = new ArrayList<AnnotationDiffer>();
			for (Map<String, AnnotationDiffer> differsByType : differsByDocThenType) {
				differs.add(differsByType.get(type));
			}
			System.out.println("calculate various stats");
			// calculate various stats
			AnnotationDiffer differ = new AnnotationDiffer(differs);
			typesNames[rowIndex] = type;
			vals[rowIndex][0] = differ.getCorrectMatches();
			vals[rowIndex][1] = differ.getMissing();
			vals[rowIndex][2] = differ.getSpurious();
			vals[rowIndex][3] = differ.getPartiallyCorrectMatches();

			double[] tempvals = getMeasureValue(differ, "f1.0-strict");
			vals[rowIndex][4] = tempvals[0];
			vals[rowIndex][5] = tempvals[1];
			vals[rowIndex][6] = tempvals[2];
			tempvals = getMeasureValue(differ, "f1.0-lenient");
			vals[rowIndex][7] = tempvals[0];
			vals[rowIndex][8] = tempvals[1];
			vals[rowIndex][9] = tempvals[2];
			tempvals = getMeasureValue(differ, "f1.0-average");
			vals[rowIndex][10] = tempvals[0];
			vals[rowIndex][11] = tempvals[1];
			vals[rowIndex][12] = tempvals[2];
		}
		System.out.println("macro summary");
		// macro summary
		int i = annotationTypes.size();
		typesNames[i] = "Macro Summary";

		for (int row = 0; row < annotationTypes.size(); row++) {
			vals[i][4] += vals[row][4];
			vals[i][5] += vals[row][5];
			vals[i][6] += vals[row][6];
			vals[i][7] += vals[row][7];
			vals[i][8] += vals[row][8];
			vals[i][9] += vals[row][9];
			vals[i][10] += vals[row][10];
			vals[i][11] += vals[row][11];
			vals[i][12] += vals[row][12];
		}
		vals[i][4] = vals[i][4] / annotationTypes.size();
		vals[i][5] = vals[i][5] / annotationTypes.size();
		vals[i][6] = vals[i][6] / annotationTypes.size();
		vals[i][7] = vals[i][7] / annotationTypes.size();
		vals[i][8] = vals[i][8] / annotationTypes.size();
		vals[i][9] = vals[i][9] / annotationTypes.size();
		vals[i][10] = vals[i][10] / annotationTypes.size();
		vals[i][11] = vals[i][11] / annotationTypes.size();
		vals[i][12] = vals[i][12] / annotationTypes.size();
		System.out.println("micro summary");
		// micro summary
		i++;
		typesNames[i] = "Micro Summary";
		for (int row = 0; row < annotationTypes.size(); row++) {
			vals[i][0] += vals[row][0];
			vals[i][1] += vals[row][1];
			vals[i][2] += vals[row][2];
			vals[i][3] += vals[row][3];
		}

		ArrayList<AnnotationDiffer> differs = new ArrayList<AnnotationDiffer>();
		for (Map<String, AnnotationDiffer> differsByType : differsByDocThenType) {
			differs.addAll(differsByType.values());
		}
		AnnotationDiffer differ = new AnnotationDiffer(differs);

		double[] tempvals = getMeasureValue(differ, "f1.0-strict");
		vals[i][4] = tempvals[0];
		vals[i][5] = tempvals[1];
		vals[i][6] = tempvals[2];
		tempvals = getMeasureValue(differ, "f1.0-lenient");
		vals[i][7] = tempvals[0];
		vals[i][8] = tempvals[1];
		vals[i][9] = tempvals[2];
		tempvals = getMeasureValue(differ, "f1.0-average");
		vals[i][10] = tempvals[0];
		vals[i][11] = tempvals[1];
		vals[i][12] = tempvals[2];

		// populate the html table with values
		return toTable(typesNames, null, vals, colnames);
	}

	private static String[] toTable(String[] types, String[] anchorsOnFirstCol,
			double[][] vals, String[] columnNames) {
		StringBuffer buffer = new StringBuffer();
		buffer.append("\n");
		buffer.append("\n");

		// add column titles
		for (String s : columnNames) {
			buffer.append(s);
			buffer.append("\t");
		}
		buffer.append("\t\n");

		// doc name followed by values as calculated earlier
		for (int i = 0; i < types.length; i++) {
			buffer.append(types[i] + "\t");
			NumberFormat f = NumberFormat.getInstance(Locale.ENGLISH);
			double[] colvals = vals[i];
			for (double v : colvals) {
				buffer.append(f.format(v));
				buffer.append("\t");
			}
			buffer.append("\t\n");
		}
		int k = types.length - 1;
		String[] results = new String[2];
		results[0] = buffer.toString();
		results[1] = Double.toString(vals[k][12]);
		return results;
	}

	private static double[] getMeasureValue(AnnotationDiffer differ,
			String measure) {
		double[] vals = new double[3];
		// recall
		if (measure.endsWith("strict")) {
			vals[0] = differ.getRecallStrict();
		} else if (measure.endsWith("lenient")) {
			vals[0] = differ.getRecallLenient();
		} else {
			vals[0] = differ.getRecallAverage();
		}
		// precision
		if (measure.endsWith("strict")) {
			vals[1] = differ.getPrecisionStrict();
		} else if (measure.endsWith("lenient")) {
			vals[1] = differ.getPrecisionLenient();
		} else {
			vals[1] = differ.getPrecisionAverage();
		}

		// f-measure
		double beta = Double
				.valueOf(measure.substring(1, measure.indexOf('-')));
		if (measure.endsWith("strict")) {
			vals[2] = differ.getFMeasureStrict(beta);
		} else if (measure.endsWith("lenient")) {
			vals[2] = differ.getFMeasureLenient(beta);
		} else {
			vals[2] = differ.getFMeasureAverage(beta);
		}
		return vals;
	}
}