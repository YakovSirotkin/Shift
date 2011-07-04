package ru.spbau.bioinf.shift;

import org.apache.log4j.Logger;
import org.jdom.Document;
import org.jdom.Element;
import ru.spbau.bioinf.shift.util.ReaderUtil;
import ru.spbau.bioinf.shift.util.XmlUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Properties;

public class ProteinFinder {


    private static final Logger log = Logger.getLogger(ProteinFinder.class);

    public static final double EPSILON = 0.01;
    private File outputDir;
    private ArrayList<Protein> proteins;
    private List<Spectrum> spectrums;
    private PrintWriter matchFile;


    public static void main(String[] args) throws Exception {
        String dataset = "data/set8";
        String proteinDatabaseName = "prot.fasta";

        if (args.length > 0) {
            dataset = args[0];
        }

        if (args.length > 1) {
            proteinDatabaseName = args[1];
        }

        ProteinFinder processor = new ProteinFinder();


        processor.processResult(new File(dataset), proteinDatabaseName);
    }

    public void processResult(File datasetDir, String proteinDatabaseFilename) throws Exception {
        log.debug("star result processing");
        File msinputDir = new File(datasetDir, "msinput");
        BufferedReader input = ReaderUtil.getBufferedReader(new File(msinputDir, "input_data"));
        spectrums = new ArrayList<Spectrum>();
        outputDir = new File(datasetDir, "xml");
        File msoutputDir = new File(datasetDir, "msouput");
        msoutputDir.mkdirs();
        matchFile = createOutputFile(msoutputDir, "match.txt");


        new File(outputDir, "spectrums").mkdirs();

        Properties properties;

        Spectrum s205 = null;
        while ((properties = ReaderUtil.readPropertiesUntil(input, "PRECURSOR_MASS")).size() > 0) {
            Spectrum spectrum = new Spectrum(properties, input);
            spectrums.add(spectrum);
            if (spectrum.getId() == 205) {
                s205 = spectrum;
            }
        }
        log.debug("spectrums data loaded");
        FastaReaderString fastaReader = new FastaReaderString(new File(msinputDir, proteinDatabaseFilename));

        proteins = fastaReader.getProteins();

        log.debug("protein database loaded");

        int[] stat = new int[50];

        int total = 0;


        for (Spectrum spectrum : spectrums) {

            total += processSpectrum(proteins, spectrum);

        }
        matchFile.close();
        /*
        for (int i = 0; i < stat.length; i++) {
            if (stat[i]>0) {
                System.out.println(i + " " +  stat[i]);
            }
        } */
        System.out.println("total = " + total);
    }

    private void testSave(int spectrumId) throws IOException {
        Spectrum sp = spectrums.get(spectrumId);
        Protein p = proteins.get(13822);
        Map<String, List<Double>> pos = getPositions(sp);



        Match m = getScore2(p, sp, pos);
        ArrayList<Match> test = new ArrayList<Match>();
        test.add(m);
        saveSpectrum(sp, pos, test);
    }

    private int processSpectrum(List<Protein> proteins, Spectrum spectrum) throws IOException {
        int ans = 0;
        int bestScore = 0;
        ArrayList<Match> best = new ArrayList<Match>();
        Map<String,List<Double>> positions = getPositions(spectrum);
        for (Protein protein : proteins) {
            Match match = getScore2(protein, spectrum, positions);
            if (match != null) {
                int score = match.getScore();
                if (score >= bestScore) {
                    if (score > bestScore) {
                        best.clear();
                    }
                    best.add(match);
                    bestScore = score;
                }
            }
        }
        if (best.size() > 0 && best.size() < 20) {
            ans = 1;
            //saveSpectrum(spectrum, positions, best);
            for (Match match : best) {
                matchFile.println(spectrum.getId() + " " + match.getProteinId() + " " + match.getScore());
                matchFile.flush();
            }
        } else {
            log.debug("Spectrum " + spectrum.getId() + " is too bad, best.size() is " +  best.size());
        }
        return ans;
    }

    private void saveSpectrum(Spectrum spectrum,  Map<String, List<Double>> positions, ArrayList<Match> best) throws IOException {
        for (Match bestMatch : best) {
            log.debug("bestScore = " + bestMatch.getScore() + " bestProtein " + bestMatch.getProteinId());
            //System.out.println(proteins.get(bestProtein));
        }
        Document doc = new Document();
        Element root = new Element("spectrum");
        doc.setRootElement(root);
        spectrum.addToXml(root, new HashMap<Integer, List<Break>>());
        Element matches = new Element("matches");
        for (Match match : best) {
            Protein protein  = proteins.get(match.getProteinId());double[] sd = spectrum.getData();
            double[] pd = protein.getSpectrum();
            List<Double> shifts = getShifts(protein, positions);

            for (double shift : shifts) {
                match.addShift(shift, getScore(sd, pd, shift));
            }

            matches.addContent(match.toXml());
        }
        root.addContent(matches);

        XmlUtil.saveXml(doc, outputDir + "/spectrums/spectrum" + spectrum.getId() + ".xml");
        log.debug("Match for spectrum " + spectrum.getId() + " saved");
    }


    private Match getScore2(Protein protein, Spectrum spectrum, Map<String, List<Double>> positions) {
        int bestScore = 0;
        Match bestMatch = null;
        double[] sd = spectrum.getData();
        double[] pd = protein.getSpectrum();
        List<Double> shifts = getShifts(protein, positions);

        for (double shift : shifts) {
            int score = getScore(sd, pd, shift);

            if (score > bestScore) {
                bestScore = score;
                bestMatch = new Match(protein, bestScore);
                //System.out.println("s = " + s);
            }


        }
        return bestMatch;
    }

    public static List<Double> getShifts(Protein protein, Map<String, List<Double>> positions) {
        double[] pd = protein.getSpectrum();
        List<Double> shifts = new ArrayList<Double>();
        String acids = protein.getSimplifiedAcids();
        for (Map.Entry<String, List<Double>> entry : positions.entrySet()) {
            String key = entry.getKey();
            List<Double> values = entry.getValue();
            int cur = acids.indexOf(key);
            while (cur >= 0) {
                for (double value : values) {
                    shifts.add(pd[cur]-value);
                }
                cur = acids.indexOf(key, cur +1);
            }
        }
        Collections.sort(shifts);
        shifts = merge(shifts);
        return shifts;
    }

    public static int getScore(double[] sd, double[] pd, double shift) {
        int score = 0;
        int j = 0;
        for (int i = 0; i < pd.length; i++) {
            while (j < sd.length && sd[j] < pd[i] - shift - 0.1) {
                j++;
            }
            if (j < sd.length) {
                if (Math.abs(sd[j] - pd[i] + shift) < 0.1) {
                    score++;
                }
            }
        }
        return score;
    }

    public static Map<String, List<Double>> getPositions(Spectrum spectrum) {
        Map<String, List<Double>> ans = new HashMap<String, List<Double>>();
        double[] values = spectrum.getData();
        for (int i = 0; i < values.length - 1; i++) {
            for (int j = i + 1; j < values.length; j++) {
                double diff = values[j] - values[i];
                if (diff > 300)
                    break;
                for (Map.Entry<Character, Double> entry : Acids.acids.entrySet()) {
                    if (Math.abs(entry.getValue() - diff) < EPSILON) {
                        for (int k = j + 1; k < values.length; k++) {
                            double diff2 = values[k] - values[j];
                            if (diff2 > 300)
                                break;
                            for (Map.Entry<Character, Double> entry2 : Acids.acids.entrySet()) {
                                if (Math.abs(entry2.getValue() - diff2) < EPSILON) {
                                    for (int m = k + 1; m < values.length; m++) {
                                        double diff3 = values[m] - values[k];
                                        if (diff3 > 300)
                                            break;
                                        for (Map.Entry<Character, Double> entry3 : Acids.acids.entrySet()) {
                                            if (Math.abs(entry3.getValue() - diff3) < EPSILON) {
                                                String key = "" + entry.getKey() + entry2.getKey() + entry3.getKey();
                                                if (!ans.containsKey(key)) {
                                                    ans.put(key, new ArrayList<Double>());
                                                }
                                                ans.get(key).add(values[i]);
                                            }
                                        }
                                    }

                                }
                            }
                        }

                    }
                }
            }
        }
        return ans;
    }

    public static List<Double> merge(List<Double> v) {
        ArrayList<Double> a = new ArrayList<Double>();
        double prev = 0;
        int count = 0;
        double sum = 0;
        for (double s : v) {
            if (s - prev < EPSILON) {
                count++;
                sum += s;
                prev = s;

            } else {
                if (count > 0) {
                    a.add(sum / count);
                }
                count = 1;
                prev = s;
                sum = s;
            }
        }
        a.add(sum / count);
        return a;
    }

    void process(List<Solid> ans, List<PeakUnion> peaks, int pos, String res, LinkedList<PeakUnion> path) {
        PeakUnion cur = peaks.get(pos);
        boolean finish = true;
        for (int i = pos + 1; i < peaks.size(); i++) {
            PeakUnion next = peaks.get(i);
            for (Map.Entry<Character, Double> acid : Acids.acids.entrySet()) {
                if (Math.abs(next.getAverage() - cur.getAverage() - acid.getValue()) < EPSILON) {
                    finish = false;
                    path.add(next);
                    process(ans, peaks, i, res += acid.getKey(), path);
                    path.removeLast();
                }
            }
        }
        if (finish && res.length() > 1) {
            ans.add(new Solid(path.get(0).getAverage(), res));
            /*
            System.out.print(res);
            for (PeakUnion peakUnion : path) {
                System.out.print(" " + peakUnion.getAverage());
            }
            System.out.println();
            */

        }
    }

    void process(List<Solid> ans, double[] peaks, int pos, String res, LinkedList<Double> path) {
        Double cur = peaks[pos];
        boolean finish = true;
        for (int i = pos + 1; i < peaks.length; i++) {
            Double next = peaks[i];
            for (Map.Entry<Character, Double> acid : Acids.acids.entrySet()) {
                if (Math.abs(next - cur - acid.getValue()) < EPSILON) {
                    finish = false;
                    path.add(next);
                    process(ans, peaks, i, res += acid.getKey(), path);
                    path.removeLast();
                }
            }
        }
        if (finish && res.length() > 1) {
            ans.add(new Solid(path.get(0), res));
            if (res.length()>=bl) {
                System.out.println(res);
                bl = res.length();
            }
            /*
            System.out.print(res);
            for (PeakUnion peakUnion : path) {
                System.out.print(" " + peakUnion.getAverage());
            }
            System.out.println();
            */

        }
    }

    int bl =10;

    public static class Solid {
        private double mass;
        private String acids;

        public Solid(double mass, String acids) {
            this.mass = mass;
            this.acids = acids;
        }

        public double getMass() {
            return mass;
        }

        public String getAcids() {
            return acids;
        }
    }

    private static class PeakUnion {
        ArrayList<Double> values = new ArrayList<Double>();


        public double getAverage() {
            double sum = 0;
            for (Double value : values) {
                sum += value;
            }
            return sum / values.size();
        }

        public void add(double value) {
            values.add(value);
        }
    }

    public static PrintWriter createOutputFile(File msoutputDir, String fileName)
            throws UnsupportedEncodingException, FileNotFoundException {
        return new PrintWriter(new OutputStreamWriter(
                new FileOutputStream(new File(msoutputDir, fileName)), "UTF-8"));
    }
}
