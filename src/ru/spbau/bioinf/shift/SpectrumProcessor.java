package ru.spbau.bioinf.shift;

import org.apache.log4j.Logger;
import org.jdom.Document;
import org.jdom.Element;
import ru.spbau.bioinf.shift.util.ReaderUtil;
import ru.spbau.bioinf.shift.util.XmlUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Properties;

public class SpectrumProcessor {


    private static final Logger log = Logger.getLogger(SpectrumProcessor.class);

    private static final double EPSILON = 0.01;
    private File outputDir;
    private ArrayList<Protein> proteins;
    private List<Spectrum> spectrums;


    public static void main(String[] args) throws Exception {
        String dataset = "data/set8";
        String proteinDatabaseName = "prot.fasta";

        if (args.length > 0) {
            dataset = args[0];
        }

        if (args.length > 1) {
            proteinDatabaseName = args[1];
        }

        SpectrumProcessor processor = new SpectrumProcessor();


        processor.processResult(new File(dataset), proteinDatabaseName);
    }

    public void processResult(File datasetDir, String proteinDatabaseFilename) throws Exception {
        log.debug("star result processing");
        File msinputDir = new File(datasetDir, "msinput");
        BufferedReader input = ReaderUtil.getBufferedReader(new File(msinputDir, "input_data"));
        spectrums = new ArrayList<Spectrum>();
        outputDir = new File(datasetDir, "xml");

        new File(outputDir, "spectrums").mkdirs();

        Properties properties;

        while ((properties = ReaderUtil.readPropertiesUntil(input, "PRECURSOR_MASS")).size() > 0) {
            Spectrum spectrum = new Spectrum(properties, input);
            spectrums.add(spectrum);
        }
        log.debug("spectrums data loaded");
        FastaReaderString fastaReader = new FastaReaderString(new File(msinputDir, proteinDatabaseFilename));

        proteins = fastaReader.getProteins();

        log.debug("protein database loaded");


        int total = 0;

        for (Spectrum spectrum : spectrums) {
            total += processSpectrum(proteins, spectrum);

        }
        log.debug("total = " + total);
    }

    private void testSave(int spectrumId) throws IOException {
        Spectrum sp = spectrums.get(spectrumId);
        Protein p = proteins.get(13822);
        Map<Character, ArrayList<Double>> pos = getPositions(sp);


        initProteinShifts(p, pos, 0);
        Match m = getScore(p, sp);
        ArrayList<Match> test = new ArrayList<Match>();
        test.add(m);
        saveSpectrum(sp, test);
    }

    private int processSpectrum(ArrayList<Protein> proteins, Spectrum spectrum) throws IOException {
        int ans = 0;
        int bestScore = 4;
        ArrayList<Match> best = new ArrayList<Match>();
        Map<Character, ArrayList<Double>> positions = getPositions(spectrum);
        int max = 3;
        for (Protein protein : proteins) {
            max = Math.max(max, initProteinShifts(protein, positions, max));
        }

        log.debug("Max shift score is " + max);
        for (Protein protein : proteins) {
            Match match = getScore(protein, spectrum);
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
            saveSpectrum(spectrum,  best);
        } else {
            log.debug("Spectrum " + spectrum.getId() + " is too bad, best.size() is " +  best.size());
        }
        return ans;
    }

    private void saveSpectrum(Spectrum spectrum, ArrayList<Match> best) throws IOException {
        for (Match bestMatch : best) {
            log.debug("bestScore = " + bestMatch.getScore() + " bestProtein " + bestMatch.getProteinId());
        }
        Document doc = new Document();
        Element root = new Element("spectrum");
        doc.setRootElement(root);
        spectrum.addToXml(root, new HashMap<Integer, List<Break>>());
        Element matches = new Element("matches");
        for (Match match : best) {
            Protein p = proteins.get(match.getProteinId());
            ArrayList<Double> shifts = new ArrayList<Double>();
            for (int i = 0; i < p.shiftScores.size(); i++) {
                 if (p.shiftScores.get(i) > 5) {
                    shifts.add(p.shifts.get(i));
                }
            }
            shifts = removeSides(shifts);
            for (Double shift : shifts) {
                match.addShift(shift, 0);
            }

            matches.addContent(match.toXml());
        }
        root.addContent(matches);

        XmlUtil.saveXml(doc, outputDir + "/spectrums/spectrum" + spectrum.getId() + ".xml");
        log.debug("Match for spectrum " + spectrum.getId() + " saved");
    }


    private Match getScore(Protein protein, Spectrum spectrum) {
        int bestScore = 0;
        Match bestMatch = null;
        double[] sd = spectrum.getData();
        double[] pd = protein.getSpectrum();

        for (int c = 0; c < protein.shiftScores.size(); c++) {
            double shift = protein.shifts.get(c);
            int score = 0;
            int j = 0;
            for (int i = 0; i < pd.length; i++) {
                while (j < sd.length && sd[j] < pd[i] - shift - EPSILON) {
                    j++;
                }
                if (j < sd.length) {
                    if (Math.abs(sd[j] - pd[i] + shift) < EPSILON) {
                        score++;
                    }
                }
            }

            if (score > bestScore) {
                bestScore = score;
                bestMatch = new Match(shift, bestScore, protein);
            }
        }
        return bestMatch;
    }

    private int initProteinShifts(Protein protein, Map<Character, ArrayList<Double>> positions, int max) {
        protein.shifts.clear();
        protein.shiftScores.clear();
        ArrayList<Double> shifts = new ArrayList<Double>();
        String acids = protein.getSimplifiedAcids();
        double cur = 0;
        for (int i = 0; i < acids.length(); i++) {
            ArrayList<Double> p = positions.get(acids.charAt(i));
            if (p != null) {
                for (double d : p) {
                    shifts.add(cur - d);
                }
            }

            cur += Acids.acids.get(acids.charAt(i));
        }
        Collections.sort(shifts);
        shifts.add((double) Integer.MAX_VALUE);
        double prev = 0;
        int count = 0;
        double sum = 0;
        for (Double s : shifts) {
            if (s - prev < EPSILON) {
                count++;
                sum += s;
                prev = s;
                continue;
            }

            double shift = sum / count;
            if (count > 0) {
                protein.shifts.add(shift);
                protein.shiftScores.add(count);
            }
            if (count > max) {
                max = count;
            }
            count = 1;
            prev = s;
            sum = s;
        }
        return max;
    }


    private Map<Character, ArrayList<Double>> getPositions(Spectrum spectrum) {
        Map<Character, ArrayList<Double>> ans = new HashMap<Character, ArrayList<Double>>();
        double[] values = merge(spectrum.getData());
        for (int i = 0; i < values.length - 1; i++) {
            for (int j = i + 1; j < values.length; j++) {
                double diff = values[j] - values[i];
                if (diff > 300)
                    break;
                for (Map.Entry<Character, Double> entry : Acids.acids.entrySet()) {
                    if (Math.abs(entry.getValue() - diff) < EPSILON) {
                        Character key = entry.getKey();
                        if (!ans.containsKey(key)) {
                            ans.put(key, new ArrayList<Double>());
                        }
                        ans.get(key).add(values[i]);
                    }
                }
            }
        }
        return ans;
    }

    private double[] merge(double[] v) {
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
                if (count > 1) {
                    a.add(sum / count);
                }
                count = 1;
                prev = s;
                sum = s;
            }
        }
        a.add(sum / count);
        return convertToArray(a);
    }

    private ArrayList<Double> removeSides(ArrayList<Double> v) {
        int n = v.size();
        boolean[] filtered = new boolean[n];
        for (int i = 0; i < n -2; i++) {
            int plus1 = -1;
            int plus2 = -1;
             int ammoniaPos = -1;
            int waterPos = -1;
            for (int j = i+1; j < n; j++) {
                double delta = v.get(j) - v.get(i);
                if (delta > 18.5) {
                    break;
                }
                if (Math.abs(delta -1) < EPSILON) {
                    plus1 = j;
                }
                if (Math.abs(delta -2) < EPSILON) {
                    plus2 = j;
                }
                if (Math.abs(delta - Consts.AMMONIA) < EPSILON) {
                    ammoniaPos= j;
                }
                if (Math.abs(delta - Consts.WATER) < EPSILON) {
                    waterPos = j;
                }
            }
            if (plus1 > 0 && plus2 > 0) {
                filtered[i] = true;
                filtered[plus2] = true;
            }
            if (ammoniaPos > 0 && waterPos > 0) {
                filtered[ammoniaPos] = true;
                filtered[waterPos] = true;
            }
        }
        for (int i = n-1; i >0; i--) {
            int ammoniaPos = -1;
            int waterPos = -1;
            for (int j = i-1; j >=0; j--) {
                double delta = v.get(i) - v.get(j);
                if (delta > 19) {
                    break;
                }
                if (Math.abs(delta - Consts.AMMONIA) < EPSILON) {
                    ammoniaPos= j;
                }
                if (Math.abs(delta - Consts.WATER) < EPSILON) {
                    waterPos = j;
                }
            }
            if (ammoniaPos > 0 && waterPos > 0) {
                filtered[ammoniaPos] = true;
                filtered[waterPos] = true;
            }
        }

        ArrayList<Double> a = new ArrayList<Double>();
        for (int i = 0; i < n; i++) {
            if (!filtered[i]) {
                a.add(v.get(i));
            }
        }
        return a;
    }

    private double[] convertToArray(ArrayList<Double> a) {
        double[] ans = new double[a.size()];
        for (int i = 0; i < ans.length; i++) {
            ans[i] = a.get(i);
        }

        return ans;
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

    public static class Match {
        private double shift;
        private int score;
        private String protein;
        private int proteinId;

        ArrayList<Double> shifts = new ArrayList<Double>();
        ArrayList<Integer> shiftsScores = new ArrayList<Integer>();

        public Match(double shift, int score, Protein protein) {
            this.shift = shift;
            this.score = score;
            this.protein = protein.getAcids();
            this.proteinId = protein.getProteinId();
        }

        public void addShift(double v, int score) {
            shifts.add(v);
            shiftsScores.add(score);
        }
        public int getProteinId() {
            return proteinId;
        }


        public int getScore() {
            return score;
        }

        public Element toXml() {
            Element match = new Element("match");
            XmlUtil.addElement(match, "shift", shift);
            XmlUtil.addElement(match, "score", score);
            XmlUtil.addElement(match, "protein", protein);
            XmlUtil.addElement(match, "protein-id", proteinId);
            Element sh = new Element("shifts");
            match.addContent(sh);
            for (int i = 0; i < shifts.size(); i++) {
                Element shift = new Element("shift");
                XmlUtil.addElement(shift, "value", shifts.get(i));
                XmlUtil.addElement(shift, "score", shiftsScores.get(i));
                sh.addContent(shift);
            }
            return match;
        }
    }
}
