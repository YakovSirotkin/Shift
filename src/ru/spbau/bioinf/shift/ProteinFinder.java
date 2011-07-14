package ru.spbau.bioinf.shift;

import org.apache.log4j.Logger;
import ru.spbau.bioinf.shift.util.ReaderUtil;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ProteinFinder {


    private static final Logger log = Logger.getLogger(ProteinFinder.class);

    public static final double EPSILON = 0.01;
    private List<Protein> proteins;
    private List<Spectrum> spectrums;
    private PrintWriter matchFile;
    private Map<String,List<ProteinPosition>> index;

    private Configuration config;

    public List<Protein> getProteins() {
        return proteins;
    }

    public List<Spectrum> getSpectrums() {
        return spectrums;
    }

    public static void main(String[] args) throws Exception {
        Configuration config = new Configuration(args);
        ProteinFinder processor = new ProteinFinder(config);

        processor.processAll();
    }

    public ProteinFinder(Configuration config) throws Exception {
        this.config = config;
        log.debug("star result processing");
        spectrums = config.getSpectrums();

        log.debug("spectrums data loaded");

        proteins = config.getProteins();

        log.debug("protein database loaded");

        index = getIndex(proteins);

        log.debug("index loaded");
    }

    public void processAll() throws Exception {
        matchFile = ReaderUtil.createOutputFile(config.getMatchFile());
        int total = 0;

        for (Spectrum spectrum : spectrums) {
            total += processSpectrum(spectrum);
            spectrum.clearData();
        }
        matchFile.close();

        log.debug("total = " + total);
    }

    public Map<String, List<ProteinPosition>> getIndex(List<Protein> proteins) {
        Map<String, List<ProteinPosition>> ans = new HashMap<String, List<ProteinPosition>>();
        for (Protein protein : proteins) {
            double[] pd = protein.getSpectrum();
            String acids = protein.getSimplifiedAcids();
            int proteinId = protein.getProteinId();
            for (int i = 0; i < acids.length()- 2; i++) {
                String key = acids.substring(i, i+3);
                if (!ans.containsKey(key)) {
                    ans.put(key, new ArrayList<ProteinPosition>());
                }
                ans.get(key).add(new ProteinPosition(proteinId, pd[i]));
            }
        }
        return ans;
    }

    private int processSpectrum(Spectrum spectrum) throws IOException {
        int ans = 0;
        ArrayList<Match> best = getSpectrumMatches(spectrum);
        if (best.size() > 0 && best.size() < 20) {
            ans = 1;
            for (Match match : best) {
                matchFile.println(spectrum.getId() + " " + match.getProteinId() + " " + match.getScore());
                matchFile.flush();
            }
            log.debug("Spectrum " + spectrum.getId() + " save with " + best.size() + " answers.");
        } else {
            log.debug("Spectrum " + spectrum.getId() + " is too bad, best.size() is " +  best.size());
        }
        return ans;
    }

    public ArrayList<Match> getSpectrumMatches(Spectrum spectrum) {
        int bestScore = 0;
        ArrayList<Match> best = new ArrayList<Match>();
        Map<String,List<Double>> positions = getPositions(spectrum);

        for (Map.Entry<String, List<Double>> entry : positions.entrySet()) {
            String key = entry.getKey();
            List<Double> values = entry.getValue();
            List<ProteinPosition> pp = index.get(key);

            if (pp!= null) {
                for (ProteinPosition pos : pp) {
                    int proteinId = pos.getProteinId();
                    for (double value : values) {
                        int score = getScore(spectrum.getData(), proteins.get(proteinId).getSpectrum(), pos.getPos() - value);
                            if (score >= bestScore) {
                                if (score > bestScore) {
                                    best.clear();
                                }
                                boolean contains = false;
                                for (Match matchOld : best) {
                                    if (matchOld.getProteinId() == proteinId) {
                                        contains = true;
                                    }
                                }
                                if (!contains) {
                                    best.add(new Match(proteins.get(proteinId), score));
                                }
                                bestScore = score;
                            }
                    }
                }
            }
        }
        return best;
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
        int i = 0;
        int j = 0;
        do {
            double diff = sd[j] - pd[i] + shift;
            if (diff < -0.1) {
                j++;
            } else if (diff > 0.1) {
                i++;
            } else {
                score++;
                i++;
                j++;
            }
        } while (i < pd.length && j < sd.length);
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
        if (v.size() == 0) {
            return v;
        }
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
}
