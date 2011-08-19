package ru.spbau.bioinf.shift;

import ru.spbau.bioinf.shift.util.ReaderUtil;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MSAlignDiff {
    public static void main(String[] args) throws Exception {
        Configuration config = new Configuration(args);

        ProteinFinder finder = new ProteinFinder(config);
        Map<Integer, Spectrum> spectrums = finder.getSpectrums();
        List<Protein> proteins = finder.getProteins();

        BufferedReader ms = ReaderUtil.getBufferedReader(new File(config.getMatchFile().getParent(), "ms_0001.txt"));
        BufferedReader sh = ReaderUtil.getBufferedReader(new File(config.getMatchFile().getParent(), "match02.txt"));


        List<MsMatch> msHits = new ArrayList<MsMatch>();
        Map<Integer, Integer> msStat = new HashMap<Integer, Integer>();
        String s;
        while ((s = ms.readLine()) != null) {
            MsMatch m = new MsMatch(s);
            msHits.add(m);
            int proteinId = m.proteinId;
            if (!msStat.containsKey(proteinId)) {
                msStat.put(proteinId, 0);
            }
            msStat.put(proteinId, msStat.get(proteinId) + 1);
        }

        Map<Integer, Integer> sHits = new HashMap<Integer, Integer>();
        while ((s = sh.readLine()) != null) {
            String[] data = s.split(" ");
            int spectrumId = Integer.parseInt(data[0]);
            int proteinId = Integer.parseInt(data[1]);
            sHits.put(spectrumId, proteinId);
        }

        System.out.println(msHits.size() + " matches are reported by MS-Align.");

        int same = 0;
        int diff = 0;

        List<MsMatch> wrong = new ArrayList<MsMatch>();
        for (MsMatch msHit : msHits) {
            if (new Integer(msHit.proteinId).equals(sHits.get(msHit.spectrumId))) {
                same++;
            } else {
                diff++;
                wrong.add(msHit);
            }
        }

        System.out.println(same + " matches are the same.");
        System.out.println(diff + " matches are different in shift.");

        ScoringFunction scoringFunction = config.getScoringFunction();

        Collections.sort(wrong);

        Map<String, List<ProteinPosition>> index = finder.getIndex(proteins);

        for (MsMatch m : wrong) {
            int spectrumId = m.spectrumId;
            int proteinId = m.proteinId;
            System.out.println(spectrumId + " " + proteinId + " " + m.eValue + " http://msalign.ucsd.edu/set8-7/html/prsms/prsm" + spectrumId + ".html");
            System.out.println(msStat.get(proteinId) + " matches for this protein were reported by MS-align.");
            Spectrum spectrum = spectrums.get(spectrumId);
            Map<String,List<Double>> positions = ProteinFinder.getPositions(spectrum);
            System.out.println(positions.size() + " tags were discovered in the spectrum.");
            if (positions.size() > 0) {
                Protein protein = proteins.get(proteinId);
                List<Double> shiftsList = ProteinFinder.getShifts(protein, positions);
                double bestScore = 0;
                for (double shift : shiftsList) {
                    double  score = scoringFunction.getScore(spectrum, protein, shift);
                    if (score > bestScore) {
                        bestScore = score;
                    }
                }
                if (bestScore == 0) {
                    System.out.println("Protein does not contains tags from this spectrum.");
                } else {
                    Map<Integer, Double> scores = new HashMap<Integer, Double>();
                    double maxScore = bestScore;
                    for (Map.Entry<String, List<Double>> entry : positions.entrySet()) {
                        String key = entry.getKey();
                        List<Double> values = entry.getValue();
                        List<ProteinPosition> pp = index.get(key);
                        if (pp != null) {
                            for (ProteinPosition pos : pp) {
                                int pId = pos.getProteinId();
                                for (double value : values) {
                                    double score = scoringFunction.getScore(spectrum, proteins.get(pId), pos.getPos() - value);
                                    if (score >= bestScore) {
                                        if (score > maxScore) {
                                            maxScore = score;
                                        }
                                        if (scores.containsKey(pId)) {
                                            if (scores.get(pId) < score) {
                                                scores.put(pId, score);
                                            }
                                        } else {
                                            scores.put(pId, score);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    List<SMatch> others = new ArrayList<SMatch>();
                    for (Map.Entry<Integer, Double> entry : scores.entrySet()) {
                        others.add(new SMatch(entry));
                    }
                    Collections.sort(others);
                    System.out.println(bestScore + " is the score in Shift. " +others.size() + " proteins aren't worse than that:");
                    int c = 0;
                    for (SMatch other : others) {
                        c++;
                        if (c > 5) {
                            break;
                        }
                        System.out.println(other.proteinId + " " + other.score);
                    }
                }
                spectrum.clearData();
            }
            System.out.println();
        }
        Collections.sort(msHits);
        int limit = 100 * 100;
        int hl = limit/2;
        int[] statB = new int[limit];
        int[] statY = new int[limit];

        for (MsMatch msHit : msHits) {
            int spectrumId = msHit.spectrumId;
            int proteinId = msHit.proteinId;
            if (msHit.eValue < 0.000001 && new Integer(proteinId).equals(sHits.get(spectrumId))) {
                Protein protein = proteins.get(proteinId);
                SpectrumProteinMatch match = new SpectrumProteinMatch(spectrums.get(spectrumId), protein, scoringFunction);
                double bestShift = match.getBestShift();
                List<Double> sharedPeaks = new ArrayList<Double>();
                Spectrum spectrum = spectrums.get(spectrumId);
                double[] pd = protein.getSpectrum();
                for (double v : pd) {
                    sharedPeaks.add(v-bestShift);
                }
                spectrum.clearData();
                List<Peak> peaks = spectrum.getPeaks();
                for (Peak peak : peaks) {
                    double m = peak.getMonoisotopicMass();
                    double bestDiff = IonStatisticsEngine.getBestDiff(sharedPeaks, m);
                    int pos = hl + (int)Math.round(bestDiff * 100);
                    if (pos > 0 && pos < limit) {
                        statB[pos]++;
                    }
                    bestDiff = IonStatisticsEngine.getBestDiff(sharedPeaks, spectrum.getPrecursorMass() - m - Consts.WATER);
                    pos = hl + (int)Math.round(bestDiff * 100);
                    if (pos > 0 && pos < limit) {
                        statY[pos]++;
                    }
                }
            }
        }
        for (int i = 0; i < limit; i++) {
            if (statB[i] > 10) {
                System.out.println(((i-hl)/100d) + " " + statB[i]);
            }
        }
    }

    public static class MsMatch implements Comparable<MsMatch> {
        int proteinId;
        int spectrumId;
        double eValue;

        public MsMatch(String s) {
            String[] data = s.split(" ");
            spectrumId = Integer.parseInt(data[0]);
            proteinId = Integer.parseInt(data[1]);
            eValue = Double.parseDouble(data[2]);
        }

        public int compareTo(MsMatch o) {
            double diff = eValue - o.eValue;
            if (diff < 0) {
                return -1;
            }
            if (diff > 0) {
                return 1;
            }
            return 0;
        }
    }

    public static class SMatch implements Comparable<SMatch> {
        int proteinId;
        double score;

        public SMatch(Map.Entry<Integer, Double> entry) {
            this.proteinId = entry.getKey();
            this.score = entry.getValue();
        }

        public int compareTo(SMatch o) {
            double diff = o.score - score;
            if (diff < 0) {
                return -1;
            }
            if (diff > 0) {
                return 1;
            }
            diff = proteinId - o.proteinId;
            if (diff < 0) {
                return -1;
            }
            if (diff > 0) {
                return 1;
            }
            return 0;
        }
    }

}
